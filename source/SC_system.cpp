#define PROFILE

//openFHE
#include "openfhe.h"
//Serialization
#include "ciphertext-ser.h"
#include "cryptocontext-ser.h"
#include "key/key-ser.h"
#include "scheme/ckksrns/ckksrns-ser.h"

//normal libs
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <algorithm>
#include <typeinfo>
#include <chrono>


const std::string DATAFOLDER = "../../data"; //foulder
std::string ccLocation       = "/cryptocontext.txt"; //context
std::string pubKeyLocation   = "/key_pub.txt";   // Pub key
std::string multKeyLocation  = "/key_mult.txt";  // relinearization key
std::string rotKeyLocation   = "/key_rot.txt";   // automorphism / rotation key
std::string cipherOneLocation = "/ciphertext/CT"; //encrypted data
std::string cipherMultLocation   = "/ciphertextMult.txt"; //result from mul
std::string cipherAddLocation    = "/ciphertextAdd.txt"; //result from ADD

std::vector<std::string> split(const std::string &s, char delim) {
    std::stringstream ss(s);
    std::string item;
    std::vector<std::string> tokens;
    while (getline(ss, item, delim)) {
        tokens.push_back(item);
    }
    return tokens;
}

double linear(double x){
    return x;
}

double relu(double x){
    double res=0;
    if(x<0){
        res=0;
    }
    else{
        res=x;
    }
    return res;
}

double sigmoid(double x){
    return 1 / (1 + exp(-x));
}

std::vector <std::complex<double>> operator/(const std::vector <std::complex<double>>& m2, const float m1){
    
    /*  Returns the product of a float and a vectors (elementwise multiplication).
     Inputs:
     m1: float
     m2: vector
     Output: vector, m1 * m2, product of two vectors m1 and m2
     */
    
    const unsigned long VECTOR_SIZE = m2.size();
    std::vector <std::complex<double>> product (VECTOR_SIZE);
    
    for (unsigned i = 0; i != VECTOR_SIZE; ++i){
        product[i] = m2[i] * std::complex<double>{1/m1};
    };
    
    return product;
}


class User{
    public:
        lbcrypto::CryptoContext<lbcrypto::DCRTPoly> UserCC; //crypto context
        User(int multDepth, int scaleModSize,int batchSize){ //constructor will generate context and keys
            lbcrypto::CCParams<lbcrypto::CryptoContextCKKSRNS> parameters;
            parameters.SetMultiplicativeDepth(multDepth);
            #if NATIVEINT == 128
                usint scalingModSize = 78;
                usint firstModSize   = 89;
            #else
                usint scalingModSize = 59;
                usint firstModSize   = 60;
            #endif
            parameters.SetScalingModSize(scalingModSize);
            parameters.SetFirstModSize(firstModSize);
            parameters.SetBatchSize(batchSize);
            parameters.SetSecurityLevel(lbcrypto::HEStd_NotSet);
            parameters.SetRingDim(32768);
            UserCC = GenCryptoContext(parameters);
            UserCC->Enable(PKE);
            UserCC->Enable(KEYSWITCH);
            UserCC->Enable(LEVELEDSHE);
            UserCC->Enable(ADVANCEDSHE);

            std::cout << "Cryptocontext generated" << std::endl;

            UserKP = UserCC->KeyGen();
            UserCC->EvalMultKeyGen(UserKP.secretKey);
            UserCC->EvalRotateKeyGen(UserKP.secretKey, {1, 2, -1, -2});
            std::cout<<"Keys generated" << std::endl;
            std::cout<<"Ring size:"<<UserCC->GetRingDimension()<<std::endl;
        }
        void encript_send(std::vector<std::complex<double>> X_train){ //will encript and serialize data
            lbcrypto::Ciphertext<lbcrypto::DCRTPoly> CipherMat[784];
            for(int i=0;i<784;i++){
                std::vector<std::complex<double>> x_send;
                for(int j=0;j<14000;j++){
                    x_send.push_back(X_train[i+j*784]);
                }
                lbcrypto::Plaintext PlainX = UserCC->MakeCKKSPackedPlaintext(x_send);
                CipherMat[i] = UserCC->Encrypt(UserKP.publicKey, PlainX);
            }
            //std::string number=std::to_string(i);
            if (!lbcrypto::Serial::SerializeToFile(DATAFOLDER + cipherOneLocation+".txt", CipherMat, lbcrypto::SerType::BINARY)) {
                std::cerr << " Error writing ciphertext 1" << std::endl;
                std::exit(1);
            }
            std::cout<<"Encripted done"<<std::endl;

        }
        void send_CKKS(){ //will send CKKS values, keys and context(dont know if the context matters)
            if (!lbcrypto::Serial::SerializeToFile(DATAFOLDER + ccLocation, UserCC, lbcrypto::SerType::BINARY)) {
                std::cerr << "Error writing serialization of the crypto context to ""cryptocontext.txt"<< std::endl;
                std::exit(1);
            }
            std::ofstream multKeyFile(DATAFOLDER + multKeyLocation, std::ios::out | std::ios::binary);
            if (multKeyFile.is_open()) {
                if (!UserCC->SerializeEvalMultKey(multKeyFile, lbcrypto::SerType::BINARY)) {
                    std::cerr << "Error writing eval mult keys" << std::endl;
                    std::exit(1);
                }
                std::cout << "EvalMult/ relinearization keys have been serialized" << std::endl;
                multKeyFile.close();
            }
            else {
                std::cerr << "Error serializing EvalMult keys" << std::endl;
                std::exit(1);
            }

            std::ofstream rotationKeyFile(DATAFOLDER + rotKeyLocation, std::ios::out | std::ios::binary);
            if (rotationKeyFile.is_open()) {
                if (!UserCC->SerializeEvalAutomorphismKey(rotationKeyFile, lbcrypto::SerType::BINARY)) {
                    std::cerr << "Error writing rotation keys" << std::endl;
                    std::exit(1);
                }
                std::cout << "Rotation keys have been serialized" << std::endl;
            }
            else {
                std::cerr << "Error serializing Rotation keys" << std::endl;
                std::exit(1);
            }
        }
        std::vector<double> decript(std::vector<double> &vec_res){ //will deserialize and decript result to vec_res
            lbcrypto::Ciphertext<lbcrypto::DCRTPoly> res_cipher_mat[10];
            lbcrypto::Ciphertext<lbcrypto::DCRTPoly> aux;
            lbcrypto::Plaintext res_plaintext_mat;
            std::vector<double> vecadd[10];
            lbcrypto::Serial::DeserializeFromFile(DATAFOLDER + "/Result.txt", res_cipher_mat, lbcrypto::SerType::BINARY);
            std::cout<<"Deserialized the inference"<<std::endl;
            for(int i=0;i<10;i++){
                aux=res_cipher_mat[i];
                UserCC->Decrypt(UserKP.secretKey, aux, &res_plaintext_mat);
                vecadd[i]=res_plaintext_mat->GetRealPackedValue();
            }
            std::cout<<"Descripted inference"<<std::endl;
            for(int j=0;j<14000;j++){
                for(int i=0;i<10;i++){
                    vec_res.push_back(vecadd[i][j]);
                    //std::cout<<vec_res[i+j*128]<<" ";
                }
            }
            std::cout << "Estimated precision in bits: " << res_plaintext_mat->GetLogPrecision() << std::endl;
            return vec_res;
        }
        ~User(){} //destructor
    private:
        lbcrypto::KeyPair<lbcrypto::DCRTPoly> UserKP; //KEYS ARE PRIVATE, NN WILL NOT HAVE ACESS TO THIS
};


void clientProcess(int sig,char pres) { //NN process
    lbcrypto::CryptoContext<lbcrypto::DCRTPoly> clientCC;

    //reset crypto context for NN
    clientCC->ClearEvalMultKeys();
    clientCC->ClearEvalAutomorphismKeys();
    lbcrypto::CryptoContextFactory<lbcrypto::DCRTPoly>::ReleaseAllContexts();

    //Deserialization of data
    if (!lbcrypto::Serial::DeserializeFromFile(DATAFOLDER + ccLocation, clientCC, lbcrypto::SerType::BINARY)) {
        std::cerr << "I cannot read serialized data from: " << DATAFOLDER << "/cryptocontext.txt" << std::endl;
        std::exit(1);
    }
    std::cout << "Client CC deserialized"<<std::endl;

    
    std::ifstream multKeyIStream(DATAFOLDER + multKeyLocation, std::ios::in | std::ios::binary);
    if (!multKeyIStream.is_open()) {
        std::cerr << "Cannot read serialization from " << DATAFOLDER + multKeyLocation << std::endl;
        std::exit(1);
    }
    if (!clientCC->DeserializeEvalMultKey(multKeyIStream, lbcrypto::SerType::BINARY)) {
        std::cerr << "Could not deserialize eval mult key file" << std::endl;
        std::exit(1);
    }

    std::cout << "Deserialized eval mult keys"<< std::endl;
    std::ifstream rotKeyIStream(DATAFOLDER + rotKeyLocation, std::ios::in | std::ios::binary);
    if (!rotKeyIStream.is_open()) {
        std::cerr << "Cannot read serialization from " << DATAFOLDER + multKeyLocation << std::endl;
        std::exit(1);
    }
    if (!clientCC->DeserializeEvalAutomorphismKey(rotKeyIStream, lbcrypto::SerType::BINARY)) {
        std::cerr << "Could not deserialize eval rot key file" << std::endl;
        std::exit(1);
    }

    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> EncData[784];

    if (!lbcrypto::Serial::DeserializeFromFile(DATAFOLDER + cipherOneLocation+".txt", EncData, lbcrypto::SerType::BINARY)) {
        std::cerr << "Cannot read serialization from " << DATAFOLDER + cipherOneLocation << std::endl;
    }
    std::cout << "Deserialized Input" << '\n' << std::endl;

    //GET WEIGHTS:
    std::ifstream in_weights;
	in_weights.open("../../data/Weights.txt"); 
    std::vector<float> W1;
    std::vector<float> W2;
    std::vector<float> W3;
    int W1_size=128;
    int W2_size=64;
    float value;
    int i,k;
    std::cout << "read W1\n";
    int n = 784*W1_size;
	for (i = 0; i < (n>>2)<<2; i+=8) {
        //in_weights.read(reinterpret_cast<char*>(&value), sizeof(value));
        in_weights >> value;
		W1.push_back(value);
        in_weights >> value;
		W1.push_back(value);
        in_weights >> value;
		W1.push_back(value);
        in_weights >> value;
		W1.push_back(value);
        in_weights >> value;
		W1.push_back(value);
        in_weights >> value;
		W1.push_back(value);
        in_weights >> value;
		W1.push_back(value);
        in_weights >> value;
		W1.push_back(value);
	}

    for (i=(n>>2)<<2; i < n; i++){
        in_weights >> value;
		W1.push_back(value);
    }
	std::cout<< "read W2\n";
    n=W1_size*W2_size;
	for (i = 0; i < (n>>2)<<2; i+=4) {
        //in_weights.read(reinterpret_cast<char*>(&value), sizeof(value));
        in_weights >> value;
		W2.push_back(value);
        in_weights >> value;
		W2.push_back(value);
        in_weights >> value;
		W2.push_back(value);
        in_weights >> value;
		W2.push_back(value);
	}
    for (i=(n>>2)<<2; i < n; i++){
        in_weights >> value;
		W2.push_back(value);
    }
	std::cout << "read W3\n";
	for (k = 0; k < (W2_size*10); k++) {
        //in_weights.read(reinterpret_cast<char*>(&value), sizeof(value));
        in_weights >> value;
		W3.push_back(value);
    }

    //-------------------------------------2-----------------------------------
    std::cout<<"Started doing input layer"<<std::endl;
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> AfterW1[128];
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> AfterW2[64];
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> AfterW3[10];
    double lowerBound1=-36;    //relu=-32;                 sigmoid=-30    calculated -28.8796; //-23
    double upperBound1=32;      //relu=32;                  sigmoid=26     calculated 24.9469; //22
    double lowerBound2=-20; //relu=-245;                sigmoid=-12     calculated -11.0627;
    double upperBound2=15;  //relu=110;                 sigmoid=11      calculated 10.6347;
    uint32_t polydegree=15;
    if(sig==1){
        if(pres=='L') polydegree=15;
        else if(pres=='M') polydegree=50;
        else  polydegree=70;
        std::cout<<"We are using the polydegree "<<polydegree<<std::endl;
    }
    for (int i=0;i<128;i++){
        AfterW1[i] = clientCC->EvalMult(EncData[0],W1[i]);
        /*
        for ( int j=1; j < (784>>2)<<2; j+=6){
            auto sum = clientCC->EvalMult(EncData[j],W1[i+j*128]);
            AfterW1[i] = clientCC->EvalAdd(AfterW1[i], sum);
            sum=clientCC->EvalMult(EncData[j],W1[i+(j+1)*128]) ;
            AfterW1[i] = clientCC->EvalAdd(AfterW1[i], sum);
            sum=clientCC->EvalMult(EncData[j],W1[i+(j+2)*128]) ;
            AfterW1[i] = clientCC->EvalAdd(AfterW1[i], sum);
            sum=clientCC->EvalMult(EncData[j],W1[i+(j+3)*128]) ;
            AfterW1[i] = clientCC->EvalAdd(AfterW1[i], sum);
            sum=clientCC->EvalMult(EncData[j],W1[i+(j+4)*128]) ;
            AfterW1[i] = clientCC->EvalAdd(AfterW1[i], sum);
            sum=clientCC->EvalMult(EncData[j],W1[i+(j+5)*128]) ;
            AfterW1[i] = clientCC->EvalAdd(AfterW1[i], sum);
        }
        for (int j=(784>>2)<<2;j<784;j++){
            auto sum = clientCC->EvalMult(EncData[j],W1[i+j*128]);
            AfterW1[i] = clientCC->EvalAdd(AfterW1[i], sum);
        }
        */
        for(int j=1;j<784;j++){
            auto sum=clientCC->EvalMult(EncData[j],W1[i+j*128]);
            AfterW1[i] = clientCC->EvalAdd(AfterW1[i], sum);
        }
        if(sig==1) AfterW1[i]=clientCC->EvalChebyshevFunction([](double x)->double{return sigmoid(x);},AfterW1[i],lowerBound1,upperBound1,polydegree);
        std::cout<<"finished neuron "<<i<<" out of "<<128<<std::endl;
    }

    for (int i=0; i < 64;i++){
        AfterW2[i] = clientCC->EvalMult(AfterW1[0],W2[i]);
        for ( int j=1; j < 128; j++){
            auto sum = clientCC->EvalMult(AfterW1[j],W2[i+j*64]);
            AfterW2[i] = clientCC->EvalAdd(AfterW2[i], sum);
        }
        if(sig==1) AfterW2[i]=clientCC->EvalChebyshevFunction([](double x)->double{return sigmoid(x);},AfterW2[i],lowerBound2,upperBound2,polydegree);
        std::cout<<"finished neuron "<<i<<" out of "<<64<<std::endl;
    }

    for (int i=0; i < 10;i++){
        AfterW3[i] = clientCC->EvalMult(AfterW2[0],W3[i]);
        for ( int j=1; j < 64; j++){
            auto sum = clientCC->EvalMult(AfterW2[j],W3[i+j*10]);
            AfterW3[i] = clientCC->EvalAdd(AfterW3[i], sum);
        }
        std::cout<<"finished neuron "<<i<<" out of "<<10<<std::endl;
    }
    
    lbcrypto::Serial::SerializeToFile(DATAFOLDER + "/Result.txt", AfterW3, lbcrypto::SerType::BINARY);

    //----------------------------------2------------------------------------

    //serialization results

    std::cout << "Serialized all ciphertexts from client" << '\n' << std::endl;
}



int main(){
    
    //user inputs
    char conf;
    char pres;
    int sig=0;
    std::cout<<"Please insert activation function:";
    std::cin>>conf;
    while(conf!='L' and conf!='S'){
        std::cout<<"\n";
        std::cout<<"Please insert a valid input, either S for sigmoid or L for linear: ";
        std::cin>>conf;
    }
    std::cout<<"\n";
    if(conf=='S'){
        sig=1;
        std::cout<<"Please insert desired presision:";
        std::cin>>pres;
        while(pres!='L' and pres!='M' and pres!='H'){
            std::cout<<"\n";
            std::cout<<"Please insert a valid input, either L for low, M for medium or H for high: ";
            std::cin>>pres;
        }
        std::cout<<"\n";
    }
    else pres='N';





    //declare CKKS variables, will likely need tinckering
    uint32_t multDepth = 4; //18
    if(pres=='L' && sig==1)multDepth=16;
    else if(pres=='M' && sig==1)multDepth=18;
    else if(pres=='H' && sig==1)multDepth=20;
    uint32_t scaleModSize = 78;
    uint32_t batchSize = 16384;
    

    //--------------------------------1----------------------------------------------------------------//
    auto start = std::chrono::high_resolution_clock::now();
    User NNclient(multDepth,scaleModSize,batchSize);
    NNclient.send_CKKS();
    std::string line;
    std::vector<std::string> line_v;

    std::cout << "Loading data Test ...\n";
    std::vector<std::complex<double>> X_train;
    std::vector<float> y_train;
    std::ifstream myfile (DATAFOLDER+"/train.txt");
    float hit=0;
    float it=0;
    if (myfile.is_open())
    {
        
        int line_count = 28000; // Posiciona o pointer em 28000 para de seguida ler restante 1/3
    
        for (int i = 0; i < ((line_count>>2)<<2); i+=4)  {
            getline(myfile, line);
            getline(myfile, line);
            getline(myfile, line);
            getline(myfile, line);
        }
        for (int i=((line_count>>2)<<2); i < line_count; i++) getline(myfile, line);
       
        while ( getline (myfile,line) )
        {
            line_v = split(line, '\t');
            int digit = strtof((line_v[0]).c_str(),0);
            for (int i = 0; i < 10; ++i) {
                if (i == digit)
                {
                    y_train.push_back(1.);
                }
                else y_train.push_back(0.);
            }
            
            int size = static_cast<int>(line_v.size());
            for (int i = 1; i < int(size); ++i) {
                X_train.push_back(strtof((line_v[i]).c_str(),0));
            }
        }
        X_train = X_train/255.0;
        myfile.close();
    }
    NNclient.encript_send(X_train);
    clientProcess(sig,pres);
    std::vector<double> vec_res;
    NNclient.decript(vec_res);


    
    //check results



    float max=-10;
    float min=10;
        for(int i=0;i<int(vec_res.size());i++){
            if (vec_res[i]<min){
                min=vec_res[i];
            }
            if (vec_res[i]>max){
                max=vec_res[i];
            }
        }
        std::cout<<"hello "<< min<<" "<<max;
    int max_ind;
    for(int i=0;i<14000;i++){
        it++;
        max=0;
        max_ind=0;
        for(int j=0;j<10;j++){
            if (vec_res[j+i*10]>max){
                max=vec_res[j+i*10];
                max_ind=j;
            }
        }
        if (y_train[max_ind+i*10]==1){
            hit++;
        }
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout<<"time: "<<double(duration.count()/60000000)<<" minutes"<<std::endl;
    std::cout<<"hit: "<<hit<<std::endl;
    std::cout<<"total: "<<it<<std::endl;
    std::cout<<"Final Accuracy:"<<float(hit/it)<<std::endl;


     
    //---------------------------------1----------------------------------------------------------//

    


    return 0;
}