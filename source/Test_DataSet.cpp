#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <algorithm>


using namespace std;

void print ( const vector <float>& m, int n_rows, int n_columns ) {
    
    /*  "Couts" the input vector as n_rows x n_columns matrix.
     Inputs:
     m: vector, matrix of size n_rows x n_columns
     n_rows: int, number of rows in the left matrix m1
     n_columns: int, number of columns in the left matrix m1
     */
    
    for( int i = 0; i != n_rows; ++i ) {
        for( int j = 0; j != n_columns; ++j ) {
            cout << m[ i * n_columns + j ] << " ";
        }
        cout << '\n';
    }
    cout << endl;
}

vector <float> sigmoid(const vector <float>& m1) {

	/*  Returns the value of the sigmoid function f(x) = 1/(1 + e^-x).
	 Input: m1, a vector.
	 Output: 1/(1 + e^-x) for every element of the input matrix m1.
	 */

	const unsigned long VECTOR_SIZE = m1.size();
	vector <float> output(VECTOR_SIZE);


	for (unsigned i = 0; i != VECTOR_SIZE; ++i) {
		output[i] = 1 / (1 + exp(-m1[i]));
	}

	return output;
}

vector <float> transpose (float *m, const int C, const int R) {
    
    /*  Returns a transpose matrix of input matrix.
     Inputs:
     m: vector, input matrix
     C: int, number of columns in the input matrix
     R: int, number of rows in the input matrix
     Output: vector, transpose matrix mT of input matrix m
     */
    
    vector <float> mT (C*R);
    
    for(unsigned n = 0; n != C*R; n++) {
        unsigned i = n/C;
        unsigned j = n%C;
        mT[n] = m[R*j + i];
    }
    
    return mT;
}
vector <float> softmax (const vector <float>& z, const int dim) {
    
    const int zsize = static_cast<int>(z.size());
    vector <float> out;
    
    for (unsigned i = 0; i != zsize; i += dim) {
        vector <float> foo;
        for (unsigned j = 0; j != dim; ++j) {
            foo.push_back(z[i + j]);
        }
        
        float max_foo = *max_element(foo.begin(), foo.end());

        for (unsigned j = 0; j != dim; ++j) {
            foo[j] = exp(foo[j] - max_foo);
        }      

        float sum_of_elems = 0.0;
        for (unsigned j = 0; j != dim; ++j) {
            sum_of_elems = sum_of_elems + foo[j];
        }
        
        for (unsigned j = 0; j != dim; ++j) {
            out.push_back(foo[j]/sum_of_elems);
        }
    }
    return out;
}
vector <float> relu(const vector <float>& z){
    int size = z.size();
    vector <float> output;
    for( int i = 0; i < size; ++i ) {
        if (z[i] < 0){
            output.push_back(0.0);
        }
        else output.push_back(z[i]);
    }
    return output;
}
vector <float> dot (const vector <float>& m1, const vector <float>& m2, const int m1_rows, const int m1_columns, const int m2_columns) {
    
    /*  Returns the product of two matrices: m1 x m2.
     Inputs:
     m1: vector, left matrix of size m1_rows x m1_columns
     m2: vector, right matrix of size m1_columns x m2_columns (the number of rows in the right matrix
     must be equal to the number of the columns in the left one)
     m1_rows: int, number of rows in the left matrix m1
     m1_columns: int, number of columns in the left matrix m1
     m2_columns: int, number of columns in the right matrix m2
     Output: vector, m1 * m2, product of two vectors m1 and m2, a matrix of size m1_rows x m2_columns
     */
    
    vector <float> output (m1_rows*m2_columns);
    
    for( int row = 0; row != m1_rows; ++row ) {
        for( int col = 0; col != m2_columns; ++col ) {
            output[ row * m2_columns + col ] = 0.f;
            for( int k = 0; k != m1_columns; ++k ) {
                output[ row * m2_columns + col ] += m1[ row * m1_columns + k ] * m2[ k * m2_columns + col ];
            }
        }
    }
    
    return output;
}
vector<string> split(const string &s, char delim) {
    stringstream ss(s);
    string item;
    vector<string> tokens;
    while (getline(ss, item, delim)) {
        tokens.push_back(item);
    }
    return tokens;
}
vector <float> operator/(const vector <float>& m2, const float m1){
    
    /*  Returns the product of a float and a vectors (elementwise multiplication).
     Inputs:
     m1: float
     m2: vector
     Output: vector, m1 * m2, product of two vectors m1 and m2
     */
    
    const unsigned long VECTOR_SIZE = m2.size();
    vector <float> product (VECTOR_SIZE);
    
    for (unsigned i = 0; i != VECTOR_SIZE; ++i){
        product[i] = m2[i] / m1;
    };
    
    return product;
}
vector <float> operator-(const vector <float>& m1, const vector <float>& m2){
    
    /*  Returns the difference between two vectors.
     Inputs:
     m1: vector
     m2: vector
     Output: vector, m1 - m2, difference between two vectors m1 and m2.
     */
    
    const unsigned long VECTOR_SIZE = m1.size();
    vector <float> difference (VECTOR_SIZE);
    
    for (unsigned i = 0; i != VECTOR_SIZE; ++i){
        difference[i] = m1[i] - m2[i];
    };
    
    return difference;
}

vector <float> operator*(const vector <float>& m1, const vector <float>& m2){
    
    /*  Returns the product of two vectors (elementwise multiplication).
     Inputs:
     m1: vector
     m2: vector
     Output: vector, m1 * m2, product of two vectors m1 and m2
     */
    
    const unsigned long VECTOR_SIZE = m1.size();
    vector <float> product (VECTOR_SIZE);
    
    for (unsigned i = 0; i != VECTOR_SIZE; ++i){
        product[i] = m1[i] * m2[i];
    };
    
    return product;
}

int main(int argc, const char * argv[]) {
    char conf;
    int sig=0;
    std::cout<<"Please insert activation function:";
    std::cin>>conf;
    int flag=0;
    while(flag==0){
        if(conf!='L' and conf!='S'){
            std::cout<<"\n";
            std::cout<<"Please insert a valid input, either S for sigmoid or L for linear: ";
            std::cin>>conf;
        }
        else flag=1;
    }
    std::cout<<"\n";
    if(conf=='S') sig=1;
    
    string line;
    vector<string> line_v;

    cout << "Loading data Test ...\n";
    vector<float> X_train;
    vector<float> y_train;
    ifstream myfile ("../data/train.txt");
    
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
            for (unsigned i = 0; i < 10; ++i) {
                if (i == digit)
                {
                    y_train.push_back(1.);
                }
                else y_train.push_back(0.);
            }
            
            int size = static_cast<int>(line_v.size());
            for (unsigned i = 1; i < size; ++i) {
                X_train.push_back(strtof((line_v[i]).c_str(),0));
            }
        }
        X_train = X_train/255.0;
        myfile.close();
    }
    
    else cout << "Unable to open file" << '\n';
    
    int xsize = static_cast<int>(X_train.size());
    int ysize = static_cast<int>(y_train.size());
    
    // Some hyperparameters for the NN
    int BATCH_SIZE = 14000;


// Leitura dos pesos da NN
cout << "Loading weights ...\n";
    ifstream in_weights;
	in_weights.open("../data/Weights.txt"); // 28000 linhas
    //in_weights.open("C:\\Users\\pedro\\Weights_1");// 42000 linhas
	if (!in_weights) {
		cout << "Error in openning the file\n";
        return 1;
	}
	vector<float> W1;
    vector<float> W2;
    vector<float> W3;
    int W1_size=128;
    int W2_size=64;
    float value;
    int i,j,k;
    cout << "read W1\n";
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
	cout << "read W2\n";
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
	cout << "read W3\n";
	for (k = 0; k < (W2_size*10); k++) {
        //in_weights.read(reinterpret_cast<char*>(&value), sizeof(value));
        in_weights >> value;
		W3.push_back(value);
    }
	cout << "\n";
    //cout << "-------------------------Writing W3-----------------------------\n";
	in_weights.close();
    cout << "Size of W3 = " << W3.size();
    cout << "\n";
    cout << "y_train Vector Size = "<< y_train.size();
    cout << "\n";
    cout << "x_train Vector Size = "<< X_train.size();
    cout << "\n";
    //cout << "Testing the model ...\n";
    
            // Building batches of input variables (X) and labels (y)
        //int randindx = rand() % (42000-BATCH_SIZE);
        int randindx = 0;
        vector<float> b_X;
        vector<float> b_y;
        for (unsigned j = randindx*784; j < (((randindx+BATCH_SIZE)*784)>>2)<<2; j+=10){
            b_X.push_back(X_train[j]);
            b_X.push_back(X_train[j+1]);
            b_X.push_back(X_train[j+2]);
            b_X.push_back(X_train[j+3]);
            b_X.push_back(X_train[j+4]);
            b_X.push_back(X_train[j+5]);
            b_X.push_back(X_train[j+6]);
            b_X.push_back(X_train[j+7]);
            b_X.push_back(X_train[j+8]);
            b_X.push_back(X_train[j+9]);
        }
        for (unsigned j = (((randindx+BATCH_SIZE)*784)>>2)<<2; j < (randindx+BATCH_SIZE)*784; j++) b_X.push_back(X_train[j]);
        
        for (unsigned k = randindx*10; k < (randindx+BATCH_SIZE)*10; ++k){
            b_y.push_back(y_train[k]);
        }

        // Feed forward
        
        vector<float> a1 = dot( b_X, W1, BATCH_SIZE, 784, W1_size );
        float max1=-100;
        float min1=100;
        for(int i=0;i<a1.size();i++){
            if(a1[i]>max1)max1=a1[i];
            if(a1[i]<min1)min1=a1[i];
        }
        if (sig==1) a1=sigmoid(a1);
        vector<float> a2 = dot( a1, W2, BATCH_SIZE, W1_size, W2_size );
        float max2=-100;
        float min2=100;
        for(int i=0;i<a2.size();i++){
            if(a2[i]>max2)max2=a2[i];
            if(a2[i]<min2)min2=a2[i];
        }
        if (sig==1) a2=sigmoid(a2);
        vector<float> res = dot( a2, W3, BATCH_SIZE, W2_size, 10 );
        vector<float> yhat=softmax(res,10);
        
        cout<<"The max and min of a1 are: "<<max1<<" "<<min1<<endl;
        cout<<"The max and min of a2 are: "<<max2<<" "<<min2<<endl;
        cout << "b_y Vector Size = "<< b_y.size();
        cout << "\n";
        cout << "b_X Vector Size = "<< b_X.size();
        cout << "\n";
        cout << "yhat Vector Size = "<< yhat.size();
        cout << "\n";
      
   // -------------------calculo do erro ---------------------------------
            float max_value;
            int indice, error;
            error=0;
            double Error_Percentage;
            for (i=0; i<BATCH_SIZE; i++){
                max_value=0;
                for (j=0; j<10; j++){
                    if (yhat[i*10+j]>max_value) {
                            max_value=yhat[i*10+j];
                            indice=i*10+j;
                    }
                }
                if (!b_y[indice]) error++;

                //cout << "b_y[indice] = " << b_y[indice] << ", indice = " << indice <<"yhat[indice] = "<< yhat[indice];
                //cout << "\n";
                            }
            Error_Percentage=(static_cast<double>(error)/static_cast<double>(BATCH_SIZE))*100;
            cout << "Error_abs " << error;
            cout << "\n";
            cout << "Taxa de Successo = " << 100-Error_Percentage <<"%";
            cout << "\n";
            //print ( yhat, 40, 10 );
            //cout << "Ground truth:" << "\n";
            //print ( b_y, 40, 10 );
            vector<float> loss_m = yhat - b_y;
            float loss = 0.0;
            for (unsigned k = 0; k < BATCH_SIZE*10; ++k){
                loss += loss_m[k]*loss_m[k];
            }
            cout << "                                            Loss " << loss/BATCH_SIZE <<"\n";
            cout << "--------------------------------------------End of Program :(------------------------------------------------" <<"\n";
        
    
    return(0);
}