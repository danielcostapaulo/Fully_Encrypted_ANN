# Fully Encrypted ANN

Artificial intelligence and Machine Learning(ML) algorithms are exciting and powerful tools for the present and future and are being utilized in more and more industries. A particularly useful industry for this technology is medicine, in which, a patient provides his personal information to a trusted agent, such as a doctor, and this returns with a diagnosis. Utilizing ML to provide a diagnosis has been proven to be a successful and viable option, however, these models directly utilize the personal information of a patient, which a maleficent agent can steal or manipulate. This project implements an Artificial Neural Network(ANN) in a completely encrypted manner, using CKKS, in which the ML model can provide a diagnosis without ever having access to personal information. This project is mostly demonstrative, however, it shows the potential of this technology and its main drawback which is the prediction time.

Some notes on the CKKS scheme are that it can perform arithmetic computations on encrypted vector data. This is far slower as it stands than plain data and has an additional limitation which is the fact that non-arithmetic functions need to be approximated. In the case of ANN, the activation functions of the neurons can be non-arithmetic, such as the sigmoid function, in which this project uses the Chebyshev approximation. This can be a more complex (and slower) approximation or a less complex (and faster) approximation.


## Architecture

The main architecture is the user (a doctor for example) that has access to the patient's information, which will encrypt the data with a key. This key will not be shared. It then generates a new key, MultEvalKey, utilized in the CKKS scheme. This key, alongside the encrypted data, is sent to the ANN, which will derive the diagnosis without ever decrypting the input, that is, it will never have direct access to the personal information of the patients. When the diagnosis is obtained, it is sent back to the user and it can then be decrypted using the same key utilized to encrypt.

## Usage

Firstly, one should execute the training file, to generate the weights and afterward one can execute the test file (to check the accuracy) or the 
system to run the encrypted version of the NN.

### Training:

Run a file called "SC_train" in the source folder.
Additionally, one can compile the source code "SC_train.cpp" with their C++ compiler.

### Testing:

Run a file called "Test_Dataset" in the source folder.
This program will ask for the pretended activation function which can be :
    L-linear
    S-sigmoid
Additionally, one can compile the source code "Test_Dataset.cpp" with their C++ compiler.

### Encrypted system:

It is important to note that the library OpenFHE is mandatory for this program
Run a file called "SC_enc" in the directory "/source/build".
This program will ask firstly for the pretended activation function which can be:
    L-linear (no activation function)
    S-sigmoid
In case the user chooses sigmoid, he/she will be asked to provide the intended level of precision, the more precise the 
slower it will be. The possible values are:
    L-low
    M-medium
    H-high
Additionally, one can compile using the command "make" in the directory "source/build". the CMake file is also available,
please refer to the OpenFHE documentation if you intend to change it.

In case the SC_enc program gives any type of error related to high error, this is due to the fact that the lower and upper
bounds in the Chebyshev approximation are incorrect. In this case, please run "Test_Dataset", and copy the maximum and minimum
indicated and replace the values in the "SC_enc".

## Results

The two metrics analyzed are accuracy and execution time on a dataset with 75MB of data.

| Activation Function  | Chebyshev PolyDegree | Accuracy | Time(m) |
| ------------- | ------------- | ------------- | ------------- |
| Linear  | -  | 32.8% | 12 |
| Sigmoid  | 15  | 96.4429%  | 37  |
| Sigmoid  | 50  | 96.9357%  | 44  |
| Sigmoid  | 70  | 96.9214%  | 48  |


    
