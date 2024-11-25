HOW TO RUN

Firstly, one should execute the training file, to generate the weights and afterwards one can execute the test file (to check the accuracy) or the 
system to run the encripted version of the NN.
In this case, the weights are already on the data folder, so training isint necessary.

How to run training:

Run file called "SC_train" in the source folder.
Aditionaly, one can compile the source code "SC_train.cpp" with their c++ compiler.

How to run testing:

Run file called "Test_Dataset" in the source folder.
This program will ask for the pretended activation function which can be :
    L-linear
    S-sigmoid
Aditionaly, one can compile the source code "Test_Dataset.cpp" with their c++ compiler.

How to run encripted system:

It is important to note that the library OpenFHE is mandatory for this program
Run file called "SC_enc" in the directory "/source/build".
This program will ask firstly for the pretended activation function which can be:
    L-linear (no activation function)
    S-sigmoid
In case the user chooses sigmoid, he/she will be asked to provide the intended level of presision, the more presise the 
slower it will be. The possible values are:
    L-low
    M-medium
    H-high
Aditionaly, one can compile using the comand "make" in the directory "source/build". the cmake file is also available,
please refer to the OpenFHE documentation if you intend on changing it.

In case the SC_enc program gives any type of error related to high error, this is due to the fact that the lower and upper
bounds in the chebyshev aproximation are incorrect. In this case, please run "Test_Dataset", and copy the maximum and minimum
indicated and replace the values in the "SC_enc".
    