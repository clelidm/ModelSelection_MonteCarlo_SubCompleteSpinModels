# README

## Requirements:

`MC_Algo1_32bits.cpp` uses the C++11 version of C++.

## Usage

### Input datafile

The input datafile must be a __list of *n*-bit binary numbers__, each bit representing a spin variable (binary variable).
Each line is a datapoint. The file must have no space between the columns (the bits of the same number).

   __Ex.__ 4 first datapoints of an input file with *n=10* spin variables:
   
          0001000000
          1000100001
          0000000100
          0000100000
          
   Each column is a spin; Each line is a datapoint.

At __maximum *n*=32__:  datapoints are treated by the program as integers encoded on 32 bits, 
the lowest bit being the last column of your datafile.
The program could easily be extended to larger values of *n*. 
However the present algorithm is efficient only for small values of *n*, and have difficulties to find the 
global maximum for *n* larger than 25. 
The search for the best sub-complete models in larger systems require an evolved version of this algorithm.

### Specification at the beginning of the `MC_Algo1_32bits.cpp` file

At the beginning of the `MC_Algo1_32bits.cpp` file, __you must specify__:
 - `const string data_file_name`: the location and name of your datafile;
 - `const int n`: the total number *n* of spin variables in your datafile, *n*<32;
 - and `const int m`: the chosen number *m* of operators on which the sub-complete model is based, *m*<*n*. 
 
The two values of *n* and *m* fully specify the class of sub-complete models the program will explore.

You can modify:
 - `const string directory`: the output directory;
 - `const string MC_fileOUT`: the name of the output file instead of `outputfile_name`, 
the details of the MC steps will be printed in this file.

You can also play with the parameters for the MC:
 - `const int N_MCsample`: total number of MC steps (including shuffling);
 - `const double beta`: the inverse of the temperature.
 
 ### Compiling
 
 For compiling, use the command:
 
     g++ -std=c++11 -O3 MC_Algo1_32bits.cpp
     
You may want to use the option `-O3` of g++ to turn on some optimisations. If `-std=c++11` doesn't work, please try `-std=c++0x` instead. This creates an executable `a.out`, that can be run using:

     ./a.out

You might want to rename the executable using the `-c` option of g++.
These two command lines are reminded at the beginning of the `MC_Algo1_32bits.cpp` file.
 
 
