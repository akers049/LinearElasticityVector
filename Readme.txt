

By Andy Akerson Aug 3 2022 -- updated Sept 24 2026

This does the linear elasticity for a 2D cantelever using deal.ii
(tested for version 9.6.1)


Calling format to compile and run this:

cmake . -DDEAL_II_DIR="dealii install location"
make
printf "inputFiles/inputFile_test.in" | ./run_stuff 


it will dump some vtu outputs into ./output/run_00/lagrangian_solution/

Files included:

  ./CMakeList.txt                   : Cmake file
  ./source/LinearElasticity.cc      : Source code for the system assembly and solver and stuff
  ./source/Constitutive.cc          : Source file for the linear elastic constitutive
  ./source/run_stuff.cc             : Source file for the int main()
  ./include/LinearElasticity.h      : Header file for the system assembly and shit
  ./include/Constitutive.h          : Header file fro the constitutive
  ./inputFiles/inputFile_test.in    : example input file. Says the file format inside it

