# iron-magma

Model for prediction of iron redox in magmatic liquids

Copyright Roberto Moretti (2005)

# Citation

# Licence

MIT licence, see Licence file.

# Contributors

Roberto Moretti (IPGP), moretti@ipgp.fr

Charles Le Losq (IPGP), lelosq@ipgp.fr

# Dependencies

A working fortran compiler. We suggest using gfortran, tested on Mac and Linux. It works well with this software!

# How to use

Download the repository, and use the provided example input file. It first requires compilation of the FORTRAN source, then running the compilated software.

## Compilation

To create the program, with gfortran, in the terminal on Linux or MacOS:

`$ gfortran iron2.for -o iron2.o`

## Running the software

The software takes an input file, iron2.in, which contains the compositions of interest.

It returns an output file, iron2.RDX

Run in the terminal, after compilation, run the command:

`$ ./iron2.o`

# Input file: INPUT.txt

This is a text files, where values are separated by commas.

First line: enters the number of lines to process, leave other columns empty.

then subsequent lines, enter your composition and conditions in the order

SiO2, TiO2, Al2O3, Fe2O3, Cr2O3, FeO, Mno, MgO, CaO, Na2O, K2O, P2O5, H2O, T°C, Pbar, logfO2, INTEGER, FLAG,

oxides are wt% (normalization done internally)

INTEGER: whatever you want, leave 1

FLAG: 0 or 2; if O takes logfO2 and calculates FeO and Fe2O3; if 2 takes FeO and Fe2O3 and calculates logfO2

# Results: OUTPUT.txt

Results are outputed in the file OUTPUT.txt. Everything is self explanatory in this file, with a header.