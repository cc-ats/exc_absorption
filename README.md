# exc_absorption
exciton absorption spectra for H2OBPC thin films

Description: 

1. Methods folder has necessary python scripts

   
2. Systems folder has QChem output files for H2OBPC molecule for functionals wPBE, PBE, PBE0

Usage: 

       1. extract.py extracts the data from Qchem .out files, saves them in JSON files 


       2. tb_ham.py builds the TB hamiltonian for different systems using the output from extract.py
