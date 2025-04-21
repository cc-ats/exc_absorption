# Exciton Absorption 

This repository provides a codebase for parsing Q-Chem output and constructing tight-binding Hamiltonians for excitonic systems. The codebase is organized into two main stages and supports flexible system dimensionality (1D, 2D, 3D) for model building and comparison with experiment.

---


---

## Workflow

### **Stage 1: Extract Data from Q-Chem Output**

Run `extract.py` to parse Q-Chem `.out` files and generate corresponding `.json` files.

#### **Command Format**

python3 extract.py <molecule> <functional>

python3 extract.py H2OBPC PBE0

### **Stage 2: Generate Tight Binding Hamiltonian**

Run `tb_ham.py` using the generated JSON to construct the model Hamiltonian. The script supports 1D, 2D, and 3D arrangements with customizable coupling schemes and spectral parameters.

#### **Command Format**

Examples for comparison with experimental data:

1D-AB model:

python3 tb_ham.py H2OBPC PBE0 "600,1,1" --zigzag "1" --coupling "J_coul" --rcut "26" --eshift 0.4 --gamma 0.45 --digitized ../systems/raw-sum-tabs.txt

2D model:

python3 tb_ham.py H2OBPC PBE0 "60,10,1" --coupling "J_coul" --rcut "17" --eshift 0.4 --gamma 0.45 --digitized ../systems/raw-sum-tabs.txt

3D model:

python3 tb_ham.py H2OBPC PBE0 "30,10,2" --coupling "J_coul" --rcut "10" --eshift 0.4 --gamma 0.45 --digitized ../systems/raw-sum-tabs.txt


This material is based on work supported by National Science Foundation under NSF OAC-2311442. YS and PH were supported by NSF grant OAC-2311442 for part of this study, specifically the code development for modeling exciton absorption spectra.









