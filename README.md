# Dataset-of-HRP-38-test-system
This repository includes the dataset of HRP-38 system in both .xlsx and .mat files. The index of the selected typical day used in the Case Study is also included.

The codes for operation simulation applied in this paper are also provided in the dataset linkage. Based on the programs, researchers can compare their customized planning scheme with the benchmarks present in the paper. The codes can also be modified to consider new elements by researchers themselves. Thus, the proposed system can facilitate comparisons and collaboration between different TEP studies worldwide.
run_all_op.m is the main file. The simulation for the four cases will run automatically one by one after running the main file.
Yalmip, matpower toolbox, and CPLEX solver are required to run the codes.

HRP-38-50 is the dataset where the renewable energy penetration is scaled up to 50%. 

The related documents and explanation can be found in the following linkage:
https://ieeexplore.ieee.org/document/8931650

errata:
The units of generator investment cost are 10^7CNY rather than 10^4CNY
