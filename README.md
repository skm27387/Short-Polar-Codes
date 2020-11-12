# Short-Polar-Codes
This repository is the result of a set of experiments carried out on Polar codes and their modifications for short block lengths.
Short codes are useful for applications in 5G and beyond which require ultra reliability and low latency communication (URLLC).

This work is inspired by Arikan's 2019 Shannon Lecture and accompanying papers.

[1] Arıkan, Erdal. "From sequential decoding to channel polarization and back again." arXiv preprint arXiv:1908.09594 (2019). (https://arxiv.org/abs/1908.09594)

[2] Yao, Hanwen, Arman Fazeli, and Alexander Vardy. "List Decoding of Arikan's PAC Codes." arXiv preprint arXiv:2005.13711 (2020). (https://arxiv.org/abs/2005.13711)

[3] Rowshan, Mohammad, Andreas Burg, and Emanuele Viterbo. "Polarization-adjusted Convolutional (PAC) Codes: Fano Decoding vs List Decoding." arXiv preprint arXiv:2002.06805 (2020). (https://arxiv.org/abs/2002.06805)

[4] Li, Bin, Hui Shen, and David Tse. "An adaptive successive cancellation list decoder for polar codes with cyclic redundancy check." IEEE Communications Letters 16.12 (2012): 2044-2047. (https://arxiv.org/abs/1208.3091)

[5] Sason, Igal, and Shlomo Shamai. Performance analysis of linear codes under maximum-likelihood decoding: A tutorial. Now Publishers Inc, 2006. (https://webee.technion.ac.il/people/sason/monograph_postprint.pdf)

Notes : 
1. The curves that indicate SP Polar Codes are selectively precoded polar codes (https://arxiv.org/pdf/2011.04930.pdf).
2. ShortPolarCodesSim.pdf shows the BLER simulation results of short polar code variants for R = 0.5, N = 128. In case of CA-SCL, 8-bit CRC is used.
3. TubShortPolarCodes.pdf shows the Truncated Union Bound analysis of short polar code variants for R = 0.5, N = 128.

This is a work in progress.
We will update the results as simulation progresses.

We have added a C++ file spp.cpp which can be compiled using g++ compiler. 
Using this code, codeword weight analysis can be done for both SPP and PAC codes.
fRM_128.txt has the information bit indices to be used for RM rate profiling. 


Disclaimer : We do not guarantee that the codes are of commercial quality or that it is designed to meet all scenarios. We just provide a basic framework. Please use at your own risk.
