# trans-C
### Identification of *trans* interacting DNA domains

Trans-C is a tool that finds in Hi-C data biologically meaningful sets of loci that interact strongly in *trans*. 

It initiates probabilistic random walks with restarts from a set of seed loci to traverse an input Hi-C contact map and returns dense subnetwork containing the seed loci.

You can find more about the theory behind trans-C by reading our [paper](link.pdf). If you use Trans-C in an academic setting please cite us.


![Alt text](fig1.png)



### Dependencies
trans-C is implemented in Python and uses the common numpy, sklearn, and scipy packages. Additionally, it relies on the iced and cooler packages, so make sure to install these dependencies (pip installed or conda install work effortlessly).



### How to run
The user needs to provide a Hi-C map, a file specifying the chromosome lengths, the resolution of the analysis and a file with the seed loci. Optionally, you may provide a path to the desired output directory and a value for the restart parameter alpha (default = 0.5). These inputs need to be specified in a config.yml or passed as a command line arguments.

Then, simply run:

\%> python bin/trans-C.py config.yml

or

\%> python bin/trans-C.py hic\_file(.cool/.hic/.npy)  chrom\_sizes.txt -100000 seed\_file.bed out\_dir (-alpha 0.5)
 

1. The Hi-C matrix needs to be in one of the three common formats: cooler (.cool), hic file (.hic), or numpy object (.npy)
2. The chromosome sizes file should follow the canonical structure chromosome\_name chromsome\_length, i.e chr1 248956422
3. The resolution of the analysis is the size of the bins of the Hi-C matrix, i.e 100kb
4. The seed loci should be provided in .bed format, i.e chr4 1350000 1362000

### Questions
Feel free to email Borislav Hristov: borislav at uw.edu 
