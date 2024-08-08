# Overview

Each folder contains results of 100 [mss_sim.py](https://github.com/ampivirotto/mss_forward_sim) runs for a specific selection coefficient. The selection coefficient is given in the folder name. 

Each run began with selecting a random bacterial gene alignment and generating a random ancestral sequence from the codon usage in that alignment.

Each run generated a fasta file with sample of 11 sequences, one from each species, on the specified phylogeny. Each run also generated an output file summarizing the simulation. The gene names for each run are used to name the fasta file and the results file. (*e.g.*, yeaG in S0_rand_yeaG.fa and s0_yeaG_results.txt)

Each fasta file was copied onto a new file with an anonomyzed name for hyphy analyses. The anonymized files are in a subfolder within each folder. 
