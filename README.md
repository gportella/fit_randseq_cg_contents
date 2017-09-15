
# Generates BED file from random regions that best matches the CG contents of a reference region

This program is intended to help generating control regions with similar CG contents. 

Given a reference bed file *with intervals of the same length* and another bed file *with intervals of the same
length*, it generates an output bed file from regions from the "control" bed file such that the CG content
profile along the sequences enoded in the reference bed files match with those of the control. The output bed file has
the same size as the refernce one.

It requires a genome file, i.e. a single fasta file that is compatible with the bed files.


### How to compile

Check out the `do_cmake_*.sh` scripts. 


### TODO

Remove unnecessary dependencies (OpenMP, Eigen), as they are not used, boilerplate code 
inherited from similar programs.

#### Requires

1. Seqan 2.2.0 (otehr versions untested, might work as well)


