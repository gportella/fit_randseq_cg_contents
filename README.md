
# Generates BED file from random regions that best matches the CG contents of a reference region

This program generates "control" regions with similar CG contents. It also outputs the CG enrichement profile
for a bunch of equal length genomic intervals (aka `bed` file), without actually generating the control regions.

The CG-enrichment is defined in windows (option -w) as 
```
(#C + #G )  / (#C + #G + #A + #T + #N) - CG-background 
```

where #N is used in the case the base is not assigned in the genome. The CG-background depends on the genome (mouse is 0.42). 

Given a reference bed file *with intervals of the same length* and another bed file *with intervals of the same
length*, it generates an output bed file from regions from the "control" bed file such that the CG content
profile along the sequences enoded in the reference bed files match with those of the control. The output bed file has
the same size as the reference one. *Obviously* the length of the "other" bed file should be larger.

The program randomly samples from the "other" bed file and proposes a CG-profile. It uses a MonteCarlo move to accept or reject the regions from the "other" bed file. It calculates the differences with the reference profile, and if their root mean square difference is below a threshold it outputs the candidate random regions and stops. If not, it continues until the maximum number of iterations is exausted.  

It requires a genome file, i.e. a single fasta file that is compatible with the bed files.
Issue `cg_contents -h ` to get a help message, with all commands available.  

Assuming the `cg_content` is in the path, a typical execution to simply generate a CG enrichment profile would be

```bash
cg_content -i genome.fa -b genomic_regions.bed -pf -cgref cg_enrichment_profile.dat 
```

or to generate the random control regions with the same CG profile, 

```bash
cg_content -i genome.fa -b reference.bed -br random_regions.bed -max_iter 5000000 -cutoff 0.00010 -bo rand_regions_fit_CG_.bed -v
```

### How to compile

Clone this repository, and move to the root directory. Then simply issue
```bash 
cmake . 
make
# or  make -j 4 to speed up

```
and you should get the compiled program.  You would need a modern version of CMake (3.5 or above), and Seqan library (2.2 or above).
If `cmake` can not find the dependencies by itself, either fix the CMakeLists.txt file or give it a little help, e.g.

```bash
cmake . -DCMAKE_C_COMPILER=/opt/local/bin/clang-mp-3.9 -DCMAKE_CXX_COMPILER=/opt/local/bin/clang++-mp-3.9 -DSEQAN_INCLUDE_PATH=/YOUR_PATH/seqan-library-2.2.0/include
```


##### Requirements

1. Compiler supporting at least C++14 (required to compile Seqan anyway)
2. CMake 3.5.0 or above
3. Seqan 2.2.0 (other versions untested, might work as well)



