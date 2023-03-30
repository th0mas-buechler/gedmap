# Preamble

GED-MAP is a prototype for pangemonic read mapping.
This software is part of the article *Efficient short read mapping to a pangenome that is represented by a graph of ED strings* by Thomas Büchler, Jannik Olbrich and Enno Ohlebusch. (Submitted to Bioinformatics)

# Requirements

A modern, C++11 ready compiler such as g++ version 4.9 or higher or clang version 3.2 or higher .

Installation of the SDSL library. (https://github.com/simongog/sdsl-lite)

# Installation

	git clone https://github.com/thomas-buechler-ulm/gedmap/
	cd gedmap
	make gedmap

First, if sdsl is not installed in your home directory, edit the PATHS file and set the correct path to sdsl-lite directory, before executing 'make'.
This produces the executable 'gedmap'.

## Short description

The binary 'gedmap' contains 4 programs:
- gedmap parse: Parses a FASTA and a VCF file into and GEDS (+ files for adjacency and transformation to FASTA positions)
- gedmap index: Calculated a minimizer index from of the GEDS (+ files for adjacency and transformation to FASTA positions) 
- gedmap align: Aligns reads of a FASTQ file to the GEDS with help of the index.
- gedmap sample: Sampled reads from a GEDS

Run 'gedmap parse -h', etc. to get the information about this program.
The parameter for all programs are also described on the bottom of this page.

## A little example
A little example is given in the directory with this name.  (It also contains a README file.)
To test your installation you can run:

	make a_little_example



# Requierements for the experiments

Installations of:
- hisat2 (version 2.2.1) from http://daehwankimlab.github.io/hisat2/manual/)
- samtools
- vg ( version v1.37.0-11-g2f6837d33) from https://github.com/vgteam/vg


# Rerun experiments

Set correct paths in the PATHS file. (Not nescessary if all software tools are installed in your home directory.)

	make data
	make all_eperiments

The first command downloads all human sequence data into the directory 'data'.
The second command runs the three mapping tools. 
The indexes and alignments will be calculated and stored in 'exp/gedmap',  'exp/hista2' and  'exp/giraffe'.
Samples will be generated and stored in 'exp/sample'.
The program 'eval_sam' and 'com_sams' (also generated by 'make all_eperiments') evaluate the sam files to obtain the results given in the paper.

At the end the directories 'exp/gedmap',  'exp/hista2' and  'exp/giraffe'  will contain following files: 'index.out', 'map.out' and 'map.eval'
The first and second contain the information about time and space consumption measured by the time program.
The files 'map.eval' contains the information about mapping rate and the accuracy

# Parameter

### gedmap parse:
                                                       
                                                       
'gedmap parse' parses a FA file and a VCF file to an EDS graph.

2 Arguments expected:

	[1] filename of FA
	[2] filename of VCF
	

Optional parameters: 

	-nosv        , do not include structural variants, i.e. copy number variation and other large variation. (Only an EDS is generated)
	-tmp tmp_dir , to set tmp direcoty (default /tmp)
	-lim l       , variants with ref or alt larger than l are handled as structural variants (default 50)
	-o           , output prefix for the geds, adj and 2fa file (default [1])


                                     
                                     
                                     
### gedmap index:                  
                                                       
'gedmap index' calculates a minimizer index of a given EDS graph.
If [1] was build with allowing structural variants, its strongly
 recommended to provide the adjacency file (parameter -a). This also accelerates the indexing process, because nodes can be indexed parallel. Poviding parameter -2fa can also acellerate the the indexing process if the GEDS contains multiple sequences.

1 Arguments expected:

	[1] filename of GEDS


Optional parameters: 
... IO

	-a   fname  , file name of the adijacency file.
	-2fa fname  , file of the 2fa file
	-o fname    , file name of minimizer index (default [1].min)

... Index parameter

	-k k        , kmer size (default 20)
	-w w        , window size (default 5)
	-n n        , maximum number of N allowed in a minimizer (default 2)
	-trim x     , trims kmer sets greater than x from index (default 1000 , 0 for no trim)

... Control

	-tmp tmp_dir, to set tmp direcoty (default /tmp)
	-t x      , maximum number of threads used (default 128)


                                                       
                      
### gedmap align:             
                                       
'gedmap align' algings reads to the given GEDS and MINIMIZER INDEX
If [2] was build with allowing structural variants, its strongly
 recommended to provide the adjacency file (parameter -a)

3 Arguments expected:

	[1] filename of FASTQ
	[2] filename of GEDS
	[3] filename of MINI


Optional parameters: 

 IO

	-2fa             , .2fa-file , if given this is used to transform GEDS-positions to FA positions
	-a fname         , file name of the adijacency file
	-o               , fname, output will be stored in file fname (default  [2].sam)
	-oa              , only aligned reads will be reported
	-mao x           , max number of alignments in output (default 1)
	-io              , output reads in the same order as in the input (may be a bit slower and with higher memory)

Hotspot finding

	-rc              , reversed complement of pattern will be searched, too
	-mc x            , minimizer count, maximum x minimizers will be looked up per read(default 80)
	-ws x            , window size of hotspot (default 500)
	-wh x            , minimum hotspot score (default 1)

DP parameter

	-d x             , max distance in alignment (default 30)
	-weights s,c,l,h , weights used for alignment-DP: s,c,l,h are costs for gap start, gap continue, minimum mismatch cost and maximum mismatch cost, respectively (default 6,1,0,4)
	-sd x            , satisfying distance, when best alignment has a distance x or smaller, don't align further reads (default 7)
	-mac x           , max number of alignments completely calculated (default 5)
	-mat x           , max number of alignments tried to calculate (default 10)

Paired end parameter

	-mp              , filename of FASTQ containing the mates (optional, presence indicates paired-end mode)
	-fragment-mean   , mean length of fragment in paired-end mode (ignored if -mp is not present, default 700)
	-fallback        , fallback for paired-end mapping (ignored if -mp is not present)
	-fmat            , maximum number of alignments tried for fallback (ignored if -mp and -fallback are not present, default 100)
	-mam x           , max number of alignments used for pairing (default 1,10)

Control

	-tmp tmp_dir    , to set tmp direcoty (default /tmp)
	-t x            , maximum number of threads used per index copy (default uses as many as avaiable)


                                                       
                                                                                         
### gedmap sample:             
   
'gedmap sample' samples reads from the given GEDS.

1 Arguments expected:

	[1] filename of GEDS


Optional parameters: 

	-c x , sample x reads (count, DEFAULT x=100)
	-l x , reads length = x (length, DEFAULT x=100)
	-e x , probability of a base to be false = 1/x (error rate, DEFAULT x=100)
	-e_d x , probability of a base indel = 1/x, (insertion rate, DEFAULT x=1000)
	-e_i x , probability of a base indel = 1/x, (deletion rate, DEFAULT x=1000)
	-rc x , probability of a read to be a reverse complement = 1/x, (rev complement rate, DEFAULT x=2)
	-s x , seed for rng=x, (seed, DEFAULT x=0)
	-o fname , output will be written to file fname (DEFAULT  fname=[1].sample.fastq)
	-fragment-mean x , enables paired-end mode, x is length of fragment from which the reads-pairs come from (second file will be {argument to "-o"}.2)
	-2fa 2fa-file , if given this is used to transform GEDS-positions to FA positions
	-a fname,  file name of the adijacency file
	-ae fname,  file name of the adijacency file, only sample reads that go over edges in the graph

for all probabilities above: prob = 1/x (f.e. if x=4 then prob=25%) and prob=0 for x=0 
