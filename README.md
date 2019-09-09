![Logo](https://github.com/bcgsc/arks/blob/master/arks-logo.png)

# Important Note
ARKS has been integrated into ARCS as of version 1.1.0: https://github.com/bcgsc/arcs

### To run ARCS in ARKS mode:
```arcs --arks [options] <linked reads fastq>```

### The full ARKS pipeline can be run with the supplied Makefile:
```Examples/arcs-make arks draft=myassembly reads=myreads```

Future code maintenance for ARKS will be done in the ARCS repository ONLY.

# ARKS

Scaffolding genome sequence assemblies using 10X Genomics GemCode/Chromium data.
This project is a new kmer-based (alignment free) implementation of
[ARCS](https://github.com/bcgsc/arcs). It provides improved runtime performance
over the original ARCS implementation by removing the requirement to perform
alignments with `bwa mem`.

### Dependencies
* Boost (tested on 1.61)
* GCC (tested on 4.4.7)
* Autotools (if cloning directly from repository) 
* LINKS (tested on 1.8)



### Installation:
#### Install ARKS using Brew
```
brew install brewsci/bio/arks
```

#### Install ARKS from the source code
If cloning directly from the repository run:
```
./autogen.sh
```
To compile ARKS run:
```
./configure && make
```
To install ARKS in a specified directory:
```
./configure --prefix=/ARKS/PATH && make install
```
If your boost library headers are not in your PATH you can specify their location:
```
./configure –-with-boost=/boost/path --prefix=/ARKS/PATH && make install
```

### ARKS+LINKS Pipeline 

ARKS requires three input files:
* Draft assembly fasta file
* Interleaved linked reads file (Barcode sequence expected in the BX tag of the read header; Run [Long Ranger basic](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/what-is-long-ranger) on raw chromium reads to produce this interleaved file)
* CSV file listing barcode multiplicities (Generate with command: bin/calcBarcodeMultiplicities.pl reads.fof > read_multiplicities.csv, where reads.fof is file with name of linked read file)

The Makefile located here: Examples/arks-make will run the full ARKS pipeline, including generating the barcode multiplicity file above, and the 3 steps below. It will also optionally run the misassembly corrector [Tigmint](https://github.com/bcgsc/tigmint) prior to scaffolding with ARKS.

There are three steps to the pipeline:

1. Run ARKS to generate a Graphviz Dot file (.gv). Nodes in the graph are the sequences to scaffold, and edges show that there is evidence to suggest nodes are linked based on the data obtained from the GemCode/Chromium reads.

2. Run the python script Examples/makeTSVfile.py to generate a file named XXX.tigpair_checkpoint file from the ARKS graph file. The XXX.tigpair_checkpoint file will be provided to LINKS in step 3.

3. Run LINKS with the XXX.tigpair_checkpoint file as input. To do this, the base name (-b) must be set to the same name as XXX.

An example bash script on how to run the ARKS+LINKS pipeline using the ARKS Makefile can be found at: Examples/pipeline_example.sh

you can test your installation by following instructions at: Examples/arks_test-demo/README.txt
and compare your output to the files provided at: Examples/arks_test-demo/output/ 

### Citing ARKS

Paper :
https://doi.org/10.1186/s12859-018-2243-x

LINKS :
http://www.bcgsc.ca/platform/bioinfo/software/links
https://github.com/warrenlr/LINKS


### License  

ARKS Copyright (c) 2016 British Columbia Cancer Agency Branch.  All rights reserved.

ARKS is released under the GNU General Public License v3

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

For commercial licensing options, please contact Patrick Rebstein <prebstein@bccancer.bc.ca>
