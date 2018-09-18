#!/bin/bash
#LEC2018
echo "Ensure the following are in your PATH: arks, LINKS, arks-make"
echo "	NOTE: arks-make Makefile is found in Examples directory of repo"
echo Running ARKS makefile....
arks-make arks draft=test_scaffolds reads=test_reads m=50-6000 a=0.9 time=1 
