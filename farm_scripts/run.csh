#!/bin/tcsh 

set ENV_FILE    = "/u/home/dmriser/analysis-projects-clas12/environment/farm.csh"
set GROOVY_FILE = "/u/home/dmriser/analysis-projects-clas12/projects/simple_examples/parallelLoop.groovy"
set OUTPUT_DIR  = "/volatile/clas12/dmriser"
set DATA_DIR    = "/work/clas12/rg-a/trains/v16_v2/skim4_inclusive"

source $ENV_FILE
run-groovy $GROOVY_FILE $DATA_DIR/* 
cp output.hipo $OUTPUT_DIR
