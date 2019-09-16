#!/usr/bin/env python 

import os 

if __name__ == '__main__':

    env_file = '/u/home/dmriser/analysis-projects-clas12/environment/farm.csh'
    groovy_file = '/u/home/dmriser/analysis-projects-clas12/projects/simple_examples/parallelLoop.groovy'
    output_dir = '/volatile/clas12/dmriser'
    data_dir = '/work/clas12/rg-a/trains/v16_v2/skim4_inclusive'

    os.system('source {}'.format(env_file))
    os.system('run-groovy {} {}/*'.format(groovy_file, data_dir))
    os.system('cp output.hipo {}'.format(output_dir))
