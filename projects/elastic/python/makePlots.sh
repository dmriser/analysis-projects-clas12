#!/bin/bash

#python data_sim.py -d=rga-in.root -s=esepp.root -o=compare-in
#python mon.py -i=rga-in.root -o=rga-in
#python mon.py -i=esepp.root -o=esepp
python monrad.py -i=rga-rad.root -o=rga-rad
python monrad.py -i=esepp-rad.root -o=esepp-rad


