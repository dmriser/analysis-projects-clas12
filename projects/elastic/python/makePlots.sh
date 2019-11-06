#!/bin/bash

#python data_sim_lite.py -d=rga-in-lite.hipo.root -s=esepp-lite.hipo.root -o=compare-in
#python pres.py -i=rga-in.root -o=rga-in

python data_sim.py -d=rga-in.root -s=esepp.root -o=compare-in
python mon.py -i=rga-in.root -o=rga-in
python mon.py -i=esepp.root -o=esepp
#python data_sim.py -d=rga-out.root -s=sim04-in.root -o=compare-out
