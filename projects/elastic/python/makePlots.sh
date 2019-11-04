#!/bin/bash

#python data_sim_lite.py -d=mon-lite-rga-in.root -s=mon-lite-sim02-in.root -o=compare-in
python mon.py -i=rga-in.root -o=rga-in
#python pres.py -i=rga-in.root -o=rga-in

#python mon.py -i=rga-out.root -o=rga-out
#python mon.py -i=sim04-in.root -o=sim04-in

#python data_sim.py -d=rga-in.root -s=sim04-in.root -o=compare-in
#python data_sim.py -d=rga-out.root -s=sim04-in.root -o=compare-out
