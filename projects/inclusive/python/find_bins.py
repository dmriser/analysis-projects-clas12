#!/usr/bin/env python

import matplotlib.pyplot as plt 
import numpy as np

if __name__ == '__main__':

    input_filename = '../groovy/q2.csv'
    with open(input_filename) as input_file:
        tokens = [float(tok.strip(',')) for tok in input_file.readlines()[0].split()]
        q2 = np.array(tokens)


    q2bins = 10
    quants = np.linspace(0,1,q2bins+1)
    edges = [np.quantile(q2,q) for q in quants]

    plt.hist(q2, bins=np.linspace(0,10,100), edgecolor='k', color='orange')
    for edge in edges:
        plt.axvline(edge, color='k', linewidth=1)
    plt.show()
    

    print(edges)
