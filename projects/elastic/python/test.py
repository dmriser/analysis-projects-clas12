#!/usr/bin/env python 


if __name__ == '__main__':

    # Try to load root 
    try:
        from ROOT import TH1F, TFile, TCanvas
        print('ROOT loaded')
    except ImportError as e:
        print('ROOT failed to load: {}'.format(e))
    
    # Try to do some work
    hist = TH1F('hist', 'hist', 100, -5, 5)
    hist.FillRandom('gaus', 1000)
    
    can = TCanvas('can', 'can', 800, 600)
    hist.Draw()
    can.Print('test.pdf')
    
    


