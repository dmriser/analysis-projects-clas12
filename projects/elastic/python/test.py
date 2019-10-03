#!/usr/bin/env python 

def load_histos(file):
    h = {} 
    for k in file.GetListOfKeys():
        h[k.GetName()] = file.Get(k.GetName())
    return h 

if __name__ == '__main__':

    # Try to load root 
    try:
        from ROOT import (TH1F, TH2F, TFile, TCanvas,
                          gPad, gStyle)

        print('ROOT loaded')
    except ImportError as e:
        print('ROOT failed to load: {}'.format(e))
    
    # Try to do some work
    input_filename = 'histos.root'
    input_file = TFile(input_filename)
    histos = load_histos(input_file)
    print(histos.keys())

    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)
    gStyle.SetPalette(57)

    can = TCanvas('can', 'can', 800, 600)

    can.Divide(3,2)
    for i in range(1,7):
        can.cd(i)
        histos['histos_dwtheta_'+str(i)].Draw('colz')
    
    can.Print('test.pdf[')
    can.Clear() 

    histos['histos_dwphi'].Draw('colz')
    can.Print('test.pdf')
    can.Clear() 

    can.Divide(3,2)
    for i in range(1,7):
        can.cd(i)
        histos['histos_dthetaele_'+str(i)].Draw('hist')
    can.Print('test.pdf')
    can.Clear() 


    can.Divide(3,2)
    for i in range(1,7):
        can.cd(i)
        histos['histos_dpele_'+str(i)].Draw('hist')
    
    can.Print('test.pdf]')
    
    


