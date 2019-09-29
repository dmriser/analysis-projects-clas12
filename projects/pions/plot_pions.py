#!/usr/bin/env python

from ROOT import TCanvas, TFile, TH1F, TH2F
from ROOT import gPad, gStyle

if __name__ == '__main__':

    input_filename = 'output.hipo.root'
    root_file = TFile(input_filename)

    histos = {}
    colors = {'pip':99, 'pim':55}
    
    for had in ['pip', 'pim']:
        for axis in ['x', 'z', 'q2', 'pt']:
            for isect in range(1,7):
                title = 'pions_{}_{}_{}'.format(axis,had,isect)
                print('Getting ' + title)
                histos[title] = root_file.Get(title + '_eid')
                histos[title].SetLineColor(colors[had])
                histos[title].SetXTitle(axis)

    for had in ['pip', 'pim']:
        for axis in ['xq2', 'zpt']:
            for isect in range(1,7):
                title = 'pions_{}_{}_{}'.format(axis,had,isect)
                histos[title] = root_file.Get(title + '_eid')

    for had in ['pip', 'pim']:
        for axis in ['beta_p']:
            for isect in range(1,7):
                title = 'pions_{}_{}_{}'.format(axis,had,isect)
                print('Getting ' + title)
                histos[title] = root_file.Get(title + '_eid')


    for had in ['pip', 'pim']:
        for axis in ['dc1_xy', 'dc2_xy', 'dc3_xy']:
            title = 'pions_{}_{}'.format(axis,had)
            print('Getting ' + title)
            histos[title] = root_file.Get(title + '_eid')

    #  pions_dc1_xy_pip_chi2
    for axis in ['dc1_xy_pip', 'dc1_xy_ele']:
        title = 'pions_{}_chi2'.format(axis)
        print('Getting ' + title)
        histos[title] = root_file.Get(title)
        histos[title].Scale(1.0 / histos[title].GetEntries())
        histos[title].SetMinimum(0.0)
        histos[title].SetMaximum(2.0)
        
    for axis in ['dc1_xy_ele']:
        title = 'pions_{}_chi2_eid'.format(axis)
        print('Getting ' + title)
        histos[title] = root_file.Get(title)
        histos[title].Scale(1.0 / histos[title].GetEntries())
        histos[title].SetMinimum(0.0)
        histos[title].SetMaximum(2.0)

    canvas = TCanvas('canvas', 'canvas', 1200, 800)
    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)

    for axis in ['x', 'z', 'q2', 'pt']:
        canvas.Clear()
        canvas.Divide(3,2)

        for isect in range(1,7):
            for had in ['pip', 'pim']:
                canvas.cd(isect)
                histos['pions_{}_{}_{}'.format(axis,had,isect)].Draw('same')

        canvas.Print('{}.pdf'.format(axis))

    for axis in ['xq2', 'zpt']:
        canvas.Clear()
        canvas.Divide(3,2)
        for had in ['pip', 'pim']:
            for isect in range(1,7):
                canvas.cd(isect)
                gPad.SetLogz()
                histos['pions_{}_{}_{}'.format(axis,had,isect)].Draw('colz')
                
            canvas.Print('{}_{}.pdf'.format(had,axis))

    for axis in ['beta_p']:
        canvas.Clear()
        canvas.Divide(3,2)
        for had in ['pip', 'pim']:
            for isect in range(1,7):
                canvas.cd(isect)
                gPad.SetLogz()
                histos['pions_{}_{}_{}'.format(axis,had,isect)].Draw('colz')
                
            canvas.Print('{}_{}.pdf'.format(had,axis))

    # Chi2 plots
    canvas.Clear()
    histos['pions_dc1_xy_ele_chi2'].Draw('colz')
    canvas.Print('dc1_xy_ele_chi2.pdf')

    canvas.Clear()
    histos['pions_dc1_xy_pip_chi2'].Draw('colz')
    canvas.Print('dc1_xy_pip_chi2.pdf')

    canvas.Clear()
    histos['pions_dc1_xy_ele_chi2_eid'].Draw('colz')
    canvas.Print('dc1_xy_ele_chi2_eid.pdf')

