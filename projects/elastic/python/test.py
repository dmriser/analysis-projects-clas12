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
                          gPad, gStyle, TLatex)

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

    can = TCanvas('can', 'can', 800, 800)
    tex = TLatex()
    tex.SetNDC()
    tex.SetTextFont(42)
    tex.SetTextSize(0.05)
    tex.SetTextColor(1)
    
    can.Print('elastic.pdf[')
    can.Divide(2,2)

    # Show w spectrum 
    can.cd(1)
    histos['histos_w_base'].SetLineColor(1)
    histos['histos_w_base'].Draw()
    histos['histos_w_ctof_proton'].SetLineColor(99)
    histos['histos_w_ctof_proton'].Draw('same')

    can.cd(2)
    histos['histos_w_ctof_proton'].SetLineColor(1)
    histos['histos_w_ctof_proton'].Draw('same')
    histos['histos_w_ctof_proton_pass_angle'].SetLineColor(55)
    histos['histos_w_ctof_proton_pass_angle'].Draw('same')
    histos['histos_w_ctof_proton_fail_angle'].SetLineColor(99)
    histos['histos_w_ctof_proton_fail_angle'].Draw('same')
    tex.DrawLatex(0.1, 0.925, 'CTOF protons #color[55]{PASS} and #color[99]{FAIL} angle cut')

    # Show angular spectrum
    can.cd(3)
    histos['histos_opening_angle_ep_ctof_proton'].GetXaxis().SetRangeUser(155,180)
    histos['histos_opening_angle_ep_ctof_proton'].SetLineColor(1)
    histos['histos_opening_angle_ep_ctof_proton'].Draw()
    histos['histos_opening_angle_ep_ctof_proton_pass_w'].GetXaxis().SetRangeUser(155,180)
    histos['histos_opening_angle_ep_ctof_proton_pass_w'].SetLineColor(99)
    histos['histos_opening_angle_ep_ctof_proton_pass_w'].Draw('same')

    can.cd(4)
    histos['histos_opening_angle_ep_tof_proton'].SetLineColor(1)
    histos['histos_opening_angle_ep_tof_proton'].Draw()
    histos['histos_opening_angle_ep_tof_proton_pass_w'].SetLineColor(99)
    histos['histos_opening_angle_ep_tof_proton_pass_w'].Draw('same')

    
    can.Print('elastic.pdf')
    can.Print('elastic.pdf]')
    


