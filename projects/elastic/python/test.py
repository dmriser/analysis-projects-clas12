#!/usr/bin/env python 

from ROOT import (TH1F, TH2F, TF1, TFile, TCanvas,
                  gPad, gStyle, TLatex)

def load_histos(file):
    ''' Use the ROOT file structure to load a dictionary of histograms. '''
    h = {} 
    for k in file.GetListOfKeys():
        h[k.GetName()] = file.Get(k.GetName())
    return h 

def plot_sector_momreso(histos, can, lab, tex, tit):
    ''' Plot the results of fitting the resolution '''
    can.Clear()
    can.Divide(3,2)
    for i in range(1,7):
        can.cd(i)
        title = tit.format(i)

        mean = histos[title].GetMean()
        std = histos[title].GetStdDev()
        fit = TF1('fit_{}'.format(i), 'gaus', mean - 1.7 * std, mean + 1.7 * std)
        fit.SetParameter(1, mean)
        fit.SetParameter(2, std)

        histos[title].Fit(fit, 'brq')
        fit.Draw('same')

        histos[title].SetLineColor(1)
        histos[title].Draw()

        par1 = '#mu = {0:6.3f}'.format(fit.GetParameter(1))
        par2 = '#sigma = {0:6.3f}'.format(fit.GetParameter(2))
        err1 = ' #pm = {0:6.3f}'.format(fit.GetParError(1))
        err2 = ' #pm = {0:6.3f}'.format(fit.GetParError(2))
        label = '#splitline{' + par1 + '}{' + par2 + '}'
        lab.DrawLatex(0.65, 0.82, label)
        lab.DrawLatex(0.45, 0.02, '#Delta p / p' if 'fraction' in title else '#Delta p')
        lab.DrawLatex(0.15, 0.82, 'Sector {}'.format(i))

    can.cd(2)
    can_tit = 'Momentum Resolution (Electron)' if 'ele' in tit else 'Momentum Resolution (Proton)'
    tex.DrawLatex(0.01, 0.95, can_tit)

    output_title = tit.split('_{}')[0].split('histos_')[-1] + '_sectors.pdf'
    can.Print(output_title)

def plot_sector_thetareso(histos, can, lab, tex, tit):
    ''' Plot the results of fitting the resolution '''
    can.Clear()
    can.Divide(3,2)
    for i in range(1,7):
        can.cd(i)
        title = tit.format(i)

        fit = TF1('fit_{}'.format(i), 'gaus')
        histos[title].Fit(fit)
        fit.Draw('same')

        histos[title].SetLineColor(1)
        histos[title].Draw()

        par1 = '#mu = {0:6.3f}'.format(fit.GetParameter(1))
        par2 = '#sigma = {0:6.3f}'.format(fit.GetParameter(2))
        err1 = ' #pm = {0:6.3f}'.format(fit.GetParError(1))
        err2 = ' #pm = {0:6.3f}'.format(fit.GetParError(2))
        label = '#splitline{' + par1 + '}{' + par2 + '}'
        lab.DrawLatex(0.65, 0.82, label)
        lab.DrawLatex(0.45, 0.02, '#Delta #theta (degrees)')
        lab.DrawLatex(0.15, 0.82, 'Sector {}'.format(i))

    can.cd(2)
    can_tit = 'Angular Resolution (Electron)' if 'ele' in tit else 'Angular Resolution (Proton)'
    tex.DrawLatex(0.01, 0.95, can_tit)

    output_title = tit.split('_{}')[0].split('histos_')[-1] + 'sectors.pdf'
    can.Print(output_title)

def plot_sector_reso2d(histos, can, lab, tex, tit):
    ''' Plot the results of fitting the resolution '''
    can.Clear()
    can.Divide(3,2)
    for i in range(1,7):
        can.cd(i)
        title = tit.format(i)

        yrange = [7.5, 11.5] if 'ele' in tit else [0.5, 4.5]
        histos[title].GetYaxis().SetRangeUser(yrange[0], yrange[1])
        histos[title].Draw('colz')
        gPad.SetLogz() 

        label = '#Delta p / p' if 'fraction' in tit else '#Delta p (GeV)' 
        lab.DrawLatex(0.45, 0.02, label)

    can.cd(2)
    can_tit = 'Momentum Resolution (Electron)' if 'ele' in tit else 'Momentum Resolution (Proton)'
    tex.DrawLatex(0.01, 0.95, can_tit)

    output_title = tit.split('_{}')[0].split('histos_')[-1] + 'sectors.pdf'
    can.Print(output_title)


if __name__ == '__main__':
    
    input_filename = 'histos.root'
    input_file = TFile(input_filename)
    histos = load_histos(input_file)
    print(histos.keys())

    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)
    gStyle.SetPalette(57)

    can = TCanvas('can', 'can', 1080, 720)
    tex = TLatex()
    tex.SetNDC()
    tex.SetTextFont(42)
    tex.SetTextSize(0.05)
    tex.SetTextColor(1)

    lab = TLatex()
    lab.SetNDC()
    lab.SetTextFont(42)
    lab.SetTextSize(0.04)
    lab.SetTextColor(1)

    """
    can.cd(1)
    histos['histos_w_ctof_proton'].SetLineColor(1)
    histos['histos_w_ctof_proton'].Draw('same')
    histos['histos_w_ctof_proton_pass_angle'].SetLineColor(55)
    histos['histos_w_ctof_proton_pass_angle'].Draw('same')
    histos['histos_w_ctof_proton_fail_angle'].SetLineColor(99)
    histos['histos_w_ctof_proton_fail_angle'].Draw('same')
    tex.DrawLatex(0.1, 0.925, 'CTOF protons #color[55]{PASS} and #color[99]{FAIL} angle cut')
    lab.DrawLatex(0.45, 0.02, 'W (GeV/c^{2})')
    can.Print('w_ctof_proton.pdf')

    # Show w spectrum 
    can.cd(1)
    histos['histos_w_base'].SetLineColor(1)
    histos['histos_w_base'].Draw()
    histos['histos_w_ctof_proton'].SetLineColor(55)
    histos['histos_w_ctof_proton'].Draw('same')
    histos['histos_w_tof_proton'].SetLineColor(99)
    histos['histos_w_tof_proton'].Draw('same')
    tex.DrawLatex(0.1, 0.925, 'Elastic electrons #color[55]{CTOF proton} and #color[99]{TOF proton}')
    lab.DrawLatex(0.45, 0.02, 'W (GeV/c^{2})')
    can.Print('w_summary.pdf')
    
    # Show angular spectrum
    can.cd(1)
    histos['histos_opening_angle_ep_ctof_proton'].SetLineColor(1)
    histos['histos_opening_angle_ep_ctof_proton'].Draw()
    histos['histos_opening_angle_ep_ctof_proton_pass_w'].SetLineColor(99)
    histos['histos_opening_angle_ep_ctof_proton_pass_w'].Draw('same')
    histos['histos_opening_angle_ep_ctof_proton_fail_w'].SetLineColor(55)
    histos['histos_opening_angle_ep_ctof_proton_fail_w'].Draw('same')
    tex.DrawLatex(0.1, 0.925, 'CTOF proton #color[99]{pass W} and #color[55]{fail W}')
    lab.DrawLatex(0.45, 0.02, '#theta_{ep} (deg)') 
    can.Print('opening_angle.pdf')
    """
    
    # Show momentum resolution 
    plot_sector_momreso(histos, can, lab, tex, 'histos_delta_p_electron_{}')
    plot_sector_momreso(histos, can, lab, tex, 'histos_delta_p_proton_{}')
    #plot_sector_momreso(histos, can, lab, tex, 'histos_p_res_pro_ctof_proton_pass_all_{}')
    #plot_sector_momreso(histos, can, lab, tex, 'histos_p_res_fraction_ele_ctof_proton_pass_all_{}')
    #plot_sector_momreso(histos, can, lab, tex, 'histos_p_res_fraction_pro_ctof_proton_pass_all_{}')

    # Show angular resolution 
    #plot_sector_thetareso(histos, can, lab, tex, 'histos_theta_res_ele_ctof_proton_pass_all_{}')
    #plot_sector_thetareso(histos, can, lab, tex, 'histos_theta_res_pro_ctof_proton_pass_all_{}')
    
    # Show 2-D resolution as a function of mom
    #plot_sector_reso2d(histos, can, lab, tex, 'histos_p_res_p_ele_ctof_proton_pass_all_{}')
    #plot_sector_reso2d(histos, can, lab, tex, 'histos_p_res_p_pro_ctof_proton_pass_all_{}')
    #plot_sector_reso2d(histos, can, lab, tex, 'histos_p_res_fraction_p_ele_ctof_proton_pass_all_{}')
    #plot_sector_reso2d(histos, can, lab, tex, 'histos_p_res_fraction_p_pro_ctof_proton_pass_all_{}')
