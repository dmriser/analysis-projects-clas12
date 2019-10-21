#!/usr/bin/env python 

from ROOT import (TH1F, TH2F, TF1, TFile, TCanvas,
                  gPad, gStyle, TLatex)

def load_histos(file):
    ''' Use the ROOT file structure to load a dictionary of histograms. '''
    h = {}
    for k in file.GetListOfKeys():
        h[k.GetName()] = file.Get(k.GetName())
    return h

def setup_global_options():
    gStyle.SetOptTitle(0)
    gStyle.SetOptStat(0)

def plot_sector_page(canvas, histos, title_formatter, label, save_name,
                     xtitle=None, ytitle=None, title=None, log=False):
    
    canvas.Clear() 
    canvas.Divide(3,2)

    for i in range(1,7):
        canvas.cd(i)
        
        if isinstance(histos[title_formatter.format(i)], TH1F):
            histos[title_formatter.format(i)].SetFillColorAlpha(55, 0.65)
            histos[title_formatter.format(i)].Draw()

        elif isinstance(histos[title_formatter.format(i)], TH2F):
            histos[title_formatter.format(i)].Draw('colz')
            if log:
                gPad.SetLogz() 
        else:
            raise NotImplementedException('plot_sector_page only supports TH1F, TH2F')
            
        if title:
            label.DrawLatex(0.1, 0.925, title)

        if xtitle:
            label.DrawLatex(0.5, 0.015, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.045, 0.5, ytitle)
            label.SetTextAngle(0)

    canvas.Print(save_name)

def plot_beam_energy(canvas, histos, label, save_name):

    canvas.Clear()
    canvas.Divide(3,2)
    
    for s in range(1,7):

        canvas.cd(s)

        # Start by styling the histograms
        histos['histos_de_beam_{}'.format(s)].SetLineColor(99)
        histos['histos_de_beam_from_angles{}'.format(s)].SetLineColor(55)

        histos['histos_de_beam_{}'.format(s)].Draw()
        histos['histos_de_beam_from_angles{}'.format(s)].Draw('same')

    canvas.Print(save_name)
        
        
        
def plot_w(canvas, histos, save_name, label,
           title=None, xtitle=None, ytitle=None):

    canvas.Clear() 

    # I need to run the code one more time
    # and add a histogram for the w after
    # enforcing CTOF contraints. 
    #
    # Style
    histos['histos_w'].SetLineColor(1)
    histos['histos_w_pass_angle_in_ctof'].SetLineColor(55)

    # Fix missing bins
    histos['histos_w'].GetXaxis().SetRangeUser(0.7, 1.3)
    histos['histos_w_pass_angle_in_ctof'].GetXaxis().SetRangeUser(0.7, 1.3)
    
    histos['histos_w'].Draw()
    histos['histos_w_pass_angle_in_ctof'].Draw('same')
            
    if title:
        label.DrawLatex(0.1, 0.925, title)

    if xtitle:
        label.DrawLatex(0.5, 0.015, xtitle)

    if ytitle:
        label.SetTextAngle(90)
        label.DrawLatex(0.04, 0.5, ytitle)
        label.SetTextAngle(0)

    canvas.Print(save_name)


def plot_theta_p(canvas, histos, save_name, label,
           title=None, xtitle=None, ytitle=None):

    canvas.Clear() 

    histos['histos_theta_p_combined'].SetLineColor(1)
    histos['histos_theta_p_tof_1'].SetLineColor(1)
    histos['histos_theta_p_tof_2'].SetLineColor(1)
    histos['histos_theta_p_tof_3'].SetLineColor(1)
    histos['histos_theta_p_ctof'].SetLineColor(1)

    histos['histos_theta_p_ctof'].SetFillColorAlpha(64, 1.0)
    histos['histos_theta_p_tof_1'].SetFillColorAlpha(77, 1.0)
    histos['histos_theta_p_tof_2'].SetFillColorAlpha(55, 1.0)
    histos['histos_theta_p_tof_3'].SetFillColorAlpha(99, 1.0)
    
    #histos['histos_theta_p_tof_1'].SetLineColor(77)
    #histos['histos_theta_p_tof_2'].SetLineColor(55)
    #histos['histos_theta_p_tof_3'].SetLineColor(99)
    #histos['histos_theta_p_ctof'].SetLineColor(64)

    histos['histos_theta_p_combined'].GetXaxis().SetRangeUser(10, 70)
    
    histos['histos_theta_p_combined'].Draw()
    histos['histos_theta_p_ctof'].Draw('same')
    #histos['histos_theta_p_tof_1'].Draw('same')
    histos['histos_theta_p_tof_2'].Draw('same')
    histos['histos_theta_p_tof_3'].Draw('same')
            
    if title:
        label.DrawLatex(0.1, 0.925, title)

    if xtitle:
        label.DrawLatex(0.5, 0.015, xtitle)

    if ytitle:
        label.SetTextAngle(90)
        label.DrawLatex(0.04, 0.5, ytitle)
        label.SetTextAngle(0)

    label.DrawLatex(0.75, 0.85, '#color[64]{CTOF}')
    label.DrawLatex(0.75, 0.80, '#color[55]{FTOF 1-B}')
    label.DrawLatex(0.75, 0.75, '#color[99]{FTOF 2}')
    canvas.Print(save_name)

    
if __name__ == '__main__':

    input_rootfile = 'histos.root'
    rootfile = TFile(input_rootfile)
    histos = load_histos(rootfile)
        
    setup_global_options() 

    can = TCanvas('can', 'can', 1100, 800)

    lab = TLatex()
    lab.SetNDC()
    lab.SetTextFont(42)
    lab.SetTextSize(0.05)
    lab.SetTextColor(1)

    plot_w(can, histos, label=lab, save_name='w_dist.pdf',
           xtitle='W (GeV/c^2)', title='Electrons (forward)')

    plot_theta_p(can, histos, label=lab, save_name='theta_p.pdf',
           xtitle='#theta_{p} (deg)', title='Angular Distribution of Elastic Protons')

    plot_sector_page(can, histos, 'histos_de_beam_de_beam_from_angles{}', label=lab,
                     save_name='de_beam_de_beam_from_angles.pdf',
                     xtitle='#Delta E_{beam} (#theta_{e}, p_{e})',
                     ytitle='#Delta E_{beam} (#theta_{e}, #theta_{p})',
                     title='#Delta E_{beam}', log=False)
    
    plot_beam_energy(can, histos, lab, 'de_beam.pdf')
