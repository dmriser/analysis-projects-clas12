#!/usr/bin/env python 

from ROOT import (TH1F, TH2F, TF1, TFile, TCanvas,
                  gPad, gStyle, TLatex, TLine)

def load_histos(file):
    ''' Use the ROOT file structure to load a dictionary of histograms. '''
    h = {}
    for k in file.GetListOfKeys():
        h[k.GetName()] = file.Get(k.GetName())
    return h

def setup_global_options():
    gStyle.SetOptTitle(0)
    gStyle.SetOptStat(0)

def plot_sector_beam(canvas, histos, title_formatter, label, save_name):
    
    canvas.Clear() 
    canvas.Divide(3,2)

    vert = TLine(0, -1.2, 0, 1.2)
    hori = TLine(-1.2, 0, 1.2, 0)

    vert.SetLineStyle(8)
    vert.SetLineStyle(1)
    vert.SetLineWidth(1)
    hori.SetLineStyle(8)
    hori.SetLineStyle(1)
    hori.SetLineWidth(1)
    
    xtitle = '#Delta E_{beam} (#theta_{e}, p_{e})'
    ytitle = '#Delta E_{beam} (#theta_{e}, #theta_{p})'

    for i in range(1,7):
        canvas.cd(i)
        
        histos[title_formatter.format(i)].GetXaxis().SetRangeUser(-1.2, 1.2)
        histos[title_formatter.format(i)].GetYaxis().SetRangeUser(-1.2, 1.2)
        histos[title_formatter.format(i)].Draw('colz')
        vert.Draw()
        hori.Draw()
        
        label.DrawLatex(0.1, 0.925, '#Delta E_{beam} = E_{beam} - E_{pred}')
        label.DrawLatex(0.5, 0.015, xtitle)
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
        
def plot_phi_vz(canvas, histos, label, save_name):
    canvas.Clear()

    # Start by styling the histograms
    histos['histos_phi_electron_delta_vz'].Draw('colz')

    line = TLine(-30, 0, 330, 0)
    line.SetLineColor(1)
    line.SetLineWidth(1)
    line.Draw()

    label.DrawLatex(0.45, 0.02, '#phi_{e}')
    #label.DrawLatex(0.35, 0.925, 'Elastic Vertex Difference')
    label.SetTextAngle(90)
    label.DrawLatex(0.04, 0.35, '#Delta v_{z} = v_{z} (e) - v_{z} (p) (cm)')
    label.SetTextAngle(0)
    
    canvas.Print(save_name)
    
        
def plot_event_selection(canvas, histos, save_name, label):
    
    canvas.Clear() 
    canvas.Divide(2,1)

    histos['histos_w_in_ctof'].SetLineColor(1)
    histos['histos_w_pass_angle_in_ctof'].SetLineColor(1)
    histos['histos_w_pass_angle_in_ctof'].SetFillColorAlpha(64,1.0)

    histos['histos_angle_ep'].SetLineColor(1)
    histos['histos_angle_ep_pass_w_in_ctof'].SetLineColor(1)
    histos['histos_angle_ep_pass_w_in_ctof'].SetFillColorAlpha(64,1.0)

    histos['histos_w_in_ctof'].SetMinimum(0)
    histos['histos_w_pass_angle_in_ctof'].SetMinimum(0)
    
    # Fix missing bins
    histos['histos_w_in_ctof'].GetXaxis().SetRangeUser(0.7, 1.3)
    histos['histos_w_pass_angle_in_ctof'].GetXaxis().SetRangeUser(0.7, 1.3)

    wleft = TLine(0.8, histos['histos_w_in_ctof'].GetMinimum(),
                  0.8, 0.8 * histos['histos_w_in_ctof'].GetMaximum())
    wright = TLine(1.08, histos['histos_w_in_ctof'].GetMinimum(),
                   1.08, 0.8 * histos['histos_w_in_ctof'].GetMaximum())

    wleft.SetLineColor(99)
    wleft.SetLineWidth(1)
    wright.SetLineColor(99)
    wright.SetLineWidth(1)
    
    canvas.cd(1)
    histos['histos_w_in_ctof'].Draw()
    histos['histos_w_pass_angle_in_ctof'].Draw('same')
    wleft.Draw()
    wright.Draw()
    
    label.DrawLatex(0.1, 0.925, 'Electron (forward) and Proton (central)')
    label.DrawLatex(0.5, 0.015, 'W (GeV/c^{2})')
    label.DrawLatex(0.67, 0.86, '#color[64]{#phi_{ep} > 177}')
    
    canvas.cd(2)
    histos['histos_angle_ep'].Draw()
    histos['histos_angle_ep_pass_w_in_ctof'].Draw('same')

    left = TLine(177, 0.0,
                 177, 0.8 * histos['histos_angle_ep'].GetMaximum())
    left.SetLineColor(99)
    left.SetLineWidth(1)
    left.Draw()
    
    label.DrawLatex(0.1, 0.925, 'Electron (forward) and Proton (central)')
    label.DrawLatex(0.5, 0.015, '#phi_{ep}')
    label.DrawLatex(0.57, 0.86, '#color[64]{W \in [0.8, 1.08]}')
    
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

    plot_event_selection(can, histos, label=lab, save_name='event_selection.pdf')

    plot_theta_p(can, histos, label=lab, save_name='theta_p.pdf',
           xtitle='#theta_{p} (deg)', title=None)

    plot_sector_beam(can, histos, 'histos_de_beam_de_beam_from_angles{}', label=lab,
                     save_name='de_beam_de_beam_from_angles.pdf')
    
    plot_beam_energy(can, histos, lab, 'de_beam.pdf')

    plot_phi_vz(can, histos, lab, 'phi_vz.pdf')
    
