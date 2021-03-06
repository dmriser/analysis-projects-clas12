#!/usr/bin/env python 

from array import array 
import argparse
from ROOT import (TH1F, TH2F, TF1, TFile, TCanvas,
                  gPad, gStyle, TLatex, TLine, TGraphErrors,
                  TVector)


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

def plot_reso_theta(canvas, histos, title_formatter,
                    label, save_name, xtitle=None, ytitle=None, title=None,
                    theta_range=None):
    
    canvas.Clear() 
    canvas.Divide(3,2)
    
    for i in range(1,7):
        canvas.cd(i)

        if theta_range:
            histos[title_formatter.format(i)].GetXaxis().SetRangeUser(
                theta_range[0], theta_range[1])
        
        histos[title_formatter.format(i)].Draw('colz')

        if title:
            label.DrawLatex(0.1, 0.925, '#Delta E_{beam} = E_{beam} - E_{pred}')

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
    label.DrawLatex(0.67, 0.86, '#color[64]{#phi_{ep} > 174}')
    
    canvas.cd(2)
    histos['histos_angle_ep'].Draw()
    histos['histos_angle_ep_pass_w_in_ctof'].Draw('same')

    left = TLine(174, 0.0,
                 174, 0.8 * histos['histos_angle_ep'].GetMaximum())
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

def fit_slices(histo, x_range, x_bin_step):
    x_start = histo.GetXaxis().FindBin(x_range[0])
    x_stop =  histo.GetXaxis().FindBin(x_range[1])

    x_values = array('d')
    slices = [] 
    fits = [] 
    for i, x_bin in enumerate(range(x_start, x_stop + 1, x_bin_step)):
        projec = histo.ProjectionY(histo.GetTitle() + '_proj{}'.format(i) , x_bin, x_bin + x_bin_step)
        fit = TF1(histo.GetTitle()+'_fit{}'.format(i), 'gaus')

        projec.Fit(fit)
        
        slices.append(projec)
        fits.append(fit)
        
        x_low = histo.GetXaxis().GetBinCenter(x_bin)
        x_high = histo.GetXaxis().GetBinCenter(x_bin + x_bin_step) 
        x_values.append(0.5 * (x_high + x_low))
        


    means = array('d')
    means_err = array('d')
    stds = array('d')
    stds_err = array('d')
    zeros = array('d')

    for f in fits:
        means.append(f.GetParameter(1))
        means_err.append(f.GetParError(1))
        stds.append(f.GetParameter(2))
        stds_err.append(f.GetParError(2))
        zeros.append(0.0)

    print(x_values)
    print(means)
    graph_mean = TGraphErrors(len(x_values), x_values, means, zeros, means_err)
    #graph_mean.SetTitle(histo.GetTitle() + '_mean')
    
    graph_std = TGraphErrors(len(x_values), x_values, stds, zeros, stds_err)
    #graph_std.SetTitle(histo.GetTitle() + '_std')
    
    return graph_mean, graph_std
        
def plot_fits(canvas, histos, x_range, x_bin_step, title_formatter,
              save_name, label, y_range=None,
              title=None, xtitle=None, ytitle=None, draw_mean=False):

    canvas.Clear() 
    canvas.Divide(3,2)


    root_is_dumb = [] 
    for i in range(1,7):
        canvas.cd(i)
        
        title = title_formatter.format(i)
        mean, std = fit_slices(histos[title], x_range, x_bin_step)
        
        if draw_mean:
            mean.SetMarkerStyle(21)
            mean.SetMarkerSize(1)

            if y_range:
                mean.GetYaxis().SetLimits(y_range[0], y_range[1])

            mean.Draw('AP')
            root_is_dumb.append(mean)

        else:
            std.SetMarkerStyle(21)
            std.SetMarkerSize(1)
            std.Draw('AP')

            if y_range:
                std.GetYaxis().SetLimits(y_range[0], y_range[1])

            root_is_dumb.append(std)

        if title:
            label.DrawLatex(0.1, 0.925, title)

        if xtitle:
            label.DrawLatex(0.5, 0.015, xtitle)

        if ytitle:
            label.SetTextAngle(90)
            label.DrawLatex(0.035, 0.5, ytitle)
            label.SetTextAngle(0)

            
    canvas.Print(save_name)
        
        
if __name__ == '__main__':

    ap = argparse.ArgumentParser()
    ap.add_argument(
        '-i',
        '--input_file',
	required=True
    )
    ap.add_argument(
        '-o',
        '--output_prefix',
        required=True
    )
    args = ap.parse_args()

    input_rootfile = args.input_file
    output_pdfname = args.output_prefix + '.pdf'
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
    
    plot_reso_theta(can, histos, 'histos_theta_electron_delta_p_electron_{}', label=lab,
                    save_name='theta_electron_delta_p_electron.pdf',
                    xtitle='#theta_{e} (deg)',
                    ytitle='#Delta p_{e}',
                    theta_range=[6.0,13.0])

    plot_reso_theta(can, histos, 'histos_theta_proton_delta_p_proton_{}', label=lab,
                    save_name='theta_proton_delta_p_proton.pdf',
                    xtitle='#theta_{p} (deg)',
                    ytitle='#Delta p_{p}',
                    theta_range=[35.0, 60.0])

    plot_fits(canvas=can, histos=histos,
              title_formatter='histos_theta_electron_delta_p_electron_{}', label=lab,
              x_range=[7.0, 11.0], x_bin_step=10, y_range=[-0.2, 0.2],
              save_name='graph_theta_electron_delta_p_electron_reso.pdf',
              draw_mean=False, xtitle='#theta_{e} (deg)', ytitle='p_{e} (resolution)')

    plot_fits(canvas=can, histos=histos,
              title_formatter='histos_theta_proton_delta_p_proton_{}', label=lab,
              x_range=[40.0, 55.0], x_bin_step=10, y_range=[-0.2, 0.2],
              save_name='graph_theta_proton_delta_p_proton_reso.pdf',
              draw_mean=False, xtitle='#theta_{p} (deg)', ytitle='p_{p} (resolution)')
