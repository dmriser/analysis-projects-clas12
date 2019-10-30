#!/usr/bin/env python 

import argparse 
import numpy as np

# Trick for docker install to run headless
# and never use plt.show() 
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt 


from array import array 
from ROOT import (TH1F, TH2F, TF1, TFile, TCanvas,
                  gPad, gStyle, TLatex, TGraphErrors)

default_histo = TH1F('default', '', 100, 0, 1)
default_histo2d = TH2F('default', '', 100, 0, 1, 100, 0, 1)

def load_histos(file):
    ''' Use the ROOT file structure to load a dictionary of histograms. '''
    h = {}
    for k in file.GetListOfKeys():
        h[k.GetName()] = file.Get(k.GetName())
    return h

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

    graph = TGraphErrors(len(x_values), x_values, means, zeros, stds)
    graph.SetName('g_' + histo.GetName())

    # return graph
    return np.array(x_values), np.array(means), np.array(stds), slices, fits 

def plot_fits_mpl(histos1, histos2, config1, config2,
                  x_range, x_bin_step, title_formatter,
                  save_name, y_range=None, title=None,
                  xtitle=None, ytitle=None, max_errorbar=0.5):
    
    fig = plt.figure(figsize=(12,16))

    # Plot options for each type. 
    opts1 = {'marker':'o', 'linestyle':'', 'color':'k'}
    opts2 = {'marker':'o', 'linestyle':'', 'color':'red'}
    
    for i in range(1,7):
        ax = fig.add_subplot(3, 2, i)

        # Get histogram slices for plotting. 
        tit = title_formatter.format(i)
        x1, mu1, sig1, slices1, fits1 = fit_slices(histos1.get(tit, default_histo2d), x_range, x_bin_step)
        x2, mu2, sig2, slices2, fits2 = fit_slices(histos2.get(tit, default_histo2d), x_range, x_bin_step)
        x1, mu1, sig1 = remove_bad_points(x1, mu1, sig1, max_errorbar)
        x2, mu2, sig2 = remove_bad_points(x2, mu2, sig2, max_errorbar)

        label1 = 'Sector {} ({})'.format(i, config1)
        label2 = 'Sector {} ({})'.format(i, config2)

        # Draw things 
        ax.errorbar(x1, mu1, sig1, label=label1, **opts1)
        ax.errorbar(x2, mu2, sig2, label=label2, **opts2)

        if y_range:
            ax.set_ylim(y_range)

        # Add center line
        ax.axhline(0.0, linestyle='--', linewidth=1, color='k', alpha=0.85)

        # Add a grid
        ax.grid(alpha=0.2)

        # Legend
        ax.legend(frameon=False)
        
        if title:
            ax.set_title(title)

        if xtitle:
            ax.set_xlabel(xtitle)

        if ytitle:
            ax.set_ylabel(ytitle)
            
    fig.tight_layout()
    fig.savefig(save_name, bbox_inches='tight')

def remove_bad_points(x, mu, sig, max_errorbar):
    condition = np.logical_and(mu != 0, sig < max_errorbar)
    idx = np.where(condition)[0]
    return x[idx], mu[idx], sig[idx]
    
if __name__ == '__main__':

    # Parse command line arguments. 
    ap = argparse.ArgumentParser()
    ap.add_argument('-d', '--data_file', required=True)
    ap.add_argument('-s', '--sim_file', required=True)
    ap.add_argument('-o', '--output_prefix', required=True)
    args = ap.parse_args()

    # Setup files
    files = {}
    files['data'] = TFile(args.data_file)
    files['sim'] = TFile(args.sim_file)

    output_pdfname = args.output_prefix + '.pdf'

    # Load histograms
    histos = {}
    for config_type, file in files.items():
        histos[config_type] = load_histos(file)
    
    # Plot fits to the resolutions
    plot_fits_mpl(histos1=histos['data'], histos2=histos['sim'], config1='Data', config2='Sim',
        x_range=[6,12], y_range=[-0.8,0.8], x_bin_step=5, title_formatter='histos_theta_electron_delta_p_electron_{}',
        save_name='theta_electron_delta_p_electron_fit_{}.pdf'.format(args.output_prefix),
        title='Electron Momentum Resolution (from $\\theta_e$)', xtitle='$\\theta_e$', ytitle='$\Delta P_{e}$',
        max_errorbar = 0.4
    ) 

    plot_fits_mpl(histos1=histos['data'], histos2=histos['sim'], config1='Data', config2='Sim',
        x_range=[40,55], y_range=[-0.8,0.8], x_bin_step=6, title_formatter='histos_theta_proton_delta_p_proton_{}',
        save_name='theta_proton_delta_p_proton_fit_{}.pdf'.format(args.output_prefix),
        title='Proton Momentum Resolution (from $\\theta_e$)', xtitle='$\\theta_p$', ytitle='$\Delta P_{p}$',
        max_errorbar = 0.8
    ) 
 
    plot_fits_mpl(histos1=histos['data'], histos2=histos['sim'], config1='Data', config2='Sim',
        x_range=[0.3,2.5], y_range=[-0.8,0.8], x_bin_step=6, title_formatter='histos_p_proton_delta_p_proton_{}',
        save_name='p_proton_delta_p_proton_fit_{}.pdf'.format(args.output_prefix),
        title='Proton Momentum Resolution (from $\\theta_e$)', xtitle='$P_p$', ytitle='$\Delta P_{p}$',
        max_errorbar = 0.8
    ) 

    #plot_fits_mpl(histos1=histos['data'], histos2=histos['sim'], config1='Data', config2='Sim',
    #    x_range=[8,11], y_range=[-0.8,0.8], x_bin_step=6, title_formatter='histos_p_electron_delta_p_electron_{}',
    #    save_name='p_electron_delta_p_electron_fit_{}.pdf'.format(args.output_prefix),
    #    title='Electron Momentum Resolution (from $\\theta_e$)', xtitle='$P_e$', ytitle='$\Delta P_{e}$',
    #    max_errorbar = 0.8
    #) 
 
    plot_fits_mpl(histos1=histos['data'], histos2=histos['sim'], config1='Data', config2='Sim',
        x_range=[44,55], y_range=[-5,5], x_bin_step=6, title_formatter='histos_theta_proton_delta_theta_proton_{}',
        save_name='theta_proton_delta_theta_proton_fit_{}.pdf'.format(args.output_prefix),
        title='Proton $\\theta$ Resolution (from $\\theta_e$)', xtitle='$\\theta_p$', ytitle='$\Delta \\theta_{p}$',
        max_errorbar = 3
    ) 

    #plot_fits_mpl(histos1=histos['data'], histos2=histos['sim'], config1='Data', config2='Sim',
    #    x_range=[40,55], y_range=[-5,5], x_bin_step=6, title_formatter='histos_theta_electron_delta_theta_electron_{}',
    #    save_name='theta_electron_delta_theta_electron_fit_{}.pdf'.format(args.output_prefix),
    #    title='Electron $\\theta$ Resolution (from $\\theta_e$)', xtitle='$\\theta_e$', ytitle='$\Delta \\theta_{e}$',
    #    max_errorbar = 3
    #) 

    plot_fits_mpl(histos1=histos['data'], histos2=histos['sim'], config1='Data', config2='Sim',
        x_range=[6.5,11.5], y_range=[-0.5,0.5], x_bin_step=6, title_formatter='histos_theta_ele_de_beam_{}',
        save_name='theta_electron_de_beam_fit_{}.pdf'.format(args.output_prefix),
        title='Beam Energy (from $\\theta_e$, $P_e$)', xtitle='$\\theta_e$', ytitle='$\Delta E_{beam}$',
        max_errorbar = 3
    ) 
