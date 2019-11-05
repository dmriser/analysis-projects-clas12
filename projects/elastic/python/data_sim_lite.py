#!/usr/bin/env python 

import argparse 
import numpy as np

# Trick for docker install to run headless
# and never use plt.show() 
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt 

from scipy.optimize import minimize 

from array import array 
from ROOT import (TH1F, TH2F, TH1D, TF1, TFile, TCanvas,
                  gPad, gStyle, TLatex, TGraphErrors)

allowed_types = (TH1F, TH1D)
default_histo = TH1F('default', '', 100, 0, 1)
default_histo2d = TH2F('default', '', 100, 0, 1, 100, 0, 1)

def load_histos(file):
    ''' Use the ROOT file structure to load a dictionary of histograms. '''
    h = {}
    for k in file.GetListOfKeys():
        h[k.GetName()] = file.Get(k.GetName())
    return h

def numpify(histo):
    """ TH1F to np.arrays. """

    if type(histo) not in allowed_types:
        raise NotImplementedError('Can not numpify type {}'.format(type(histo)))

    nbins = histo.GetNbinsX()

    # Setup output 
    x_lows = np.zeros(nbins)
    x_highs = np.zeros(nbins)
    values = np.zeros(nbins)
    errors = np.zeros(nbins)
    
    for i in range(1, nbins + 1):
        x_lows[i-1] = histo.GetBinLowEdge(i)
        x_highs[i-1] = histo.GetBinLowEdge(i) + histo.GetBinWidth(i)
        values[i-1] = histo.GetBinContent(i)
        errors[i-1] = histo.GetBinError(i)

    return x_lows, x_highs, values, errors

def chi2(data, theory, err):
    return np.sum((data-theory)**2 / (0.00001 + err**2)) 

def model(x,p):
    return p[0] * np.exp( -0.5 * (x - p[1])**2 / p[2]**2 )

def scipy_fit_slice(x, y, err, bounds):        
    
    p0 = np.random.uniform(-1, 1, 3)    
    p0[0] = np.max(y)
    p0[1] = np.mean(x*y)
    p0[2] = np.std(x*y)
    
    res = minimize(lambda p: chi2(y, model(x,p), err), x0=p0, bounds=bounds, method='Nelder-Mead')
    return res.x

def fit_slices(histo, x_range, x_bin_step):

    x_start = histo.GetXaxis().FindBin(x_range[0])
    x_stop =  histo.GetXaxis().FindBin(x_range[1])

    x_values = array('d')
    slices = []
    fits = []
    for i, x_bin in enumerate(range(x_start, x_stop + 1, x_bin_step)):
        projec = histo.ProjectionY(histo.GetTitle() + '_proj{}'.format(i) , x_bin, x_bin + x_bin_step)
        
        fmin = projec.GetMean() - 3 * projec.GetStdDev()
        fmax = projec.GetMean() + 3 * projec.GetStdDev() 
        fit = TF1(histo.GetTitle() + '_fit{}'.format(i), 'gaus', fmin, fmax)
        fit.SetTitle(histo.GetTitle() + '_fit{}'.format(i))
        fit.SetName(histo.GetTitle() + '_fit{}'.format(i))
        projec.Fit(fit, 'R', '',  fmin, fmax)
        print('Fitting in {},{}'.format(fmin,fmax))
        
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

def plot_slices(slices, fits, output_name):
    """ Plot slices with fit values. """

    ncols = 3
    nrows = int(np.ceil(len(slices) / ncols) + 1)
    
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                             figsize=(ncols * 4, nrows * 3))

    i,j = 0, 0
    for current_slice, current_fit in zip(slices, fits):
        x_low, x_high, vals, errs = numpify(current_slice)
        x_cent = 0.5 * (x_low + x_high)
        width = x_high[0] - x_low[0]
        #y_fit = np.array([current_fit.Eval(xc) for xc in x_cent])

        mu = np.mean(x_cent * vals)
        std = np.std(x_cent * vals)

        bounds = [] 
        bounds.append([0, 2*np.max(vals)])
        bounds.append([mu - 2 * std, mu + 2 * std])
        bounds.append([0, 5 * std])

        pars = scipy_fit_slice(x_cent, vals, errs, bounds)
        y_fit = model(x_cent, pars)
        print(pars)
         
        axes[i,j].bar(x_cent, vals, width=width, edgecolor='')
        axes[i,j].plot(x_cent, y_fit, linestyle='-', linewidth=1, color='red')
        axes[i,j].set_xlim([mu-1*std, mu+1*std])
        axes[i,j].grid(alpha=0.2)
        axes[i,j].set_title(current_fit.GetTitle())
        
        j += 1
        if j >= ncols:
            j = 0
            i += 1 

    fig.tight_layout()
    fig.savefig(output_name, bbox_inches='tight')
        
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
        print('Fitting: ', tit)
        x1, mu1, sig1, slices1, fits1 = fit_slices(histos1.get(tit, default_histo2d), x_range, x_bin_step)
        x2, mu2, sig2, slices2, fits2 = fit_slices(histos2.get(tit, default_histo2d), x_range, x_bin_step)
        x1, mu1, sig1 = remove_bad_points(x1, mu1, sig1, max_errorbar)
        x2, mu2, sig2 = remove_bad_points(x2, mu2, sig2, max_errorbar)
        
        # Experimental
        plot_slices(slices1, fits1, tit + '_data_slices.pdf')

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
    """
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
    """
    #plot_fits_mpl(histos1=histos['data'], histos2=histos['sim'], config1='Data', config2='Sim',
    #    x_range=[8,11], y_range=[-0.8,0.8], x_bin_step=6, title_formatter='histos_p_electron_delta_p_electron_{}',
    #    save_name='p_electron_delta_p_electron_fit_{}.pdf'.format(args.output_prefix),
    #    title='Electron Momentum Resolution (from $\\theta_e$)', xtitle='$P_e$', ytitle='$\Delta P_{e}$',
    #    max_errorbar = 0.8
    #) 
 
    #plot_fits_mpl(histos1=histos['data'], histos2=histos['sim'], config1='Data', config2='Sim',
    #    x_range=[44,55], y_range=[-5,5], x_bin_step=6, title_formatter='histos_theta_proton_delta_theta_proton_{}',
    #    save_name='theta_proton_delta_theta_proton_fit_{}.pdf'.format(args.output_prefix),
    #    title='Proton $\\theta$ Resolution (from $\\theta_e$)', xtitle='$\\theta_p$', ytitle='$\Delta \\theta_{p}$',
    #    max_errorbar = 3
    #) 

    #plot_fits_mpl(histos1=histos['data'], histos2=histos['sim'], config1='Data', config2='Sim',
    #    x_range=[40,55], y_range=[-5,5], x_bin_step=6, title_formatter='histos_theta_electron_delta_theta_electron_{}',
    #    save_name='theta_electron_delta_theta_electron_fit_{}.pdf'.format(args.output_prefix),
    #    title='Electron $\\theta$ Resolution (from $\\theta_e$)', xtitle='$\\theta_e$', ytitle='$\Delta \\theta_{e}$',
    #    max_errorbar = 3
    #) 

    #plot_fits_mpl(histos1=histos['data'], histos2=histos['sim'], config1='Data', config2='Sim',
    #    x_range=[6.5,11.5], y_range=[-0.5,0.5], x_bin_step=6, title_formatter='histos_theta_ele_de_beam_{}',
    #    save_name='theta_electron_de_beam_fit_{}.pdf'.format(args.output_prefix),
    #    title='Beam Energy (from $\\theta_e$, $P_e$)', xtitle='$\\theta_e$', ytitle='$\Delta E_{beam}$',
    #    max_errorbar = 3
    #) 

    print(50 * '-')
    print(histos['data'].keys())
    print(50 * '-')
    print(histos['sim'].keys())
    
    plt.figure(figsize=(16,9))
    for i in range(1,7):
        title = 'histos_p_ele_ctof_{}'.format(i)

        x_low_dat, x_high_dat, val_dat, err_dat = numpify(histos['data'][title])
        x_low_sim, x_high_sim, val_sim, err_sim = numpify(histos['sim'][title])
        x_cent_dat = 0.5 * (x_high_dat + x_low_dat)
        x_cent_sim = 0.5 * (x_high_sim + x_low_sim)
        width_dat = x_high_dat[1] - x_high_dat[0]
        width_sim = x_high_sim[1] - x_high_sim[0]
        
        val_dat /= np.max(val_dat)
        val_sim /= np.max(val_sim)

        plt.subplot(2,3,i)
        plt.step(x_low_dat, val_dat, where='post', label='Data', color='k')
        plt.step(x_low_sim, val_sim, where='post', label='Sim', color='red')
        plt.xlabel('$p_e$')
        plt.grid(alpha=0.2)
        plt.legend(frameon=False)
         
    plt.tight_layout()
    plt.savefig('p_ele_compare-in.pdf')

    plt.close()
    plt.figure(figsize=(16,9))
    for i in range(1,7):
        title = 'histos_theta_electron_ctof_{}'.format(i)

        x_low_dat, x_high_dat, val_dat, err_dat = numpify(histos['data'][title])
        x_low_sim, x_high_sim, val_sim, err_sim = numpify(histos['sim'][title])
        x_cent_dat = 0.5 * (x_high_dat + x_low_dat)
        x_cent_sim = 0.5 * (x_high_sim + x_low_sim)
        width_dat = x_high_dat[1] - x_high_dat[0]
        width_sim = x_high_sim[1] - x_high_sim[0]
        
        val_dat /= np.max(val_dat)
        val_sim /= np.max(val_sim)

        plt.subplot(2,3,i)
        plt.step(x_low_dat, val_dat, where='post', label='Data', color='k')
        plt.step(x_low_sim, val_sim, where='post', label='Sim', color='red')
        plt.xlabel('$\\theta_e$')
        plt.grid(alpha=0.2)
        plt.legend(frameon=False)
        
    plt.tight_layout()
    plt.savefig('theta_ele_compare-in.pdf')

    
    plt.close()
    plt.figure(figsize=(16,9))
    for i in range(1,7):
        title = 'histos_theta_proton_ctof_{}'.format(i)

        x_low_dat, x_high_dat, val_dat, err_dat = numpify(histos['data'][title])
        x_low_sim, x_high_sim, val_sim, err_sim = numpify(histos['sim'][title])
        x_cent_dat = 0.5 * (x_high_dat + x_low_dat)
        x_cent_sim = 0.5 * (x_high_sim + x_low_sim)
        width_dat = x_high_dat[1] - x_high_dat[0]
        width_sim = x_high_sim[1] - x_high_sim[0]
        
        val_dat /= np.max(val_dat)
        val_sim /= np.max(val_sim)

        plt.subplot(2,3,i)
        plt.step(x_low_dat, val_dat, where='post', label='Data', color='k')
        plt.step(x_low_sim, val_sim, where='post', label='Sim', color='red')
        plt.xlabel('$\\theta_p$')
        plt.grid(alpha=0.2)
        plt.legend(frameon=False)
        
    plt.tight_layout()
    plt.savefig('theta_pro_compare-in.pdf')

    
    plt.close()
    plt.figure(figsize=(16,9))
    for i in range(1,7):
        title = 'histos_p_pro_ctof_{}'.format(i)

        x_low_dat, x_high_dat, val_dat, err_dat = numpify(histos['data'][title])
        x_low_sim, x_high_sim, val_sim, err_sim = numpify(histos['sim'][title])
        x_cent_dat = 0.5 * (x_high_dat + x_low_dat)
        x_cent_sim = 0.5 * (x_high_sim + x_low_sim)
        width_dat = x_high_dat[1] - x_high_dat[0]
        width_sim = x_high_sim[1] - x_high_sim[0]
        
        val_dat /= np.max(val_dat)
        val_sim /= np.max(val_sim)

        plt.subplot(2,3,i)
        plt.step(x_low_dat, val_dat, where='post', label='Data', color='k')
        plt.step(x_low_sim, val_sim, where='post', label='Sim', color='red')
        plt.xlabel('$p_p$')
        plt.grid(alpha=0.2)
        plt.legend(frameon=False)
        
    plt.tight_layout()
    plt.savefig('p_pro_compare-in.pdf')
