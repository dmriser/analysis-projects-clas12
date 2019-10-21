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
    canvas.Divide(2,3)

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
            label.DrawLatex(0.04, 0.5, ytitle)
            label.SetTextAngle(0)

    canvas.Print(save_name)

def plot_page(canvas, histos, histo_title, label, save_name,
                     xtitle=None, ytitle=None, title=None, log=False):
    
    canvas.Clear() 
        
    if isinstance(histos[histo_title], TH1F):
        histos[histo_title].Draw()
    elif isinstance(histos[histo_title], TH2F):
        histos[histo_title].Draw('colz')
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
        label.DrawLatex(0.04, 0.5, ytitle)
        label.SetTextAngle(0)

    canvas.Print(save_name)

    
if __name__ == '__main__':

    input_rootfile = 'histos.root'
    output_pdfname = 'elastic_10_22_19.pdf'
    rootfile = TFile(input_rootfile)
    histos = load_histos(rootfile)

    for k,v in histos.items():
        print(k, v)
    
        
    setup_global_options() 

    can = TCanvas('can', 'can', 800, 1100)

    lab = TLatex()
    lab.SetNDC()
    lab.SetTextFont(42)
    lab.SetTextSize(0.05)
    lab.SetTextColor(1)

    can.Print('{}['.format(output_pdfname))

    plot_sector_page(can, histos, 'histos_w_{}', lab, save_name=output_pdfname,
                     title='Electron (Forward) and Positive (CTOF)', xtitle='W')

    plot_sector_page(can, histos, 'histos_w_q2_{}', lab, save_name=output_pdfname,
                     title='Electron (Forward) and Positive (CTOF)', xtitle='W', ytitle='Q^{2}', log=False)

    plot_page(can, histos, 'histos_phi_electron_w', lab, save_name=output_pdfname,
                     title='W vs. #phi_{e}', xtitle='#phi_{e}', ytitle='W', log=False)

    plot_sector_page(can, histos, 'histos_theta_electron_vz_electron_{}', lab, save_name=output_pdfname,
                     title='v_{z} (e) vs. #theta_{e}', xtitle='#theta_{e}', ytitle='v_{z} (e)', log=False)

    plot_sector_page(can, histos, 'histos_w_inclusive_{}', lab, save_name=output_pdfname,
                     title='Electron (Forward)', xtitle='W')

    plot_sector_page(can, histos, 'histos_w_q2_inclusive_{}', lab, save_name=output_pdfname,
                     title='Electron (Forward)', xtitle='W', ytitle='Q^{2}', log=False)

    plot_sector_page(can, histos, 'histos_delta_p_electron_{}', lab, save_name=output_pdfname,
                     title='#Delta P_{e} from #theta_{e}', xtitle='#Delta P_{e}', log=False)

    plot_sector_page(can, histos, 'histos_theta_electron_delta_p_electron_{}', lab, save_name=output_pdfname,
                     title='#Delta P_{e} vs #theta_{e} from #theta_{e}', xtitle='#theta_{e}', ytitle='#Delta P_{e}', log=False)

    plot_page(can, histos, 'histos_phi_electron_delta_p_electron', lab, save_name=output_pdfname,
                     title='#Delta P_{e} vs. #phi_{e} from #theta_{e}', xtitle='#phi_{e}', ytitle='#Delta P_{e}', log=False)

    plot_sector_page(can, histos, 'histos_delta_p_proton_{}', lab, save_name=output_pdfname,
                     title='#Delta P_{p} from #theta_{e}', xtitle='#Delta P_{p}', log=False)

    plot_sector_page(can, histos, 'histos_theta_proton_delta_p_proton_{}', lab, save_name=output_pdfname,
                     title='#Delta P_{p} vs #theta_{p} from #theta_{e}', xtitle='#theta_{p}', ytitle='#Delta P_{p}', log=False)

    plot_sector_page(can, histos, 'histos_p_proton_delta_p_proton_{}', lab, save_name=output_pdfname,
                     title='#Delta P_{p} vs P_{p} from #theta_{e}', xtitle='P_{p}', ytitle='#Delta P_{p}', log=False)
    
    plot_sector_page(can, histos, 'histos_delta_theta_proton_{}', lab, save_name=output_pdfname,
                     title='#Delta #theta_{p} from #theta_{e}', xtitle='#Delta #theta_{p}', log=False)

    plot_sector_page(can, histos, 'histos_theta_electron_delta_theta_proton_{}', lab, save_name=output_pdfname,
                     title='#Delta #theta_{p} vs #theta_{e} from #theta_{e}', xtitle='#theta_{e}', ytitle='#Delta #theta_{p}', log=False)

    plot_sector_page(can, histos, 'histos_theta_proton_delta_theta_proton_{}', lab, save_name=output_pdfname,
                     title='#Delta #theta_{p} vs #theta_{p} from #theta_{e}', xtitle='#theta_{p}', ytitle='#Delta #theta_{p}', log=False)

    plot_sector_page(can, histos, 'histos_theta_electron_vz_electron_{}', lab, save_name=output_pdfname,
                     title='v_{z} (e) vs #theta_{e}', xtitle='#theta_{e}', ytitle='v_{z} (e)', log=False)

    plot_page(can, histos, 'histos_phi_electron_vz_electron', lab, save_name=output_pdfname,
                     title='v_{z} (e) vs. #phi_{e}', xtitle='#phi_{e}', ytitle='v_{z} (e)', log=False)

    plot_sector_page(can, histos, 'histos_vz_electron_{}', lab, save_name=output_pdfname,
                     title='v_{z} (e)', xtitle='v_{z} (e)', log=False)

    plot_sector_page(can, histos, 'histos_theta_proton_vz_proton_{}', lab, save_name=output_pdfname,
                     title='v_{z} (p) vs #theta_{p}', xtitle='#theta_{p}', ytitle='v_{z} (p)', log=False)

    plot_page(can, histos, 'histos_phi_proton_vz_proton', lab, save_name=output_pdfname,
                     title='v_{z} (p) vs. #phi_{p}', xtitle='#phi_{p}', ytitle='v_{z} (p)', log=False)

    plot_sector_page(can, histos, 'histos_vz_proton_{}', lab, save_name=output_pdfname,
                     title='v_{z} (p)', xtitle='v_{z} (p)', log=False)

    plot_page(can, histos, 'histos_phi_electron_delta_vz', lab, save_name=output_pdfname,
                     title='#Delta v_{z} vs. #phi_{e}', xtitle='#phi_{e}', ytitle='#Delta v_{z}', log=False)

    plot_page(can, histos, 'histos_phi_proton_delta_vz', lab, save_name=output_pdfname,
                     title='#Delta v_{z} vs. #phi_{p}', xtitle='#phi_{p}', ytitle='#Delta v_{z}', log=False)

    plot_sector_page(can, histos, 'histos_delta_vz_{}', lab, save_name=output_pdfname,
                     title='#Delta v_{z}', xtitle='#Delta v_{z}', log=False)

    plot_sector_page(can, histos, 'histos_de_beam_{}', lab, save_name=output_pdfname,
                     title='#Delta E_{beam} from (#theta_{e}, P_{e})', xtitle='#Delta E_{beam}', log=False)

    plot_sector_page(can, histos, 'histos_de_beam_from_angles{}', lab, save_name=output_pdfname,
                     title='#Delta E_{beam} from (#theta_{e}, #theta_{p})', xtitle='#Delta E_{beam}', log=False)

    plot_sector_page(can, histos, 'histos_de_beam_de_beam_from_angles{}', lab, save_name=output_pdfname,
                     title='#Delta E_{beam}', xtitle='#Delta E (#theta_{e}, P_{e})',
                     ytitle='#Delta E (#theta_{e}, #theta_{p})', log=False)

    plot_sector_page(can, histos, 'histos_theta_ele_de_beam_{}', lab, save_name=output_pdfname,
                     title='#Delta E_{beam} vs #theta_{e}', xtitle='#theta_{e}',
                     ytitle='#Delta E (#theta_{e}, P_{e})', log=False)

    # Close the sucker 
    can.Print('{}]'.format(output_pdfname))
