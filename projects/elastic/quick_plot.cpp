void quick_plot(){
  TFile *inputFile = TFile::Open("mon-lite.hipo.root");

  TH1F *h = (TH1F*) inputFile->Get("histos_angle_ep_gen");
  
  TCanvas *can = new TCanvas("can", "can", 800, 600);
  h->SetAxisRange(0, 1400, "Y");
  h->Draw();

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.03);
  latex->DrawLatex(0.4, 0.02, "#phi_{ep} (generated)");

  can->Print("angle_ep_gen.pdf");
}
