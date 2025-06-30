void make_all_hists()
{
  TFile *f[5];
  f[0] = new TFile("lime_auau_hists_nominal.root", "r");
  f[1] = new TFile("lime_auau_hists_39mb.root", "r");
  f[2] = new TFile("lime_auau_hists_45mb.root", "r");
  f[3] = new TFile("lime_auau_hists_A1.root", "r");
  f[4] = new TFile("lime_auau_hists_A2.root", "r");

  TH1D *h_system[5];
  for (int i = 0; i < 5; i++)
    {
      h_system[i] = (TH1D*) f[i]->Get("hNpart");
      h_system[i]->SetName(Form("hNpart_%d", i));
      
    }
  h_system[0]->SetLineColor(kBlack);
  h_system[1]->SetLineColor(kBlue);
  h_system[2]->SetLineColor(kRed);
  h_system[3]->SetLineColor(kOrange);
  h_system[4]->SetLineColor(kGreen);
  
  TCanvas *c = new TCanvas("c","c", 500, 500);

  c->SetLogy();
  c->SetTicks(1,1);

  std::string words[5] = {"#sigma_{NN} = 42 mb",
			  "#sigma_{NN} = 39 mb",
			  "#sigma_{NN} = 45 mb",
			  "r = 6.65 fm",
			  "r = 6.25 fm"};
  h_system[0]->Draw("hist");
  for (int i = 1; i < 5; i++)
    {
      h_system[i]->Draw("hist same");
    }
  TLegend *l = new TLegend(0.5, 0.63, 0.7, 0.88);
  l->SetTextSize(0.04);
  l->SetTextFont(42);
  l->SetLineWidth(0);
  for (int i = 0; i < 5; i++)
    {
      l->AddEntry(h_system[i], words[i].c_str());
    }
  l->Draw("same");
}
