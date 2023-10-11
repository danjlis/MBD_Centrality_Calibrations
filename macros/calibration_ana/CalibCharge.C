#include "dlUtility.h"
void CalibCharge(int runnumber)
{

  TFile* f = new TFile(Form("../rootfiles/plots_%d.root", runnumber),"r");


  TF1 *f1 = new TF1("f1","gaus(0) + expo(3)",20, 700);
  f1->SetParameter(0, 100);
  f1->SetParameter(1, 150);
  f1->SetParameter(2, 50);

  f1->SetParLimits(0, 1, 10000);
  f1->SetParLimits(1, 30, 500);
  f1->SetParLimits(2, 10, 100);

  TF1 *f2 = new TF1("f2","[0]*TMath::Landau(x,[1],[2]) + expo(3)",20, 1000);
  
  f2->SetParLimits(0, 100, 50000);
  f2->SetParLimits(1, 30, 500);
  f2->SetParLimits(2, 10, 100);

  float ch[128];
  float peak_gaus[128];
  float peak_landau[128];
  float peak_gaus_err[128];
  float peak_landau_err[128];

  float cherr[128];
  TH1D *h[128];
  TH1D *h1 = new TH1D("h1","",32,0, 500);
  TH1D *h2 = new TH1D("h2","",32,0, 500);  

  TFile *fout = new TFile(Form("calib_mbd_%d.root", runnumber),"RECREATE");
  TNtuple *tn = new TNtuple("mbd_calib","mbd_calib","channel:peak:width");
  
  for (int i = 0; i < 128; i++)
    {

      h[i] = (TH1D*)f->Get(Form("h_mbd_charge_raw_ch%d", i));    

      float local_min = 0;
      float local_max = 0;
      int good = 0;
      int low_bin = 1;
      int high_bin = 1;

      for (int ic = 2 ; ic < h[i]->GetNbinsX()-2 ; ic++)
	{
	  good = 0;
	  for (int j = -2; j < 3; j++)
	    {
	      if (j == 0) continue;
	      if (h[i]->GetBinContent(ic+1) <= h[i]->GetBinContent(ic+1 + j))
		good++;

	    }
	  if (good == 4)
	    {
	      local_min = h[i]->GetBinLowEdge(ic+1);
	      low_bin = ic;
	      break;
	      
	    }
	    
	}

      for (int ic = low_bin ; ic < h[i]->GetNbinsX()-2 ; ic++)
	{
	  if (h[i]->GetBinContent(ic + 1) > local_max)
	    {
	      local_max = h[i]->GetBinContent(ic+1);
	      high_bin = ic+1;
	    }
	}

      local_max = h[i]->GetBinLowEdge(high_bin);

      f1->SetParameter(0, 100);
      f1->SetParameter(1, local_max);
      f1->SetParameter(2, 20);
      f2->SetParameter(0, 10000);
      f2->SetParameter(1, local_max);
      f2->SetParameter(2, 20);

      float upper = local_min  + 3* (h[i]->GetBinCenter(high_bin) - local_min);
      h[i]->GetXaxis()->SetRangeUser(0, 1500);
      
      h[i]->Fit("f1","Q+", "",local_min, upper);
      h[i]->Fit("f2","Q+", "",local_min, upper);
      h[i]->GetFunction("f2")->SetLineColor(kBlue);
      ch[i] = static_cast<float>(i);
      peak_gaus[i] = f1->GetParameter(1);
      peak_landau[i] = f2->GetParameter(1);
      cherr[i] = 0.5;
      peak_gaus_err[i] = f1->GetParError(1);
      peak_landau_err[i] = f2->GetParError(1);
      std::cout << i<<" : "<<f2->GetParameter(1)<<std::endl;
      std::cout << i<<" : "<<f1->GetParameter(1)<<std::endl;
      h1->Fill(f1->GetParameter(1));
      h2->Fill(f2->GetParameter(1));
      tn->Fill(i, 1./f2->GetParameter(1), f2->GetParameter(2));
    }

  TGraphErrors *g_gaus = new TGraphErrors(128, ch, peak_gaus,cherr, peak_gaus_err);
  TGraphErrors *g_landau = new TGraphErrors(128, ch, peak_landau,cherr, peak_landau_err);

  SetLineAtt(g_gaus, kRed+2, 1, 2);
  SetMarkerAtt(g_gaus, kRed+2, 3, 1);
  SetLineAtt(g_landau, kBlue+2, 1, 2);
  SetMarkerAtt(g_landau, kBlue+2, 3, 1);

  TCanvas *c1 = new TCanvas("c1","c1");
  g_gaus->Draw("AP");
  g_landau->Draw("P");

  fout->cd();
  for (int i = 0; i < 128; i++)
    {
      h[i]->Write();
    }
  h1->Write();
  h2->Write();
  fout->Write();
  fout->Close();

  
}
