#include "dlUtility.h"
#include "sPhenixStyle.C"
#include <string>
#include <string.h>
#include <iostream>
#include <filesystem>
#include <vector>
#include "TMath.h"
#include "TF1.h"
#include "TPad.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TH2.h"

void DrawMBDChargeSum(const int runnumber);
void DrawMBDChannels(const int runnumber);
void DrawMBDVertex(const int runnumber);
void DrawMBDCentralityCalibrations(const int runnumber, bool use_shifted = false, bool use_balanced = false, bool flag = false);
void DrawMBDCentralityCheck(const int runnumber);
void DrawZDCCheck(const int runnumber);

void Draw_QA_Centrality_ref(const int runnumber){
  //  DrawMBDChannels(runnumber);
  DrawZDCCheck(runnumber);
  DrawMBDChargeSum(runnumber);
  DrawMBDVertex(runnumber);
  DrawMBDCentralityCalibrations(runnumber);
  //DrawMBDCentralityCalibrations(runnumber,1);
  //DrawMBDCentralityCalibrations(runnumber,1,0,1);
  DrawMBDCentralityCheck(runnumber);
  return;
}

void Draw_QA_Centrality(const int runnumber){
  //DrawMBDChannels(runnumber);
  //DrawZDCCheck(runnumber);
  DrawMBDChargeSum(runnumber);
  DrawMBDVertex(runnumber);
  //DrawMBDCentralityCalibrations(runnumber);
  DrawMBDCentralityCalibrations(runnumber,1);
  DrawMBDCentralityCalibrations(runnumber,1,0,1);
  DrawMBDCentralityCheck(runnumber);
  return;
}


void DrawMBDCentralityCalibrations(const int runnumber, bool use_shifted = false, bool use_balanced = false, bool flag = false)
{

  const char *env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }

  gStyle->SetOptStat(0);
  SetyjPadStyle();

  TString extra = Form("_%s%s%s", (use_shifted ? "sca" :""), (use_balanced ? "bal":""), (flag ? "forc":""));
  if (!use_shifted && !use_balanced && !flag) extra = "";
  TFile *file = new TFile(Form("%s/output/plots/mbdana_centrality_trigeff_%d.root", env_p, runnumber), "r");

  if (!file)
    {
      cout<< "NOFILE" <<endl;
      return;
    }

  int maxcentbins = 99;
  TH1D *hSimBBCwTrig = (TH1D*) file->Get("hSimBBCwTrig");
  if (!hSimBBCwTrig) return;
  TH1D *hSimBBC = (TH1D*) file->Get("hSimBBC");
  if (!hSimBBC) return;
  TH1D *hRealBBC;

  TH1D *h_npart_cent[99];
  TH1D *h_npart_total = (TH1D*) file->Get("h_npart_total");

  for (int i =0;i<91;i++)
    {
      h_npart_cent[i] = (TH1D*) file->Get(Form("h_npart_cent_%d", i));
      if (!h_npart_cent[i]) 
	{
	  std::cout << "no histogram for npart "<<i<<std::endl;
	}
    }
  if (use_balanced && use_shifted) hRealBBC = (TH1D*) file->Get("h_charge_sum_min_bias_w_vertex_cut_balanced_scaled");
  else if (use_balanced) hRealBBC = (TH1D*) file->Get("h_charge_sum_min_bias_w_vertex_cut_balanced");
  else if (use_shifted) hRealBBC = (TH1D*) file->Get("h_charge_sum_min_bias_w_vertex_cut_scaled");
  hRealBBC = (TH1D*) file->Get("h_charge_sum_min_bias_w_vertex_cut");
  if (!hRealBBC) return;
  TH1D *hRatio = (TH1D*) file->Get("hRatio");
  if (!hRatio) return;

  TProfile *hpRatio = new TProfile("hpRatio","",35,0, 700);
  for (int i = 0; i < hRatio->GetNbinsX(); i++)
    {
      hpRatio->Fill(hRatio->GetBinCenter(i), hRatio->GetBinContent(i));
    }

  TF1 *trigeffcurve = (TF1*) file->Get("trigeffcurve");
  if (!trigeffcurve) return;
  TString calib_file_name = Form("%s/calib/mbdana_centrality_%d.root", env_p, runnumber);

  if (flag && use_balanced && use_shifted) calib_file_name = Form("%s/calib/mbdana_sca_bal_centrality_%d.root", env_p, 20869);
  else if (flag && use_balanced) calib_file_name = Form("%s/calib/mbdana_bal_centrality_%d.root", env_p, 20869);
  else if (flag && use_shifted) calib_file_name = Form("%s/calib/mbdana_centrality_%d.root", env_p, 20869); 
  else if (use_balanced && use_shifted) calib_file_name = Form("%s/calib/mbdana_sca_bal_centrality_%d.root", env_p, runnumber);
  else if (use_balanced) calib_file_name = Form("%s/calib/mbdana_bal_centrality_%d.root", env_p, runnumber);
  else if (use_shifted) calib_file_name = Form("%s/calib/mbdana_sca_centrality_%d.root", env_p, runnumber); 
  else calib_file_name = Form("%s/calib/mbdana_centrality_%d.root", env_p, runnumber); 

  calib_file_name = Form("%s/calib/mbdana_centrality_%d.root", env_p, runnumber); 
  TFile *calibfile = new TFile(calib_file_name, "r");

  if (!calibfile)
    {

      cout<< "NOFILE" <<endl;

      return;
     
    }

  float bin, low, high, npart;
  TNtuple *tn = (TNtuple*) calibfile->Get("tn_centrality");

  if (!tn)
    {

      cout<< "NO TNTuple" <<endl;

      return;
    }
  tn->SetBranchAddress("bin",&bin);
  tn->SetBranchAddress("low",&low);
  tn->SetBranchAddress("high",&high);

  float centrality_low[100];
  float centrality_high[100];
  float npart_cent[100];
  for (int i = 0; i < tn->GetEntries(); i++)
  {
    tn->GetEntry(i);
    centrality_low[i] = low;
    centrality_high[i] = high;
    std::cout << centrality_low[i] << ", "<<centrality_high[i] <<std::endl;
  }  

  TString calib_file_name2 = Form("%s/calib/mbdana_npart_%d.root", env_p, runnumber); 
  TFile *calibfile2 = new TFile(calib_file_name2, "r");

  TNtuple *tn2 = (TNtuple*) calibfile2->Get("tn_npart");
  tn2->SetBranchAddress("npart",&npart);

  for (int i = 0; i < tn2->GetEntries(); i++)
  {
    tn2->GetEntry(i);
    npart_cent[i] = npart;
  }  

  SetsPhenixStyle();

  TCanvas *c1 = new TCanvas("canvas_bbc_realdataandsim","canvas_bbc_realdataandsim", 700, 700);
  c1->Divide(1,2);
  c1->cd(1);
  gPad->SetTicks(1);	
  hRealBBC->SetLineWidth(2);
  hRealBBC->SetMarkerStyle(24);
  hRealBBC->SetMaximum(hSimBBC->GetBinContent(10)*3);
  hRealBBC->GetXaxis()->SetRangeUser(0.0,2500);
  gPad->SetTopMargin(.13);
  
  string foo = "MBD Charge N+S (Glauber+NBD fit)";

  hRealBBC->SetTitle(Form(";%s;N_{events}", foo.c_str()));
  hRealBBC->SetTitleFont(42);
  hRealBBC->SetTitleSize(.07);
  hRealBBC->SetTitleOffset(1);
  
  hRealBBC->SetLabelSize(.05);

  hRealBBC->GetYaxis()->SetTitleFont(42);
  hRealBBC->GetYaxis()->SetTitleOffset(1);
  hRealBBC->GetYaxis()->SetTitleSize(.07);
  hRealBBC->GetYaxis()->SetLabelSize(.05);
  hRealBBC->DrawCopy("p,e,l");

  hSimBBC->SetLineWidth(2);
  hSimBBC->SetLineColor(kYellow);
  
  hSimBBC->DrawCopy("l,same");

  int nhistbins = hRealBBC->GetNbinsX();
  double err_real; double err_sim;
  double trigeff = hRealBBC->IntegralAndError(1, nhistbins, err_real)/hSimBBC->IntegralAndError(1, nhistbins, err_sim);
  double trigeff_err = trigeff*sqrt(TMath::Power(err_real/hRealBBC->Integral(1, nhistbins), 2) + TMath::Power(err_sim/hSimBBC->Integral(1, nhistbins), 2));
  //  drawText("8/30/2023", 0.8, 0.92, 0, kBlack, 0.07);
  TLatex l;
  l.SetNDC();
  l.SetTextSize(0.04);
  
  DrawSPHENIX(0.65, 0.8, 1, 0, 0.06);
  drawText("Prod 2023p015 - ana399", 0.65, 0.64, 0, kBlack, 0.06);
  drawText(Form("Run %d", runnumber), 0.75, 0.57, 0, kBlack, 0.06);

  // now do the trigger efficiency fitting
  c1->cd(2);
  gPad->SetTicks(1);


  hpRatio->SetLineWidth(2);
  hpRatio->SetLineColor(kYellow);
  hRatio->SetMarkerStyle(24);
  hRatio->SetMaximum(1.5);
  hRatio->SetMinimum(0);
  hRatio->GetXaxis()->SetRangeUser(0.0,700);
  hRatio->SetXTitle(foo.c_str());

  hRatio->SetYTitle("Ratio Real Data / Monte Carlo");
  hRatio->SetTitleFont(42);
  hRatio->SetTitleSize(.07);
  hRatio->SetTitleOffset(1);
  hRatio->SetLabelSize(.05);
  hRatio->GetYaxis()->SetTitleFont(42);
  hRatio->GetYaxis()->SetTitleSize(.07);
  hRatio->GetYaxis()->SetTitleOffset(1);
  hRatio->GetYaxis()->SetLabelSize(.05);

  hRatio->Draw("p,e");
  hpRatio->Draw("hist same");
  drawText(Form("Trigger Efficiency = %2.1f #pm %1.1f%%", trigeff*100, trigeff_err*100), 0.37, 0.8, 0, kBlack, 0.07);
  TLine *tl = new TLine(0.0,1.0,700,1.0);
  SetLineAtt(tl, kRed, 3, 4);
  tl->Draw();
  
  c1->cd(1);
  gPad->SetTicks(1);	
  gPad->SetLogy(1);

  hSimBBCwTrig->DrawCopy("same");
  float maxrange = 2500.;
  nhistbins = hRealBBC->GetNbinsX();
  TH1D *hslice = new TH1D("hslice","hslice",nhistbins,-0.5,maxrange - 0.5); 

  float in_data[100] = {};
  float in_sim[100] = {};
  for (int j = 0; j < 100; j++)
    {
      in_data[j] = 0;
      in_sim[j] = 0;
    }
  
  for (int j=0; j<91; j++) {
    hslice->Reset();
    
    for (int i = 1+(floor(centrality_high[j])); i< 1+(floor(centrality_high[j + 1]));i++)
      {
	hslice->SetBinContent(i,hSimBBCwTrig->GetBinContent(i));
	in_data[j] += hRealBBC->GetBinContent(i);
	in_sim[j] += hSimBBC->GetBinContent(i);
      }
    cout << "Centbin = " << j << " Lowcut = " << 
      hslice->GetBinCenter(1+(floor(centrality_low[j]))) << 
      " Highcut = " << 
      hslice->GetBinCenter(1+(floor(centrality_low[j+1]))) << " -- "
	 << in_data[j] <<" / " << in_sim[j] << endl;
    
    hslice->SetFillColorAlpha(2+j,0.3);
    hslice->DrawCopy("same");
  }

  // redraw real data to be on top
  hRealBBC->DrawCopy("p,e,l,same");
  hSimBBC->DrawCopy("hist,same");  
  c1->SaveAs(Form("%s/output/centplots/run%d_glauber%s.pdf", env_p, runnumber, extra.Data()));
  c1->SaveAs(Form("%s/output/centplots/run%d_glauber%s.png", env_p, runnumber, extra.Data()));
  //===============================================================================


  TCanvas *c2 = new TCanvas("c2","c2", 800, 600);
  gPad->SetLogy();
  SetyjPadStyle();
  gPad->SetRightMargin(0.25);
 
  SetLineAtt(h_npart_total, kBlack, 3, 1);
  h_npart_total->SetMarkerStyle(21);
  h_npart_total->SetMarkerSize(1);
      
  h_npart_total->SetTitle(";N_{part};counts");
  h_npart_total->Draw();
  gStyle->SetPalette(kRainBow);
  
  TLegend *ll = new TLegend(0.76,0.1, 0.99, 0.9);
  ll->SetHeader("Cent % - <N_{part}>");
  for (int i = 0; i < 91; i++)
    {

      h_npart_cent[i]->SetMarkerStyle(21);
      h_npart_cent[i]->SetMarkerSize(1);
      h_npart_cent[i]->Draw("same PMC PLC");
      ll->AddEntry(h_npart_cent[i], Form("%d-%d %% - %3.1f", i*5, i*5+5, h_npart_cent[i]->GetMean()),"PMC PLC"); 
    }
  SetLegendStyle(ll);
  ll->Draw("same");
  DrawSPHENIX(0.45, 0.83, 1, 0, 0.03);

  
  c2->SaveAs(Form("%s/output/centplots/run%d_npart_glauber%s.pdf", env_p, runnumber, extra.Data()));
  c2->SaveAs(Form("%s/output/centplots/run%d_npart_glauber%s.png", env_p, runnumber, extra.Data()));
  return;
}

void DrawMBDTimeChannels(const int runnumber)
{

  
  const char *env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }

  gStyle->SetOptStat(0);
  SetyjPadStyle();

  TFile *file = new TFile(Form("%s/output/plots/histout_%d.root", env_p, runnumber), "r");

  if (!file)
    {

      cout<< "NOFILE" <<endl;

      return;
    }

  TH1D *h_time[128];
  TH1D *h_time_raw[128];

  for (int i = 0; i < 128; i++)
    {
      h_time_raw[i] = (TH1D*) file->Get(Form("h_tdc_channel_%d", i));
      if (!h_time_raw[i]) return;
      h_time_raw[i]->SetTitle(";TDC_{raw};Counts");
    }

  // ALL 128 Cahnnels;

  TCanvas *c = new TCanvas("c","c", 1000, 1200);

  c->Print(Form("time_%d.pdf[", runnumber));
  int colorline = kViolet + 1;
  int colorfill = kViolet - 2;

  TPad *pads[64];
  float xi = 0.05;
  float yi = .9 - xi;
  float dpx = (1. - 2*xi)/8.;
  float dpy = (.9 - 2*xi)/8.;
  DrawSPHENIX(0.05, 0.95, 1, 1, 0.03, 0);
  drawText(Form("Run - %d", runnumber), 0.05, 0.91, 0, kBlack, 0.03);
  drawText("South MBD Time Distributions", 0.95, 0.91, 1, kBlack, 0.03);
  for (int i = 0; i < 64; i++)
    {
      float x = (i%8);
      float y = (i/8);
      pads[i] = new TPad(Form("pad%d", i), Form("pad%d", i),xi + x*dpx, yi - (y+1)*dpy, xi + (x+1)*dpx, yi - (y)*dpy);
      pads[i]->SetTopMargin(0);
      pads[i]->SetRightMargin(0);
      pads[i]->SetBottomMargin(0);
      pads[i]->SetLeftMargin(0);
      pads[i]->SetTicks(1,1);
      pads[i]->Draw();
    }

  for (int i = 0 ; i < 64; i++)
    {
      pads[i]->cd();
      h_time_raw[i]->GetXaxis()->SetRangeUser(-10, 10);
      SetLineAtt(h_time_raw[i], colorline, 2, 1);
      h_time_raw[i]->SetFillColorAlpha(colorfill, 0.3);
      h_time_raw[i]->Draw();
      
      drawText(Form("S Ch. %d", i),0.1, 0.88, 0, kBlack, 0.07);
    }

  c->Print(Form("%s/output/centplots/run%d_time.pdf(", env_p, runnumber));
  c->SaveAs(Form("%s/output/centplots/run%d_time_south.png", env_p, runnumber));

  c = new TCanvas("c1","c1", 1000, 1200);
  DrawSPHENIX(0.05, 0.95, 1, 1, 0.03, 0);
  drawText(Form("Run - %d", runnumber), 0.05, 0.91, 0, kBlack, 0.03);
  drawText("North MBD Time Distributions", 0.95, 0.91, 1, kBlack, 0.03);

  TPad *pads1[64];
  for (int i = 0; i < 64; i++)
    {
      float x = (i%8);
      float y = (i/8);
      pads1[i] = new TPad(Form("pad%d", i), Form("pad%d", i),xi + x*dpx, yi - (y+1)*dpy, xi + (x+1)*dpx, yi - (y)*dpy);
      pads1[i]->SetTopMargin(0);
      pads1[i]->SetRightMargin(0);
      pads1[i]->SetBottomMargin(0);
      pads1[i]->SetLeftMargin(0);
      pads1[i]->SetTicks(1,1);
      pads1[i]->Draw();
    }

  for (int i = 0 ; i < 64; i++)
    {
      pads1[i]->cd();
      h_time_raw[i+64]->GetXaxis()->SetRangeUser(-10, 10);
      SetLineAtt(h_time_raw[i+64], colorline, 1, 1);
      h_time_raw[i+64]->SetFillColorAlpha(colorfill, 0.3);
      h_time_raw[i+64]->Draw();
      drawText(Form("N Ch. %d", i),0.1, 0.88, 0, kBlack, 0.07);
    }

  c->SaveAs(Form("%s/output/centplots/run%d_time_north.png", env_p, runnumber));
  c->Print(Form("%s/output/centplots/run%d_time.pdf)", env_p, runnumber));

}
 
void DrawMBDVertex(const int runnumber)
{
  gStyle->SetOptStat(0);
  SetyjPadStyle();

  const char *env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }


  TFile *file = new TFile(Form("%s/output/plots/mbdana_zdc_check_%d.root", env_p, runnumber), "r");

  TH1D *h_vertex = (TH1D*) file->Get("h_vertex");

  if (!h_vertex) return;

  TH1D *h_time_0 = (TH1D*) file->Get("h_time_zero");

  if (!h_time_0) return;

  SetLineAtt(h_time_0, kBlack, 2, 1);
  SetLineAtt(h_vertex, kRed+3, 2, 1);

  float max = 1.3 * h_time_0->GetBinContent(h_time_0->GetMaximumBin());
  TCanvas *c_time = new TCanvas("time","time", 600, 600);
  h_time_0->SetTitle(";Time [ns];Counts");      
  h_time_0->SetMaximum(max);
  h_time_0->Draw();
  
  DrawSPHENIX(0.2, 0.87, 1, 0, 0.04);
  drawText(Form("Run %d", runnumber), 0.2, 0.76, 0, kBlack, 0.04);
  c_time->SaveAs(Form("%s/output/centplots/run%d_time_0.png", env_p, runnumber));

  max = 1.3 * h_vertex->GetBinContent(h_vertex->GetMaximumBin());
  TCanvas *c_vertex = new TCanvas("vertex","vertex", 600, 300);

  h_vertex->SetTitle(";Vertex [cm];Counts");      
  h_vertex->SetMaximum(max);
  h_vertex->Fit("gaus", "Q","", -30, 30);

  float n, mean, std, chi2ndf;
  n = h_vertex->GetFunction("gaus")->GetParameter(0);
  mean = h_vertex->GetFunction("gaus")->GetParameter(1);
  std  = h_vertex->GetFunction("gaus")->GetParameter(2);
  chi2ndf  = h_vertex->GetFunction("gaus")->GetChisquare()/h_vertex->GetFunction("gaus")->GetNDF();    
  h_vertex->GetFunction("gaus")->SetLineColor(kBlue -3);
  h_vertex->GetFunction("gaus")->SetLineWidth(3);

  h_vertex->Draw();
  h_vertex->Draw("same");
  
  
  DrawSPHENIX(0.2, 0.87, 1, 0, 0.04);
  drawText(Form("Run %d", runnumber), 0.2, 0.76, 0, kBlack, 0.04);
  drawText(Form("<z_{vtx}> = %2.2f", mean), 0.87, 0.87, 1, kBlack, 0.04);
  drawText(Form("#sigma(z_{vtx}) = %2.2f", std), 0.87, 0.82, 1, kBlack, 0.04);
  drawText(Form("#Chi^{2}/NDF = %2.2f", chi2ndf), 0.87, 0.77, 1, kBlack, 0.04);
  c_vertex->SaveAs(Form("%s/output/centplots/run%d_zvtx.png", env_p, runnumber));
  gPad->SetLogy();
  c_vertex->SaveAs(Form("%s/output/centplots/run%d_zvtx_log.png", env_p, runnumber));
  
}

void DrawMBDChargeSum(const int runnumber)
{
  gStyle->SetOptStat(0);
  SetyjPadStyle();

  const char *env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }


  TFile *file = new TFile(Form("%s/output/plots/mbdana_charge_sum_%d.root", env_p, runnumber), "r");
  if (!file) return;
  TH1D *h_charge_sum = (TH1D*) file->Get("h_charge_sum");
  if (!h_charge_sum) return;
  TH1D *h_charge_sum_vtx = (TH1D*) file->Get("h_charge_sum_vtx");
  if (!h_charge_sum_vtx) return;
  TH1D *h_charge_sum_mb_vtx = (TH1D*) file->Get("h_charge_sum_min_bias_w_vertex_cut");
  if (!h_charge_sum_mb_vtx) return;  
  TH1D *h_charge_sum_sca = (TH1D*) file->Get("h_charge_sum_scaled");
  TH1D *h_charge_sum_sca_vtx = (TH1D*) file->Get("h_charge_sum_vtx_scaled");
  TH1D *h_charge_sum_sca_mb_vtx = (TH1D*) file->Get("h_charge_sum_min_bias_w_vertex_cut_scaled");
  
  TH2D *h_charge_sum_v_zdc =  (TH2D*) file->Get("h_charge_sum_v_zdc");
  TH2D *h_charge_sum_v_zdc_mb =  (TH2D*) file->Get("h_charge_sum_v_zdc_mb");

  TProfile *hp_charge_sum_v_zdc =  (TProfile*) file->Get("hp_charge_sum_v_zdc");
  TProfile *hp_charge_sum_v_zdc_mb =  (TProfile*) file->Get("hp_charge_sum_v_zdc_mb");

  TFile *file2 = new TFile(Form("%s/output/plots/mbdana_charge_sum_20869.root", env_p), "r");
  TProfile *hp_charge_sum_v_zdc_mbref =  (TProfile*) file2->Get("hp_charge_sum_v_zdc_mb");
  hp_charge_sum_v_zdc_mbref->GetXaxis()->SetRangeUser(10, 1600);
  hp_charge_sum_v_zdc_mb->GetXaxis()->SetRangeUser(10, 1600);

  SetLineAtt(h_charge_sum, kBlack, 2, 1);
  SetLineAtt(h_charge_sum_vtx, kBlue, 2, 1);
  SetLineAtt(h_charge_sum_mb_vtx, kSpring + 2, 2, 1);
  SetLineAtt(h_charge_sum_sca, kBlack, 2, 1);
  SetLineAtt(h_charge_sum_sca_vtx, kBlue, 2, 1);
  SetLineAtt(h_charge_sum_sca_mb_vtx, kViolet + 2, 2, 1);


  TCanvas *c = new TCanvas("c1","c1", 500, 500);

  SetyjPadStyle();
  gPad->SetLogy();
  h_charge_sum->SetTitle(";MBD Charge Sum; Counts");
  h_charge_sum->SetMaximum(10*h_charge_sum->GetBinContent(h_charge_sum->GetMaximumBin()));
  h_charge_sum->Draw("");
  h_charge_sum_vtx->Draw("same");
  h_charge_sum_mb_vtx->Draw("same");
  DrawSPHENIX(0.22, 0.85, 1, 0, 0.04);
  drawText(Form("Run %d", runnumber), 0.22, 0.74, 0, kBlack, 0.04);

  TLegend *l = new TLegend(0.6, 0.7, 0.8, 0.9);
  SetLegendStyle(l);
  l->SetHeader("MBD Charge Dist.");
  l->AddEntry(h_charge_sum, "No cuts");
  l->AddEntry(h_charge_sum_vtx, "|z_{mbd} < 60|");
  l->AddEntry(h_charge_sum_mb_vtx, "ZDC and MBD sum cut");
  l->Draw("same");

  c->SaveAs(Form("%s/output/centplots/run%d_charge_sum.png", env_p, runnumber));


  c = new TCanvas("c3","c3", 500, 500);
  float scale =   h_charge_sum_sca_mb_vtx->GetMean()/h_charge_sum_mb_vtx->GetMean();
  gPad->SetLogy();
  h_charge_sum_mb_vtx->SetTitle(";MBD Charge Sum; Counts");
  h_charge_sum_mb_vtx->SetMaximum(10*h_charge_sum_mb_vtx->GetBinContent(h_charge_sum_mb_vtx->GetMaximumBin()));
  h_charge_sum_mb_vtx->Draw("");
  h_charge_sum_sca_mb_vtx->Draw("same");
  DrawSPHENIX(0.22, 0.85, 1, 0, 0.04);
  drawText(Form("Run %d", runnumber), 0.22, 0.74, 0, kBlack, 0.04);
  drawText("Min Bias Events", 0.86, 0.74, 1, kBlack, 0.04);
  drawText(Form("Scale Factor = %1.2f", scale), 0.86, 0.69, 1, kBlack, 0.04);
  
  l = new TLegend(0.6, 0.8, 0.8, 0.9);
  SetLegendStyle(l);
  l->AddEntry(h_charge_sum_mb_vtx, "not scaled");
  l->AddEntry(h_charge_sum_sca_mb_vtx, "Scaled to run 20869");
  l->Draw("same");

  c->SaveAs(Form("%s/output/centplots/run%d_charge_sum_scaled.png", env_p, runnumber));
  if (h_charge_sum_v_zdc_mb)
    {
      c = new TCanvas("c2","c2", 500, 500);
      gPad->SetLogx();
      h_charge_sum_v_zdc_mb->SetTitle(";MBD #SigmaQ; ZDC #Sigma E");
      h_charge_sum_v_zdc_mb->Draw("");
      SetLineAtt(hp_charge_sum_v_zdc_mb, kRed, 1,3);
      SetLineAtt(hp_charge_sum_v_zdc_mbref, kBlue - 2, 1,3);
      SetMarkerAtt(hp_charge_sum_v_zdc_mb, kRed, 1,75);
      SetMarkerAtt(hp_charge_sum_v_zdc_mbref, kBlue - 2, 1,75);
      hp_charge_sum_v_zdc_mb->GetXaxis()->SetRangeUser(10, 2000);
      hp_charge_sum_v_zdc_mbref->GetXaxis()->SetRangeUser(10, 2000);

      hp_charge_sum_v_zdc_mb->Draw("same");
      hp_charge_sum_v_zdc_mbref->Draw("same");
      DrawSPHENIX(0.22, 0.87, 1, 0, 0.04);
      drawText(Form("Run %d", runnumber), 0.22, 0.77, 0, kBlack, 0.04);
      l = new TLegend(0.22, 0.2, 0.5, 0.32);
      l->AddEntry(hp_charge_sum_v_zdc_mbref,"Ref. run - 20869","P");
      l->AddEntry(hp_charge_sum_v_zdc_mb,"This run","P");
      SetLegendStyle(l);
      l->Draw("same");

      c->SaveAs(Form("%s/output/centplots/run%d_charge_sum_zdc.png", env_p, runnumber));
    }
}

void DrawMBDChargeSum(const int runnumber, const int runnumber2)
{
  gStyle->SetOptStat(0);
  SetyjPadStyle();

  const char *env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }


  TFile *file = new TFile(Form("%s/output/plots/mbdana_charge_sum_%d.root", env_p, runnumber), "r");
  if (!file) return;

  TFile *file2 = new TFile(Form("%s/output/plots/mbdana_charge_sum_%d.root", env_p, runnumber2), "r");
  if (!file2) return;

  TH1D *h_charge_sum_mb_vtx = (TH1D*) file->Get("h_charge_sum_min_bias_w_vertex_cut");
  if (!h_charge_sum_mb_vtx) return;  

  TH1D *h_charge_sum_mb_vtx2 = (TH1D*) file2->Get("h_charge_sum_min_bias_w_vertex_cut");
  if (!h_charge_sum_mb_vtx2) return;  
  double sum1 = 0;
  double sum2 = 0;
  for (int i = 1470; i <2000;i++)
    {
      sum1 += h_charge_sum_mb_vtx->GetBinContent(i);
      sum2 += h_charge_sum_mb_vtx2->GetBinContent(i);
    } 
  h_charge_sum_mb_vtx2->Scale(h_charge_sum_mb_vtx->Integral()/h_charge_sum_mb_vtx2->Integral());//sum1/sum2);
  SetLineAtt(h_charge_sum_mb_vtx, kSpring + 2, 2, 1);
  SetLineAtt(h_charge_sum_mb_vtx2, kViolet + 2, 2, 1);


  TCanvas *c = new TCanvas("c1","c1", 700, 700);

  c->Divide(1,2);
  
  c->cd(1);
  gPad->SetLogy();
  h_charge_sum_mb_vtx2->SetTitle(";MBD Charge Sum; Counts");
  h_charge_sum_mb_vtx2->SetMaximum(10*h_charge_sum_mb_vtx->GetBinContent(h_charge_sum_mb_vtx->GetMaximumBin()));
  h_charge_sum_mb_vtx2->Draw("");
  h_charge_sum_mb_vtx->Draw("same");
  DrawSPHENIX(0.23, 0.8, 1, 0, 0.05);

  TLegend *l = new TLegend(0.6, 0.65, 0.8, 0.85);
  //  SetLegendStyle(l);
  l->SetHeader("MBD Charge Dist.");
  l->AddEntry(h_charge_sum_mb_vtx, Form("Run %d - MB", runnumber));
  l->AddEntry(h_charge_sum_mb_vtx2, Form("Run %d - Central", runnumber2));
  l->Draw("same");

  // p2->cd();
  
  TH1D *hRatio = (TH1D*) h_charge_sum_mb_vtx2->Clone();

  hRatio->Reset();
  hRatio->SetName("hRatio");

  for (int i = 1; i <= h_charge_sum_mb_vtx2->GetNbinsX(); i++)
    {
      hRatio->SetBinContent(i, h_charge_sum_mb_vtx2->GetBinContent(i)/h_charge_sum_mb_vtx->GetBinContent(i));
      hRatio->SetBinError(i, sqrt(TMath::Power(h_charge_sum_mb_vtx2->GetBinError(i)/h_charge_sum_mb_vtx2->GetBinContent(i), 2) + 
				  TMath::Power(h_charge_sum_mb_vtx->GetBinError(i)/h_charge_sum_mb_vtx->GetBinContent(i), 2)) * hRatio->GetBinContent(i));
      if (hRatio->GetBinContent(i) < 0.01) hRatio->SetBinError(i, 0.01);

    }
  hRatio->SetMinimum(0);
  hRatio->SetMaximum(2.0);
  hRatio->SetTitle(";MBD Charge N+S; Central/Minimum Bias Trigger");
  hRatio->GetXaxis()->SetRangeUser(0., 1000.);
  SetLineAtt(hRatio, kBlack, 1, 1);
  SetMarkerAtt(hRatio, kBlack, 1, 24);

  TProfile *hpRatio = new TProfile("hpr","", 100, 0, 1000);
  for (int i = 1; i <= 1001; i++) hpRatio->Fill(hRatio->GetBinCenter(i), hRatio->GetBinContent(i));
  c->cd(2);
  hRatio->DrawCopy("p, e");
  SetLineAtt(hpRatio, kYellow +1, 2, 1);
  TLine *line = new TLine(0., 1., 1000., 1.);
  SetLineAtt(line, kRed, 2, 4);
  line->Draw("same");
  hpRatio->DrawCopy("hist same");
  c->SaveAs(Form("%s/output/centplots/run%d_run%d_charge_sum.png", env_p, runnumber, runnumber2));



}


void DrawMBDCentralityCheck(const int runnumber)
{
  gStyle->SetOptStat(0);
  SetyjPadStyle();
  const char *env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }
  TFile *file = new TFile(Form("%s/output/plots/mbdana_centrality_check_%d.root", env_p, runnumber), "r");

  if (!file)
    {
      return;
    }

  TCanvas *c = new TCanvas("c","c");
  TH1D *h_cent_bin_shifted;
  TH1D *h_cent_bin;

  h_cent_bin = (TH1D*) file->Get("hcent_bins");

  if (!h_cent_bin) return;

  h_cent_bin_shifted = (TH1D*) file->Get("hcent_bins_sca");
  if (!h_cent_bin_shifted) return;

  bool use_shifted = true;

  if (!h_cent_bin_shifted) use_shifted = false;


  h_cent_bin->SetTitle(";Centrality Bin; Fraction of Events");
  if (use_shifted)  SetMarkerAtt(h_cent_bin_shifted, kViolet - 2, 1, 89);
  SetMarkerAtt(h_cent_bin, kSpring + 2, 1, 90);
  if (use_shifted)  SetLineAtt(h_cent_bin_shifted, kViolet + 2, 1, 1);
  SetLineAtt(h_cent_bin, kSpring-1, 1, 1);
  
  h_cent_bin->SetMaximum(0.02);
  h_cent_bin->SetMinimum(0.005);
  h_cent_bin->Draw();
  if (use_shifted) h_cent_bin_shifted->Draw("same");
  TLine *l = new TLine(-0.5, 0.05, 99.5, 0.05);
  SetLineAtt(l, kBlack, 2, 7);
  l->Draw("same");

  TF1 *flatline = new TF1("flatline","[0]",-0.5, 99.5);
  flatline->SetParameter(0,1./.91);

  h_cent_bin_shifted->Fit(flatline,"NDORQ","");

  float chi2_shifted = 0;
  float yy_shifted = 0;
  float chi2 = flatline->GetChisquare()/flatline->GetNDF();
  float yy = flatline->GetParameter(0);
  
  if (use_shifted)
    {
      //      h_cent_bin_shifted->Scale(1./h_cent_bin_shifted->Integral());
      h_cent_bin_shifted->Fit(flatline,"QNDOR","");

      chi2_shifted = flatline->GetChisquare()/flatline->GetNDF();;
      yy_shifted = flatline->GetParameter(0);
    }
  
  DrawSPHENIX(0.19, 0.85, 1, 0);
  drawText(Form("Run %d", runnumber), 0.19, 0.73, 0, kBlack);
  flatline->SetLineColor(kRed);
  flatline->SetLineWidth(2);
  flatline->Draw("same");
  TLegend *lg = new TLegend(0.6, 0.65, 0.88, 0.89);
  SetLegendStyle(lg);
  lg->SetHeader("#chi^2 ; <Frac. of Events>");
  lg->AddEntry(h_cent_bin,Form("Not Scaled: %2.2f ; %1.4f", chi2, yy));
  if (use_shifted)   lg->AddEntry(h_cent_bin_shifted,Form("Scaled: %2.2f ; %1.4f", chi2_shifted, yy_shifted)); 
  lg->Draw("same");
  c->SaveAs(Form("%s/output/centplots/run%d_centrality.png", env_p, runnumber));
  
}

void DrawZDCCheck(const int runnumber)
{
  gStyle->SetOptStat(0);
  SetyjPadStyle();

  const char *env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }


  TFile *file = new TFile(Form("%s/output/plots/mbdana_zdc_check_%d.root", env_p, runnumber), "r");

  if (!file)
    {
      return;
    }
  
  TH1D *h_zdc_sum_n;

  TH1D *h_zdc_sum_s;

  TH1D *h_zdc_energy[6];

  for (int i = 0; i < 6; i++)
    {
      h_zdc_energy[i] = (TH1D*) file->Get(Form("h_zdc_energy_%d",i));
      if (!h_zdc_energy[i]) return;
    }
  SetLineAtt(h_zdc_energy[0], kBlack, 3, 1);
  SetLineAtt(h_zdc_energy[1], kRed, 3, 1);
  SetLineAtt(h_zdc_energy[2], kBlue, 3, 1);
  SetLineAtt(h_zdc_energy[3], kBlack, 3, 4);
  SetLineAtt(h_zdc_energy[4], kRed, 3, 4);
  SetLineAtt(h_zdc_energy[5], kBlue, 3, 4);

  h_zdc_sum_s = (TH1D*) file->Get("h_zdc_sum_s");
  if (!h_zdc_sum_s) return;
  h_zdc_sum_n = (TH1D*) file->Get("h_zdc_sum_n");
  if (!h_zdc_sum_n) return;
  h_zdc_sum_s->Scale(1./h_zdc_sum_s->Integral(5, 1000));
  h_zdc_sum_n->Scale(1./h_zdc_sum_n->Integral(5, 1000));
 
  TH1D *h_clone= (TH1D*) h_zdc_sum_s->Clone();
  h_clone->Reset();
  h_clone->Fill(100);
  float max = h_zdc_sum_s->GetBinContent(h_clone->GetMaximumBin());


  SetLineAtt(h_zdc_sum_n, kRed, 2, 1);
  SetLineAtt(h_zdc_sum_s, kBlue, 2, 1);
  h_zdc_sum_s->GetXaxis()->SetRangeUser(40, 10000);
  h_zdc_sum_s->SetMinimum(0.05 * max);
  h_zdc_sum_s->SetMaximum(100 * max);

  TCanvas *c = new TCanvas("c","c");
  c->SetLogx(1);
  c->SetLogy(1);

  h_zdc_sum_s->Draw("hist");
  h_zdc_sum_n->Draw("same hist");

  TLine *l = new TLine(100, 0.08*max, 100, max);
  SetLineAtt(l, kBlack, 3, 7);

  l->Draw("same");
  DrawSPHENIX(0.19, 0.85, 1, 0);
  drawText(Form("Run %d", runnumber), 0.19, 0.73, 0, kBlack);
  drawText("|z_{vtx}| < 30 cm", 0.19, 0.67, 0, kBlack);

  TLegend *lg = new TLegend(0.69, 0.7, 0.87, 0.87);
  lg->AddEntry(h_zdc_sum_n, "North","L");
  lg->AddEntry(h_zdc_sum_s, "South","L");
  lg->AddEntry(l, "Single Neutron","L");
  SetLegendStyle(lg);
  lg->Draw("same");

  c->SaveAs(Form("%s/output/centplots/run%d_zdc.png", env_p, runnumber));
  c->SaveAs(Form("%s/output/centplots/run%d_zdc.pdf", env_p, runnumber));

  c = new TCanvas("c","c");
  c->SetLogx(1);
  c->SetLogy(1);

  h_zdc_energy[2]->Draw("hist");
  for (int i = 0; i < 6;i++)  h_zdc_energy[i]->Draw("same hist");

  DrawSPHENIX(0.19, 0.85, 1, 0);
  drawText(Form("Run %d", runnumber), 0.19, 0.73, 0, kBlack);
  drawText("|z_{vtx}| < 30 cm", 0.19, 0.67, 0, kBlack);

  TString names[6] = {"South 1", "South 2", "South 3", "North 1", "North 2", "North 3"};
  lg = new TLegend(0.69, 0.7, 0.87, 0.87);
  for (int i = 0; i < 6 ; i++)lg->AddEntry(h_zdc_energy[i], names[i],"L");
  SetLegendStyle(lg);
  lg->Draw("same");

  c->SaveAs(Form("%s/output/centplots/run%d_zdc_singles.png", env_p, runnumber));
  c->SaveAs(Form("%s/output/centplots/run%d_zdc_singles.pdf", env_p, runnumber));


}

void DrawMBDChannels(const int runnumber)
{
  
  const char *env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }

  gStyle->SetOptStat(0);
  TFile *fout = new TFile(Form("%s/output/plots/mbdana_channels_%d.root", env_p, runnumber), "r");
  TH1D *hpeaks_fit = new TH1D("hpeaks_fit",";MPV_{calib};Tubes", 41, 0.795, 1.205);

  TH1D *h_charge[128];
  TH1D *h_time[128];
  TH1D *h_charge_mb[128];
  TH1D *h_time_mb[128];

  for (int i = 0; i < 128; i++)
    {
      h_charge_mb[i] = (TH1D*) fout->Get(Form("h_mbd_charge_ch%d", i));
      h_time[i] = (TH1D*) fout->Get(Form("h_mbd_time_ch%d", i));
      h_charge[i] = (TH1D*) fout->Get(Form("h_mbd_charge_min_bias_ch%d", i));
      h_time_mb[i] = (TH1D*) fout->Get(Form("h_mbd_time_min_bias_ch%d", i));
      if (!h_charge_mb[i]) return;
      if (!h_time_mb[i]) return;
      if (!h_charge[i]) return;
      if (!h_time[i]) return;

    }

  TCanvas *c = new TCanvas("c","c", 1000, 1200);

  c->Print(Form("fits_%d.pdf[", runnumber));
  int colorline_mb = kSpring + 4;
  int colorfill_mb = kSpring + 2;
  int colorline = kBlack;

  TPad *pads[64];
  float xi = 0.05;
  float yi = .9 - xi;
  float dpx = (1. - 2*xi)/8.;
  float dpy = (.9 - 2*xi)/8.;

  DrawSPHENIX(0.05, 0.95, 1, 1, 0.03, 0);
  drawText(Form("Run - %d", runnumber), 0.05, 0.91, 0, kBlack, 0.03);
  drawText("Dotted Line at 1", 0.05, 0.87, 0, kBlack, 0.03);
  drawText("South MBD Charge Distributions", 0.95, 0.91, 1, kBlack, 0.03);
  for (int i = 0; i < 64; i++)
    {
      float x = (i%8);
      float y = (i/8);
      pads[i] = new TPad(Form("pad%d", i), Form("pad%d", i),xi + x*dpx, yi - (y+1)*dpy, xi + (x+1)*dpx, yi - (y)*dpy);
      pads[i]->SetTopMargin(0);
      pads[i]->SetRightMargin(0);
      pads[i]->SetBottomMargin(0);
      pads[i]->SetLeftMargin(0);
      pads[i]->SetTicks(1,1);
      pads[i]->Draw();
    }

  TH1D *hh = (TH1D*) h_charge[0]->Clone();
  hh->Reset();
  hh->Fill(1);
  int binpeak = hh->GetMaximumBin();
  for (int i = 0 ; i < 64; i++)
    {
      pads[i]->cd();

      TF1 *f_lan_w_gausexp = new TF1("lan_w_gausexp","[0]*TMath::Landau(x,[1],[2],3) + [3]*TMath::Exp(-1*(x - [4])/[5])+gaus(6)", 0.5, 3);
      f_lan_w_gausexp->SetParameters(500, 1, .1, 2000, .07, 4, 80000, .05, .04);
      f_lan_w_gausexp->SetParLimits(8, 0, .2);
      f_lan_w_gausexp->SetParLimits(7, -.1, .2);
      f_lan_w_gausexp->SetParLimits(6, 0, 1000000);
      h_charge[i]->Fit("lan_w_gausexp","Q","", .2, 3);
      if (! h_charge[i]->GetFunction("lan_w_gausexp")) continue;
      h_charge[i]->GetFunction("lan_w_gausexp")->SetLineColor(kBlack);
      float peakfit =h_charge[i]->GetBinContent(binpeak);
      h_charge[i]->SetMaximum(peakfit * 1.3);

      h_charge[i]->GetXaxis()->SetRangeUser(.1, 3);
      SetLineAtt(h_charge[i], colorline, 2, 1);
      h_charge[i]->SetFillColorAlpha(colorfill_mb, 0.3);
      
      h_charge[i]->Draw();
      TLine *tl = new TLine(1, 0,1, peakfit*1.3);
      tl->SetLineColor(kBlue + 2);
      tl->SetLineStyle(2);
      tl->Draw();
      drawText(Form("S Ch. %d", i),0.5, 0.88, 0, kBlack, 0.07);
      drawText(Form("#Chi^2 = %4.4f", f_lan_w_gausexp->GetChisquare()/f_lan_w_gausexp->GetNDF()), 0.5, 0.8, 0, kBlack, 0.07);
      drawText(Form("MPV = %4.2f", f_lan_w_gausexp->GetParameter(1)), 0.5, 0.72, 0, kBlack, 0.07);
      hpeaks_fit->Fill(f_lan_w_gausexp->GetParameter(1));
    }

  c->Print(Form("%s/output/centplots/run%d_charge.pdf(", env_p, runnumber));
  c->SaveAs(Form("%s/output/centplots/run%d_charge_south.png", env_p, runnumber));

  c = new TCanvas("c1","c1", 1000, 1200);
  DrawSPHENIX(0.05, 0.95, 1, 1, 0.03, 0);
  drawText(Form("Run - %d", runnumber), 0.05, 0.91, 0, kBlack, 0.03);
  drawText("Dotted Line at 1", 0.05, 0.87, 0, kBlack, 0.03);
  drawText("North MBD Charge Distributions", 0.95, 0.91, 1, kBlack, 0.03);
  TPad *pads1[64];
  for (int i = 0; i < 64; i++)
    {
      float x = (i%8);
      float y = (i/8);
      pads1[i] = new TPad(Form("pad%d", i), Form("pad%d", i),xi + x*dpx, yi - (y+1)*dpy, xi + (x+1)*dpx, yi - (y)*dpy);
      pads1[i]->SetTopMargin(0);
      pads1[i]->SetRightMargin(0);
      pads1[i]->SetBottomMargin(0);
      pads1[i]->SetLeftMargin(0);
      pads1[i]->SetTicks(1,1);
      pads1[i]->Draw();
    }
  for (int i = 0 ; i < 64; i++)
    {
      pads1[i]->cd();

      TF1 *f_lan_w_gausexp = new TF1("lan_w_gausexp","[0]*TMath::Landau(x,[1],[2],3) + [3]*TMath::Exp(-1*(x - [4])/[5])+gaus(6)", 0.5, 3);
      f_lan_w_gausexp->SetParameters(500, 1, .1, 2000, .07, 4, 80000, .05, .04);
      f_lan_w_gausexp->SetParLimits(8, 0, .2);
      f_lan_w_gausexp->SetParLimits(7, -.1, .2);
      f_lan_w_gausexp->SetParLimits(6, 0, 1000000);
      h_charge[i+64]->Fit("lan_w_gausexp","Q","", .2, 3);
      if (!h_charge[i+64]->GetFunction("lan_w_gausexp")) continue;
      h_charge[i+64]->GetFunction("lan_w_gausexp")->SetLineColor(kBlack);
      float peakfit =h_charge[i+64]->GetBinContent(binpeak);
      h_charge[i+64]->SetMaximum(peakfit * 1.3);

      h_charge[i+64]->GetXaxis()->SetRangeUser(.1, 3);
      SetLineAtt(h_charge[i+64], colorline, 2, 1);
      h_charge[i+64]->SetFillColorAlpha(colorfill_mb, 0.3);
      
      h_charge[i+64]->Draw();

      TLine *tl = new TLine(1, 0,1, peakfit*1.3);
      tl->SetLineColor(kBlue + 2);
      tl->SetLineStyle(2);
      tl->Draw();
      drawText(Form("N Ch. %d", i),0.5, 0.88, 0, kBlack, 0.07);
      drawText(Form("#Chi^2 = %4.4f", f_lan_w_gausexp->GetChisquare()/f_lan_w_gausexp->GetNDF()), 0.5, 0.8, 0, kBlack, 0.07);
      drawText(Form("MPV = %4.2f", f_lan_w_gausexp->GetParameter(1)), 0.5, 0.72, 0, kBlack, 0.07);
      hpeaks_fit->Fill(f_lan_w_gausexp->GetParameter(1));
    }

  c->Print(Form("%s/output/centplots/run%d_charge.pdf)", env_p, runnumber));
  c->SaveAs(Form("%s/output/centplots/run%d_charge_north.png", env_p, runnumber));

  TCanvas *cc = new TCanvas("cc","cc", 500, 500);
  SetyjPadStyle();
  hpeaks_fit->Draw();
  cc->SetTicks(1,1);
  DrawSPHENIX(0.5, 0.8, 1, 0, 0.04);
  drawText(Form("Run %d", runnumber), 0.5, 0.7, 0, kBlack, 0.04);
  drawText("Official MBD Calibrations", 0.5, 0.65, 0, kBlack, 0.04);
  cc->SaveAs(Form("%s/output/centplots/run%d_charge_fits.png", env_p, runnumber));

  c = new TCanvas("c2","c2", 1000, 1200);
  colorline_mb = kViolet+1;
  colorfill_mb = kViolet - 1;

  DrawSPHENIX(0.05, 0.95, 1, 1, 0.03, 0);
  drawText(Form("Run - %d", runnumber), 0.05, 0.91, 0, kBlack, 0.03);
  drawText("Dotted Line at 0", 0.05, 0.87, 0, kBlack, 0.03);
  drawText("South MBD Time Distributions", 0.95, 0.91, 1, kBlack, 0.03);
  TPad *pads2[64];
  for (int i = 0; i < 64; i++) {
      float x = (i%8);
      float y = (i/8);
      pads2[i] = new TPad(Form("pad%d", i), Form("pad%d", i),xi + x*dpx, yi - (y+1)*dpy, xi + (x+1)*dpx, yi - (y)*dpy);
      pads2[i]->SetTopMargin(0);
      pads2[i]->SetRightMargin(0);
      pads2[i]->SetBottomMargin(0);
      pads2[i]->SetLeftMargin(0);
      pads2[i]->SetTicks(1,1);
      pads2[i]->Draw();
    }

  for (int i = 0 ; i < 64; i++)
    {
      pads2[i]->cd();
      h_time_mb[i]->GetXaxis()->SetRangeUser(-12, 12);
      h_time[i]->Scale(h_time_mb[i]->Integral()/h_time[i]->Integral());
      SetLineAtt(h_time[i], colorline, 2, 1);
      SetLineAtt(h_time_mb[i], colorline_mb, 2, 1);
      h_time_mb[i]->SetFillColorAlpha(colorfill_mb, 0.3);
      
      h_time_mb[i]->Draw();
      h_time[i]->Draw("same");
      TLine *tl = new TLine(0.4, 0,0.4, h_time_mb[i]->GetBinContent(h_time_mb[i]->GetMaximumBin()));
      tl->SetLineColor(kBlue + 2);
      tl->SetLineStyle(2);
      tl->Draw();
      drawText(Form("S Ch. %d", i),0.5, 0.88, 0, kBlack, 0.07);
    }

  c->Print(Form("%s/output/centplots/run%d_time.pdf(", env_p, runnumber));
  c->SaveAs(Form("%s/output/centplots/run%d_time_south.png", env_p, runnumber));

  c = new TCanvas("c3","c3", 1000, 1200);
  DrawSPHENIX(0.05, 0.95, 1, 1, 0.03, 0);
  drawText(Form("Run - %d", runnumber), 0.05, 0.91, 0, kBlack, 0.03);
  drawText("Dotted Line at 0", 0.05, 0.87, 0, kBlack, 0.03);
  drawText("North MBD Time Distributions", 0.95, 0.91, 1, kBlack, 0.03);
  TPad *pads3[64];
  for (int i = 0; i < 64; i++)
    {
      float x = (i%8);
      float y = (i/8);
      pads3[i] = new TPad(Form("pad%d", i), Form("pad%d", i),xi + x*dpx, yi - (y+1)*dpy, xi + (x+1)*dpx, yi - (y)*dpy);
      pads3[i]->SetTopMargin(0);
      pads3[i]->SetRightMargin(0);
      pads3[i]->SetBottomMargin(0);
      pads3[i]->SetLeftMargin(0);
      pads3[i]->SetTicks(1,1);
      pads3[i]->Draw();
    }
  for (int i = 0 ; i < 64; i++)
    {
      pads3[i]->cd();
      h_time_mb[i+64]->GetXaxis()->SetRangeUser(-12, 12);
      h_time[i+64]->Scale(h_time_mb[i+64]->Integral()/h_time[i+64]->Integral());
      SetLineAtt(h_time[i+64], colorline, 2, 1);
      SetLineAtt(h_time_mb[i+64], colorline_mb, 2, 1);
      h_time_mb[i+64]->SetFillColorAlpha(colorfill_mb, 0.3);
      
      h_time_mb[i+64]->Draw();
      h_time[i+64]->Draw("same");
      TLine *tl = new TLine(0.4, 0,0.4, h_time_mb[i+64]->GetBinContent(h_time_mb[i+64]->GetMaximumBin()));
      tl->SetLineColor(kBlue + 2);
      tl->SetLineStyle(2);
      tl->Draw();
      drawText(Form("N Ch. %d", i),0.5, 0.88, 0, kBlack, 0.07);
    }
  c->Print(Form("%s/output/centplots/run%d_time.pdf)", env_p, runnumber));
  c->SaveAs(Form("%s/output/centplots/run%d_time_north.png", env_p, runnumber));

}
