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

void Draw_QA_Centrality(const int runnumber){
  DrawMBDChannels(runnumber);
  DrawZDCCheck(runnumber);
  DrawMBDChargeSum(runnumber);
  DrawMBDVertex(runnumber);
  DrawMBDCentralityCalibrations(runnumber);
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

  TString extra = Form("%s%s%s_", (use_shifted ? "sca" :""), (use_balanced ? "bal":""), (flag ? "forc":""));
  if (!use_shifted && !use_balanced && !flag) extra = "";
  TFile *file = new TFile(Form("%s/output/plots/mbdana_centrality_trigeff_%s%d.root", env_p, extra.Data(), runnumber), "r");

  if (!file)
    {
      cout<< "NOFILE" <<endl;
      return;
    }

  int maxcentbins = 19;
  TH1D *hSimBBCwTrig = (TH1D*) file->Get("hSimBBCwTrig");
  TH1D *hSimBBC = (TH1D*) file->Get("hSimBBC");
  TH1D *hRealBBC;

  if (use_balanced && use_shifted) hRealBBC = (TH1D*) file->Get("h_charge_sum_min_bias_w_vertex_cut_balanced_scaled");
  else if (use_balanced) hRealBBC = (TH1D*) file->Get("h_charge_sum_min_bias_w_vertex_cut_balanced");
  else if (use_shifted) hRealBBC = (TH1D*) file->Get("h_charge_sum_min_bias_w_vertex_cut_scaled");
  else  hRealBBC = (TH1D*) file->Get("h_charge_sum_min_bias_w_vertex_cut");
  TH1D *hRatio = (TH1D*) file->Get("hRatio");
  TF1 *trigeffcurve = (TF1*) file->Get("trigeffcurve");

  TString calib_file_name = Form("%s/calib/mbdana_centrality_%d.root", env_p, runnumber);

  if (flag && use_balanced && use_shifted) calib_file_name = Form("%s/calib/mbdana_sca_bal_centrality_%d.root", env_p, 23696);
  else if (flag && use_balanced) calib_file_name = Form("%s/calib/mbdana_bal_centrality_%d.root", env_p, 23696);
  else if (flag && use_shifted) calib_file_name = Form("%s/calib/mbdana_sca_centrality_%d.root", env_p, 23696); 
  else if (use_balanced && use_shifted) calib_file_name = Form("%s/calib/mbdana_sca_bal_centrality_%d.root", env_p, runnumber);
  else if (use_balanced) calib_file_name = Form("%s/calib/mbdana_bal_centrality_%d.root", env_p, runnumber);
  else if (use_shifted) calib_file_name = Form("%s/calib/mbdana_sca_centrality_%d.root", env_p, runnumber); 
  TFile *calibfile = new TFile(calib_file_name, "r");

  if (!calibfile)
    {

      cout<< "NOFILE" <<endl;

      return;
     
    }

  float bin, low, high;
  TNtuple *tn = (TNtuple*) calibfile->Get("tn_centrality");
  if (!tn)
    {

      cout<< "NO TNTuple" <<endl;

      return;
    }
  tn->SetBranchAddress("bin",&bin);
  tn->SetBranchAddress("low",&low);
  tn->SetBranchAddress("high",&high);

  float centrality_low[19];
  float centrality_high[19];
  for (int i = 0; i < 19; i++)
  {
    tn->GetEntry(i);
    centrality_low[i] = low;
    centrality_high[i] = high;
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

  hSimBBC->SetLineWidth(1);
  hSimBBC->SetLineColor(1);
  
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
  drawText("|z_{MBD}| < 30 + Sum Cut", 0.65, 0.64, 0, kBlack, 0.06);

  // now do the trigger efficiency fitting
  c1->cd(2);
  gPad->SetTicks(1);

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
  drawText(Form("Trigger Efficiency = %2.1f #pm %1.1f%%", trigeff*100, trigeff_err*100), 0.37, 0.8, 0, kBlack, 0.07);
  TLine *tl = new TLine(0.0,1.0,700,1.0);
  SetLineAtt(tl, kBlue, 2, 4);
  tl->Draw();
  
  c1->cd(1);
  gPad->SetTicks(1);	
  gPad->SetLogy(1);

  hSimBBCwTrig->DrawCopy("same");
  float maxrange = 2500.;
  nhistbins = hRealBBC->GetNbinsX();
  TH1D *hslice = new TH1D("hslice","hslice",nhistbins,-0.5,maxrange - 0.5); 
  for (int j=0; j<maxcentbins; j++) {
    hslice->Reset();
    for (int i = 1+(int)centrality_low[j]; i< 1+(int)centrality_high[j];i++)
      hslice->SetBinContent(i,hSimBBCwTrig->GetBinContent(i));

      cout << "Centbin = " << j << " Lowcut = " << 
	hslice->GetBinCenter(1+(int)centrality_low[j]) << 
        " Highcut = " << 
        hslice->GetBinCenter(1+(int)centrality_high[j]) << endl;

    hslice->SetFillColorAlpha(2+j,0.3);
    hslice->DrawCopy("same");
  }

  // redraw real data to be on top
  hRealBBC->DrawCopy("p,e,l,same");
  
  c1->SaveAs(Form("%s/output/centplots/run%d_%scharge_sum.pdf", env_p, runnumber, extra.Data()));
  c1->SaveAs(Form("%s/output/centplots/run%d_%scharge_sum.png", env_p, runnumber, extra.Data()));
  //===============================================================================

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


  TFile *file = new TFile(Form("%s/output/plots/mbd_charge_sum_%d.root", env_p, runnumber), "r");

  TH1D *h_vertex = (TH1D*) file->Get("h_vertex");
  TH1D *h_vertex_c[20];
  for (int i = 0 ; i < 20; i++)
    {
      h_vertex_c[i] = (TH1D*) file->Get(Form("h_vertex_c%d", i));
      h_vertex_c[i]->Rebin(5);
      h_vertex_c[i]->GetXaxis()->SetRangeUser(-50, 50);
      h_vertex_c[i]->Fit("gaus","Q","", -30, 30);
      h_vertex_c[i]->Scale(1./h_vertex_c[i]->Integral());
    }
  TH1D *h_vertex_w_t = (TH1D*) file->Get("h_vertex_w_t");
  TH1D *h_time_0 = (TH1D*) file->Get("h_time_0");
  TEfficiency *h_voccuppancy = (TEfficiency*) file->Get("h_voccuppancy");
  TEfficiency *h_occuppancy = (TEfficiency*) file->Get("h_occuppancy");
  TEfficiency *h_ring_voccuppancy[3];
  TEfficiency *h_ring_occuppancy[3];

  for (int i = 0 ; i < 3; i++)
    {

      h_ring_voccuppancy[i] = (TEfficiency*) file->Get(Form("h_voccuppancy_ring_%d",i));
      h_ring_occuppancy[i] = (TEfficiency*) file->Get(Form("h_occuppancy_ring_%d",i));

    }
  
  // ALL 128 Cahnnels;

  

  SetLineAtt(h_time_0, kBlack, 2, 1);
  SetLineAtt(h_vertex, kRed+3, 2, 1);
  SetLineAtt(h_vertex_w_t, kBlue - 1, 2, 1);


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
  
  TCanvas *c3 = new TCanvas("c3","c3");

  c3->Divide(3,1);
  for (int i = 0; i < 3; i++)
    {
      c3->cd(i+1);
      h_ring_voccuppancy[i]->SetMarkerColor(kRed);
      h_ring_voccuppancy[i]->SetMarkerStyle(8);
      h_ring_voccuppancy[i]->SetMarkerSize(2);
      h_ring_occuppancy[i]->SetMarkerColor(kBlue);
      h_ring_occuppancy[i]->SetMarkerStyle(8);
      h_ring_occuppancy[i]->SetMarkerSize(2);

      h_ring_voccuppancy[i]->Draw();
      gPad->Update();
      auto graph = h_ring_voccuppancy[i]->GetPaintedGraph();
      graph->SetMinimum(0);
      graph->SetMaximum(1.2);
      gPad->Update();
      h_ring_occuppancy[i]->Draw("same");
    }

  TCanvas *c_vtx_c = new TCanvas("c_vtx_c","c_vtx_c", 950, 700);
  SetyjPadStyle();
  gPad->SetLogy();
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.35);
  gStyle->SetPalette(kCMYK);
  h_vertex_c[1]->SetMaximum(1);
  TLegend *l_v = new TLegend(0.66,0.3, 0.99, 0.9);
  SetLegendStyle(l_v);
  l_v->SetHeader("Cent % - #sigma , #Chi^{2}/NDF");
  for (int i = 0; i <= 16; i++)
    {
      h_vertex_c[i+1]->SetLineWidth(2);
      if (i == 0){
	h_vertex_c[i+1]->SetLineColor(kSpring - 2);
	h_vertex_c[i+1]->SetTitle(";Z_{vtx} [cm]; Arb. Units [Normalized by Integral]"); 
	h_vertex_c[i+1]->Draw("hist");
      }
      else if (i != 16) h_vertex_c[1 + i]->Draw("same hist PLC PMC");
      else
	{
	  h_vertex_c[1 + i]->SetLineColor(kRed - 3 );
	  h_vertex_c[1 + i]->Draw("same hist");
	}

      l_v->AddEntry(h_vertex_c[1 + i], Form("%d - %d %% - %2.2f cm , %2.2f", (i)*5, (i == 16 ? 100 :(i+1)*5), h_vertex_c[1+i]->GetFunction("gaus")->GetParameter(2), h_vertex_c[1+i]->GetFunction("gaus")->GetChisquare()/h_vertex_c[1+i]->GetFunction("gaus")->GetNDF()));
    }
  h_vertex_c[1]->Draw("hist same");
  DrawSPHENIX(0.15, 0.87, 1, 0, 0.04);
  drawText(Form("Run %d", runnumber), 0.15, 0.76, 0, kBlack, 0.04);
  l_v->Draw("same");
  c_vtx_c->SaveAs(Form("%s/output/centplots/run%d_zvtx_cent.png", env_p, runnumber));
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

  TH1D *h_charge_sum = (TH1D*) file->Get("h_charge_sum");
  TH1D *h_charge_sum_vtx = (TH1D*) file->Get("h_charge_sum_vtx");
  TH1D *h_charge_sum_mb_vtx = (TH1D*) file->Get("h_charge_sum_min_bias_w_vertex_cut");

  TH1D *h_charge_sum_sca = (TH1D*) file->Get("h_charge_sum_scaled");
  TH1D *h_charge_sum_sca_vtx = (TH1D*) file->Get("h_charge_sum_vtx_scaled");
  TH1D *h_charge_sum_sca_mb_vtx = (TH1D*) file->Get("h_charge_sum_min_bias_w_vertex_cut_scaled");
  

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

  c->SaveAs(Form("%s/output/centplots/run%d_chargesum.png", env_p, runnumber));


  c = new TCanvas("c2","c2", 500, 500);

  gPad->SetLogy();
  h_charge_sum->SetTitle(";MBD Charge Sum; Counts");
  h_charge_sum_mb_vtx->SetMaximum(10*h_charge_sum_mb_vtx->GetBinContent(h_charge_sum_mb_vtx->GetMaximumBin()));
  h_charge_sum_mb_vtx->Draw("");
  h_charge_sum_sca_mb_vtx->Draw("same");
  DrawSPHENIX(0.22, 0.85, 1, 0, 0.04);
  drawText(Form("Run %d", runnumber), 0.22, 0.74, 0, kBlack, 0.04);
  drawText("ZDC_{N,S} >= 40 GeV", 0.88, 0.77, 1, kBlack, 0.02);
  drawText("!(Q_{MBD,N} < 10 && Q_{MBD,S} > 150)", 0.88, 0.73, 1, kBlack, 0.02);
  l = new TLegend(0.6, 0.8, 0.8, 0.9);
  SetLegendStyle(l);
  l->AddEntry(h_charge_sum_mb_vtx, "not scaled");
  l->AddEntry(h_charge_sum_sca_mb_vtx, "Scaled to run 23696");
  l->Draw("same");

  c->SaveAs(Form("%s/output/centplots/run%d_chargesum_scaled.png", env_p, runnumber));

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
  h_cent_bin_shifted = (TH1D*) file->Get("hcent_bins_sca");

  bool use_shifted = true;
  if (!h_cent_bin_shifted) use_shifted = false;
  h_cent_bin->SetTitle(";Centrality Bin; Fraction of Events");
  if (use_shifted)  SetMarkerAtt(h_cent_bin_shifted, kViolet - 2, 2, 89);
  SetMarkerAtt(h_cent_bin, kSpring + 2, 2, 90);
  if (use_shifted)  SetLineAtt(h_cent_bin_shifted, kViolet + 2, 2, 1);
  SetLineAtt(h_cent_bin, kSpring-1, 2, 1);
  const char* hist_labels[18] = {"0-5%","5-10%", "10-15%", "15-20%", "20-25%","25-30%", "30-35%", "35-40%","40-45%", "45-50%", "50-55%", "55-60%", "60-65%", "65-70%", "70-75", "75-80%", "80-85%","85-100%"}; 

  
  int nb = h_cent_bin->GetNbinsX();
  for (int i = 0; i < nb;i++)
    {
      h_cent_bin->GetXaxis()->SetBinLabel(i+1, hist_labels[i]);
    }
  
  h_cent_bin->SetMaximum(0.08);
  h_cent_bin->SetMinimum(0.04);
  h_cent_bin->Draw();
  if (use_shifted) h_cent_bin_shifted->Draw("same");
  TLine *l = new TLine(-0.5, 0.05, 17.5, 0.05);
  SetLineAtt(l, kBlack, 2, 7);
  l->Draw("same");

  TF1 *flatline = new TF1("flatline","[0]",-0.5, 17.5);
  flatline->SetParameter(0,0.05);

  h_cent_bin->Fit(flatline,"NDORQ","");
  float chi2_shifted = 0;
  float yy_shifted = 0;
  float chi2 = flatline->GetChisquare()/flatline->GetNDF();
  float yy = flatline->GetParameter(0);
  
  flatline->SetParameter(0,0.05);

  if (use_shifted)
    {
      //      h_cent_bin_shifted->Scale(1./h_cent_bin_shifted->Integral());
      h_cent_bin_shifted->Fit(flatline,"QNDOR","");

      chi2_shifted = flatline->GetChisquare()/flatline->GetNDF();;
      yy_shifted = flatline->GetParameter(0);
    }
  
  DrawSPHENIX(0.19, 0.85, 1, 0);
  drawText(Form("Run %d", runnumber), 0.19, 0.73, 0, kBlack);

  TLegend *lg = new TLegend(0.6, 0.65, 0.88, 0.89);
  SetLegendStyle(lg);
  lg->SetHeader("#chi^2 ; <Frac. of Events>");
  lg->AddEntry(h_cent_bin,Form("Not Shifted: %2.2f ; %1.4f", chi2, yy));
  if (use_shifted)   lg->AddEntry(h_cent_bin_shifted,Form("Shifted: %2.2f ; %1.4f", chi2_shifted, yy_shifted)); 
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


  TFile *file = new TFile(Form("%s/output/plots/mbd_charge_sum_%d.root", env_p, runnumber), "r");

  if (!file)
    {
      return;
    }
  
  TH1D *h_zdc_sum_n;
  TH1D *h_zdc_sum_s;
  TH1D *h_zdc_sum_n_prime;
  TH1D *h_zdc_sum_s_prime;

  TH1D *h_zdc_sum_n_w_vtx;
  TH1D *h_zdc_sum_s_w_vtx;
  TH1D *h_zdc_sum_n_prime_w_vtx;
  TH1D *h_zdc_sum_s_prime_w_vtx;

  h_zdc_sum_s = (TH1D*) file->Get("h_zdc_sum_s");
  h_zdc_sum_n = (TH1D*) file->Get("h_zdc_sum_n");
  h_zdc_sum_s_prime = (TH1D*) file->Get("h_zdc_sum_s_prime");
  h_zdc_sum_n_prime = (TH1D*) file->Get("h_zdc_sum_n_prime");
  h_zdc_sum_s_w_vtx = (TH1D*) file->Get("h_zdc_sum_s_w_vertex_30");
  h_zdc_sum_n_w_vtx = (TH1D*) file->Get("h_zdc_sum_n_w_vertex_30");
  h_zdc_sum_s_prime_w_vtx = (TH1D*) file->Get("h_zdc_sum_s_prime_w_vertex_30");
  h_zdc_sum_n_prime_w_vtx = (TH1D*) file->Get("h_zdc_sum_n_prime_w_vertex_30");

  h_zdc_sum_s_w_vtx->Scale(1./h_zdc_sum_s_w_vtx->Integral(5, 1000));
  h_zdc_sum_n_w_vtx->Scale(1./h_zdc_sum_n_w_vtx->Integral(5, 1000));
  h_zdc_sum_s_prime_w_vtx->Scale(1./h_zdc_sum_s_prime_w_vtx->Integral(5, 1000));
  h_zdc_sum_n_prime_w_vtx->Scale(1./h_zdc_sum_n_prime_w_vtx->Integral(5, 1000));
  /* h_zdc_sum_s_w_vtx->Scale(1./1020000. ); */
  /* h_zdc_sum_n_w_vtx->Scale(1./1020000. ); */
  /* h_zdc_sum_s_prime_w_vtx->Scale(1./1020000.); */
  /* h_zdc_sum_n_prime_w_vtx->Scale(1./1020000.); */
 
  TH1D *h_clone= (TH1D*) h_zdc_sum_s_w_vtx->Clone();
  h_clone->Reset();
  h_clone->Fill(100);
  float max = h_zdc_sum_s_w_vtx->GetBinContent(h_clone->GetMaximumBin());

  SetLineAtt(h_zdc_sum_n_w_vtx, kRed, 2, 1);
  SetLineAtt(h_zdc_sum_s_w_vtx, kBlue, 2, 1);
  SetLineAtt(h_zdc_sum_n_prime_w_vtx, kBlack, 2, 1);
  SetLineAtt(h_zdc_sum_s_prime_w_vtx, kGreen+2, 2, 1);
  h_zdc_sum_s_w_vtx->GetXaxis()->SetRangeUser(40, 10000);
  h_zdc_sum_s_w_vtx->SetMinimum(0.05 * max);
  h_zdc_sum_s_w_vtx->SetMaximum(100 * max);

  TCanvas *c = new TCanvas("c","c");
  c->SetLogx(1);
  c->SetLogy(1);

  h_zdc_sum_s_w_vtx->Draw("hist");
  h_zdc_sum_n_w_vtx->Draw("same hist");

  TLine *l = new TLine(100, 0.08*max, 100, max);
  SetLineAtt(l, kBlack, 3, 7);

  l->Draw("same");
  DrawSPHENIX(0.19, 0.85, 1, 0);
  drawText(Form("Run %d", runnumber), 0.19, 0.73, 0, kBlack);
  drawText("|z_{vtx}| < 30 cm", 0.19, 0.67, 0, kBlack);

  TLegend *lg = new TLegend(0.69, 0.7, 0.87, 0.87);
  lg->AddEntry(h_zdc_sum_n_w_vtx, "North","L");
  lg->AddEntry(h_zdc_sum_s_w_vtx, "South","L");
  lg->AddEntry(l, "Single Neutron","L");
  SetLegendStyle(lg);
  lg->Draw("same");

  c->SaveAs(Form("%s/output/centplots/run%d_zdc.png", env_p, runnumber));
  c->SaveAs(Form("%s/output/centplots/run%d_zdc.pdf", env_p, runnumber));

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
      h_charge[i]->Fit("lan_w_gausexp","","", .2, 3);
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
