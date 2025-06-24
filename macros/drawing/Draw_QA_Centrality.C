
#include "dlUtility.h"
#include "sPhenixStyle.C"
#include <string>
#include <map>
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
int ndivs = 93;
int isSim = 0;
const int mbd_ring_index[64] =                                                                                                                                                                                                                                                                                                                                                                              
   {2, 2, 2, 1, 1, 2, 1, 0,   
    0, 2, 1, 0, 2, 1, 0, 1,
    0, 2, 1, 0, 2, 1, 0, 2,
    1, 0, 0, 2, 1, 1, 2, 2,
    2, 2, 2, 1, 1, 2, 1, 0,
    0, 2, 1, 0, 2, 1, 0, 1,
    0, 2, 1, 0, 2, 1, 0, 2,
    1, 0, 0, 2, 1, 1, 2, 2};

const int ringcolors[3] = {kRed - 7, kSpring - 2, kBlue - 3};

void DrawMBDChargeSum(const int runnumber);
void DrawMBDChargeSum(const int runnumber, const int runnumber2);
void DrawMBDChannels(const int runnumber);
void DrawMBDVertex(const int runnumber);
void DrawMBDCentralityCalibrations(const int runnumber, bool use_shifted = false, bool use_balanced = false, bool flag = false);
void DrawMBDCentralityCheck(const int runnumber);
void DrawMCMBDCentralityCheck(const int runnumber);
void DrawZDCCheck(const int runnumber);
void DrawRunCentralityChecks(const std::string runlist);
void Draw_QA_Centrality(const std::string runlist)
{
  std::ifstream inputFile(runlist.c_str());
   if (!inputFile.is_open()) {
       std::cerr << "Error opening file!" << std::endl;
       return; // Exit if file cannot be opened
   }


   std::string line;

   while (std::getline(inputFile, line)) {
     std::stringstream ss(line);
     int runnumber;
     int runnumber_ref = 0;
     while (ss >> runnumber){// >> runnumber_ref) {
       std::cout << runnumber << " -- > " << runnumber_ref << std::endl;
       //DrawMBDChannels(runnumber);
       //// DrawZDCCheck(runnumber);
       //  DrawMBDChannels(runnumber);
       //DrawMBDChargeSum(runnumber);
       DrawMBDChargeSum(runnumber);
       //DrawMBDVertex(runnumber);
       //DrawMBDCentralityCalibrations(runnumber, 0, 0);
       //DrawMBDCentralityCalibrations(runnumber, 0, 1);
       DrawMBDCentralityCalibrations(runnumber, 0, 1);
       //DrawMBDCentralityCalibrations(runnumber,1);
       //DrawMBDCentralityCalibrations(runnumber,1,0,1);
       DrawMBDCentralityCheck(runnumber);

     }
   }
  return;
}

void Draw_QA_Centrality_2(const int runnumber){
  //DrawMBDChannels(runnumber);
  //DrawZDCCheck(runnumber);
  DrawMBDChargeSum(runnumber);
  DrawMBDCentralityCalibrations(runnumber,0, 1);
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
  TString filename = Form("%s/output_2024/plots/mbdana_centrality_trigeff%s_%d.root", env_p, extra.Data(), runnumber);
  if (runnumber < 10)
    {
      filename = Form("%s/output_2024/plots/mbdana_centrality_trigeff%s_hijing.root", env_p, extra.Data());
    }
  TFile *file = new TFile(filename.Data(), "r");
  std::cout << "opening :"<<Form("%s/output_2024/plots/mbdana_centrality_trigeff%s_%d.root", env_p, extra.Data(), runnumber)<<std::endl;
  if (!file)
    {
      cout<< "NOFILE" <<endl;
      return;
    }

  int maxcentbins = 99;
  TH1D *hSimMBDwTrig = (TH1D*) file->Get("hSimMBDwTrig");
  if (!hSimMBDwTrig) return;
  TH1D *hSimMBD = (TH1D*) file->Get("hSimMBD");
  if (!hSimMBD) return;
  TH1D *hRealMBD;
  TH1D *hRealMBDfine;

  TH1D *h_npart_cent100[99];
  TH1D *h_ecc2_cent100[99];
  TH1D *h_ecc3_cent100[99];
  TH1D *h_b_cent100[99];
  TH1D *h_npart_total = (TH1D*) file->Get("h_npart_total");
  TH1D *h_ecc2_total = (TH1D*) file->Get("h_ecc2_all");
  TH1D *h_ecc3_total = (TH1D*) file->Get("h_ecc3_all");
  TH1D *h_b_total = (TH1D*) file->Get("h_b_all");

  h_ecc2_total->Rebin(2);
  h_ecc3_total->Rebin(2);
  h_b_total->Rebin(2);

  TH1D *h_npart_cent[20];
  TH1D *h_ecc2_cent[20];
  TH1D *h_ecc3_cent[20];
  TH1D *h_b_cent[20];
  for (int i =0;i<100;i++)
    {
      h_npart_cent100[i] =(TH1D*) file->Get(Form("h_npart_cent_%d", i)); 

      if (!h_npart_cent100[i]) 
	{
	  std::cout << "no histogram for npart "<<i<<std::endl;
	  break;
	}
      h_ecc2_cent100[i] =(TH1D*) file->Get(Form("h_ecc2_%d", i)); 
      if (!h_ecc2_cent100[i]) 
	{
	  std::cout << "no histogram for ecc2 "<<i<<std::endl;
	  break;
	}
      h_ecc2_cent100[i]->Rebin(2);
      h_ecc3_cent100[i] =(TH1D*) file->Get(Form("h_ecc3_%d", i)); 
      if (!h_ecc3_cent100[i]) 
	{
	  std::cout << "no histogram for ecc3 "<<i<<std::endl;
	  break;
	}
      h_ecc3_cent100[i]->Rebin(2);
      h_b_cent100[i] =(TH1D*) file->Get(Form("h_b_%d", i)); 
      if (!h_b_cent100[i]) 
	{
	  std::cout << "no histogram for b "<<i<<std::endl;
	  break;
	}
      h_b_cent100[i]->Rebin(2);
      if (i%5 == 0)
	{
	  h_npart_cent[i/5] = (TH1D*) h_npart_cent100[i]->Clone();
	  h_npart_cent[i/5]->SetName(Form("h_npart_cent_%d_%d", i, i+5));

	  h_ecc2_cent[i/5] = (TH1D*) h_ecc2_cent100[i]->Clone();
	  h_ecc2_cent[i/5]->SetName(Form("h_ecc2_cent_%d_%d", i, i+5));

	  h_ecc3_cent[i/5] = (TH1D*) h_ecc3_cent100[i]->Clone();
	  h_ecc3_cent[i/5]->SetName(Form("h_ecc3_cent_%d_%d", i, i+5));

	  h_b_cent[i/5] = (TH1D*) h_b_cent100[i]->Clone();
	  h_b_cent[i/5]->SetName(Form("h_b_cent_%d_%d", i, i+5));

	}
      else
	{
	  h_npart_cent[i/5]->Add(h_npart_cent100[i]);
	  h_ecc2_cent[i/5]->Add(h_ecc2_cent100[i]);
	  h_ecc3_cent[i/5]->Add(h_ecc3_cent100[i]);
	  h_b_cent[i/5]->Add(h_b_cent100[i]);
	}
    }


  if (use_balanced && use_shifted) hRealMBD = (TH1D*) file->Get("h_charge_sum_min_bias_w_vertex_cut_balanced_scaled");
  else if (use_balanced) hRealMBD = (TH1D*) file->Get("h_charge_sum_min_bias_w_vertex_cut_balanced");
  else if (use_shifted) hRealMBD = (TH1D*) file->Get("h_charge_sum_min_bias_w_vertex_cut_scaled");
  else hRealMBD = (TH1D*) file->Get("h_charge_sum_min_bias_w_vertex_cut");
  if (!hRealMBD) return;
  if (use_balanced && use_shifted) hRealMBDfine = (TH1D*) file->Get("h_charge_sum_fine_min_bias_w_vertex_cut_balanced_scaled");
  else if (use_balanced) hRealMBDfine = (TH1D*) file->Get("h_charge_sum_fine_min_bias_w_vertex_cut_balanced");
  else if (use_shifted) hRealMBDfine = (TH1D*) file->Get("h_charge_sum_fine_min_bias_w_vertex_cut_scaled");
  else hRealMBDfine = (TH1D*) file->Get("h_charge_sum_fine_min_bias_w_vertex_cut");
  if (!hRealMBDfine) return;

  TH1D *hRatio = (TH1D*) file->Get("hRatiofine");
  if (!hRatio) return;

  hRatio->Rebin(10);
  hRatio->Scale(1./10.);
  
  TProfile *hpRatio = new TProfile("hpRatio","",35,0, 1500);
  for (int i = 0; i < hRatio->GetNbinsX(); i++)
    {
      hpRatio->Fill(hRatio->GetBinCenter(i), hRatio->GetBinContent(i));
    }

  TF1 *trigeffcurve = (TF1*) file->Get("trigeffcurve");
  if (!trigeffcurve) return;
  TString calib_file_name = Form("%s/calib_2024/mbdana_centrality_%d.root", env_p, runnumber);
  if (runnumber < 10)
    {
      calib_file_name = Form("%s/calib_2024/mbdana_centrality_hijing.root", env_p);
    }

  if (flag && use_balanced && use_shifted) calib_file_name = Form("%s/calib_2024/mbdana__centralitysca_bal_%d.root", env_p, 20869);
  else if (flag && use_balanced) calib_file_name = Form("%s/calib_2024/mbdana_centrality_bal_%d.root", env_p, 20869);
  else if (flag && use_shifted) calib_file_name = Form("%s/calib_2024/mbdana_centrality_%d.root", env_p, 20869); 
  else if (use_balanced && use_shifted) calib_file_name = Form("%s/calib_2024/mbdana_centrality_sca_bal_%d.root", env_p, runnumber);
  else if (use_balanced) calib_file_name = Form("%s/calib_2024/mbdana_centrality_bal_%d.root", env_p, runnumber);
  else if (use_shifted) calib_file_name = Form("%s/calib_2024/mbdana_centrality_sca_%d.root", env_p, runnumber);
  if (isSim)
    {

      if (flag && use_balanced && use_shifted) calib_file_name = Form("%s/calib_2024/mbdana_centralitysca_bal_%s.root", env_p, "hijing");
      else if (flag && use_balanced) calib_file_name = Form("%s/calib_2024/mbdana_centrality_bal_%s.root", env_p, "hijing");
      else if (flag && use_shifted) calib_file_name = Form("%s/calib_2024/mbdana_centrality_%s.root", env_p, "hijing"); 
      else if (use_balanced && use_shifted) calib_file_name = Form("%s/calib_2024/mbdana_centrality_sca_bal_%s.root", env_p, "hijing");
      else if (use_balanced) calib_file_name = Form("%s/calib_2024/mbdana_centrality_bal_%s.root", env_p, "hijing");
      else if (use_shifted) calib_file_name = Form("%s/calib_2024/mbdana_centrality_sca_%s.root", env_p, "hijing");

    }

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

  TString calib_file_name2 = Form("%s/calib_2024/mbdana_npart_%d.root", env_p, runnumber); 
  if (isSim)
    {
      calib_file_name2 = Form("%s/calib_2024/mbdana_npart_hijing.root", env_p); 
    }
  TFile *calibfile2 = new TFile(calib_file_name2, "r");

  TNtuple *tn2 = (TNtuple*) calibfile2->Get("tn_npart");
  tn2->SetBranchAddress("npart",&npart);

  for (int i = 0; i < tn2->GetEntries(); i++)
  {
    tn2->GetEntry(i);
    npart_cent[i] = npart;
  }  

  SetsPhenixStyle();


  // LOOK HERE
  TCanvas *c1 = new TCanvas("canvas_bbc_realdataandsim","canvas_bbc_realdataandsim", 700, 700);
  ratioPanelCanvas(c1, 0.4, 0, 0.15);
  
  c1->cd(1);
  gPad->SetTicks(1);	
  hRealMBD->SetLineWidth(2);
  hRealMBD->SetMarkerStyle(24);
  hRealMBD->SetMaximum(hSimMBD->GetBinContent(10)*3);
  hRealMBD->GetXaxis()->SetRangeUser(0.0,2500);
  gPad->SetTopMargin(.13);
  
  string foo = "MBD Charge";

  hRealMBD->SetTitle(Form(";%s;N_{events}", foo.c_str()));
  hRealMBD->GetXaxis()->CenterTitle(true);
  hRealMBD->GetYaxis()->CenterTitle(true);
  hRealMBD->SetTitleFont(42);
  hRealMBD->SetTitleSize(.07);
  hRealMBD->SetTitleOffset(1);
  
  hRealMBD->SetLabelSize(.05);

  hRealMBD->GetYaxis()->SetTitleFont(42);
  hRealMBD->GetYaxis()->SetTitleOffset(1);
  hRealMBD->GetYaxis()->SetTitleSize(.07);
  hRealMBD->GetYaxis()->SetLabelSize(.05);
  hRealMBD->DrawCopy("p,e,l");

  hSimMBD->SetLineWidth(2);
  hSimMBD->SetLineColor(kRed);
  
  hSimMBD->DrawCopy("l,same");

  int nhistbins = hRealMBD->GetNbinsX();
  double err_real; double err_sim;
  double trigeff = hRealMBD->IntegralAndError(1, nhistbins, err_real)/hSimMBD->IntegralAndError(1, nhistbins, err_sim);
  double trigeff_err = trigeff*sqrt(TMath::Power(err_real/hRealMBD->Integral(1, nhistbins), 2) + TMath::Power(err_sim/hSimMBD->Integral(1, nhistbins), 2));
  //  drawText("8/30/2023", 0.8, 0.92, 0, kBlack, 0.07);
  TLatex l;
  l.SetNDC();
  l.SetTextSize(0.04);
  
  DrawSPHENIXemma(0.65, 0.8, 1, 0, 0.06);
  drawText(Form("Run %d", runnumber), 0.65, 0.67, 0, kBlack, 0.06);
  TLegend *l1 = new TLegend(0.25, 0.67, 0.5, 0.8);
  l1->SetTextSize(0.06);
  l1->SetLineWidth(0);
  l1->AddEntry(hRealMBD, "Data", "p, e, l");
  l1->AddEntry(hSimMBD, "Glauber + NBD Fit", "p, e, l");
  l1->Draw("same");
  // now do the trigger efficiency fitting
  c1->cd(2);
  gPad->SetTicks(1);

  hpRatio->SetLineWidth(2);
  hpRatio->SetLineColor(kRed);
  hRatio->SetMarkerStyle(24);
  hRatio->SetMaximum(1.5);
  hRatio->SetMinimum(0);
  hRatio->GetXaxis()->SetRangeUser(0.0,400);
  hRatio->SetXTitle(foo.c_str());

  hRatio->SetYTitle("Data / MC");
  hRatio->SetTitleFont(42);
  hRatio->SetTitleSize(.09);
  hRatio->SetTitleOffset(0.8);
  hRatio->SetLabelSize(.07);
  hRatio->GetYaxis()->SetTitleFont(42);
  hRatio->GetYaxis()->SetTitleSize(.09);
  hRatio->GetYaxis()->SetTitleOffset(0.8);
  hRatio->GetYaxis()->SetLabelSize(.07);
  hRatio->GetYaxis()->CenterTitle(true);
  hRatio->GetXaxis()->CenterTitle(true);

  hRatio->Draw("p,e");
  //hpRatio->Draw("hist same");
  //  drawText(Form("Trigger Efficiency = %2.1f #pm %1.1f%%", trigeff*100, trigeff_err*100), 0.37, 0.8, 0, kBlack, 0.07);
  TLine *tl = new TLine(0.0,1.0,400,1.0);
  SetLineAtt(tl, kRed, 3, 4);
  tl->Draw();
  
  c1->cd(1);
  gPad->SetTicks(1);	
  gPad->SetLogy(1);

  hSimMBDwTrig->DrawCopy("same");
  float maxrange = 2500.;
  nhistbins = hRealMBD->GetNbinsX();
  TH1D *hslice = new TH1D("hslice","hslice",nhistbins,-0.5,maxrange - 0.5); 

  float in_data[100] = {};
  float in_sim[100] = {};
  for (int j = 0; j < 100; j++)
    {
      in_data[j] = 0;
      in_sim[j] = 0;
    }
  
  for (int j=0; j< ndivs; j++) {
    if (j%5 == 0)
      {
	hslice->Reset();
      }
    for (int i = 1+(floor(centrality_high[j+1])); i< 1+(floor(centrality_high[j]));i++)
      {
	hslice->SetBinContent(i,hSimMBDwTrig->GetBinContent(i));
	in_data[j] += hRealMBD->GetBinContent(i);
	in_sim[j] += hSimMBD->GetBinContent(i);
      }
    cout << "Centbin = " << j << " Lowcut = " << 
      hslice->GetBinCenter(1+(floor(centrality_low[j]))) << 
      " Highcut = " << 
      hslice->GetBinCenter(1+(floor(centrality_low[j+1]))) << " -- "
	 << in_data[j] <<" / " << in_sim[j] << endl;
    
    if ((j + 1)%5 == 0)
      {
	hslice->SetFillColorAlpha(100 - j,0.3);
	hslice->DrawCopy("same");
      }
  }

  // redraw real data to be on top
  hRealMBD->DrawCopy("p,e,l,same");
  hSimMBD->DrawCopy("hist,same");  
  c1->SaveAs(Form("%s/output_2024/centrality_2024/run%d_glauber%s.pdf", env_p, runnumber, extra.Data()));
  c1->SaveAs(Form("%s/output_2024/centrality_2024/run%d_glauber%s.png", env_p, runnumber, extra.Data()));
  //===============================================================================


  float npart_cents[16] = {
    350.0,
    301.4,
    256.8,
    218.3,
    185.0,
    155.2,
    129.2,
    106.3,
    86.3,
    68.7,
    53.7,
    40.9,
    30.3,
    21.7,
    15.2,
    10.2 
  };

  TCanvas *c2 = new TCanvas("c2","c2", 800, 600);
  gPad->SetLogy();
  SetyjPadStyle();
  gPad->SetRightMargin(0.28);
 
  SetLineAtt(h_npart_total, kBlack, 3, 1);
  h_npart_total->SetMarkerStyle(21);
  h_npart_total->SetMarkerSize(1);
      
  h_npart_total->SetTitle(";N_{part};counts");
  h_npart_total->Draw();
  gStyle->SetPalette(kRainBow);
  
  TLegend *ll = new TLegend(0.73,0.15, 0.99, 0.95);
  ll->SetHeader("Cent % - <N_{part}> #pm #sigma(N_{part})");
  for (int i = 0; i < 16; i++)
    {

      h_npart_cent[i]->SetMarkerStyle(21);
      h_npart_cent[i]->SetMarkerSize(1);
      h_npart_cent[i]->Draw("same PMC PLC");
      ll->AddEntry(h_npart_cent[i], Form("%d-%d %% - %3.1f #pm %3.1f", i*5, i*5+5, npart_cents[i], h_npart_cent[i]->GetRMS()),"PMC PLC"); 
    }
  SetLegendStyle(ll);
  ll->Draw("same");
  DrawSPHENIX(0.43, 0.88, 1, 0, 0.04);
  
  c2->SaveAs(Form("%s/output_2024/centrality_2024/run%d_npart_glauber%s.pdf", env_p, runnumber, extra.Data()));
  c2->SaveAs(Form("%s/output_2024/centrality_2024/run%d_npart_glauber%s.png", env_p, runnumber, extra.Data()));


  // ecc2
 
  SetLineAtt(h_ecc2_total, kBlack, 3, 1);
  h_ecc2_total->SetMarkerStyle(21);
  h_ecc2_total->SetMarkerSize(1);
      
  h_ecc2_total->SetTitle(";N_{part};counts");
  h_ecc2_total->Draw();
  gStyle->SetPalette(kRainBow);
  gPad->SetLogy(0);
  ll = new TLegend(0.73,0.15, 0.99, 0.95);
  ll->SetHeader("Cent % - <#epsilon_{2}> #pm #sigma(#epsilon_{2})");
  for (int i = 0; i < 18; i++)
    {

      h_ecc2_cent[i]->SetMarkerStyle(21);
      h_ecc2_cent[i]->SetMarkerSize(1);
      h_ecc2_cent[i]->Draw("same PMC PLC");
      ll->AddEntry(h_ecc2_cent[i], Form("%d-%d %% - %1.2f #pm %1.2f", i*5, i*5+5, h_ecc2_cent[i]->GetMean(), h_ecc2_cent[i]->GetRMS()),"PMC PLC"); 
    }
  SetLegendStyle(ll);
  ll->Draw("same");
  DrawSPHENIX(0.43, 0.88, 1, 0, 0.04);
  
  c2->SaveAs(Form("%s/output_2024/centrality_2024/run%d_ecc2_glauber%s.pdf", env_p, runnumber, extra.Data()));
  c2->SaveAs(Form("%s/output_2024/centrality_2024/run%d_ecc2_glauber%s.png", env_p, runnumber, extra.Data()));

  // ecc3
 
  SetLineAtt(h_ecc3_total, kBlack, 3, 1);
  h_ecc3_total->SetMarkerStyle(21);
  h_ecc3_total->SetMarkerSize(1);
      
  h_ecc3_total->SetTitle(";N_{part};counts");
  h_ecc3_total->Draw();
  gStyle->SetPalette(kRainBow);
  
  ll = new TLegend(0.73,0.15, 0.99, 0.95);
  ll->SetHeader("Cent % - <#epsilon_{3}> #pm #sigma(#epsilon_{3})");
  for (int i = 0; i < 18; i++)
    {

      h_ecc3_cent[i]->SetMarkerStyle(21);
      h_ecc3_cent[i]->SetMarkerSize(1);
      h_ecc3_cent[i]->Draw("same PMC PLC");
      ll->AddEntry(h_ecc3_cent[i], Form("%d-%d %% - %1.2f #pm %1.2f", i*5, i*5+5, h_ecc3_cent[i]->GetMean(), h_ecc3_cent[i]->GetRMS()),"PMC PLC"); 
    }
  SetLegendStyle(ll);
  ll->Draw("same");
  DrawSPHENIX(0.43, 0.88, 1, 0, 0.04);
  
  c2->SaveAs(Form("%s/output_2024/centrality_2024/run%d_ecc3_glauber%s.pdf", env_p, runnumber, extra.Data()));
  c2->SaveAs(Form("%s/output_2024/centrality_2024/run%d_ecc3_glauber%s.png", env_p, runnumber, extra.Data()));


  // ecc2
 
  SetLineAtt(h_b_total, kBlack, 3, 1);
  h_b_total->SetMarkerStyle(21);
  h_b_total->SetMarkerSize(1);
      
  h_b_total->SetTitle(";b [fm];counts");
  h_b_total->Draw();
  gStyle->SetPalette(kRainBow);
  
  ll = new TLegend(0.73,0.15, 0.99, 0.95);
  ll->SetHeader("Cent % - <b> #pm #sigma(b)");
  for (int i = 0; i < 18; i++)
    {

      h_b_cent[i]->SetMarkerStyle(21);
      h_b_cent[i]->SetMarkerSize(1);
      h_b_cent[i]->Draw("same PMC PLC");
      ll->AddEntry(h_b_cent[i], Form("%d-%d %% - %2.2f #pm %2.2f", i*5, i*5+5, h_b_cent[i]->GetMean(), h_b_cent[i]->GetRMS()),"PMC PLC"); 
    }
  SetLegendStyle(ll);
  ll->Draw("same");
  DrawSPHENIX(0.43, 0.88, 1, 0, 0.04);
  
  c2->SaveAs(Form("%s/output_2024/centrality_2024/run%d_b_glauber%s.pdf", env_p, runnumber, extra.Data()));
  c2->SaveAs(Form("%s/output_2024/centrality_2024/run%d_b_glauber%s.png", env_p, runnumber, extra.Data()));


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

  TFile *file = new TFile(Form("%s/output_2024/plots/histout_%d.root", env_p, runnumber), "r");

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
      h_time_raw[i]->GetXaxis()->SetRangeUser(-25, 25);
      SetLineAtt(h_time_raw[i], colorline, 2, 1);
      h_time_raw[i]->SetFillColorAlpha(colorfill, 0.3);
      h_time_raw[i]->Draw();
      
      drawText(Form("S Ch. %d", i),0.1, 0.88, 0, kBlack, 0.07);
    }

  c->Print(Form("%s/output_2024/centrality_2024/run%d_time.pdf(", env_p, runnumber));
  c->SaveAs(Form("%s/output_2024/centrality_2024/run%d_time_south.png", env_p, runnumber));

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
      h_time_raw[i+64]->GetXaxis()->SetRangeUser(-25, 25);
      SetLineAtt(h_time_raw[i+64], colorline, 1, 1);
      h_time_raw[i+64]->SetFillColorAlpha(colorfill, 0.3);
      h_time_raw[i+64]->Draw();
      drawText(Form("N Ch. %d", i),0.1, 0.88, 0, kBlack, 0.07);
    }

  c->SaveAs(Form("%s/output_2024/centrality_2024/run%d_time_north.png", env_p, runnumber));
  c->Print(Form("%s/output_2024/centrality_2024/run%d_time.pdf)", env_p, runnumber));

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


  TFile *file = new TFile(Form("%s/output_2024/plots/mbdana_zdc_check_%d.root", env_p, runnumber), "r");

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
  c_time->SaveAs(Form("%s/output_2024/centrality_2024/run%d_time_0.pdf", env_p, runnumber));

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
  c_vertex->SaveAs(Form("%s/output_2024/centrality_2024/run%d_zvtx.pdf", env_p, runnumber));
  gPad->SetLogy();
  c_vertex->SaveAs(Form("%s/output_2024/centrality_2024/run%d_zvtx_log.pdf", env_p, runnumber));
  
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

  TString filename = Form("%s/output_2024/plots/mbdana_charge_sum_%d.root", env_p, runnumber);
  if (isSim)
    {
      filename = Form("%s/output_2024/plots/mbdana_charge_sum_hijing.root", env_p);
    }
  TFile *file = new TFile(filename.Data(), "r");
  if (!file) return;
  TH1D *h_charge_sum = (TH1D*) file->Get("h_charge_sum");
  if (!h_charge_sum) return;
  TH1D *h_charge_sum_vtx = (TH1D*) file->Get("h_charge_sum_vtx");
  if (!h_charge_sum_vtx) return;
  TH2D *h2_charge_sum_vtx_ns = (TH2D*) file->Get("h2_charge_sum_vtx_ns");
  if (!h2_charge_sum_vtx_ns) 
    {
      std::cout << "AHJHHHHH" << std::endl;
      return;
    }
  TH1D *h_charge_sum_mb = (TH1D*) file->Get("h_charge_sum_min_bias");
  if (!h_charge_sum_mb) return;
  TH1D *h_charge_sum_mb_vtx = (TH1D*) file->Get("h_charge_sum_min_bias_w_vertex_cut");
  if (!h_charge_sum_mb_vtx) return;  

  TH2D *h_charge_sum_v_vtx_mb_vtx = (TH2D*) file->Get("h_charge_sum_v_vtx_min_bias_w_vertex_cut");
  if (!h_charge_sum_v_vtx_mb_vtx) return;  

  TH2D *h_charge_sum_balanced_v_vtx_mb_vtx = (TH2D*) file->Get("h_charge_sum_balanced_v_vtx_min_bias_w_vertex_cut");
  if (!h_charge_sum_balanced_v_vtx_mb_vtx) return;  

  TH1D *h_charge_sum_tv[3];
  TH1D *h_charge_sum_vtx_tv[3];
  TH1D *h_charge_sum_mb_tv[3];
  TH1D *h_charge_sum_mb_vtx_tv[3];

  for (int i = 0; i < 3; i++)
    {
      h_charge_sum_tv[i] = (TH1D*) file->Get(Form("h_charge_sum_trigvtx_%d", i));
      if (!h_charge_sum_tv[i]) return;
      h_charge_sum_vtx_tv[i] = (TH1D*) file->Get(Form("h_charge_sum_vtx_trigvtx_%d", i));
      if (!h_charge_sum_vtx_tv[i]) return;
      h_charge_sum_mb_tv[i] = (TH1D*) file->Get(Form("h_charge_sum_min_bias_trigvtx_%d", i));
      if (!h_charge_sum_mb_tv[i]) return;
      h_charge_sum_mb_vtx_tv[i] = (TH1D*) file->Get(Form("h_charge_sum_min_bias_w_vertex_cut_trigvtx_%d", i));
      if (!h_charge_sum_mb_vtx_tv[i]) return;
    }


  TH1D *h_charge_sum_sca = (TH1D*) file->Get("h_charge_sum_scaled");
  TH1D *h_charge_sum_sca_vtx = (TH1D*) file->Get("h_charge_sum_vtx_scaled");
  TH1D *h_charge_sum_sca_mb_vtx = (TH1D*) file->Get("h_charge_sum_min_bias_w_vertex_cut_scaled");
  
  TH2D *h_charge_sum_v_zdc =  (TH2D*) file->Get("h_charge_sum_v_zdc");
  TH2D *h_charge_sum_v_zdc_mb =  (TH2D*) file->Get("h_charge_sum_v_zdc_mb");

  TProfile *hp_charge_sum_v_zdc =  (TProfile*) file->Get("hp_charge_sum_v_zdc");
  TProfile *hp_charge_sum_v_zdc_mb =  (TProfile*) file->Get("hp_charge_sum_v_zdc_mb");

  TFile *file2 = new TFile(Form("%s/output_2024/plots/mbdana_charge_sum_54912.root", env_p), "r");
  TProfile *hp_charge_sum_v_zdc_mbref =  (TProfile*) file2->Get("hp_charge_sum_v_zdc_mb");
  hp_charge_sum_v_zdc_mbref->GetXaxis()->SetRangeUser(10, 1500);
  hp_charge_sum_v_zdc_mb->GetXaxis()->SetRangeUser(10, 1500);

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

  c->SaveAs(Form("%s/output_2024/centrality_2024/run%d_charge_sum.pdf", env_p, runnumber));

  c->Clear();

  SetyjPadStyle();
  gPad->SetLogz();
  gPad->SetRightMargin(0.2);

  h2_charge_sum_vtx_ns->SetTitle(";MBD Charge Sum South; MBD Charge Sum North");
  h2_charge_sum_vtx_ns->Rebin2D(4,4);
  h2_charge_sum_vtx_ns->GetXaxis()->SetRangeUser(0, 1500);
  h2_charge_sum_vtx_ns->GetYaxis()->SetRangeUser(0, 1500);

  h2_charge_sum_vtx_ns->Draw("colz");

  DrawSPHENIX(0.22, 0.85, 1, 0, 0.04);
  drawText(Form("Run %d", runnumber), 0.22, 0.74, 0, kBlack, 0.04);
  drawText("MBD NS #geq 2", 0.22, 0.69, 0, kBlack, 0.04);
  drawText("|z_{vtx}| < 60 cm", 0.22, 0.64, 0, kBlack, 0.04);
  c->SaveAs(Form("%s/output_2024/centrality_2024/run%d_charge_sum_2d_bkg.pdf", env_p, runnumber));


  c->Clear();
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
  l->AddEntry(h_charge_sum_sca_mb_vtx, "Scaled to run 54912");
  l->Draw("same");

  c->SaveAs(Form("%s/output_2024/centrality_2024/run%d_charge_sum_scaled.pdf", env_p, runnumber));
  if (h_charge_sum_v_zdc_mb)
    {
      c->Clear();
      gPad->SetLogx();
      gPad->SetLogy(0);
      gPad->SetRightMargin(0.16);

      h_charge_sum_v_zdc_mb->Rebin2D(1, 4);
      h_charge_sum_v_zdc_mb->GetZaxis()->SetTitleOffset(1.5);
      h_charge_sum_v_zdc_mb->GetXaxis()->SetTitleOffset(1.5);
      h_charge_sum_v_zdc_mb->GetYaxis()->SetMaxDigits(2);
      h_charge_sum_v_zdc_mb->SetTitle(";MBD Charge Sum; ZDC Energy Sum [GeV]; N_{events}");
      h_charge_sum_v_zdc_mb->Rebin2D(2,2);//Draw("same");
      h_charge_sum_v_zdc_mb->Draw("colz");
      SetLineAtt(hp_charge_sum_v_zdc_mb, kRed, 1,1);
      SetLineAtt(hp_charge_sum_v_zdc_mbref, kBlue - 2, 1,1);
      SetMarkerAtt(hp_charge_sum_v_zdc_mb, kRed, 1,75);
      SetMarkerAtt(hp_charge_sum_v_zdc_mbref, kBlue - 2, 1,75);
      hp_charge_sum_v_zdc_mb->GetXaxis()->SetRangeUser(10, 1700);
      hp_charge_sum_v_zdc_mbref->GetXaxis()->SetRangeUser(10, 1700);

      hp_charge_sum_v_zdc_mb->Draw("same");
      //hp_charge_sum_v_zdc_mbref->Draw("same");
      DrawSPHENIXemma(0.45, 0.95, 1, 1, 0.035);
      drawText(Form("Run %d", runnumber), 0.25, 0.95, 0, kBlack, 0.035);
      //drawText(Form("Run %d", runnumber), 0.22, 0.77, 0, kBlack, 0.04);
      l = new TLegend(0.22, 0.2, 0.5, 0.32);
      //l->AddEntry(hp_charge_sum_v_zdc_mbref,"Ref. run - 54912","P");
      //l->AddEntry(hp_charge_sum_v_zdc_mb,"This run","P");
      SetLegendStyle(l);
      //l->Draw("same");

      c->SaveAs(Form("%s/output_2024/centrality_2024/run%d_charge_sum_zdc.png", env_p, runnumber));
      c->SaveAs(Form("%s/output_2024/centrality_2024/run%d_charge_sum_zdc.pdf", env_p, runnumber));
    }

  TCanvas *c_charge = new TCanvas("c_charge","c_charge", 500, 700);

  ratioPanelCanvas(c_charge);

  int colors[3] = {kRed - 7 ,  kViolet - 4, kAzure + 2};
  SetLineAtt(h_charge_sum, kBlack, 1, 1);
  SetLineAtt(h_charge_sum_mb, kBlack, 1, 1);
  SetLineAtt(h_charge_sum_mb_vtx, kBlack, 1, 1);
  SetLineAtt(h_charge_sum_vtx, kBlack, 1, 1);
  TH1D *h_charge_sum_coarse = (TH1D*) h_charge_sum->Clone();
  h_charge_sum_coarse->Rebin(10);
  TH1D *h_charge_sum_mb_coarse = (TH1D*) h_charge_sum_mb->Clone();
  h_charge_sum_mb_coarse->Rebin(10);
  TH1D *h_charge_sum_vtx_coarse = (TH1D*) h_charge_sum_vtx->Clone();
  h_charge_sum_vtx_coarse->Rebin(10);
  TH1D *h_charge_sum_mb_vtx_coarse = (TH1D*) h_charge_sum_mb_vtx->Clone();
  h_charge_sum_mb_vtx_coarse->Rebin(10);
  h_charge_sum->SetMaximum(1000000);
  h_charge_sum_vtx->SetMaximum(1000000);
  h_charge_sum_mb->SetMaximum(1000000);
  h_charge_sum_mb_vtx->SetMaximum(1000000);
  TH1D *h_charge_sum_tv_coarse[3];
  TH1D *h_charge_sum_vtx_tv_coarse[3];
  TH1D *h_charge_sum_mb_tv_coarse[3];
  TH1D *h_charge_sum_mb_vtx_tv_coarse[3];
  for (int i = 0; i < 3; i++) 
    {
      SetLineAtt(h_charge_sum_tv[i], colors[i], 1, 1);
      SetLineAtt(h_charge_sum_vtx_tv[i], colors[i], 1, 1);
      SetLineAtt(h_charge_sum_mb_tv[i], colors[i], 1, 1);
      SetLineAtt(h_charge_sum_mb_vtx_tv[i], colors[i], 1, 1);
      h_charge_sum_tv_coarse[i] = (TH1D*) h_charge_sum_tv[i]->Clone();
      h_charge_sum_tv_coarse[i]->Rebin(10);

      h_charge_sum_vtx_tv_coarse[i] = (TH1D*) h_charge_sum_vtx_tv[i]->Clone();
      h_charge_sum_vtx_tv_coarse[i]->Rebin(10);

      h_charge_sum_mb_tv_coarse[i] = (TH1D*) h_charge_sum_mb_tv[i]->Clone();
      h_charge_sum_mb_tv_coarse[i]->Rebin(10);

      h_charge_sum_mb_vtx_tv_coarse[i] = (TH1D*) h_charge_sum_mb_vtx_tv[i]->Clone();
      h_charge_sum_mb_vtx_tv_coarse[i]->Rebin(10);


    }

  TGraphAsymmErrors *h_vertex_ratios[3];
  TGraphAsymmErrors *h_vertex_ratios_mb[3];
  TGraphAsymmErrors *h_vertex_ratios_vtx[3];
  TGraphAsymmErrors *h_vertex_ratios_mb_vtx[3];

  for (int i = 0; i < 3; i++)
    {
      h_vertex_ratios[i] = new TGraphAsymmErrors(h_charge_sum_tv_coarse[i],h_charge_sum_coarse,"cl=0.683 b(1,1) mode");      
      h_vertex_ratios_vtx[i] = new TGraphAsymmErrors(h_charge_sum_vtx_tv_coarse[i],h_charge_sum_vtx_coarse,"cl=0.683 b(1,1) mode");      
      h_vertex_ratios_mb[i] = new TGraphAsymmErrors(h_charge_sum_mb_tv_coarse[i],h_charge_sum_mb_coarse,"cl=0.683 b(1,1) mode");      
      h_vertex_ratios_mb_vtx[i] = new TGraphAsymmErrors(h_charge_sum_mb_vtx_tv_coarse[i],h_charge_sum_mb_vtx_coarse,"cl=0.683 b(1,1) mode");      
      h_vertex_ratios[i]->SetLineColor(colors[i]);
      h_vertex_ratios_mb[i]->SetLineColor(colors[i]);
      h_vertex_ratios_vtx[i]->SetLineColor(colors[i]);
      h_vertex_ratios_mb_vtx[i]->SetLineColor(colors[i]);
    }
  
  c_charge->cd(1);
  gPad->SetLogy();
  h_charge_sum->Draw();
  for (int i = 0; i < 3; i++)
    {
      h_charge_sum_tv[i]->Draw("same");
    }
  DrawSPHENIX(0.6, 0.85, 1, 0, 0.04);
  drawText(Form("Run %d", runnumber), 0.6, 0.74, 0, kBlack, 0.04);
  drawText("N_{tube}^{N,S} >= 2", 0.6, 0.69, 0, kBlack, 0.04);
  l = new TLegend(0.22, 0.7, 0.55, 0.88);
  l->SetLineWidth(0);
  l->AddEntry(h_charge_sum, "MBD NS >= 2");
  l->AddEntry(h_charge_sum_tv[0], "MBD NS >= 2 w/ |z_{vtx}| < 10 cm");
  l->AddEntry(h_charge_sum_tv[1], "MBD NS >= 2 w/ |z_{vtx}| < 60 cm");
  l->AddEntry(h_charge_sum_tv[2], "MBD NS >= 2 w/ |z_{vtx}| < 150 cm");
  l->Draw("same");
  c_charge->cd(2);
  TH1D *h = (TH1D*) h_charge_sum->Clone();
  h->GetXaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleSize(0.06);
  h->Reset();
  h->SetMaximum(1.3);
  h->SetMinimum(0.0);
  h->Draw();
  for (int i = 0; i < 3; i++)
    {
      h_vertex_ratios[i]->Draw("same");
    }

  c_charge->SaveAs(Form("%s/output_2024/centrality_2024/run%d_charge_sum_vertex.pdf", env_p, runnumber));
  c_charge->cd(1);
  gPad->SetLogy();
  h_charge_sum_mb->Draw();
  for (int i = 0; i < 3; i++)
    {
      h_charge_sum_mb_tv[i]->Draw("same");
    }
  DrawSPHENIX(0.6, 0.85, 1, 0, 0.04);
  drawText(Form("Run %d", runnumber), 0.6, 0.74, 0, kBlack, 0.04);
  drawText("N_{tube}^{N,S} >= 2", 0.6, 0.69, 0, kBlack, 0.04);
  drawText("#Sigma E_{ZDC}^{N,S} > 40 GeV", 0.6, 0.64, 0, kBlack, 0.04);
  drawText("MBD Bkg Cut", 0.6, 0.59, 0, kBlack, 0.04);

  l->Draw("same");
  c_charge->cd(2);
  h->Draw();
  for (int i = 0; i < 3; i++)
    {
      h_vertex_ratios_mb[i]->Draw("same");
    }

  c_charge->SaveAs(Form("%s/output_2024/centrality_2024/run%d_charge_sum_vertex_mb.pdf", env_p, runnumber));

  c_charge->cd(1);
  gPad->SetLogy();

  h_charge_sum_vtx->Draw();
  for (int i = 0; i < 3; i++)
    {
      h_charge_sum_vtx_tv[i]->Draw("same");
    }
  DrawSPHENIX(0.6, 0.85, 1, 0, 0.04);
  drawText(Form("Run %d", runnumber), 0.6, 0.74, 0, kBlack, 0.04);
  drawText("N_{tube}^{N,S} >= 2", 0.6, 0.69, 0, kBlack, 0.04);
  drawText("MBD Bkg Cut", 0.6, 0.64, 0, kBlack, 0.04);

  l->Draw("same");
  c_charge->cd(2);
  h->Draw();
  for (int i = 0; i < 3; i++)
    {
      h_vertex_ratios_vtx[i]->Draw("same");
    }

  c_charge->SaveAs(Form("%s/output_2024/centrality_2024/run%d_charge_sum_vertex_vtx.pdf", env_p, runnumber));


  c_charge->cd(1);
  gPad->SetLogy();
  h_charge_sum_mb_vtx->Draw();
  for (int i = 0; i < 3; i++)
    {
      h_charge_sum_mb_vtx_tv[i]->Draw("same");
    }
  DrawSPHENIX(0.6, 0.85, 1, 0, 0.04);
  drawText(Form("Run %d", runnumber), 0.6, 0.74, 0, kBlack, 0.04);
  drawText("N_{tube}^{N,S} >= 2", 0.6, 0.69, 0, kBlack, 0.04);
  drawText("#Sigma E_{ZDC}^{N,S} > 40 GeV", 0.6, 0.64, 0, kBlack, 0.04);
  drawText("MBD Bkg Cut", 0.6, 0.59, 0, kBlack, 0.04);
  drawText("|z_{vtx}| < 60 cm", 0.6, 0.54, 0, kBlack, 0.04);
  l->Draw("same");
  c_charge->cd(2);

  h->Draw();
  for (int i = 0; i < 3; i++)
    {
      h_vertex_ratios_mb_vtx[i]->Draw("same");
    }

  c_charge->SaveAs(Form("%s/output_2024/centrality_2024/run%d_charge_sum_vertex_mb_vtx.pdf", env_p, runnumber));

  // Integrate and normalize each coloumn.
  
  TCanvas *cc1 = new TCanvas("cc1","cc1", 500, 500);
  cc1->SetLogz();
  gPad->SetRightMargin(0.18);
  h_charge_sum_v_vtx_mb_vtx->SetTitle(";Vertex [cm]; MBD Charge Sum");
  TProfile *h_charge_p = (TProfile*) h_charge_sum_v_vtx_mb_vtx->ProfileX();
  h_charge_sum_v_vtx_mb_vtx->RebinY(10);
  for (int ib = 1; ib <= h_charge_sum_v_vtx_mb_vtx->GetXaxis()->GetNbins(); ib++)
    {
      double scale = h_charge_sum_v_vtx_mb_vtx->Integral(ib, ib, 1, -1);
      for (int iby = 1; iby <= h_charge_sum_v_vtx_mb_vtx->GetYaxis()->GetNbins(); iby++)
	{
	  int bin = h_charge_sum_v_vtx_mb_vtx->GetBin(ib, iby);
	  h_charge_sum_v_vtx_mb_vtx->SetBinContent(ib, iby, h_charge_sum_v_vtx_mb_vtx->GetBinContent(bin)/scale);
	}
    }

  h_charge_sum_v_vtx_mb_vtx->Draw("colz");
   
  h_charge_p->SetLineColor(kRed);
  h_charge_p->Draw("same");

  DrawSPHENIX(0.2, 0.95, 1, 1, 0.03);
  drawText(Form("Run %d", runnumber), 0.2, 0.85, 0, kBlack, 0.03);
  drawText("MBD NS >=2,#Sigma E_{ZDC} > 40 GeV", 0.2, 0.81, 0, kBlack, 0.03);
  TLegend *lred = new TLegend(0.6, 0.81, 0.8, 0.88);
  lred->SetLineWidth(0);
  lred->AddEntry(h_charge_p, "<#Sigma Q_{MBD}>");
  lred->Draw("same");
  cc1->SaveAs(Form("%s/output_2024/centrality_2024/run%d_charge_sum_v_vtx.pdf", env_p, runnumber));

  TCanvas *cc2 = new TCanvas("cc2","cc2", 500, 500);
  cc2->SetLogz();
  gPad->SetRightMargin(0.18);
  h_charge_sum_balanced_v_vtx_mb_vtx->SetTitle(";Vertex [cm]; MBD Charge Sum Vertex Corrected");
  TProfile *h_charge_bal_p = (TProfile*) h_charge_sum_balanced_v_vtx_mb_vtx->ProfileX();
  h_charge_sum_balanced_v_vtx_mb_vtx->RebinY(10);
  const int vertexbins = h_charge_sum_balanced_v_vtx_mb_vtx->GetXaxis()->GetNbins();
  TH1D *vertex_charge_sum_balanced[vertexbins];  
  for (int ib = 1; ib <= h_charge_sum_balanced_v_vtx_mb_vtx->GetXaxis()->GetNbins(); ib++)
    {

      double scale = h_charge_sum_balanced_v_vtx_mb_vtx->Integral(ib, ib, 1, -1);
      for (int iby = 1; iby <= h_charge_sum_balanced_v_vtx_mb_vtx->GetYaxis()->GetNbins(); iby++)
	{
	  int bin = h_charge_sum_balanced_v_vtx_mb_vtx->GetBin(ib, iby);
	  h_charge_sum_balanced_v_vtx_mb_vtx->SetBinContent(ib, iby, h_charge_sum_balanced_v_vtx_mb_vtx->GetBinContent(bin)/scale);
	}
      vertex_charge_sum_balanced[ib - 1] = (TH1D*) h_charge_sum_balanced_v_vtx_mb_vtx->ProjectionY(Form("h_proj_charge_%d", ib), ib, ib);
    }

  h_charge_sum_balanced_v_vtx_mb_vtx->Draw("colz");
   
  h_charge_bal_p->SetLineColor(kRed);
  h_charge_bal_p->Draw("same");

  DrawSPHENIX(0.2, 0.95, 1, 1, 0.03);
  drawText(Form("Run %d", runnumber), 0.2, 0.85, 0, kBlack, 0.03);
  drawText("MBD NS >=2,#Sigma E_{ZDC} > 40 GeV", 0.2, 0.81, 0, kBlack, 0.03);
  lred = new TLegend(0.6, 0.81, 0.8, 0.88);
  lred->SetLineWidth(0);
  lred->AddEntry(h_charge_bal_p, "<#Sigma Q_{MBD}>");
  lred->Draw("same");
  cc2->SaveAs(Form("%s/output_2024/centrality_2024/run%d_charge_sum_balanced_v_vtx.pdf", env_p, runnumber));

  TCanvas *cc3 = new TCanvas("cc3","cc3", 500, 700);
  ratioPanelCanvas(cc3);
  cc3->cd(1);
  gPad->SetLogy();
  vertex_charge_sum_balanced[7]->Rebin(5);
  vertex_charge_sum_balanced[8]->Rebin(5);
  vertex_charge_sum_balanced[9]->Rebin(5);
  vertex_charge_sum_balanced[10]->Rebin(5);
  SetLineAtt(vertex_charge_sum_balanced[7], kBlack, 1, 1);
  SetLineAtt(vertex_charge_sum_balanced[8], colors[0], 1, 1);
  SetLineAtt(vertex_charge_sum_balanced[9], colors[1], 1, 1);
  SetLineAtt(vertex_charge_sum_balanced[10], colors[2], 1, 1);
  vertex_charge_sum_balanced[7]->SetMaximum(vertex_charge_sum_balanced[7]->GetBinContent(vertex_charge_sum_balanced[7]->GetMaximumBin())*50);
  vertex_charge_sum_balanced[7]->SetTitle(";MBD Charge Sum;Normalized by Vertex region");
  vertex_charge_sum_balanced[7]->Draw();
  vertex_charge_sum_balanced[8]->Draw("same");
  vertex_charge_sum_balanced[9]->Draw("same");
  vertex_charge_sum_balanced[10]->Draw("same");

  lred = new TLegend(0.2, 0.69, 0.5, 0.88);
  lred->SetLineWidth(0);

  lred->AddEntry(vertex_charge_sum_balanced[7], "-1 < z_{vtx} < 1");
  lred->AddEntry(vertex_charge_sum_balanced[8], "1 < z_{vtx} < 5");
  lred->AddEntry(vertex_charge_sum_balanced[9], "5 < z_{vtx} < 10");
  lred->AddEntry(vertex_charge_sum_balanced[10], "10 < z_{vtx} < 15");
  lred->SetTextSize(0.03);
  lred->Draw("same");
  DrawSPHENIX(0.6, 0.85, 1, 0, 0.04);
  drawText(Form("Run %d", runnumber), 0.6, 0.74, 0, kBlack, 0.04);
  drawText("N_{tube}^{N,S} >= 2", 0.6, 0.69, 0, kBlack, 0.04);
  drawText("#Sigma E_{ZDC}^{N,S} > 40 GeV", 0.6, 0.64, 0, kBlack, 0.04);
  drawText("MBD Bkg Cut", 0.6, 0.59, 0, kBlack, 0.04);
  drawText("|z_{vtx}| < 60 cm", 0.6, 0.54, 0, kBlack, 0.04);

  cc3->cd(2);
  TH1D *h_02 = (TH1D*) vertex_charge_sum_balanced[8]->Clone();
  TH1D *h_01 = (TH1D*) vertex_charge_sum_balanced[9]->Clone();
  TH1D *h_12 = (TH1D*) vertex_charge_sum_balanced[10]->Clone();

  h_02->Divide(vertex_charge_sum_balanced[7]);
  h_01->Divide(vertex_charge_sum_balanced[7]);
  h_12->Divide(vertex_charge_sum_balanced[7]);

  TH1D *hg = (TH1D*)  h_02->Clone();
  hg->SetTitle(";MBD Charge Sum; Ratio");

  hg->GetXaxis()->SetTitleFont(42);
  hg->GetXaxis()->SetTitleSize(0.06);
  hg->GetYaxis()->SetTitleFont(42);
  hg->GetYaxis()->SetTitleSize(0.06);

  hg->SetMinimum(0);
  hg->SetMaximum(2.0);
  hg->Reset();
  hg->Draw();
  SetLineAtt(h_02, colors[0], 1, 1);
  SetLineAtt(h_01, colors[1], 1, 1);
  SetLineAtt(h_12, colors[2], 1, 1);
  h_02->Draw("same");
  h_01->Draw("same");
  h_12->Draw("same");
  lred = new TLegend(0.2, 0.7, 0.5, 0.95);
  lred->SetLineWidth(0);
  lred->SetTextSize(0.04);
  lred->AddEntry(h_02, "1 < z_{vtx} < 5 / -1 < z_{vtx} < 1");
  lred->AddEntry(h_01, "5 < z_{vtx} < 10 / -1 < z_{vtx} < 1");
  lred->AddEntry(h_12, "10 < z_{vtx} < 15 / -1 < z_{vtx} < 1");

  lred->Draw("same");

  cc3->SaveAs(Form("%s/output_2024/centrality_2024/run%d_charge_sum_balanced_vertex.pdf", env_p, runnumber));

  file->Close();
  delete file;
  
  c->Close();
  cc3->Close();
  cc2->Close();
  cc1->Close();
  c_charge->Close();
  
  delete c;
  delete cc3;
  delete cc2;
  delete cc1;
  delete c_charge;
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


  TFile *file = new TFile(Form("%s/output_2024/plots/mbdana_charge_sum_%d.root", env_p, runnumber), "r");
  if (!file) return;

  TFile *file2 = new TFile(Form("%s/output_2024/plots/mbdana_charge_sum_%d.root", env_p, runnumber2), "r");
  if (!file2) return;

  TH1D *h_charge_sum_mb_vtx = (TH1D*) file->Get("h_charge_sum_min_bias_w_vertex_cut_balanced_scaled");
  if (!h_charge_sum_mb_vtx) return;  

  TH1D *h_charge_sum_mb_vtx2 = (TH1D*) file2->Get("h_charge_sum_min_bias_w_vertex_cut_balanced");
  if (!h_charge_sum_mb_vtx2) return;  

  double sum1 = 0;
  double sum2 = 0;
  for (int i = 1300; i <1600;i++)
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
  l->AddEntry(h_charge_sum_mb_vtx, Form("Run %d", runnumber));
  l->AddEntry(h_charge_sum_mb_vtx2, Form("Run %d - Ref", runnumber2));
  l->Draw("same");

  // p2->cd();

  c->cd(2);
  
  TH1D *hRatio = (TH1D*) h_charge_sum_mb_vtx->Clone();
  TH1D *hRatiod = (TH1D*) h_charge_sum_mb_vtx2->Clone();
  hRatio->Rebin(10);
  hRatiod->Rebin(10);
  
  hRatio->SetName("hRatio");
  hRatio->Divide(hRatiod);
  hRatio->SetTitle(";MBD Charge N+S; Run/Ref");
  hRatio->SetMinimum(0);
  hRatio->SetMaximum(2.0);

  SetLineAtt(hRatio, kBlack, 1, 1);
  SetMarkerAtt(hRatio, kBlack, 1, 24);

  hRatio->Draw("P");
  TLine *line = new TLine(0., 1., 2500., 1.);
  SetLineAtt(line, kRed, 2, 4);
  line->Draw("same");
  c->SaveAs(Form("%s/output_2024/centrality_2024/run%d_run%d_charge_sum.pdf", env_p, runnumber, runnumber2));



}


void DrawMBDCentralityCheck(const int runnumber)
{
  gStyle->SetOptStat(0);
  SetyjPadStyle();



  std::string datestring = "new";
  const char *env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");



  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }
  TString filename = Form("%s/output_2024/plots/mbdana_centrality_check_vtx_%d.root", env_p, runnumber);
  if (isSim)
    {
      filename = Form("%s/output_2024/plots/mbdana_centrality_check_vtx_hijing.root", env_p);
    }
  // TFile *file_vtx = new TFile(filename.Data(), "r");

  // if (!file_vtx)
  //   {
  //     return;
  //   }
  // std::cout << __LINE__ << std::endl;
  // TH1D *h_cent_bin_vtx_shifted;
  // TH1D *h_cent_bin_vtx;
  // std::cout << __LINE__ << std::endl;
  // h_cent_bin_vtx = (TH1D*) file_vtx->Get("hcent_bins");
  // if (!h_cent_bin_vtx) return;
  // std::cout << __LINE__ << std::endl;
  // TH1D *h_vertex_vtx = (TH1D*) file_vtx->Get("h_vertex");
  // if (!h_vertex_vtx) return;
  // std::cout << __LINE__ << std::endl;
  // h_cent_bin_vtx_shifted = (TH1D*) file_vtx->Get("hcent_bins_sca");

  // TH1D *h_cent_vtx_vtx[20];
  // TH1D *h_cent_vertex_vtx[20];
  // TH1D *h_cent_vtx_sca_vtx[20];
  // for (int i = 0; i < 20; i++)
  //   {
  //     h_cent_vertex_vtx[i] = (TH1D*) file_vtx->Get(Form("h_cent_vertex_%d", i));
  //     if (!h_cent_vertex_vtx[i]) return;
  //     h_cent_vtx_vtx[i] = (TH1D*) h_cent_vertex_vtx[i]->Clone();
  //     for (int ib = 0; ib < h_cent_vtx_vtx[i]->GetNbinsX(); ib++)
  // 	{
  // 	  h_cent_vtx_vtx[i]->SetBinContent(ib+1, h_cent_vertex_vtx[i]->GetBinContent(ib+1)/ h_vertex_vtx->GetBinContent(ib+1));
  // 	  h_cent_vtx_vtx[i]->SetBinError(ib+1, h_cent_vtx_vtx[i]->GetBinContent(ib+1)*TMath::Power(h_cent_vertex_vtx[i]->GetBinError(ib+1)/h_cent_vertex_vtx[i]->GetBinContent(ib+1), 2) + TMath::Power( h_vertex_vtx->GetBinError(ib+1)/h_vertex_vtx->GetBinContent(ib+1), 2));
  // 	}
  //     h_cent_vtx_sca_vtx[i] = (TH1D*) file_vtx->Get(Form("h_cent_vtx_sca_%d", i));
  //     if (!h_cent_vtx_sca_vtx[i]) return;
  //   }
  // std::cout << __LINE__ << std::endl;
  TString filename2 = Form("%s/output_2024/plots/mbdana_centrality_check_%d.root", env_p, runnumber);
  if (isSim)
    {
      filename2 = Form("%s/output_2024/plots/mbdana_centrality_check_hijing.root", env_p);
    }
  TFile *filer = new TFile(filename2.Data(), "r");

  if (!filer)
    {
      return;
    }
  std::cout << __LINE__ << std::endl;
  TH1D *h_cent_bin_r_shifted;
  TH1D *h_cent_bin_r;
  h_cent_bin_r = (TH1D*) filer->Get("hcent_bins");
  TH1D *h_vertex_r = (TH1D*) filer->Get("h_vertex");
  if (!h_vertex_r) return;
  std::cout << __LINE__ << std::endl;
  if (!h_cent_bin_r) return;
  std::cout << __LINE__ << std::endl;
  h_cent_bin_r_shifted = (TH1D*) filer->Get("hcent_bins_sca");

  TH1D *h_cent_vtx_r[20];
  TH1D *h_cent_vertex_r[20];
  TH1D *h_cent_vtx_sca_r[20];
  for (int i = 0; i < 20; i++)
    {
      h_cent_vertex_r[i] = (TH1D*) filer->Get(Form("h_cent_vertex_%d", i));
      if (!h_cent_vertex_r[i]) return;
      std::cout << __LINE__ << std::endl;
      h_cent_vtx_r[i] = (TH1D*) filer->Get(Form("h_cent_vtx_%d", i));h_cent_vertex_r[i]->Clone();
      h_cent_vtx_sca_r[i] = (TH1D*) filer->Get(Form("h_cent_vtx_sca_%d", i));
      if (!h_cent_vtx_sca_r[i]) return;
    }

  std::cout << __LINE__ << std::endl;
  TString filename3 = Form("%s/output_2024/plots/mbdana_centrality_check_bal_%d.root", env_p, runnumber);
  if (isSim)
    {
      filename3 = Form("%s/output_2024/plots/mbdana_centrality_check_bal_hijing.root", env_p);
    }
  TFile *file = new TFile(filename3.Data(), "r");

  if (!file)
    {
      return;
    }
  std::cout << __LINE__ << std::endl;
  TCanvas *c = new TCanvas("c","c", 500, 500);
  TH1D *h_cent_bin_shifted;
  TH1D *h_cent_bin;

  TH1D *h_cent_vtx[20];
  TH1D *h_cent_vertex[20];
  TH1D *h_cent_vtx_sca[20];
  TH1D *h_vertex = (TH1D*) filer->Get("h_vertex");
  for (int i = 0; i < 20; i++)
    {
      h_cent_vertex[i] = (TH1D*) file->Get(Form("h_cent_vertex_%d", i));
      if (!h_cent_vertex[i]) return;
      h_cent_vtx[i] = (TH1D*) file->Get(Form("h_cent_vtx_%d", i));
      h_cent_vtx_sca[i] = (TH1D*) file->Get(Form("h_cent_vtx_sca_%d", i));
      if (!h_cent_vtx_sca[i]) return;
    }

  h_cent_bin = (TH1D*) file->Get("hcent_bins");

  if (!h_vertex) return;

  if (!h_cent_bin) return;

  h_cent_bin->Scale(100);
  h_cent_bin_r->Scale(100);
  for (int ib = 1; ib <= h_cent_bin_r->GetNbinsX(); ib++)
    {
      if (h_cent_bin->GetBinError(ib) < 0.01) h_cent_bin->SetBinError(ib, 0.01);
      if (h_cent_bin_r->GetBinError(ib) < 0.01) h_cent_bin_r->SetBinError(ib, 0.01);
    }

  h_cent_bin_shifted = (TH1D*) file->Get("hcent_bins_sca");

  if (!h_cent_bin_shifted) return;
  h_cent_bin_shifted->Scale(100);
  h_cent_bin_r_shifted->Scale(100);

  h_cent_bin_r->SetTitle(";Centrality Bin; Percentage of Events");
  h_cent_bin->SetTitle(";Centrality Bin; Percentage of Events");

  SetMarkerAtt(h_cent_bin_shifted, kViolet - 2, 1, 89);
  SetMarkerAtt(h_cent_bin, kSpring + 2, 1, 90);
  SetMarkerAtt(h_cent_bin_r, kAzure - 6, 1, 90);
  SetLineAtt(h_cent_bin_shifted, kViolet + 2, 1, 1);
  h_cent_bin_r->GetXaxis()->SetTitleFont(42);  
  h_cent_bin_r->GetXaxis()->SetTitleSize(0.04);  
  h_cent_bin_r->GetYaxis()->SetTitleOffset(1.9);  
  h_cent_bin_r->GetYaxis()->SetTitleFont(42);  
  h_cent_bin_r->GetYaxis()->SetTitleSize(0.04);  
  h_cent_bin_r->GetXaxis()->SetLabelSize(0.04);  
  h_cent_bin_r->GetYaxis()->SetLabelSize(0.04);  
  h_cent_bin->GetXaxis()->SetTitleFont(42);  
  h_cent_bin->GetXaxis()->SetTitleSize(0.04);  
  h_cent_bin->GetYaxis()->SetTitleFont(42);  
  h_cent_bin->GetYaxis()->SetTitleOffset(1.9);  
  h_cent_bin->GetYaxis()->SetTitleSize(0.04);  
  h_cent_bin->GetXaxis()->SetLabelSize(0.04);  
  h_cent_bin->GetYaxis()->SetLabelSize(0.04);  

  h_cent_bin->SetMaximum(1.4);
  h_cent_bin->SetMinimum(0.9);
  h_cent_bin->Draw();
  //  h_cent_bin_r->Draw("same");
  //h_cent_bin_vtx->Draw("same");

  TF1 *flatline = new TF1("flatline","[0]",-0.5, 91.5);
  flatline->SetParameter(0,100./(float)  ndivs);
  h_cent_bin->Fit(flatline,"NDORQ","");

  float chi2_r = 0;
  float yy_r = 0;
  float chi2_vtx = 0;
  float yy_vtx = 0;
  float chi2 = flatline->GetChisquare()/flatline->GetNDF();
  float yy = flatline->GetParameter(0);
  

  h_cent_bin_r->Fit(flatline,"QNDOR","");
  
  chi2_r = flatline->GetChisquare()/flatline->GetNDF();;
  yy_r = flatline->GetParameter(0);
    
  DrawSPHENIX(0.19, 0.85, 1, 0, 0.04);
  drawText(Form("Run %d %s", runnumber, datestring.c_str()), 0.19, 0.74, 0, kBlack, 0.04);
  flatline->SetLineColor(kRed);
  flatline->SetLineWidth(2);
  flatline->SetParameter(0, 100./(float)ndivs);
  flatline->Draw("same");
  TLegend *lg = new TLegend(0.6, 0.65, 0.82, 0.89);
  SetLegendStyle(lg);
  lg->SetHeader("Flatline #chi^2/NDF");
  lg->AddEntry(h_cent_bin_r,Form("Uncorr. : %2.3f", chi2_r));
  lg->AddEntry(h_cent_bin,Form("Vtx Scaled: %2.3f", chi2)); 
  lg->Draw("same");
  c->SaveAs(Form("%s/output_2024/centrality_2024/run%d_centrality.pdf", env_p, runnumber));

  // h_cent_bin_r->SetTitle(";Centrality Bin; Fraction of Events");
  // SetMarkerAtt(h_cent_bin_r_shifted, kViolet - 2, 1, 89);
  // SetMarkerAtt(h_cent_bin_r, kSpring + 2, 1, 90);
  // SetLineAtt(h_cent_bin_r_shifted, kViolet + 2, 1, 1);
  // SetLineAtt(h_cent_bin_r, kSpring-1, 1, 1);
  
  // h_cent_bin_r->SetMaximum(1.4);
  // h_cent_bin_r->SetMinimum(0.9);
  // h_cent_bin_r->Draw();
  // h_cent_bin_r_shifted->Draw("same");

  // flatline->SetParameter(0,100./(float) ndivs);

  // h_cent_bin_r->Fit(flatline,"NDORQ","");

  // float chi2_shifted = 0;
  // float yy_shifted = 0;
  // chi2 = flatline->GetChisquare()/flatline->GetNDF();
  // yy = flatline->GetParameter(0);
  

  // h_cent_bin_r_shifted->Fit(flatline,"QNDOR","");
  
  // chi2_shifted = flatline->GetChisquare()/flatline->GetNDF();;
  // yy_shifted = flatline->GetParameter(0);
  
  
  // DrawSPHENIX(0.19, 0.85, 1, 0, 0.03);
  // drawText(Form("Run %d", runnumber), 0.19, 0.76, 0, kBlack, 0.03);
  // flatline->SetLineColor(kRed);
  // flatline->SetLineWidth(2);
  // flatline->Draw("same");
  // lg = new TLegend(0.7, 0.65, 0.88, 0.89);
  // SetLegendStyle(lg);
  // lg->SetHeader("Flatline #chi^2/NDF");
  // lg->AddEntry(h_cent_bin_r,Form("Not Scaled: %2.3f", chi2));
  // lg->AddEntry(h_cent_bin_r_shifted,Form("Scaled: %2.3f", chi2_shifted)); 
  // lg->Draw("same");
  // c->SaveAs(Form("%s/output_2024/centrality_2024/run%d_centrality_r.pdf", env_p, runnumber));

  TCanvas *c1 = new TCanvas("c1","c1", 500, 500);
  int colors[4] = {kBlack, kRed - 7, kAzure + 6, kViolet - 4};
  int classes[4] = {0, 1, 9, 17};
  std::cout << __LINE__ << std::endl;
  for (int i = 0; i < 4; i++)
    {
      SetLineAtt(h_cent_vtx[classes[i]], colors[i], 1, 1);
      SetMarkerAtt(h_cent_vtx[classes[i]], colors[i], 0.6, 8);
      SetLineAtt(h_cent_vtx_r[classes[i]], colors[i], 1, 4);
      SetMarkerAtt(h_cent_vtx_r[classes[i]], colors[i], 0.8, 8);

    }
  std::cout << __LINE__ << std::endl;
  h_cent_vtx[classes[0]]->SetMinimum(0.00);
  h_cent_vtx[classes[0]]->SetMaximum(0.18);
  h_cent_vtx[classes[0]]->Draw();
  h_cent_vtx[classes[1]]->Draw("same");
  h_cent_vtx[classes[2]]->Draw("same"); 
  h_cent_vtx[classes[3]]->Draw("same");
  std::cout << __LINE__ << std::endl;  
  DrawSPHENIX(0.19, 0.85, 1, 0);
  drawText(Form("Run %d", runnumber), 0.19, 0.73, 0, kBlack);
  TLegend *l2 = new TLegend(0.65, 0.69, 0.86, 0.88);
  l2->SetLineWidth(0);
  l2->SetTextSize(0.03);
  l2->SetHeader("Vertex Scale Factor");
  l2->AddEntry(h_cent_vtx[classes[0]], "0-5%");
  l2->AddEntry(h_cent_vtx[classes[1]], "5-10%");
  l2->AddEntry(h_cent_vtx[classes[2]], "45-50%");
  l2->AddEntry(h_cent_vtx[classes[3]], "85-90%");
  l2->Draw("same");

  c1->SaveAs(Form("%s/output_2024/centrality_2024/run%d_vtx_centrality.pdf", env_p, runnumber));

  std::cout << __LINE__ << std::endl;
  h_cent_vtx_r[classes[0]]->SetMinimum(0.00);
  h_cent_vtx_r[classes[0]]->SetMaximum(0.18);
  h_cent_vtx_r[classes[0]]->Draw();
  h_cent_vtx_r[classes[1]]->Draw("same");
  h_cent_vtx_r[classes[2]]->Draw("same"); 
  h_cent_vtx_r[classes[3]]->Draw("same");
  DrawSPHENIX(0.19, 0.85, 1, 0);
  drawText(Form("Run %d", runnumber), 0.19, 0.73, 0, kBlack);
  l2 = new TLegend(0.65, 0.69, 0.86, 0.88);
  l2->SetLineWidth(0);
  l2->SetTextSize(0.03);
  l2->SetHeader("No Vertex Scale Factor");
  l2->AddEntry(h_cent_vtx[classes[0]], "0-5%");
  l2->AddEntry(h_cent_vtx[classes[1]], "5-10%");
  l2->AddEntry(h_cent_vtx[classes[2]], "45-50%");
  l2->AddEntry(h_cent_vtx[classes[3]], "85-90%");
  l2->Draw("same");
  
  c1->SaveAs(Form("%s/output_2024/centrality_2024/run%d_vtx_r_centrality.pdf", env_p, runnumber));
  std::cout << __LINE__ << std::endl;

  for (int i = 0; i < 4; i++)
    {
      SetLineAtt(h_cent_vertex[classes[i]], colors[i], 1, 1);
      SetMarkerAtt(h_cent_vertex[classes[i]], colors[i], 0.6, 1);
    }
  h_cent_vertex[classes[0]]->SetMaximum(  h_cent_vertex[classes[0]]->GetBinContent(  h_cent_vertex[classes[0]]->GetMaximumBin())*1.3);
  h_cent_vertex[classes[0]]->Draw();
  h_cent_vertex[classes[1]]->Draw("same");
  h_cent_vertex[classes[2]]->Draw("same"); 
  h_cent_vertex[classes[3]]->Draw("same");
  
  DrawSPHENIX(0.19, 0.85, 1, 0);
  drawText(Form("Run %d", runnumber), 0.19, 0.73, 0, kBlack);
  l2 = new TLegend(0.65, 0.69, 0.86, 0.88);
  l2->SetLineWidth(0);
  l2->SetTextSize(0.03);
  l2->SetHeader("Centrality Range");
  l2->AddEntry(h_cent_vertex[classes[0]], "0-5%");
  l2->AddEntry(h_cent_vertex[classes[1]], "5-10%");
  l2->AddEntry(h_cent_vertex[classes[2]], "45-50%");
  l2->AddEntry(h_cent_vertex[classes[3]], "85-90%");
  l2->Draw("same");
  c1->SaveAs(Form("%s/output_2024/centrality_2024/run%d_vertex_centrality.pdf", env_p, runnumber));

  for (int i = 0; i < 4; i++)
    {
      SetLineAtt(h_cent_vertex_r[classes[i]], colors[i], 1, 1);
      SetMarkerAtt(h_cent_vertex_r[classes[i]], colors[i], 0.6, 1);
    }
  h_cent_vertex_r[classes[0]]->SetMaximum(  h_cent_vertex_r[classes[0]]->GetBinContent(  h_cent_vertex_r[classes[0]]->GetMaximumBin())*1.3);
  h_cent_vertex_r[classes[0]]->Draw();
  h_cent_vertex_r[classes[1]]->Draw("same");
  h_cent_vertex_r[classes[2]]->Draw("same"); 
  h_cent_vertex_r[classes[3]]->Draw("same");
  
  DrawSPHENIX(0.19, 0.85, 1, 0);
  drawText(Form("Run %d", runnumber), 0.19, 0.73, 0, kBlack);
  l2 = new TLegend(0.65, 0.69, 0.86, 0.88);
  l2->SetLineWidth(0);
  l2->SetTextSize(0.03);
  l2->SetHeader("Centrality Range");
  l2->AddEntry(h_cent_vertex_r[classes[0]], "0-5%");
  l2->AddEntry(h_cent_vertex_r[classes[1]], "5-10%");
  l2->AddEntry(h_cent_vertex_r[classes[2]], "45-50%");
  l2->AddEntry(h_cent_vertex_r[classes[3]], "85-90%");
  l2->Draw("same");
  c1->SaveAs(Form("%s/output_2024/centrality_2024/run%d_vertex_r_centrality.pdf", env_p, runnumber));


  c1->Close();
  delete c1;

  c->Close();
  delete c;

  filer->Close();
  file->Close();
}
void DrawMCMBDCentralityCheck(const int runnumber)
{
  gStyle->SetOptStat(0);
  SetyjPadStyle();



  std::string datestring = "new";
  const char *env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");


  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }
  TString filename2 = Form("%s/output_2024/plots/mbdana_centrality_check_%d.root", env_p, runnumber);
  if (isSim)
    {
      filename2 = Form("%s/output_2024/plots/mbdana_centrality_check_hijing.root", env_p);
    }
  TFile *filer = new TFile(filename2.Data(), "r");

  if (!filer)
    {
      return;
    }
  TH1D *h_cent_bin_r_shifted;
  TH1D *h_cent_bin_r;
  h_cent_bin_r = (TH1D*) filer->Get("hcent_bins");

  if (!h_cent_bin_r) return;

  h_cent_bin_r_shifted = (TH1D*) filer->Get("hcent_bins_sca");

  TH1D *h_cent_vtx_r[20];
  TH1D *h_cent_vertex_r[20];
  TH1D *h_cent_vtx_sca_r[20];
  for (int i = 0; i < 20; i++)
    {
      h_cent_vertex_r[i] = (TH1D*) filer->Get(Form("h_cent_vertex_%d", i));
      if (!h_cent_vertex_r[i]) return;
      h_cent_vtx_r[i] = (TH1D*) filer->Get(Form("h_cent_vtx_%d", i));
      if (!h_cent_vtx_r[i]) return;
      h_cent_vtx_sca_r[i] = (TH1D*) filer->Get(Form("h_cent_vtx_sca_%d", i));
      if (!h_cent_vtx_sca_r[i]) return;
    }

  TCanvas *c = new TCanvas("c","c", 500, 500);
  h_cent_bin_r->Scale(100);
  for (int ib = 1; ib <= h_cent_bin_r->GetNbinsX(); ib++)
    {
      if (h_cent_bin_r->GetBinError(ib) < 0.01) h_cent_bin_r->SetBinError(ib, 0.01);
    }
  h_cent_bin_r->SetTitle(";Centrality Bin; Percentage of Events");

  SetMarkerAtt(h_cent_bin_r, kAzure - 6, 1, 90);

  h_cent_bin_r->GetXaxis()->SetTitleFont(42);  
  h_cent_bin_r->GetXaxis()->SetTitleSize(0.04);  
  h_cent_bin_r->GetYaxis()->SetTitleOffset(1.9);  
  h_cent_bin_r->GetYaxis()->SetTitleFont(42);  
  h_cent_bin_r->GetYaxis()->SetTitleSize(0.04);  
  h_cent_bin_r->GetXaxis()->SetLabelSize(0.04);  
  h_cent_bin_r->GetYaxis()->SetLabelSize(0.04);  
  h_cent_bin_r->Draw();

  TF1 *flatline = new TF1("flatline","[0]",-0.5, 91.5);
  flatline->SetParameter(0,100./(float)  ndivs);

  float chi2_r = 0;
  float yy_r = 0;
  float chi2_vtx = 0;
  float yy_vtx = 0;
  float chi2 = 0;
  float yy = 0;
  
  h_cent_bin_r->Fit(flatline,"QNDOR","");
  
  chi2_r = flatline->GetChisquare()/flatline->GetNDF();;
  yy_r = flatline->GetParameter(0);
  
  DrawSPHENIX(0.19, 0.85, 1, 0, 0.04);
  drawText(Form("Run %d %s", runnumber, datestring.c_str()), 0.19, 0.74, 0, kBlack, 0.04);
  flatline->SetLineColor(kRed);
  flatline->SetLineWidth(2);
  flatline->SetParameter(0, 100./(float)94.);
  flatline->Draw("same");

  TLegend *lg = new TLegend(0.6, 0.65, 0.82, 0.89);
  SetLegendStyle(lg);
  lg->SetHeader("Flatline #chi^2/NDF");
  lg->AddEntry(h_cent_bin_r,Form("Uncorr. : %2.3f", chi2_r));
  lg->Draw("same");
  c->SaveAs(Form("%s/output_2024/centrality_2024/run%d_centrality.pdf", env_p, runnumber));

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


  TFile *file = new TFile(Form("%s/output_2024/plots/mbdana_zdc_check_%d.root", env_p, runnumber), "r");

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
  TCanvas *c = new TCanvas("c","c", 500, 500);
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

  c->SaveAs(Form("%s/output_2024/centrality_2024/run%d_zdc.png", env_p, runnumber));
  c->SaveAs(Form("%s/output_2024/centrality_2024/run%d_zdc.pdf", env_p, runnumber));

  c->Clear();
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

  c->SaveAs(Form("%s/output_2024/centrality_2024/run%d_zdc_singles.png", env_p, runnumber));
  c->SaveAs(Form("%s/output_2024/centrality_2024/run%d_zdc_singles.pdf", env_p, runnumber));


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
  std::string filepath = Form("%s/output_2024/plots/mbdana_channels_%d.root", env_p, runnumber);
  if (runnumber == 0)
    {
      filepath = Form("%s/output_2024/plots/mbdana_channels_hijing.root", env_p);
    }
  else if (runnumber == 1)
    {
      filepath = Form("%s/output_2024/plots/mbdana_channels_ampt.root", env_p);
    }
  else if (runnumber == 2)
    {
      filepath = Form("%s/output_2024/plots/mbdana_channels_epos.root", env_p);
    }
  TFile *fout = new TFile(filepath.c_str(), "r");
  TH1D *hpeaks_fit = new TH1D("hpeaks_fit",";MPV_{calib};Tubes", 41, 0.795, 1.205);

  TH1D *h_charge[128];
  TH1D *h_time[128];
  TH1D *h_charge_mb[128];
  TH1D *h_time_mb[128];

  TProfile *hp_charge_vertex[128];
  TProfile *hp_time_vertex[128];
  for (int i = 0; i < 128; i++)
    {
      h_charge[i] = (TH1D*) fout->Get(Form("h_mbd_charge_ch%d", i));
      h_time[i] = (TH1D*) fout->Get(Form("h_mbd_time_ch%d", i));
      h_charge_mb[i] = (TH1D*) fout->Get(Form("h_mbd_charge_min_bias_ch%d", i));
      h_time_mb[i] = (TH1D*) fout->Get(Form("h_mbd_time_min_bias_ch%d", i));
      hp_charge_vertex[i] = (TProfile*) fout->Get(Form("hp_charge_vertex_ch%d", i));
      //      hp_time_vertex[i] = (TProfile*) fout->Get(Form("hp_time_vertex_ch%d", i));
      if (!h_charge_mb[i]) return;
      if (!h_time_mb[i]) return;
      if (!h_charge[i]) return;
      if (!h_time[i]) return;

      if (!hp_charge_vertex[i]) return;
      //if (!hp_time_vertex[i]) return;

    }

  for (int i = 0; i < 128; i++)
    {
      SetLineAtt(hp_charge_vertex[i], ringcolors[mbd_ring_index[i%64]], 1,1 );
      SetMarkerAtt(hp_charge_vertex[i], ringcolors[mbd_ring_index[i%64]], 1, 8 );
      //      SetLineAtt(hp_time_vertex[i], ringcolors[mbd_ring_index[i]], 1,1 );
      //SetMarkerAtt(hp_time_vertex[i], ringcolors[mbd_ring_index[i]], 1, 8 );
    }

  // TCanvas *ctv = new TCanvas("ctv","ctv", 1000, 700);
  // ctv->Divide(2, 1);
  // ctv->cd(1);
  // hp_time_vertex[0]->SetMinimum(-10);
  // hp_time_vertex[0]->SetMaximum(10);
  // hp_time_vertex[0]->Draw();
  // for (int i = 1; i < 64; i++)
  //   {
  //     hp_time_vertex[i]->Draw("same");
  //   }
  // ctv->cd(2);
  // hp_time_vertex[64]->SetMinimum(-10);
  // hp_time_vertex[64]->SetMaximum(10);
  // hp_time_vertex[64]->Draw();
  // for (int i = 1; i < 64; i++)
  //   {
  //     hp_time_vertex[i+64]->Draw("same");
  //   }

  // ctv->Print(Form("%s/output_2024/centrality_2024/run%d_time_vtx.pdf", env_p, runnumber));
  // ctv->SaveAs(Form("%s/output_2024/centrality_2024/run%d_time_vtx.pdf", env_p, runnumber));

  TCanvas *ccv = new TCanvas("ccv","ccv", 1000, 700);
  ccv->Divide(2, 1);
  ccv->cd(1);
  hp_charge_vertex[0]->SetMinimum(0);
  hp_charge_vertex[0]->SetMaximum(14);
  hp_charge_vertex[0]->SetTitle(";vertex [cm]; <q_{i}>");
  hp_charge_vertex[64]->SetTitle(";vertex [cm]; <q_{i}>");
  hp_charge_vertex[0]->Draw();
  DrawSPHENIX(0.2, 0.86, 1, 0, 0.04);
  drawText(Form("run %d", runnumber), 0.2,0.76, 0, kBlack, 0.04);
  drawText("South", 0.8,0.86, 1, kBlack, 0.04);
  for (int i = 1; i < 64; i++)
    {
      hp_charge_vertex[i]->Draw("same");
    }
  ccv->cd(2);
  hp_charge_vertex[64]->SetMinimum(0);
  hp_charge_vertex[64]->SetMaximum(14);
  hp_charge_vertex[64]->Draw();
  drawText("North", 0.8,0.86, 1, kBlack, 0.04);
  for (int i = 1; i < 64; i++)
    {
      hp_charge_vertex[i+64]->Draw("same");
    }
  TLegend *lgr = new TLegend(0.2, 0.75, 0.4, 0.88);
  lgr->SetLineWidth(0);
  lgr->AddEntry(hp_charge_vertex[7], "Inner Ring","p");
  lgr->AddEntry(hp_charge_vertex[6], "Middle Ring","p");
  lgr->AddEntry(hp_charge_vertex[0], "Outer Ring","p");
  lgr->Draw("same");
  ccv->Print(Form("%s/output_2024/centrality_2024/run%d_charge_vtx.pdf", env_p, runnumber));
  ccv->SaveAs(Form("%s/output_2024/centrality_2024/run%d_charge_vtx.png", env_p, runnumber));

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
      h_charge[i]->Fit("lan_w_gausexp","0RQ","", .2, 3);
      if (! h_charge[i]->GetFunction("lan_w_gausexp")) continue;
      h_charge[i]->GetFunction("lan_w_gausexp")->SetLineColor(kBlack);
      h_charge[i]->GetXaxis()->SetRangeUser(0.1, 3);
      float peakfit =h_charge[i]->GetBinContent(h_charge[i]->GetMaximumBin());
      h_charge[i]->SetMaximum(peakfit * 1.3);

      //      h_charge[i]->GetXaxis()->SetRangeUser(.1, 3);
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

  c->Print(Form("%s/output_2024/centrality_2024/run%d_charge.pdf(", env_p, runnumber));
  c->SaveAs(Form("%s/output_2024/centrality_2024/run%d_charge_south.png", env_p, runnumber));

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
      h_charge[i+64]->Fit("lan_w_gausexp","0RQ","", .2, 3);
      if (!h_charge[i+64]->GetFunction("lan_w_gausexp")) continue;
      h_charge[i+64]->GetFunction("lan_w_gausexp")->SetLineColor(kBlack);

      h_charge[i+64]->GetXaxis()->SetRangeUser(0.1, 3);
      float peakfit =h_charge[i+64]->GetBinContent(h_charge[i+64]->GetMaximumBin());
      h_charge[i+64]->SetMaximum(peakfit * 1.3);

      //      h_charge[i+64]->GetXaxis()->SetRangeUser(.1, 3);
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

  c->Print(Form("%s/output_2024/centrality_2024/run%d_charge.pdf)", env_p, runnumber));
  c->SaveAs(Form("%s/output_2024/centrality_2024/run%d_charge_north.png", env_p, runnumber));

  TCanvas *cc = new TCanvas("cc","cc", 500, 500);
  SetyjPadStyle();
  hpeaks_fit->Draw();
  cc->SetTicks(1,1);
  DrawSPHENIX(0.5, 0.8, 1, 0, 0.04);
  drawText(Form("Run %d", runnumber), 0.5, 0.7, 0, kBlack, 0.04);
  drawText("Official MBD Calibrations", 0.5, 0.65, 0, kBlack, 0.04);
  cc->SaveAs(Form("%s/output_2024/centrality_2024/run%d_charge_fits.png", env_p, runnumber));

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
      h_time_mb[i]->GetXaxis()->SetRangeUser(-25, 25);
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

  c->Print(Form("%s/output_2024/centrality_2024/run%d_time.pdf(", env_p, runnumber));
  c->SaveAs(Form("%s/output_2024/centrality_2024/run%d_time_south.png", env_p, runnumber));

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
      h_time_mb[i+64]->GetXaxis()->SetRangeUser(-25, 25);
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
  c->Print(Form("%s/output_2024/centrality_2024/run%d_time.pdf)", env_p, runnumber));
  c->SaveAs(Form("%s/output_2024/centrality_2024/run%d_time_north.png", env_p, runnumber));


  c = new TCanvas("cc1","cc1", 500, 500);
  h_charge[0]->SetTitle(";Charge;Counts");
  float maximum =   1.3*h_charge[0]->GetBinContent(  h_charge[0]->GetMaximumBin());
  h_charge[0]->SetMaximum(maximum);
  h_charge[0]->Draw();
  DrawSPHENIX(0.5, 0.85, 1, 0, 0.04);
  drawText(Form("Run - %d", runnumber), 0.5, 0.75, 0, kBlack, 0.04);
  drawText("South - Tube 0", 0.5, 0.7, 0, kBlack, 0.04);
  TLine *l = new TLine(0.5, 0, 0.5, maximum);
  l->SetLineColor(kBlack);
  l->SetLineWidth(2);
  l->Draw("same");

  c->SaveAs(Form("%s/output_2024/centrality_2024/run%d_one_charge.png", env_p, runnumber));
  h_time_mb[0]->SetTitle(";Time [ns];Counts");
  maximum =   1.3*h_time_mb[0]->GetBinContent(  h_time_mb[0]->GetMaximumBin());
  h_time_mb[0]->SetMaximum(maximum);

  h_time_mb[0]->Draw();
  DrawSPHENIX(0.2, 0.85, 1, 0, 0.04);
  drawText(Form("Run - %d", runnumber), 0.2, 0.75, 0, kBlack, 0.04);
  drawText("South - Tube 0", 0.2, 0.70, 0, kBlack, 0.04);

  TLine *l1 = new TLine(-20, 0, -20, maximum);
  l1->SetLineColor(kBlack);
  l1->SetLineWidth(2);
  l1->Draw("same");
  TLine *l2 = new TLine(20, 0, 20, maximum);
  l2->SetLineColor(kBlack);
  l2->SetLineWidth(2);
  l2->Draw("same");
  c->SaveAs(Form("%s/output_2024/centrality_2024/run%d_one_time.pdf", env_p, runnumber));


}
void DrawRunCentralityChecks(const std::string runlist)
{

  const char *env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }
  std::ifstream file("../rundb/time_table.csv");
  std::string runline;
  std::map<int, std::string> runnumbertimes;

  if (file.is_open()) {
    while (std::getline(file, runline)) {
      std::stringstream ss(runline);
      std::string token;
      std::vector<std::string> values;

      while (std::getline(ss, token, ',')) {
        values.push_back(token);
	std::cout << token << std::endl;
      }
      runnumbertimes[std::stoi(values[0])] = std::string(values[1]);
    }
    file.close();
  }
  else {
    std::cerr << "Error opening file." << std::endl;
  }

  std::ifstream inputFile(runlist.c_str());
  std::map<int, TH1D*> runhistmap;
  std::map<int, TGraphErrors*> runvertexmap;
  std::map<int, TGraphErrors*> runMBmap;
  std::map<int, TGraphErrors*> runZDCmap;
  if (inputFile.is_open()) { // Check if the file opened successfully
    string line;

    while (getline(inputFile, line)) { // Read each line from the file

      int runnumber = std::stoi(line);
      TFile *file_bal = new TFile(Form("%s/output_2024/plots/mbdana_centrality_check_bal_%d.root", env_p, runnumber), "r");

      if (!file_bal)
	{
	  return;
	}

        TH1D *h_cent_bin  = (TH1D*) file_bal->Get("hcent_bins");
	h_cent_bin->Scale(100);
	h_cent_bin->Rebin(5);
	h_cent_bin->SetName(Form("h_cent_bin_%d", runnumber));
	for (int ib = 1; ib <= h_cent_bin->GetNbinsX(); ib++)
	  {
	    if (h_cent_bin->GetBinError(ib) < 0.5) h_cent_bin->SetBinError(ib, 0.05);
	  }
	runhistmap[runnumber] = h_cent_bin;


	TFile *file_zdc = new TFile(Form("%s/output_2024/plots/mbdana_zdc_check_%d.root", env_p, runnumber), "r");

	if (!file_zdc)
	  {
	    return;
	  }

        TGraphErrors *gv  = (TGraphErrors*) file_zdc->Get("g_running_average_vertex");
	gv->SetName(Form("g_running_vertex_%d", runnumber));

	runvertexmap[runnumber] = gv;
        TGraphErrors *gmb  = (TGraphErrors*) file_zdc->Get("g_running_average_MB");
	gmb->SetName(Form("g_running_MB_%d", runnumber));

	runMBmap[runnumber] = gmb;
        TGraphErrors *gz  = (TGraphErrors*) file_zdc->Get("g_running_average_ZDC");
	gz->SetName(Form("g_running_ZDC_%d", runnumber));

	runZDCmap[runnumber] = gz;



    }

    inputFile.close(); // Close the file when done
  } else {
    cout << "Error opening file." << endl;
  }


  const int nruns = runhistmap.size();
  TGraphErrors *g_run_centbin[20];
  TGraphErrors *g_run_centbin_date[20];
  for (int ic = 0 ; ic < 20; ic++)
    {
      g_run_centbin[ic] = new TGraphErrors(nruns);
      g_run_centbin_date[ic] = new TGraphErrors(nruns);

    }
  int ir = 0;
  for (auto runhist : runhistmap)
    {
      int runnumber = runhist.first;
      TH1D *hcentbin = runhist.second;

      TDatime da(runnumbertimes[runnumber].c_str());
      for (int ic = 0; ic < 20; ic++)
	{
	  g_run_centbin[ic]->SetPoint(ir, runnumber, hcentbin->GetBinContent(ic+1));
	  g_run_centbin[ic]->SetPointError(ir, 0, hcentbin->GetBinError(ic+1));
	  
	  g_run_centbin_date[ic]->SetPoint(ir, da.Convert(), hcentbin->GetBinContent(ic+1));
	  g_run_centbin_date[ic]->SetPointError(ir, 0, hcentbin->GetBinError(ic+1));

	}
      ir++;
    }
  for (int ic = 0 ; ic < 20; ic++)
    {
      SetLineAtt(g_run_centbin_date[ic], kBlack, 1, 1);

      g_run_centbin_date[ic]->SetMarkerColor(kRed);
      g_run_centbin_date[ic]->SetMarkerStyle(24);
      g_run_centbin_date[ic]->SetTitle(";Run Start Time; Percent in bin");
      SetAxisSizes(g_run_centbin_date[ic], 42, 0.08);

      g_run_centbin_date[ic]->GetXaxis()->SetTimeDisplay(1);
      g_run_centbin_date[ic]->GetXaxis()->SetNdivisions(405, kFALSE);
      g_run_centbin_date[ic]->GetXaxis()->SetTimeFormat("%h-%d %H");
      g_run_centbin_date[ic]->GetXaxis()->SetTimeOffset(0,"gmt-7");

      SetLineAtt(g_run_centbin[ic], kBlack, 1, 1);

      g_run_centbin[ic]->SetMarkerColor(kRed);
      g_run_centbin[ic]->SetMarkerStyle(24);
      g_run_centbin[ic]->SetTitle(";Run number; Percent in bin");
      SetAxisSizes(g_run_centbin[ic], 42, 0.08);

    }

  SetyjPadStyle();
  gStyle->SetOptStat(0);
  TLine *l  = new TLine(54530, 500./(float) ndivs,54974, 500./(float)ndivs);
  l->SetLineColor(kGreen);
  TF1 *flatline = new TF1("flatline","[0]");
  TCanvas  *c = new TCanvas("c","c", 1000, 2000);
  c->Divide(2,10);
  for (int ic = 0; ic < 19; ic++)
    {
      c->cd(ic+2);
      g_run_centbin[ic]->Draw("AP");
      l->Draw("same");
      g_run_centbin[ic]->Fit("flatline","NDOQ");
      float chi2 = flatline->GetChisquare()/flatline->GetNDF();

      drawText(Form("%d - %d %%", ic*5, (ic+1)*5), 0.2, 0.8, 0, kBlack, 0.05);
      drawText(Form("#Chi^{2} / NDF = %2.2f", chi2), 0.8, 0.8, 1, kBlack, 0.05);
    }
  c->cd(1);
  DrawSPHENIX(0.2, 0.8, 1, 0, 0.1);

  l = new TLine(g_run_centbin_date[0]->GetPointX(0), 500./(float)ndivs, g_run_centbin_date[0]->GetPointX(nruns - 1), 500./(float)ndivs);
  l->SetLineColor(kGreen);
  TCanvas  *c1 = new TCanvas("c1","c1", 1000, 2000);
  c1->Divide(2,10);
  for (int ic = 0; ic < 19; ic++)
    {
      c1->cd(ic+2);
      g_run_centbin_date[ic]->Draw("AP");
      if (ic == 18)
	{  
	  l = new TLine(g_run_centbin_date[0]->GetPointX(0), 400./(float)ndivs, g_run_centbin_date[0]->GetPointX(nruns - 1), 400./(float)ndivs);
	  l->SetLineColor(kGreen);
	}
      l->Draw("same");
      g_run_centbin_date[ic]->Fit("flatline","NDOQ");
      float chi2 = flatline->GetChisquare()/flatline->GetNDF();
      drawText(Form("%d - %d %%", ic*5, (ic+1)*5), 0.2, 0.8, 0, kBlack, 0.05);
      drawText(Form("#Chi^{2} / NDF = %2.2f", chi2), 0.8, 0.8, 1, kBlack, 0.05);
    }
  c1->cd(1);
  DrawSPHENIX(0.2, 0.8, 1, 0, 0.1);

  c->SaveAs(Form("%s/output_2024/centralitycheck_run24.pdf", env_p));
  c1->SaveAs(Form("%s/output_2024/centralitycheckdate_run24.pdf", env_p));

  for (auto gp : runvertexmap)
    {

      for (int ic = 0; ic < gp.second->GetN(); ic++)
	{	  
	  float err = gp.second->GetErrorY(ic);
	  if (err > 0)
	    {	  
	      gp.second->SetPointError(ic, 0, err);
	    }
	}

      gp.second->SetTitle(";Event / 10000; MBD Z Vertex");
      gp.second->SetMarkerColor(kRed);
      gp.second->SetLineColor(kBlack);
      gp.second->SetMarkerStyle(20);
      gp.second->SetMarkerSize(0.8);
    }

  for (auto gp : runMBmap)
    {

      for (int ic = 0; ic < gp.second->GetN(); ic++)
	{	  
	  float err = gp.second->GetErrorY(ic);
	  if (err > 0)
	    {	  
	      gp.second->SetPointError(ic, 0, err);
	    }
	}

      gp.second->SetTitle(";Event / 10000; Minbias Fraction");
      gp.second->SetMarkerColor(kSpring);
      gp.second->SetLineColor(kBlack);
      gp.second->SetMarkerStyle(20);
      gp.second->SetMarkerSize(0.8);
    }
  for (auto gp : runZDCmap)
    {
      for (int ic = 0; ic < gp.second->GetN(); ic++)
	{	  
	  float err = gp.second->GetErrorY(ic);
	  if (err > 0)
	    {	  
	      gp.second->SetPointError(ic, 0, err);
	    }
	}

      gp.second->SetTitle(";Event / 10000; ZDC Fraction");
      gp.second->SetMarkerColor(kCyan);
      gp.second->SetLineColor(kBlack);
      gp.second->SetMarkerStyle(20);
      gp.second->SetMarkerSize(0.8);
    }

  TCanvas *c3 = new TCanvas("c3","c3", 300, 600);
  c3->Divide(1, 3);
  c3->Print(Form("%s/output_2024/rundrama.pdf[", env_p),"pdf");
  for (auto gp : runvertexmap)
    {
      int runnumber = gp.first;

      c3->cd(1);
      runvertexmap[runnumber]->Draw("AP");
      drawText(Form("Run %d", runnumber), 0.2, 0.8, 0, kBlack, 0.1);
      c3->cd(2);
      runMBmap[runnumber]->Draw("AP");
      c3->cd(3);
      runZDCmap[runnumber]->Draw("AP");
      c3->Print(Form("%s/output_2024/rundrama.pdf",env_p), "pdf");
    }
  c3->Print(Form("%s/output_2024/rundrama.pdf]", env_p), "pdf");
}
