#include "dlUtility.h"

void DrawMBDTimeChannels(const int runnumber)
{
  gStyle->SetOptStat(0);
  SetyjPadStyle();

  TFile *file = new TFile(Form("../plots/histout_%d.root", runnumber), "r");

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

  c->Print(Form("time_%d.pdf(", runnumber));

  c->SaveAs(Form("time_%d_%s.png", runnumber,"South"));

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

  c->SaveAs(Form("time_%d_%s.png", runnumber,"North"));


  c->Print(Form("time_%d.pdf)", runnumber));


  ///

}

void DrawMBDTOAChannels(const int runnumber)
{
  gStyle->SetOptStat(0);
  SetyjPadStyle();

  TFile *file = new TFile(Form("../plots/histout_%d.root", runnumber), "r");

  TH1D *h_time[128];
  TH1D *h_time_raw[128];
  TH1D *h_tdc_raw[128];

  for (int i = 0; i < 128; i++)
    {
      h_tdc_raw[i] = (TH1D*) file->Get(Form("h_tdc_channel_%d", i));
      h_time_raw[i] = (TH1D*) file->Get(Form("h_toa_channel_%d", i));
      h_time_raw[i]->SetTitle(";Time (ns); Arb. Units (Normalized by peak height)");
      h_time_raw[i]->Rebin(20);
      h_tdc_raw[i]->Rebin(20);
      float center = h_time_raw[i]->GetBinCenter(h_time_raw[i]->GetMaximumBin());
      h_time_raw[i]->GetXaxis()->SetRangeUser(center - 10, center + 10);
      h_tdc_raw[i]->GetXaxis()->SetRangeUser(center - 10, center + 10);
      float max = h_tdc_raw[i]->GetBinContent(h_tdc_raw[i]->GetMaximumBin());
      h_tdc_raw[i]->Scale(1./max);

    }

  // ALL 128 Cahnnels;

  TCanvas *c = new TCanvas("c","c", 1000, 1200);

  c->Print(Form("time_%d.pdf[", runnumber));
  int colorline = kBlue + 1;
  int colorfill = kBlue - 2;

  TPad *pads[64];
  float xi = 0.05;
  float yi = .9 - xi;
  float dpx = (1. - 2*xi)/8.;
  float dpy = (.9 - 2*xi)/8.;
  DrawSPHENIX(0.05, 0.95, 1, 1, 0.03, 0);
  drawText(Form("Run - %d", runnumber), 0.05, 0.91, 0, kBlack, 0.03);
  drawText("South MBD Charge Arrival Distributions", 0.95, 0.91, 1, kBlack, 0.03);
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
      SetLineAtt(h_time_raw[i], colorline, 2, 1);
      h_time_raw[i]->SetFillColorAlpha(colorfill, 0.3);
      h_time_raw[i]->Draw();
      
      drawText(Form("S Ch. %d", i),0.1, 0.88, 0, kBlack, 0.07);
    }

  c->Print(Form("time_%d.pdf(", runnumber));

  c->SaveAs(Form("time_%d_%s.png", runnumber,"South"));

  c = new TCanvas("c1","c1", 1000, 1200);
  DrawSPHENIX(0.05, 0.95, 1, 1, 0.03, 0);
  drawText(Form("Run - %d", runnumber), 0.05, 0.91, 0, kBlack, 0.03);
  drawText("North MBD Charge Arrival Time Distributions", 0.95, 0.91, 1, kBlack, 0.03);

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
      SetLineAtt(h_time_raw[i+64], colorline, 2, 1);
      h_time_raw[i+64]->SetFillColorAlpha(colorfill, 0.3);
      h_time_raw[i+64]->Draw();
      drawText(Form("N Ch. %d", i),0.1, 0.88, 0, kBlack, 0.07);
    }

  c->SaveAs(Form("charge_time_%d_%s.png", runnumber,"North"));


  c->Print(Form("charge_time_%d.pdf)", runnumber));

  TCanvas *c2 = new TCanvas("c2","c2", 600, 600);
  SetyjPadStyle();
  SetLineAtt(h_tdc_raw[0], kRed+3, 2, 1);
  h_tdc_raw[0]->SetFillColorAlpha(kRed - 2, .3);
  float max = h_time_raw[0]->GetBinContent(h_time_raw[0]->GetMaximumBin());
  h_time_raw[0]->Scale(1./max);

  h_time_raw[0]->SetMaximum(1.3);
  h_time_raw[0]->Draw("hist");
  h_tdc_raw[0]->Draw("same hist");
  DrawSPHENIX(0.2, 0.87, 1, 0, 0.04);

  TLegend *l_south = new TLegend(0.6, 0.75, 0.85, 0.87);
  l_south->AddEntry(h_time_raw[0], "Charge Time of Arrival");
  l_south->AddEntry(h_tdc_raw[0], "Timing Channel");
  SetLegendStyle(l_south);
  l_south->Draw("same");

  ///

}

void DrawMBDChannels(const  int runnumber = 21813)
{
  gStyle->SetOptStat(0);
  SetyjPadStyle();


  TFile *file = new TFile(Form("../plots/mbd_calib_plots_%d.root", runnumber), "r");

  TH1D *h_time[128];
  TH1D *h_time_raw[128];

  for (int i = 0; i < 128; i++)
    {
      h_time[i] = (TH1D*) file->Get(Form("h_mbd_time_ch%d", i));
      h_time_raw[i] = (TH1D*) file->Get(Form("h_mbd_time_raw_ch%d", i));
      h_time_raw[i]->SetTitle(";TDC_{raw};Counts");
    }

  // ALL 128 Cahnnels;

  

  SetLineAtt(h_time_raw[0], kBlack, 2, 1);
  SetLineAtt(h_time_raw[16], kRed, 2, 1);
  SetLineAtt(h_time_raw[32], kBlue, 2, 1);
  SetLineAtt(h_time_raw[48], kSpring + 2, 2, 1);
  SetLineAtt(h_time_raw[64], kBlack, 2, 3);
  SetLineAtt(h_time_raw[80], kRed, 2, 3);
  SetLineAtt(h_time_raw[96], kBlue, 2, 3);
  SetLineAtt(h_time_raw[112], kSpring + 2, 2, 3);
  
  TCanvas *c_south = new TCanvas("south","south", 600, 600);
  
    
  for (int i = 0; i < 4; i++)
    {
      h_time_raw[(3 - i)*16]->Rebin(10);
      h_time_raw[(3 - i)*16]->GetXaxis()->SetRangeUser(7000, 15000);
      h_time_raw[(3 - i)*16]->SetMaximum(h_time_raw[(3 - i)*16]->GetBinContent(h_time_raw[(3 - i)*16]->GetMaximumBin())*1.4);
      h_time_raw[(3 - i)*16]->Draw(i? "same":"");
    }

  TLegend *l_south = new TLegend(0.2, 0.5, 0.4, 0.7);
  l_south->AddEntry(h_time_raw[0], "Channel 0");
  l_south->AddEntry(h_time_raw[16], "Channel 16");
  l_south->AddEntry(h_time_raw[32], "Channel 32");
  l_south->AddEntry(h_time_raw[48], "Channel 48");
  l_south->Draw("same");

  DrawSPHENIX(0.2, 0.87, 1, 0, 0.04);
  drawText(Form("Run %d", runnumber), 0.2, 0.76, 0, kBlack, 0.04);
  c_south->SaveAs(Form("South_time_%d.png", runnumber));
  TCanvas *c_north_south = new TCanvas("north_south","north_south", 600, 600);

  h_time_raw[0]->Draw();
  h_time_raw[64]->Rebin(10);
  h_time_raw[64]->Draw("same");
  TLegend *l_ns = new TLegend(0.2, 0.6, 0.4, 0.7);
  l_ns->AddEntry(h_time_raw[0], "South - Ch 0");
  l_ns->AddEntry(h_time_raw[64], "North - Ch 0");
  l_ns->Draw("same");

  DrawSPHENIX(0.2, 0.87, 1, 0, 0.04);
  drawText(Form("Run %d", runnumber), 0.2, 0.76, 0, kBlack, 0.04);
  c_north_south->SaveAs(Form("North_South_time_%d.png", runnumber));
  TCanvas *c_good_bad = new TCanvas("good_bad","good_bad", 600, 600);

  h_time_raw[54]->Rebin(10);
  h_time_raw[55]->Rebin(10);
  h_time_raw[56]->Rebin(10);
  SetLineAtt(h_time_raw[54], kViolet - 2, 2, 1);
  SetLineAtt(h_time_raw[55], kBlack, 2, 1);
  SetLineAtt(h_time_raw[56], kRed - 4, 2, 1);
  h_time_raw[55]->GetXaxis()->SetRangeUser(2000, 15000);
  h_time_raw[55]->SetMaximum(h_time_raw[55]->GetBinContent(h_time_raw[55]->GetMaximumBin())*1.4);
  h_time_raw[55]->Draw();
  h_time_raw[54]->Draw("same");
  h_time_raw[56]->Draw("same");

  TLegend *l_gb = new TLegend(0.2, 0.6, 0.4, 0.7);
  l_gb->AddEntry(h_time_raw[55], "South Ch. 55");
  l_gb->AddEntry(h_time_raw[54], "South Ch. 54");
  l_gb->AddEntry(h_time_raw[56], "South Ch. 56");
  l_gb->Draw("same");

  DrawSPHENIX(0.2, 0.87, 1, 0, 0.04);
  drawText(Form("Run %d", runnumber), 0.2, 0.76, 0, kBlack, 0.04);
  c_good_bad->SaveAs(Form("good_bad_time_%d.png", runnumber));

  TCanvas *c_zoom = new TCanvas("zoom","zoom", 600, 600);

  makeMultiPanelCanvas(c_zoom, 1, 2);

  c_zoom->cd(2);
  SetLineAtt(h_time_raw[1], kBlue - 2, 2, 1);
  h_time_raw[1]->SetFillColorAlpha(kBlue - 3, .3);
  h_time_raw[1]->SetTitle(";TDC_{raw};Counts");
  h_time_raw[1]->GetXaxis()->SetTitleOffset(1);
  h_time_raw[1]->GetXaxis()->SetRangeUser( 11000, 14000);
  //h_time_raw[1]->SetMaximum(h_time_raw[1]->GetBinContent(h_time_raw[1]->GetMaximumBin())*1.7);
  h_time_raw[1]->Draw();

  c_zoom->cd(1);

  TH1D *hclone = (TH1D*) h_time_raw[1]->Clone(); 
  hclone->GetXaxis()->SetRangeUser( 12000, 13000);
  hclone->SetMaximum(h_time_raw[1]->GetBinContent(h_time_raw[1]->GetMaximumBin())*1.6);

  hclone->Draw();

  DrawSPHENIX(0.18, 0.7, 1, 0, 0.07);
  drawText(Form("Run %d", runnumber), 0.18, 0.54, 0, kBlack, 0.07);
  c_zoom->SaveAs(Form("zoom_time_%d.png", runnumber));


}
void DrawTwoRunMBDChannels()
{
  gStyle->SetOptStat(0);
  SetyjPadStyle();
  int runnumber1 = 21813;
  int runnumber2 = 24768;
  TFile *file1 = new TFile(Form("../plots/mbd_calib_plots_%d.root", runnumber1), "r");
  TFile *file2 = new TFile(Form("../plots/mbd_calib_plots_%d.root", runnumber2), "r");


  TH1D *h_time_raw1[128];
  TH1D *h_time_raw2[128];

  for (int i = 0; i < 128; i++)
    {

      h_time_raw1[i] = (TH1D*) file1->Get(Form("h_mbd_time_raw_ch%d", i));
      h_time_raw1[i]->SetTitle(";TDC_{raw};Arb. Units [Normalized by Peak]");

      h_time_raw2[i] = (TH1D*) file2->Get(Form("h_mbd_time_raw_ch%d", i));
      h_time_raw2[i]->SetTitle(";TDC_{raw};Arb. Units [Normalized by Peak]");

    }

  // ALL 128 Cahnnels;

  

  SetLineAtt(h_time_raw1[0], kRed, 2, 1);
  SetLineAtt(h_time_raw2[0], kBlue, 2, 1);
  
  TCanvas *c_south = new TCanvas("south","south", 600, 600);
  
  h_time_raw1[0]->Rebin(10);
  h_time_raw1[0]->GetXaxis()->SetRangeUser(7000, 15000);
  h_time_raw1[0]->Scale(1./h_time_raw1[0]->GetBinContent(h_time_raw1[0]->GetMaximumBin()));
  h_time_raw1[0]->SetMaximum(1.3);
  h_time_raw1[0]->Draw("hist");

  h_time_raw2[0]->Rebin(10);
  h_time_raw2[0]->GetXaxis()->SetRangeUser(7000, 15000);
  h_time_raw2[0]->Scale(1./h_time_raw2[0]->GetBinContent(h_time_raw2[0]->GetMaximumBin()));
  h_time_raw2[0]->Draw("hist same");
      

  TLegend *l_south = new TLegend(0.2, 0.55, 0.4, 0.7);
  l_south->SetHeader("South Channel 0");
  l_south->AddEntry(h_time_raw1[0], "Run 21813");
  l_south->AddEntry(h_time_raw2[0], "Run 24768");
  l_south->Draw("same");

  DrawSPHENIX(0.2, 0.87, 1, 0, 0.04);

  c_south->SaveAs("time_dist_run21813_24768.png");
}

void DrawVertexComparison(int runnumber)
{
  gStyle->SetOptStat(0);
  SetyjPadStyle();


  TFile *file = new TFile(Form("../plots/mbd_charge_sum_%d.root", runnumber), "r");

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
  c_time->SaveAs(Form("Time_time_%d.png", runnumber));

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
  c_vertex->SaveAs(Form("Vertex_comp_%d.png", runnumber));
  gPad->SetLogy();
  c_vertex->SaveAs(Form("Vertex_comp_log_%d.png", runnumber));
  
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
}


void DrawEvent(int runnumber)
{
  gStyle->SetOptStat(0);
  SetyjPadStyle();


  TFile *file = new TFile(Form("../plots/mbd_charge_sum_%d.root", runnumber), "r");

  TH2Poly *h_hitmap[100][2];
  TH2Poly *h_timemap[100][2];
  TH1D *h_event[100][2];
  for (int i = 0 ; i < 100; i++)
    {

      for (int j = 0 ; j < 2; j++)
	{
	  h_event[i][j] = (TH1D*) file->Get(Form("event_%d_with_weird_time_%s", i, (j?"n":"s")));
	  h_hitmap[i][j] = (TH2Poly*) file->Get(Form("h_hitmap_%s_tot_%d", (j?"n":"s"), i));
	  h_timemap[i][j] = (TH2Poly*) file->Get(Form("h_timemap_%s_tot_%d", (j?"n":"s"), i));
	}

    }
  
  TCanvas *c = new TCanvas("c","c", 1000, 750);
  TPad *padsquare[4];
  for (int i = 0; i < 4; i++)
    {
      
      padsquare[i] = new TPad(Form("pad_%d", i), "", (i%2)*0.5, .3 + (1 - (i/2))*0.35, (i%2+1)*.5, 0.3 + (2 - (i/2))*0.35);
      padsquare[i]->SetTicks(1,1);
      padsquare[i]->Draw();
    }

  for (int i = 0; i < 4; i++)
    {
      padsquare[i]->cd();
      if (i < 2)
	{
	  //	  h_hitmap[0][i%2]->GetZaxis()->SetRangeUser(0, 1000);
	  h_hitmap[0][i%2]->Draw("colz0");
	}
      else
	{
	  h_timemap[0][i%2]->GetZaxis()->SetRangeUser(-5, 25);
	  h_timemap[0][i%2]->Draw("colz0");
	}
      
    }
}
