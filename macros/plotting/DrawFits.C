#include "dlUtility.h"
void DrawFits(int runnumber = 21813)
{
  gStyle->SetOptStat(0);
  TFile *fout = new TFile(Form("../plots/mbd_calib_plots_%d.root", runnumber), "r");

  TH1D *h_charge_raw[128];
  for (int i = 0; i < 128; i++)
    {
      h_charge_raw[i] = (TH1D*) fout->Get(Form("h_mbd_charge_raw_ch%d", i));
    }

  TCanvas *c = new TCanvas("c","c", 1000, 1200);

  c->Print(Form("fits_%d.pdf[", runnumber));
  int colorline = kSpring + 4;
  int colorfill = kSpring + 2;
  int colorfit = kBlack;

  TPad *pads[64];
  float xi = 0.05;
  float yi = .9 - xi;
  float dpx = (1. - 2*xi)/8.;
  float dpy = (.9 - 2*xi)/8.;
  DrawSPHENIX(0.05, 0.95, 1, 1, 0.03, 0);
  drawText(Form("Run - %d", runnumber), 0.05, 0.91, 0, kBlack, 0.03);
  drawText("Dotted Line at 0.4 #times MPV_{fit}", 0.05, 0.87, 0, kBlack, 0.03);
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

  for (int i = 0 ; i < 64; i++)
    {
      pads[i]->cd();
      std::cout << i << std::endl;
      TF1 *fitfunc = h_charge_raw[i]->GetFunction("lan_w_gausexp");
      if (!fitfunc)
	{
	  fitfunc = h_charge_raw[i]->GetFunction("lan_w_exp");
	}
      if (!fitfunc) continue;
      float peakfit = fitfunc->GetMaximum();
      h_charge_raw[i]->SetMaximum(peakfit * 1.3);
      h_charge_raw[i]->GetXaxis()->SetRangeUser(0, fitfunc->GetParameter(1)*2);
      SetLineAtt(h_charge_raw[i], colorline, 2, 1);
      h_charge_raw[i]->SetFillColorAlpha(colorfill, 0.3);
      SetLineAtt(fitfunc, colorfit, 3, 1);
      h_charge_raw[i]->Draw();
      TLine *tl = new TLine(fitfunc->GetParameter(1)*0.4, 0,fitfunc->GetParameter(1)*0.4, peakfit*1.3);
      tl->SetLineColor(kBlue + 2);
      tl->SetLineStyle(2);
      tl->Draw();
      drawText(Form("S Ch. %d", i),0.5, 0.88, 0, kBlack, 0.07);
      drawText(Form("#Chi^2 = %4.4f", fitfunc->GetChisquare()/fitfunc->GetNDF()), 0.5, 0.8, 0, kBlack, 0.07);
      drawText(Form("MPV = %4.2f", fitfunc->GetParameter(1)), 0.5, 0.72, 0, kBlack, 0.07);
    }
  c->Print(Form("fits_%d.pdf(", runnumber));

  c->SaveAs(Form("fits_%d_%s.png", runnumber,"South"));
  c = new TCanvas("c1","c1", 1000, 1200);
  DrawSPHENIX(0.05, 0.95, 1, 1, 0.03, 0);
  drawText(Form("Run - %d", runnumber), 0.05, 0.91, 0, kBlack, 0.03);
  drawText("Dotted Line at 0.4 #times MPV_{fit}", 0.05, 0.87, 0, kBlack, 0.03);
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
      TF1 *fitfunc = h_charge_raw[i+64]->GetFunction("lan_w_gausexp");
      if (!fitfunc)
	{
	  fitfunc = h_charge_raw[i+64]->GetFunction("lan_w_exp");
	}
      if (!fitfunc) continue;
      float peakfit = fitfunc->GetMaximum();
      h_charge_raw[i+64]->SetMaximum(peakfit * 1.3);
      h_charge_raw[i+64]->GetXaxis()->SetRangeUser(0, fitfunc->GetParameter(1)*2);
      SetLineAtt(h_charge_raw[i+64], colorline, 1, 1);
      h_charge_raw[i+64]->SetFillColorAlpha(colorfill, 0.3);
      SetLineAtt(fitfunc, colorfit, 3, 1);
      h_charge_raw[i+64]->Draw();
      TLine *tl = new TLine(fitfunc->GetParameter(1)*0.4, 0,fitfunc->GetParameter(1)*0.4, peakfit*1.3);
      tl->SetLineColor(kBlue + 2);
      tl->SetLineStyle(2);
      tl->Draw();
      drawText(Form("N Ch. %d", i),0.5, 0.88, 0, kBlack, 0.07);
      drawText(Form("#Chi^2 = %4.4f", fitfunc->GetChisquare()/fitfunc->GetNDF()), 0.5, 0.8, 0, kBlack, 0.07);
      drawText(Form("MPV = %4.2f", fitfunc->GetParameter(1)), 0.5, 0.72, 0, kBlack, 0.07);
    }

  c->SaveAs(Form("fits_%d_%s.png", runnumber,"North"));


  c->Print(Form("fits_%d.pdf)", runnumber));
}

void DrawTOAFits(int runnumber = 21813)
{
  gStyle->SetOptStat(0);
  TFile *fout = new TFile(Form("fitsout_%d.root", runnumber), "r");

  TH1D *h_toa_channel[128];
  for (int i = 0; i < 128; i++)
    {
      h_toa_channel[i] = (TH1D*) fout->Get(Form("h_toa_channel_%d", i));
    }

  TCanvas *c = new TCanvas("c","c", 1000, 1200);

  c->Print(Form("fits_%d.pdf[", runnumber));
  int colorline;
  int colorfill;
  int colorfit;
  colorfill = kRed - 2;
  colorfit = kBlack;
  colorline = kRed + 2;
  TPad *pads[64];
  float xi = 0.05;
  float yi = .9 - xi;
  float dpx = (1. - 2*xi)/8.;
  float dpy = (.9 - 2*xi)/8.;
  DrawSPHENIX(0.05, 0.95, 1, 1, 0.03, 0);
  drawText(Form("Run - %d", runnumber), 0.05, 0.91, 0, kBlack, 0.03);
  //drawText("Dotted Line at 0.4 #times MPV_{fit}", 0.05, 0.87, 0, kBlack, 0.03);
  drawText("South MBD Timing Distributions", 0.95, 0.91, 1, kBlack, 0.03);
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
      std::cout << i << std::endl;
      TF1 *fitfunc = h_toa_channel[i]->GetFunction("triple_gaussian");
      if (!fitfunc)
	{
	  fitfunc = h_toa_channel[i]->GetFunction("double_gaussian");
	}
      if (!fitfunc) continue;
      float peakfit = fitfunc->GetMaximum();
      h_toa_channel[i]->SetMaximum(peakfit * 1.3);
      double x_low = -10;
      double x_high = 10;
      
      h_toa_channel[i]->GetXaxis()->SetRangeUser(x_low, x_high);
      SetLineAtt(h_toa_channel[i], colorline, 2, 1);
      h_toa_channel[i]->SetFillColorAlpha(colorfill, 0.3);
      SetLineAtt(fitfunc, colorfit, 3, 1);
      h_toa_channel[i]->Draw();

      drawText(Form("S Ch. %d", i),0.5, 0.88, 0, kBlack, 0.07);
      drawText(Form("#Chi^2 = %4.4f", fitfunc->GetChisquare()/fitfunc->GetNDF()), 0.5, 0.8, 0, kBlack, 0.07);
      //      drawText(Form("MPV = %4.2f", fitfunc->GetParameter(1)), 0.5, 0.72, 0, kBlack, 0.07);
    }
  c->Print(Form("fitstoa_%d.pdf(", runnumber));

  c->SaveAs(Form("fitstoa_%d_%s.png", runnumber,"South"));

  c = new TCanvas("c1","c1", 1000, 1200);
  DrawSPHENIX(0.05, 0.95, 1, 1, 0.03, 0);
  drawText(Form("Run - %d", runnumber), 0.05, 0.91, 0, kBlack, 0.03);
  //  drawText("Dotted Line at 0.4 #times MPV_{fit}", 0.05, 0.87, 0, kBlack, 0.03);
  drawText("North MBD Timing Distributions", 0.95, 0.91, 1, kBlack, 0.03);

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
      
      TF1 *fitfunc = h_toa_channel[i]->GetFunction("triple_gaussian");
      if (!fitfunc)
	{
	  fitfunc = h_toa_channel[i]->GetFunction("double_gaussian");
	}
      if (!fitfunc) continue;
      float peakfit = fitfunc->GetMaximum();
      h_toa_channel[i]->SetMaximum(peakfit * 1.3);
      double x_low, x_high;
      fitfunc->GetRange(x_low, x_high);
      
      h_toa_channel[i]->GetXaxis()->SetRangeUser(x_low, x_high);
      SetLineAtt(h_toa_channel[i], colorline, 2, 1);
      h_toa_channel[i]->SetFillColorAlpha(colorfill, 0.3);
      SetLineAtt(fitfunc, colorfit, 3, 1);
      h_toa_channel[i]->Draw();

      drawText(Form("S Ch. %d", i),0.5, 0.88, 0, kBlack, 0.07);
      drawText(Form("#Chi^2 = %4.4f", fitfunc->GetChisquare()/fitfunc->GetNDF()), 0.5, 0.8, 0, kBlack, 0.07);
      //      drawText(Form("MPV = %4.2f", fitfunc->GetParameter(1)), 0.5, 0.72, 0, kBlack, 0.07);
    }

  c->SaveAs(Form("fitstoa_%d_%s.png", runnumber,"North"));


  c->Print(Form("fitstoa_%d.pdf)", runnumber));
}

void DrawTimeFits(int runnumber = 21813)
{
  gStyle->SetOptStat(0);
  TFile *fout = new TFile(Form("fitsout_%d.root", runnumber), "r");

  TH1D *h_toa_channel[128];
  TH1D *h_tdc_channel[128];
  for (int i = 0; i < 128; i++)
    {
      h_toa_channel[i] = (TH1D*) fout->Get(Form("h_toa_channel_%d", i));
      h_tdc_channel[i] = (TH1D*) fout->Get(Form("h_tdc_channel_%d", i));
    }

  TCanvas *c = new TCanvas("c","c", 1000, 1200);

  c->Print(Form("fits_%d.pdf[", runnumber));
  int colorline;
  int colorfill;
  int colorfit;
  colorfill = kRed - 2;
  colorfit = kBlack;
  colorline = kRed + 2;
  TPad *pads[64];
  float xi = 0.05;
  float yi = .9 - xi;
  float dpx = (1. - 2*xi)/8.;
  float dpy = (.9 - 2*xi)/8.;
  DrawSPHENIX(0.05, 0.95, 1, 1, 0.03, 0);
  drawText(Form("Run - %d", runnumber), 0.05, 0.91, 0, kBlack, 0.03);
  //drawText("Dotted Line at 0.4 #times MPV_{fit}", 0.05, 0.87, 0, kBlack, 0.03);
  drawText("South MBD Timing Distributions", 0.95, 0.91, 1, kBlack, 0.03);
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
      std::cout << i << std::endl;
      TF1 *fitfunc = h_tdc_channel[i]->GetFunction("triple_gaussian");
      if (!fitfunc)
	{
	  fitfunc = h_tdc_channel[i]->GetFunction("double_gaussian");
	}
      if (!fitfunc) continue;
      float peakfit = fitfunc->GetMaximum();
      double x_low = -10;
      double x_high = 10;
      
      h_tdc_channel[i]->GetXaxis()->SetRangeUser(x_low, x_high);
      SetLineAtt(h_tdc_channel[i], colorline, 2, 1);
      h_tdc_channel[i]->SetFillColorAlpha(colorfill, 0.3);
      SetLineAtt(fitfunc, colorfit, 3, 1);
      h_tdc_channel[i]->Draw();

      drawText(Form("S Ch. %d", i),0.5, 0.88, 0, kBlack, 0.07);
      drawText(Form("#Chi^2 = %4.4f", fitfunc->GetChisquare()/fitfunc->GetNDF()), 0.5, 0.8, 0, kBlack, 0.07);
      //      drawText(Form("MPV = %4.2f", fitfunc->GetParameter(1)), 0.5, 0.72, 0, kBlack, 0.07);
    }
  c->Print(Form("fitstdc_%d.pdf(", runnumber));

  c->SaveAs(Form("fitstdc_%d_%s.png", runnumber,"South"));

  c = new TCanvas("c1","c1", 1000, 1200);
  DrawSPHENIX(0.05, 0.95, 1, 1, 0.03, 0);
  drawText(Form("Run - %d", runnumber), 0.05, 0.91, 0, kBlack, 0.03);
  //  drawText("Dotted Line at 0.4 #times MPV_{fit}", 0.05, 0.87, 0, kBlack, 0.03);
  drawText("North MBD Timing Distributions", 0.95, 0.91, 1, kBlack, 0.03);

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
      
      TF1 *fitfunc = h_tdc_channel[i]->GetFunction("triple_gaussian");
      if (!fitfunc)
	{
	  fitfunc = h_tdc_channel[i]->GetFunction("double_gaussian");
	}
      if (fitfunc)
	{
	  float peakfit = fitfunc->GetMaximum();
	  double x_low = -10;
	  double x_high = 10;
	  
	  h_tdc_channel[i]->GetXaxis()->SetRangeUser(x_low, x_high);
	  SetLineAtt(h_tdc_channel[i], colorline, 2, 1);
	  h_tdc_channel[i]->SetFillColorAlpha(colorfill, 0.3);
	  SetLineAtt(fitfunc, colorfit, 3, 1);
	}
      h_tdc_channel[i]->Draw();
      if (!fitfunc) continue;
      drawText(Form("S Ch. %d", i),0.5, 0.88, 0, kBlack, 0.07);
      drawText(Form("#Chi^2 = %4.4f", fitfunc->GetChisquare()/fitfunc->GetNDF()), 0.5, 0.8, 0, kBlack, 0.07);
      //      drawText(Form("MPV = %4.2f", fitfunc->GetParameter(1)), 0.5, 0.72, 0, kBlack, 0.07);
    }

  c->SaveAs(Form("fitstdc_%d_%s.png", runnumber,"North"));


  c->Print(Form("fitstdc_%d.pdf)", runnumber));
}
