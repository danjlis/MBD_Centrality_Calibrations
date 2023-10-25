#include <iostream>
#include <fstream>
#include <cctype>
#include <string>
#include <math.h>

#include "TStyle.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TRandom.h"
#include "TMath.h"
#include "TLine.h"
#include "TLegend.h"
#include "sPhenixStyle.C"
#include "dlUtility.h"
using namespace std;

// simple routine (j.nagle - june 04, 2009)
// 1. read in histogram of Glauber Ncoll distribution
// 2. calculate with NBD the expected BBC South distribution
// 3. compare to real data and fit trigger turn on curve

// one can also calculate the shifted SimBBC spectra if one
// of the binary collisions is biased by a hard process (for example)
// and in this case we assume that the process scales with nbinary !!!

// add code to have a different number of centrality bins (e.g. 4 or 8)
// we will also need to calculate the real constraints on mu, k from data
// we will also need to run multiple cases with +/- uncertainties on
// some of the input parameters
//
// --> mu +/- , k +/-, biased_mu +/-, biased_k +/-
//     trigger curve for biased events +/- (to still agree with ~ 75% value)
//
//     different Glauber input distributions (read in root file)-> that is not too hard
//     however, for each set of glauber inputs, need to fit to find best mu,k
//     and then fit regular trigger curve, and then how to find the biased trigger curve
//     step through parameter values and find best 75% match (done)

double NBD_getValue(int n, double mu, double k);

// here is the main routine
void nbd_check(int runnumber = 21813,
	       int maxcentbins = 19,
	       bool flag_determineNBDparameters = true, 
	       bool flag_determineBiasFactors = false, string name = "./lemon_glauber_auau200_100k_histo.root", double alpha = 1.0, int vtxflag = 1, double particlealpha = 1.00, bool npartscaling = true, bool centraltrigger = false,
	       bool npartminusone = false, double forcetrigfrac = 93., bool sphenix = true, int which=0) {

  //===============================================================================================
  double mu = 4.4;  //defaults that are used if not determining these parameters again
  double k  = 2.14;

  double biasfactor = 1.55;
  double biased_mu = mu * biasfactor; // larger number of hits
  double biased_k  =  k * biasfactor;

  int lowfit = 700; // lowest end of standard fit range, was 30 in PHENIX
  // utyy
  //===============================================================================================

  //===============================================================================================
  // the below file contains the "hNcoll" histogram to be sampled from (step #1)
  TFile *fglauber = new TFile(name.c_str()); 
  TH1D  *hNcoll;
  if (!npartscaling) {
    hNcoll = (TH1D *) fglauber->Get("hNcoll"); // this is the default case...
  } else {
    hNcoll = (TH1D * ) fglauber->Get("hNpart"); // note that this is only npart in the Au nucleus
    // this creates a problem because there is never an Npart = 1, always Npart >= 2
    // so bias correction procedure has a flaw; so in this case shift down by one
    for (int i=1;i<=hNcoll->GetNbinsX()-1;i++) {
      if (npartminusone) hNcoll->SetBinContent(i, hNcoll->GetBinContent(i+1)); // shift down by one...
    }
  }

  // the below file contains the "hbbcQs6" histogram of BBC data distribution (from z-vertex range)
  //=======================================================
  char *fdataname = new char[100];
  if (!sphenix) sprintf(fdataname,"./fout_run14auau200gev_dist.root");
  if (sphenix) sprintf(fdataname,"../..//output/run%d/plots/mbd_charge_sum_%d.root", runnumber, runnumber);  
//  if (sphenix) sprintf(fdataname,"../output/run%d/trees_%d.root", runnumber, runnumber);
  //if (sphenix) fdataname = "./fout_sphenix2023_auau200_23722_both_z20.root";
  //  if (sphenix) fdataname = "./fout_sphenix2023_auau200_22979_both_z20.root";

  char runnum[100];
  sprintf(runnum, "000%d", runnumber);
  char zselect[100] = "|z_{MBD}|< 30 cm";
  //=======================================================  
  
  TFile *fdata = new TFile(fdataname);
  TH1D *hRealBBC;
  //if (which==0) hRealBBC = (TH1D *) fdata->Get("h_mbd_charge_ns_w_zdc_cut_w_mbd_cut_and_vertex_10"); // a particular z-vertex range
  if (which==0) hRealBBC = (TH1D *) fdata->Get("h_charge_sum_min_bias_w_vertex_30"); // a particular z-vertex range
  if (which==1) hRealBBC = (TH1D *) fdata->Get("hmbdchargesouth"); // a particular z-vertex range
  if (which==2) hRealBBC = (TH1D *) fdata->Get("hmbdchargenorth"); // a particular z-vertex range  
  if (!hRealBBC)   hRealBBC = (TH1D *) fdata->Get("hbbcQs6"); // a particular z-vertex range (old PHENIX style)
  //===============================================================================================

  //hRealBBC->Scale(1./hRealBBC->Integral());
  int nhistbins = hRealBBC->GetNbinsX();
  // j.nagle - also change 199.5 to be maxrange = ((float)nhistbins) - 0.5
  float maxrange = ((float) nhistbins) - 0.5;
  // might be good to have max Ncoll / Npart as well... (CuAu = 197 + 63 = 270)
  // Au+Au case 197 + 197 = 396
  int nncollmax = 400; // really npart here
  float maxrangencoll = ((float) nncollmax) - 0.5;

  TH1D *hSimBBCHardUnbiased = new TH1D("hSimBBCHardUnbiased","hSimBBCHardUnbiased",nhistbins,-0.5,maxrange);
  TH1D *hSimBBCHardBiased = new TH1D("hSimBBCHardBiased","hSimBBCHardBiased",nhistbins,-0.5,maxrange);

  // for these two, the trigger is not required to have fired ... (BBC vs Ncoll)
  TH2D *hSimBBCNcoll = new TH2D("hSimBBCNcoll","hSimBBCNcoll",nhistbins,-0.5,maxrange,nncollmax,-0.5,maxrangencoll);
  TH2D *hSimBBCNcollHard = new TH2D("hSimBBCNcollHard","hSimBBCNcollHard",nhistbins,-0.5,maxrange,nncollmax,-0.5,maxrangencoll);

  // this now includes the trigger efficiency turn-on and can be used to calculate <Ncoll> for each centrality bin!!!
  TH2D *hSimBBCNcoll_wtrig = new TH2D("hSimBBCNcoll_wtrig","hSimBBCNcoll_wtrig",nhistbins,-0.5,maxrange,nncollmax,-0.5,maxrangencoll);

  // step #2
  TH1D *hSimBBC = new TH1D("hSimBBC","hSimBBC",nhistbins,-0.5,maxrange);

  TF1 *flatline = new TF1("flatline","[0]",lowfit,400.0);
  flatline->SetParameter(0,1.0);

  //-----------------------------------------------------------

  double bestmu = mu; double bestk = k;
  if (flag_determineNBDparameters) {

    double lowmu = 3.0;
    double highmu = 5.0;
    double lowk = 0.3;
    double highk = 2.0;

    if (which > 0) {
      // only one side
      lowmu = lowmu/2.0;
      highmu = highmu/2.0;
      lowk = lowk/2.0;
      highk = highk/2.0;
    }
    
    int nmusteps = 20;
    int nksteps = 20;

    double binsizemu = (highmu-lowmu)/((double)nmusteps);
    double binsizek = (highk-lowk)/((double)nksteps);

    TH2D *hmuk_chi2 = new TH2D("hmuk_chi2","hmuk_chi2",nmusteps,lowmu-0.5*binsizemu,highmu+0.5*binsizemu,
			       nksteps,lowk-0.5*binsizek,highk+0.5*binsizek);

    for (int imu=0;imu<nmusteps;imu++) {
      for (int ik=0;ik<nksteps;ik++) {

	double mutemp = hmuk_chi2->GetXaxis()->GetBinCenter(imu+1);
	double ktemp = hmuk_chi2->GetYaxis()->GetBinCenter(ik+1);

	hSimBBC->Reset();
	// loop over nNcoll values
	for (int ib=1; ib<=hNcoll->GetNbinsX(); ib++) {

	  int ncoll = (int) hNcoll->GetBinCenter(ib);
	  double event_weightfactor = hNcoll->GetBinContent(ib);
	  if (event_weightfactor <= 0) continue;

	  for (int ihit=0; ihit< (20 * ncoll * (int) mutemp); ihit++) {

	    // TRY THIS FOR AUAU 
	    if (ncoll > 10 && ihit > (1.5 * ncoll * (int) mu)) continue;

	    double nbdvalue = NBD_getValue(ihit, mutemp * (double) (pow(ncoll,alpha)), ktemp * (double) (pow(ncoll,alpha)));
	    hSimBBC->Fill((double) ihit, nbdvalue * event_weightfactor);

	  }
	} 
	// now calculate the chi2 above some minimum value - staying away from the trigger turn on
	hSimBBC->Scale(hRealBBC->Integral(lowfit,nhistbins)/hSimBBC->Integral(lowfit,nhistbins));	
	TH1D *hRatio = new TH1D("hRatio","hRatio",nhistbins,-0.5,maxrange);
	hRatio->Reset();
	for (int i=1;i<=nhistbins;i++) {
	  if (hSimBBC->GetBinContent(i) > 0 && hRealBBC->GetBinContent(i) > 0) {
	    hRatio->SetBinContent(i, hRealBBC->GetBinContent(i)/hSimBBC->GetBinContent(i));
	    hRatio->SetBinError(i, (sqrt(hRealBBC->GetBinContent(i))/hSimBBC->GetBinContent(i)));
	  }
	}	
	flatline->SetParameter(0,1.0);
	hRatio->Fit(flatline,"NDORQ","",(double)lowfit,maxrange);
	
	cout << "Parameters mu = " << mutemp << " k = " << ktemp << " chi^2 = " << 
	  flatline->GetChisquare() << " and prob = " << TMath::Prob(flatline->GetChisquare(),100) << endl;
	
	
	hRatio->Delete();

	hmuk_chi2->SetBinContent(imu+1,ik+1,flatline->GetChisquare());

      } // end loop over k steps
    } // end loop over mu steps

    // calculate the best value and +/- 1 sigma ranges
    double lowestchi2 = 99999999.;
    for (int imu=0;imu<nmusteps;imu++) {
      for (int ik=0;ik<nksteps;ik++) {
	if (hmuk_chi2->GetBinContent(imu+1,ik+1) < lowestchi2 && hmuk_chi2->GetBinContent(imu+1,ik+1) != 0.0) {
	  lowestchi2 = hmuk_chi2->GetBinContent(imu+1,ik+1);
	  bestmu = hmuk_chi2->GetXaxis()->GetBinCenter(imu+1);
	  bestk  = hmuk_chi2->GetYaxis()->GetBinCenter(ik+1);
	}
      }
    }
    double bestmu_minusonesigma = 99999.; double bestmu_plusonesigma = -99999.;
    double bestk_minusonesigma = 99999.; double bestk_plusonesigma = -99999.;
    for (int imu=0;imu<nmusteps;imu++) {
      for (int ik=0;ik<nksteps;ik++) {
	if (hmuk_chi2->GetBinContent(imu+1,ik+1) < lowestchi2+1.0 && hmuk_chi2->GetBinContent(imu+1,ik+1) != 0.0) {

	  if (hmuk_chi2->GetXaxis()->GetBinCenter(imu+1) < bestmu_minusonesigma) 
	    bestmu_minusonesigma = hmuk_chi2->GetXaxis()->GetBinCenter(imu+1);
	  if (hmuk_chi2->GetXaxis()->GetBinCenter(imu+1) > bestmu_plusonesigma) 
	    bestmu_plusonesigma = hmuk_chi2->GetXaxis()->GetBinCenter(imu+1);

	  if (hmuk_chi2->GetYaxis()->GetBinCenter(ik+1) < bestk_minusonesigma) 
	    bestk_minusonesigma = hmuk_chi2->GetYaxis()->GetBinCenter(ik+1);
	  if (hmuk_chi2->GetYaxis()->GetBinCenter(ik+1) > bestk_plusonesigma) 
	    bestk_plusonesigma = hmuk_chi2->GetYaxis()->GetBinCenter(ik+1);

	}
      }
    }
    // end loop over parameters

    cout << "=====================================================" << endl;
    cout << "Best chi2 = " << lowestchi2 << " with mu = " << 
      bestmu << " +/- " << bestmu-bestmu_minusonesigma << " / " << bestmu_plusonesigma-bestmu << 
      " and k = " <<
      bestk << " +/- " << bestk-bestk_minusonesigma << " / " << bestk_plusonesigma-bestk << endl;
    cout << "=====================================================" << endl;
	
    TCanvas *cmuk = new TCanvas("canvas_mu_k_constraint","canvas_mu_k_constraint");
    cmuk->cd();
    gPad->SetTicks(1);	
    hmuk_chi2->SetXTitle("NBD mu parameter");
    hmuk_chi2->SetYTitle("NBD k parameter");
    hmuk_chi2->DrawCopy("text");
    hmuk_chi2->DrawCopy("colz");    

    // use the best values in the rest of the calculation
    mu = bestmu;
    k  = bestk;
    biased_mu = mu * biasfactor;
    biased_k  =  k * biasfactor;

  } // end section for determining NBD values

  hSimBBC->Reset();

  //-------------------------------------------------------------------------------------------
  // now use the best values (or default values) and re-calculate things, including bias factors...
  // loop over nNcoll values
  for (int ib=1; ib<=hNcoll->GetNbinsX(); ib++) {

    int    ncoll = (int) hNcoll->GetBinCenter(ib);

    double event_weightfactor = hNcoll->GetBinContent(ib);
    if (event_weightfactor <= 0) continue;

    double hardprocess_weightfactor = pow(ncoll,particlealpha); // assume for now that hard process scales with ncoll

    // scale parameters by mu' = my * ncoll and k' = k * ncoll
    for (int ihit=0; ihit< (20 * ncoll * (int) mu); ihit++) {

      double nbdvalue = NBD_getValue(ihit, mu * (double) (pow(ncoll,alpha)), k * (double) (pow(ncoll,alpha)));
      hSimBBC->Fill((double) ihit, nbdvalue * event_weightfactor);
      // below value really needs to be ncoll to work correctly!!!!!
      hSimBBCNcoll->Fill((double) ihit, (double) ncoll, nbdvalue * event_weightfactor);
      hSimBBCHardUnbiased->Fill((double) ihit, nbdvalue * event_weightfactor * hardprocess_weightfactor);

    }
    
    if (flag_determineBiasFactors) {
      // now try biased case

      // IHIT REPRESENTS NUMBER OF HITS IN THE BBC DETECTOR - SO NEED TO LOOP UP TO HIGH ENOUGH VALUES TO GET THE TAIL
      for (int ihit=0; ihit< (20 * ncoll * (int) mu); ihit++) {

	// TRY THIS FOR AUAU 
	if (ncoll > 10 && ihit > (1.5 * ncoll * (int) mu)) continue;

	// now the problem is that you can get the total ihit by the sum of two values
	for (int ihitbackground=0; ihitbackground<=ihit; ihitbackground++) {

	  // for ncoll==1, then no background hits are possible

	  int keyncollindex = 1; // Ncoll=1
	  if (npartscaling) keyncollindex = 2;
	  if (npartscaling && npartminusone) keyncollindex = 1;

	  if (ncoll <= keyncollindex && ihitbackground != 0) continue;
	  
	  double nbdvalue = NBD_getValue(ihitbackground, mu * (double) (pow((ncoll-1),alpha)), k * (double) (pow((ncoll-1),alpha)));
	  if (ncoll <= keyncollindex) nbdvalue = 1.0;

	  // if the one biased ncoll needs to account for too many of the remaining hits, skip calc
	  int ihithard = ihit - ihitbackground;
	  double hardnbdvalue = NBD_getValue(ihithard, biased_mu, biased_k);
	  if (nbdvalue < 0.000000001 || hardnbdvalue < 0.000000001) continue;
	  hSimBBCHardBiased->Fill((double) ihit, 
				  nbdvalue * hardnbdvalue * event_weightfactor * hardprocess_weightfactor);
	  hSimBBCNcollHard->Fill((double) ihit, (double) ncoll, 
				 nbdvalue * hardnbdvalue * event_weightfactor * hardprocess_weightfactor);
	}
      }
    }

  } // end loop over possible Ncoll values

  cout << "Ended loop over creating fake data histograms" << endl;

  // step #3 

  SetsPhenixStyle();

  TCanvas *c1 = new TCanvas("canvas_bbc_realdataandsim","canvas_bbc_realdataandsim", 700, 700);
  c1->Divide(1,2);
  c1->cd(1);
  gPad->SetTicks(1);	
  // plot the comparison of the Real Data BBC and the Simulation
  hSimBBC->Scale(hRealBBC->Integral(lowfit,nhistbins)/hSimBBC->Integral(lowfit,nhistbins));
  hRealBBC->SetLineWidth(2);
  hRealBBC->SetMarkerStyle(24);
  hRealBBC->SetMaximum(hSimBBC->GetBinContent(10)*3);
  hRealBBC->GetXaxis()->SetRangeUser(0.0,2500);
  gPad->SetTopMargin(.13);
  
  char foo[100];
  if (sphenix) sprintf(foo,"MBD Charge N+S (Glauber+NBD fit #mu=%3.2f, k=%3.2f)",bestmu,bestk);
  if (sphenix && which==1) sprintf(foo,"MBD Charge South (Glauber+NBD fit #mu=%3.2f, k=%3.2f)",bestmu,bestk);
  if (sphenix && which==2) sprintf(foo,"MBD Charge North (Glauber+NBD fit #mu=%3.2f, k=%3.2f)",bestmu,bestk);  

  hRealBBC->SetTitle(Form(";%s;N_{events}", foo));
  if (!sphenix) hRealBBC->SetXTitle("BBC Charge Total (North + South)");  
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
  //  drawText("8/30/2023", 0.8, 0.92, 0, kBlack, 0.07);
  TLatex l;
  l.SetNDC();
  l.SetTextSize(0.04);
  
  DrawSPHENIX(0.65, 0.8, 1, 0, 0.06);
  drawText("|z_{MBD}| < 30 cm", 0.65, 0.64, 0, kBlack, 0.06);

  // now do the trigger efficiency fitting
  c1->cd(2);
  gPad->SetTicks(1);
  TH1D *hRatio = new TH1D("hRatio","hRatio",nhistbins,-0.5, maxrange);
  for (int i=1;i<=nhistbins;i++) {
    if (hSimBBC->GetBinContent(i) > 0 && hRealBBC->GetBinContent(i) > 0) {
      hRatio->SetBinContent(i, hRealBBC->GetBinContent(i)/hSimBBC->GetBinContent(i));
      hRatio->SetBinError(i, (sqrt(hRealBBC->GetBinContent(i))/hSimBBC->GetBinContent(i)));
    }
  }
  hRatio->SetMarkerStyle(24);
  hRatio->SetMaximum(1.5);
  hRatio->SetMinimum(0);
  hRatio->GetXaxis()->SetRangeUser(0.0,700);
  hRatio->SetXTitle(foo);

  if (!sphenix) hRatio->SetXTitle("BBC Charge Total (North + South)");
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
  
  TF1 *trigeffcurve = new TF1("trigeffcurve","1.0-TMath::Exp(-pow((x/[0]), [1]))",0.0,maxrange);
  trigeffcurve->SetParameters(0.732,0.532); // just initial value guesses
  trigeffcurve->SetParLimits(0,0.2,20.0);    
  trigeffcurve->SetParLimits(1,0.2,5.0);
  hRatio->Fit(trigeffcurve,"NDOR","",1.0,80.0);
  trigeffcurve->SetLineColor(4);
  trigeffcurve->SetLineWidth(3);
  trigeffcurve->Draw("same");
  cout << "The trigeff curve (0:120) has chi^2 = " << 
    trigeffcurve->GetChisquare() << " and prob = " << TMath::Prob(trigeffcurve->GetChisquare(),120-2) << endl;

  // also fit a flat line above say 20 -> 140 ==> but just constraint it exactly at 1.0
  flatline->SetParameter(0,1.0);
  hRatio->Fit(flatline,"NDOR","",(double)lowfit,maxrange);
  cout << "The line at one (lowfit:120) has chi^2 = " << 
    flatline->GetChisquare() << " and prob = " << TMath::Prob(flatline->GetChisquare(),100) << endl;

  TLine *tl = new TLine(0.0,1.0,maxrange,1.0);
  tl->Draw("same");

  
  
  // trigger efficiency from integral comparison
  double trigeffintegral = hRealBBC->Integral()/hSimBBC->Integral();
  cout << "Trigger Efficiency from Integrals = " << trigeffintegral << endl;


  c1->cd(1);
  gPad->SetTicks(1);	
  gPad->SetLogy(1);
  TH1D *hSimBBCwTrig = new TH1D("hSimBBCwTrig","hSimBBCwTrig",nhistbins,-0.5,maxrange);
  for (int i=1;i<=nhistbins;i++) {
    if (i==1) {
      hSimBBCwTrig->SetBinContent(i,0.0); // no chance to fire the trigger if no BBC south hits
    } else {
      hSimBBCwTrig->SetBinContent(i,hSimBBC->GetBinContent(i)*trigeffcurve->Eval(hSimBBC->GetBinCenter(i)));
    }
  }
  hSimBBCwTrig->DrawCopy("same");

  for (int i=1;i<=nhistbins;i++) { // loop over BBC hits
    for (int j=1;j<=nncollmax;j++) { // loop over Ncoll values
      if (i==1) {
	hSimBBCNcoll_wtrig->SetBinContent(i,j,0.0);
      } else {
	hSimBBCNcoll_wtrig->SetBinContent(i,j,
		 hSimBBCNcoll->GetBinContent(i,j)*trigeffcurve->Eval(hSimBBC->GetBinCenter(i)));
      }
    }
  }

  // ---------------------------------------------------------------------
  // Divide the Sim distribution into 20/forcetrigfrac, 40/forcetrigfrac, 60/forcetrigfrac regions by cut values
  double centrality_low[20] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double centrality_high[20] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  centrality_low[maxcentbins-1] = 0.0;
  int ihigh;
  for (ihigh = nhistbins; ihigh > 0; ihigh--)
    {
      if (hRealBBC->GetBinContent(ihigh))
	{
	  centrality_high[0] = hRealBBC->GetBinCenter(ihigh);
	  break;
	}
    }
  for (int i=ihigh;i>=1;i--) {
    if (maxcentbins == 4) {
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (20./forcetrigfrac) && centrality_low[0] == 0.0) {
	centrality_low[0] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[1] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (40./forcetrigfrac) && centrality_low[1] == 0.0) {
	centrality_low[1] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[2] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (60./forcetrigfrac) && centrality_low[2] == 0.0) {
	centrality_low[2] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[3] = hSimBBCwTrig->GetBinCenter(i);
      }
    } 
    if (maxcentbins == 8) {
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (10./forcetrigfrac) && centrality_low[0] == 0.0) {
	centrality_low[0] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[1] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (20./forcetrigfrac) && centrality_low[1] == 0.0) {
	centrality_low[1] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[2] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (30./forcetrigfrac) && centrality_low[2] == 0.0) {
	centrality_low[2] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[3] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (40./forcetrigfrac) && centrality_low[3] == 0.0) {
	centrality_low[3] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[4] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (50./forcetrigfrac) && centrality_low[4] == 0.0) {
	centrality_low[4] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[5] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (60./forcetrigfrac) && centrality_low[5] == 0.0) {
	centrality_low[5] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[6] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (70./forcetrigfrac) && centrality_low[6] == 0.0) {
	centrality_low[6] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[7] = hSimBBCwTrig->GetBinCenter(i);
      }
    }
    if (maxcentbins == 9) {
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (10./forcetrigfrac) && centrality_low[0] == 0.0) {
	centrality_low[0] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[1] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (20./forcetrigfrac) && centrality_low[1] == 0.0) {
	centrality_low[1] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[2] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (30./forcetrigfrac) && centrality_low[2] == 0.0) {
	centrality_low[2] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[3] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (40./forcetrigfrac) && centrality_low[3] == 0.0) {
	centrality_low[3] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[4] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (50./forcetrigfrac) && centrality_low[4] == 0.0) {
	centrality_low[4] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[5] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (60./forcetrigfrac) && centrality_low[5] == 0.0) {
	centrality_low[5] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[6] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (70./forcetrigfrac) && centrality_low[6] == 0.0) {
	centrality_low[6] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[7] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (80./forcetrigfrac) && centrality_low[7] == 0.0) {
	centrality_low[7] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[8] = hSimBBCwTrig->GetBinCenter(i);
      }
    }

    if (maxcentbins ==19) {
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (5./forcetrigfrac) && centrality_low[0] == 0.0) {
	centrality_low[0] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[1] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (10./forcetrigfrac) && centrality_low[1] == 0.0) {
	centrality_low[1] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[2] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (15./forcetrigfrac) && centrality_low[2] == 0.0) {
	centrality_low[2] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[3] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (20./forcetrigfrac) && centrality_low[3] == 0.0) {
	centrality_low[3] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[4] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (25./forcetrigfrac) && centrality_low[4] == 0.0) {
	centrality_low[4] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[5] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (30./forcetrigfrac) && centrality_low[5] == 0.0) {
	centrality_low[5] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[6] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (35./forcetrigfrac) && centrality_low[6] == 0.0) {
	centrality_low[6] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[7] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (40./forcetrigfrac) && centrality_low[7] == 0.0) {
	centrality_low[7] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[8] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (45./forcetrigfrac) && centrality_low[8] == 0.0) {
	centrality_low[8] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[9] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (50./forcetrigfrac) && centrality_low[9] == 0.0) {
	centrality_low[9] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[10] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (55./forcetrigfrac) && centrality_low[10] == 0.0) {
	centrality_low[10] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[11] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (60./forcetrigfrac) && centrality_low[11] == 0.0) {
	centrality_low[11] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[12] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (65./forcetrigfrac) && centrality_low[12] == 0.0) {
	centrality_low[12] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[13] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (70./forcetrigfrac) && centrality_low[13] == 0.0) {
	centrality_low[13] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[14] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (75./forcetrigfrac) && centrality_low[14] == 0.0) {
	centrality_low[14] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[15] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (80./forcetrigfrac) && centrality_low[15] == 0.0) {
	centrality_low[15] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[16] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (85./forcetrigfrac) && centrality_low[16] == 0.0) {
	centrality_low[17] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[18] = hSimBBCwTrig->GetBinCenter(i);
      }
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > (90./forcetrigfrac) && centrality_low[17] == 0.0) {
	centrality_low[18] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[19] = hSimBBCwTrig->GetBinCenter(i);
      }
    }

  }

/*
  cout << "Centrality  " << centrality_low[0] << " " << centrality_low[1] << " " <<
    centrality_low[2] << " " << centrality_low[3] << endl;
  cout << "Centrality  " << centrality_high[0] << " " << centrality_high[1] << " " <<
    centrality_high[2] << " " << centrality_high[3] << endl;
*/

  TH1D *hslice = new TH1D("hslice","hslice",nhistbins,-0.5,maxrange); 
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
  
  c1->SaveAs(Form("glauber_fit_%d.pdf", runnumber));
  c1->SaveAs(Form("glauber_fit_%d.png", runnumber));
  //===============================================================================
  
  TCanvas *c3 = new TCanvas("canvas_bbchits_ncoll","canvas_bbchits_ncoll");
  c3->Divide(2,1);
  c3->cd(1);
  gPad->SetTicks(1);	
  c3->SetLogz(1);
  gStyle->SetPalette(1,0);
  hSimBBCNcoll->SetXTitle("Simulated BBC Hits");
  hSimBBCNcoll->SetYTitle("Number of Binary Collisions");
  hSimBBCNcoll->GetYaxis()->SetRangeUser(0.0,maxrange);
  hSimBBCNcoll->Draw("colz");

  c3->cd(2);
  c3->SetLogz(1);
  hSimBBCNcollHard->SetXTitle("Simulated BBC Hits");
  hSimBBCNcollHard->SetYTitle("Number of Binary Collisions");
  hSimBBCNcollHard->GetYaxis()->SetRangeUser(0.0,40.0);
  hSimBBCNcollHard->Draw("colz");

  // step through the Ncoll = 1 case and calculate the trigger efficiency
  double trigger_eff_ncoll1 = 0.0;
  double trigger_eff_ncoll1_weight = 0.0;
  int keyindex = 2; // index = 2 --> Ncoll=1
  if (npartscaling) keyindex = 3;
  if (npartscaling && npartminusone) keyindex = 2;
  for (int i=1;i<=nhistbins;i++) trigger_eff_ncoll1_weight += hSimBBCNcoll->GetBinContent(i,keyindex);
  for (int i=1;i<nhistbins;i++) {
    if (i<keyindex) {
      trigger_eff_ncoll1 += 0.0; // no chance to fire the trigger with BBCsouth hits = 0
    } else {
      trigger_eff_ncoll1 += hSimBBCNcoll->GetBinContent(i,keyindex) * trigeffcurve->Eval(hSimBBC->GetBinCenter(i));
    }
  }
  cout << "Trigger Efficiency for Ncoll=1 Case = " << trigger_eff_ncoll1/trigger_eff_ncoll1_weight << endl;
  cout << "The above should have eff curve match to 50 percent" << endl;

  // now what is the trigger efficiency for Ncoll=1 Case for Hard Collisions 
  TF1 *trigeffcurvehard = new TF1("trigeffcurvehard","1.0-TMath::Exp(-pow((x/[0]), [1]))",0.0,maxrange);
  trigeffcurvehard->SetParameters(0.235,trigeffcurve->GetParameter(1)); 

  // the above is just set by hand for now !!!
  // one could vary the first parameter (for example) until one gets a 75% match... by
  // the code below for checking.......
  double pizerotrigeff = 0.75;
  pizerotrigeff = 0.59; // new value for auau with BBCLL1>1
  
  double bestpar1 = 0.;
  double besteff = 99999.;
  for (int ipar1=0;ipar1<300;ipar1++) {

    trigeffcurvehard->SetParameter(0, 0.0 + 0.02*(double) ipar1);
    trigger_eff_ncoll1 = 0.0;
    trigger_eff_ncoll1_weight = 0.0;

    for (int i=1;i<=nhistbins;i++) {
      trigger_eff_ncoll1_weight += hSimBBCNcoll->GetBinContent(i,keyindex); 
    }

    for (int i=1;i<=nhistbins;i++) {
      if (i<keyindex) {
	trigger_eff_ncoll1 += 0.0; // no chance to fire trigger with no BBC south hit
      } else {
	trigger_eff_ncoll1 += hSimBBCNcollHard->GetBinContent(i,keyindex) * trigeffcurvehard->Eval(hSimBBC->GetBinCenter(i));
      }
    }

    double trigeff = (trigger_eff_ncoll1/trigger_eff_ncoll1_weight);
    if (fabs((float) (trigeff-pizerotrigeff)) < fabs((float) (besteff-pizerotrigeff))) {
      bestpar1 = trigeffcurvehard->GetParameter(0);
      besteff = trigeff;
    }

  }

  cout << "Trigger Efficiency for Ncoll=1 [Hard] Best Case = " << besteff << " with par: " << bestpar1 << endl;
  cout << "The above should have eff curve match to " << pizerotrigeff << " percent" << endl;

  trigeffcurvehard->SetParameter(0,bestpar1);
  c1->cd(2);
  trigeffcurvehard->SetLineColor(2);
  trigeffcurvehard->SetLineStyle(2);
  trigeffcurvehard->Draw("same");

  //============================================================================
  // bias results below  

  double biascorrections[16] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double minbiascorrection = 0.0;

  if (flag_determineBiasFactors) {
    TCanvas *c2 = new TCanvas("canvas_biasratio","canvas_biasratio");
    c2->Divide(1,2);
    c2->cd(1);
    hSimBBCHardUnbiased->DrawCopy();
    hSimBBCHardBiased->SetLineColor(4);
    hSimBBCHardBiased->DrawCopy("same");
    
    TH1D *hSimBBCHardRatio = new TH1D("hSimBBCHardRatio","hSimBBCHardRatio",nhistbins,-0.5,maxrange);
    for (int i=1;i<=nhistbins;i++) {
      if (hSimBBCHardBiased->GetBinContent(i) > 0 && hSimBBCHardUnbiased->GetBinContent(i) > 0) {
	hSimBBCHardRatio->SetBinContent(i, hSimBBCHardBiased->GetBinContent(i)/hSimBBCHardUnbiased->GetBinContent(i));
      }
    }
    c2->cd(2);
    hSimBBCHardRatio->DrawCopy();

    // take the integral within each nominal centrality category of counts in
    // hSimBBCHardBiased and hSimBBCHardUnbiased
    // and also fold the trigger efficiencies...
    
    double counts_unbiased[16]= {0.,0.,0.,0.,0.,0.,0.,0.};
    double counts_biased[16]= {0.,0.,0.,0.,0.,0.,0.,0.};
    double counts_unbiased_wtrig[16]= {0.,0.,0.,0.,0.,0.,0.,0.};
    double counts_biased_wtrig[16] = {0.,0.,0.,0.,0.,0.,0.,0.};
    
    for (int icent = 0; icent<maxcentbins; icent++) {

      // loop over BBC hits
      for (int i = 1+(int)centrality_low[icent]; i< 1+(int)centrality_high[icent];i++) {

	counts_unbiased[icent] += hSimBBCHardUnbiased->GetBinContent(i);
	counts_biased[icent] += hSimBBCHardBiased->GetBinContent(i);

	if (i==1) {
	  counts_unbiased_wtrig[icent] += 0.0;  // trig eff is zero if no BBC south hits
	} else {
	  counts_unbiased_wtrig[icent] += hSimBBCHardUnbiased->GetBinContent(i) * 
	    trigeffcurve->Eval(hSimBBCHardUnbiased->GetBinCenter(i));
	}
	if (i==1) {
	  counts_biased_wtrig[icent] += 0.0; // trig eff is zero if no BBC south hits
	} else {
	  counts_biased_wtrig[icent] += hSimBBCHardBiased->GetBinContent(i) * 
	    trigeffcurvehard->Eval(hSimBBCHardBiased->GetBinCenter(i));
	}
      }
      
      // print out bias factors
      cout << "Bias correction factors:::  " << counts_unbiased[icent]/counts_biased[icent] << " and w/trig " << 
	counts_unbiased_wtrig[icent]/counts_biased_wtrig[icent] << " and just trig " <<
	(counts_unbiased_wtrig[icent]/counts_biased_wtrig[icent]) / (counts_unbiased[icent]/counts_biased[icent]) << endl;

      biascorrections[icent] = counts_unbiased_wtrig[icent]/counts_biased_wtrig[icent];
      
      // plot the distribution of Ncoll for each centrality bin !!!!

    } // end loop over centrality bins

    // now calculate the bias factor for correcting to 0-100%
    double counts_mb_unbiased = 0.;
    double counts_mb_biased = 0.;
    for (int i = 1; i<= nhistbins; i++) {
      counts_mb_unbiased += hSimBBCHardUnbiased->GetBinContent(i);
      if (i==1) {
	counts_mb_biased += 0.0; // if no hits in BBCsouth then no possibility to fire the trigger
      } else {
	counts_mb_biased += hSimBBCHardBiased->GetBinContent(i) * 
	  trigeffcurvehard->Eval(hSimBBCHardBiased->GetBinCenter(i));
      }
    }
    
    cout << "Min. Bias 0-100% correction factor = " << 
      (1.0/(counts_mb_biased/counts_mb_unbiased))/(1.0/0.93) << endl;

    minbiascorrection = (1.0/(counts_mb_biased/counts_mb_unbiased))/(1.0/0.93);

  } // end flag on determining the bias factors

  double ncollstore[100];

  // now calculate the plot ncoll distribution
  TCanvas *c4 = new TCanvas("canvas_ncolldist","canvas_ncolldist");
  c4->cd();
  // now calculate the <nbinary> for each bin and the distribution
  TH1D *hCentNcollDist = new TH1D("hCentNcollDist","hCentNcollDist",nncollmax,-0.5,maxrangencoll);
  hCentNcollDist->SetMaximum(1.0);
  hCentNcollDist->GetXaxis()->SetRangeUser(0.0,nncollmax);
  hCentNcollDist->SetXTitle("Number of Binary Collisions");
  if (npartscaling) hCentNcollDist->SetXTitle("Number of Participants");
  hCentNcollDist->DrawCopy();
  for (int icent=0;icent<maxcentbins;icent++) {
    hCentNcollDist->Reset();
    // only check the ibbc bins for this centrality class
    for (int ibbc = 1+(int)centrality_low[icent]; ibbc< 1+(int)centrality_high[icent];ibbc++) {
      for (int inc=1; inc<=nncollmax; inc++) {
	hCentNcollDist->Fill(-1.0+(double)inc,hSimBBCNcoll_wtrig->GetBinContent(ibbc,inc));
      }
    }
    cout << "Centrality bin: " << icent << " <Ncoll> = " << hCentNcollDist->GetMean() << endl;
    ncollstore[icent+1] = hCentNcollDist->GetMean();
    hCentNcollDist->Scale(1.0/hCentNcollDist->Integral());
    hCentNcollDist->SetLineColor(1+icent);
    hCentNcollDist->SetLineWidth(3);
    hCentNcollDist->DrawCopy("same");
  }

  // now do -100% case
  hCentNcollDist->Reset();
  // only check the ibbc bins for this centrality class
  for (int ibbc = 1; ibbc<=nhistbins;ibbc++) {
    for (int inc=1; inc<=nncollmax; inc++) {
      hCentNcollDist->Fill(-1.0+(double)inc,hSimBBCNcoll->GetBinContent(ibbc,inc));
    }
  }

  hCentNcollDist->Scale(1.0/hCentNcollDist->Integral());
  hCentNcollDist->SetLineColor(6);
  hCentNcollDist->SetLineStyle(2);
  // do not draw the minimum bias distribution (0-100% for now)... j.nagle 1/27/2013
  hCentNcollDist->DrawCopy("same");

  cout << "Centrality bin:  0-100% has <Ncoll> = " << hCentNcollDist->GetMean() << endl;
  ncollstore[0] = hCentNcollDist->GetMean();

  // some kind of TLegend on the previous plot
  char ncolltext[150];
  sprintf(ncolltext,"<Ncoll>= %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f",
	  ncollstore[9],ncollstore[8],ncollstore[7],ncollstore[6],ncollstore[5],
	  ncollstore[4],ncollstore[3],ncollstore[2],ncollstore[1]);
  TLegend *tlegncoll = new TLegend(0.1,0.7,0.9,0.9,ncolltext,"brNDC");
  tlegncoll->Draw("same");

  cout << "Final Summary Output - Best mu,k Trigeff, NcollMB, NcollCent..." << endl;
  cout << "FINAL " << bestmu << " " << bestk << " " << trigeffintegral << " " << ncollstore[0];
  for (int icent=0;icent<maxcentbins;icent++) cout << " " << ncollstore[icent+1];
  cout << " " << minbiascorrection;
  for (int icent=0;icent<maxcentbins;icent++) cout << " " << biascorrections[icent];
  cout << " " << endl;

  // j.nagle - 01/21/2013
  // add parsing over the gl TNtuple in the fglauber file
  // then use NBD parameterss and trigger curve to re-calculate <Ncoll> and also <Npart> and <Npart,targ,proj>

  // for a given Ncoll what is the best way to pick the BBC from NBD
  // dirt dumb - just make a histogram for NBD and then GetRandom for each ncoll (slow but simple)
  // double nbdvalue = NBD_getValue(ihit, mu * (double) (pow(ncoll,alpha)), k * (double) (pow(ncoll,alpha)));

  bool EXTRANPARTTEST = false;

  if (EXTRANPARTTEST) {

    TH1F *hncollperiph = new TH1F("hncollperiph","hncollperiph",100,-0.5,99.5);
    TH1F *hncollperipht = new TH1F("hncollperipht","hncollperipht",100,-0.5,99.5);

    TH1F *hnbd = new TH1F("hnbd","hnbd",nhistbins,-0.5,maxrange);
    for (int ibbc=0;ibbc<nhistbins;ibbc++) {
      float nbdreturn = NBD_getValue(ibbc, bestmu, bestk);
      //      cout << "hnbdcheck ibbc = " << ibbc << " bestmu,k = " << bestmu << " " << bestk << " nbd = " << nbdreturn << endl;
      hnbd->SetBinContent(ibbc+1, nbdreturn);
    }
    // this does not include the alpha variation (maybe not such a big deal???)
    TF1 *fflat = new TF1("fflat","1.0",0.0,1.0);
    TH1F *hbbc2 = new TH1F("hbbc2","hbbc2",nhistbins,-0.5,maxrange);
    
    int ecount[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double ncollcount[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double npartcount[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double npartTcount[16] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double npartPcount[16] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double e2count[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double e2bcount[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double e2sqrcount[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double e2bsqrcount[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double avefakebbc[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    // e2gaus
    // e2disk
    // e2disknbd
    // overlappoint
    // overlapgaus
    // overlapdisk
    // overlapdisknbd
    // rbarlappoint
    // rbarlapgaus
    // rbarlapdisk
    // rbarlapdisknbd
    double e2diskcount[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double e2disknbdcount[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double overlapcount[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double overlapgauscount[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double overlapdiskcount[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double overlapdisknbdcount[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double rbarcount[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double rbargauscount[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double rbardiskcount[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double rbardisknbdcount[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

    double e3count[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double e3bcount[16]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

    // histogram plotting option for N centrality bins (distributions of quantities)
    TH1F *haction[16];
    for (int ih=0;ih<=maxcentbins;ih++) {
      char fooaction[100];
      sprintf(fooaction,"haction%d",ih);
      haction[ih] = new TH1F(fooaction,fooaction,100,-0.5,99.5);
    }
    
    TTree *gl = (TTree *) fglauber->Get("gl"); 
    Double_t        b;
    Int_t           npartproj;
    Int_t           nparttarg;
    Int_t           ncoll;
    Double_t        newe2;
    Double_t        newe2b;
    Double_t        newe3;
    Double_t        newe3b;
    Double_t        newe2disk;
    Double_t        newe2disknbd;
    Double_t        newoverlap;
    Double_t        newoverlapb;
    Double_t        newoverlapdisk;
    Double_t        newoverlapdisknbd;
    Double_t        newrbar;
    Double_t        newrbarb;
    Double_t        newrbardisk;
    Double_t        newrbardisknbd;

    gl->SetBranchAddress("b",&b);
    gl->SetBranchAddress("npartproj",&npartproj);
    gl->SetBranchAddress("nparttarg",&nparttarg);
    gl->SetBranchAddress("ncoll",&ncoll);
    gl->SetBranchAddress("newe2",&newe2);
    gl->SetBranchAddress("newe2b",&newe2b);
    gl->SetBranchAddress("newe3",&newe3);
    gl->SetBranchAddress("newe3b",&newe3b);
    gl->SetBranchAddress("newe2disk",&newe2disk);
    gl->SetBranchAddress("newe2disknbd",&newe2disknbd);
    gl->SetBranchAddress("newoverlap",&newoverlap);
    gl->SetBranchAddress("newoverlapb",&newoverlapb);
    gl->SetBranchAddress("newoverlapdisk",&newoverlapdisk);
    gl->SetBranchAddress("newoverlapdisknbd",&newoverlapdisknbd);
    gl->SetBranchAddress("newrbar",&newrbar);
    gl->SetBranchAddress("newrbarb",&newrbarb);
    gl->SetBranchAddress("newrbardisk",&newrbardisk);
    gl->SetBranchAddress("newrbardisknbd",&newrbardisknbd);

    Long64_t nentries = gl->GetEntries();
    Long64_t nbytes = 0;
    // loop over 10000 Glauber entries maximum
    if (nentries > 9000) nentries = 9000;
    for (Long64_t i=0; i<nentries;i++) {
      nbytes += gl->GetEntry(i);
      // determine fake BBC charge
      double fakebbc = 0.0;

      // IF WE WANT INSTEAD OF NCOLL (+) NBD, NPARTTARG (+) NBD NEED TO CHANGE THE PART BELOW...
      if (!npartscaling) {
	for (int in=0;in<ncoll; in++) fakebbc += hnbd->GetRandom();
      } else {
	// Cu+Au need both nparttarg+npartproj, also in Au+Au case
	for (int in=0;in<nparttarg+npartproj; in++) fakebbc += hnbd->GetRandom();
      }
      // did the trigger fire
      bool   trigfired = false;
      if (fflat->GetRandom() < trigeffcurve->Eval(fakebbc)) trigfired = true;
      
      hbbc2->Fill(fakebbc);
      
      // what centrality bin is this...
      int centbin = 8;

      if (centraltrigger) {
	//======================================================================================================
	// J.NAGLE - 11/20/2014 - COULD HACK SOMETHING IN HERE TO SIMULATION CENTRAL TRIGGER HE3AU TURNON...
	//======================================================================================================
	TF1 *erf = new TF1("erf","[2]*(0.5+0.5*ROOT::Math::erf((x-[0])/[1]))",0,200);
	//	erf->SetParameters(1.0168e2,1.9512e1,9.2934e-1);
	// alternate set
	//	erf->SetParameters(1.0337e2,2.06233e1,9.9e-1);

	// p+Au values from one low-luminosity run 434824
	erf->SetParameters(5.197e1,1.315e1,0.99);
	TF1 *erfflat = new TF1("erfflat","1.0",0.0,1.0);
	
	if (erfflat->GetRandom() < erf->Eval(fakebbc)) // then central event !!!
	  centbin = 0;
	
	erf->Delete();
	erfflat->Delete();

      } else {
	if (fakebbc >= centrality_low[0] && fakebbc <= centrality_high[0]) centbin = 0;
      }
      if (fakebbc >= centrality_low[1] && fakebbc <= centrality_high[1]) centbin = 1;
      if (fakebbc >= centrality_low[2] && fakebbc <= centrality_high[2]) centbin = 2;
      if (fakebbc >= centrality_low[3] && fakebbc <= centrality_high[3]) centbin = 3;
      if (fakebbc >= centrality_low[4] && fakebbc <= centrality_high[4]) centbin = 4;
      if (fakebbc >= centrality_low[5] && fakebbc <= centrality_high[5]) centbin = 5;
      if (fakebbc >= centrality_low[6] && fakebbc <= centrality_high[6]) centbin = 6;
      if (fakebbc >= centrality_low[7] && fakebbc <= centrality_high[7]) centbin = 7;
      if (fakebbc >= centrality_low[8] && fakebbc <= centrality_high[8]) centbin = 8;
      
      if (centbin == 8) hncollperiph->Fill(ncoll);

      ecount[maxcentbins]++;  // minbias event 0-100%
      ncollcount[maxcentbins] += (double) ncoll;
      npartcount[maxcentbins] += (double) nparttarg + (double) npartproj;
      npartTcount[maxcentbins] += (double) nparttarg;
      npartPcount[maxcentbins] += (double) npartproj;
      e2count[maxcentbins] += (double) newe2;
      e2bcount[maxcentbins] += (double) newe2b;
      e2sqrcount[maxcentbins] += (double) (newe2*newe2);
      e2bsqrcount[maxcentbins] += (double) (newe2b*newe2b);
      avefakebbc[maxcentbins] += (double) (fakebbc);
      e2diskcount[maxcentbins] += (double) newe2disk;
      if (newe2disknbd>0.0 && newe2disknbd < 1000.0) e2disknbdcount[maxcentbins] += (double) newe2disknbd;
      overlapcount[maxcentbins] += (double) newoverlap;
      overlapgauscount[maxcentbins] += (double) newoverlapb;
      overlapdiskcount[maxcentbins] += (double) newoverlapdisk;
      if (newoverlapdisknbd>0.0 && newoverlapdisknbd<1000.0) overlapdisknbdcount[maxcentbins] += (double) newoverlapdisknbd;
      rbarcount[maxcentbins] += (double) newrbar;
      rbargauscount[maxcentbins] += (double) newrbarb;
      rbardiskcount[maxcentbins] += (double) newrbardisk;
      if (newrbardisknbd>0.0 && newrbardisknbd<1000.0)       rbardisknbdcount[maxcentbins] += (double) newrbardisknbd;

      e3count[maxcentbins] += (double) newe3;
      e3bcount[maxcentbins] += (double) newe3b;


      haction[maxcentbins]->Fill(nparttarg);
 
      if (trigfired) {

	if (centbin == 8) 	hncollperipht->Fill(ncoll);

	ecount[centbin]++; 
	ncollcount[centbin] += (double) ncoll;
	npartcount[centbin] += (double) nparttarg + (double) npartproj;
	npartTcount[centbin] += (double) nparttarg;
	npartPcount[centbin] += (double) npartproj;
	e2count[centbin] += (double) newe2;
	e2bcount[centbin] += (double) newe2b;
	e2sqrcount[centbin] += (double) (newe2*newe2);
	e2bsqrcount[centbin] += (double) (newe2b*newe2b);
	avefakebbc[centbin] += (double) (fakebbc);

	e2diskcount[centbin] += (double) newe2disk;
	if (newe2disknbd>0.0 && newe2disknbd < 1000.0) 	e2disknbdcount[centbin] += (double) newe2disknbd;
	overlapcount[centbin] += (double) newoverlap;
	overlapgauscount[centbin] += (double) newoverlapb;
	overlapdiskcount[centbin] += (double) newoverlapdisk;
	if (newoverlapdisknbd>0.0 && newoverlapdisknbd<1000.0) overlapdisknbdcount[centbin] += (double) newoverlapdisknbd;
	rbarcount[centbin] += (double) newrbar;
	rbargauscount[centbin] += (double) newrbarb;
	rbardiskcount[centbin] += (double) newrbardisk;
	if (newrbardisknbd>0.0 && newrbardisknbd<1000.0) rbardisknbdcount[centbin] += (double) newrbardisknbd;

	e3count[centbin] += (double) newe3;
	e3bcount[centbin] += (double) newe3b;

	haction[centbin]->Fill(nparttarg);
      }
      
    } // end gl loop

    TCanvas *foobar = new TCanvas("foobar","foobar");
    hncollperiph->Draw();
    hncollperipht->SetLineColor(2);
    hncollperipht->Draw("same");

    cout << "Checker = " << hncollperiph->GetMean() << " " << hncollperipht->GetMean() << endl;
    
    // extra checks
    cout << "Trigger fraction = " << ((double) (ecount[0]+ecount[1]+ecount[2]+ecount[3]+ecount[4]+ecount[5]+ecount[6]+ecount[7]+ecount[8])) / ((double) ecount[maxcentbins]) << endl;
    
    TCanvas *caction = new TCanvas();
    caction->cd();
    haction[0]->DrawCopy();
    for (int ih=0;ih<=maxcentbins;ih++) {
      haction[ih]->SetLineColor(ih+1);
      haction[ih]->DrawCopy("same");
    }

    TCanvas *ccc = new TCanvas("ccc","ccc");
    ccc->Divide(1,2);
    ccc->cd(1);
    hSimBBC->Draw();
    hbbc2->Scale(hSimBBC->Integral() / hbbc2->Integral());
    hbbc2->SetLineColor(3);
    hbbc2->SetLineWidth(4);
    hbbc2->DrawCopy("same");
    ccc->cd(2);
    hnbd->DrawCopy();
    cout << "hnbd mean = " << hnbd->GetMean() << endl;
    double nbdmean = 0.0;
    for (int i=1;i<=nhistbins;i++) nbdmean += hnbd->GetBinCenter(i) * hnbd->GetBinContent(i);
    cout << "hnbd mean full calc = " << nbdmean / (hnbd->Integral()) << endl;

    // now cout the final information to parse....
    //------------------
    cout << "FINALPARSENCOLL " << bestmu << " " << bestk << " " << trigeffintegral << " " << ncollcount[maxcentbins]/(double)ecount[maxcentbins];
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << ncollcount[icent]/(double)ecount[icent];
    cout << " " << minbiascorrection;
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << biascorrections[icent];
    cout << " " << endl;

    cout << "FINALPARSE_FAKEBBC " << bestmu << " " << bestk << " " << trigeffintegral << " " << (avefakebbc[maxcentbins]/(double)ecount[maxcentbins]);
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << (avefakebbc[icent]/(double)ecount[icent]);
    cout << " " << minbiascorrection;
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << biascorrections[icent];
    cout << " " << endl;
    
    cout << "FINALPARSENPART " << bestmu << " " << bestk << " " << trigeffintegral << " " << npartcount[maxcentbins]/(double)ecount[maxcentbins];
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << npartcount[icent]/(double)ecount[icent];
    cout << " " << minbiascorrection;
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << biascorrections[icent];
    cout << " " << endl;
    
    cout << "FINALPARSENPARTP " << bestmu << " " << bestk << " " << trigeffintegral << " " << npartPcount[maxcentbins]/(double)ecount[maxcentbins];
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << npartPcount[icent]/(double)ecount[icent];
    cout << " " << minbiascorrection;
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << biascorrections[icent];
    cout << " " << endl;
    
    cout << "FINALPARSENPARTT " << bestmu << " " << bestk << " " << trigeffintegral << " " << npartTcount[maxcentbins]/(double)ecount[maxcentbins];
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << npartTcount[icent]/(double)ecount[icent];
    cout << " " << minbiascorrection;
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << biascorrections[icent];
    cout << " " << endl;

    //------------------
   
    cout << "FINALPARSE_E2 " << bestmu << " " << bestk << " " << trigeffintegral << " " << e2count[maxcentbins]/(double)ecount[maxcentbins];
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << e2count[icent]/(double)ecount[icent];
    cout << " " << minbiascorrection;
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << biascorrections[icent];
    cout << " " << endl;

    cout << "FINALPARSE_E2RMS " << bestmu << " " << bestk << " " << trigeffintegral << " " << sqrt(e2sqrcount[maxcentbins]/(double)ecount[maxcentbins]);
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << sqrt(e2sqrcount[icent]/(double)ecount[icent]);
    cout << " " << minbiascorrection;
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << biascorrections[icent];
    cout << " " << endl;

    cout << "FINALPARSE_E2B " << bestmu << " " << bestk << " " << trigeffintegral << " " << e2bcount[maxcentbins]/(double)ecount[maxcentbins];
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << e2bcount[icent]/(double)ecount[icent];
    cout << " " << minbiascorrection;
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << biascorrections[icent];
    cout << " " << endl;

    cout << "FINALPARSE_E2DISK " << bestmu << " " << bestk << " " << trigeffintegral << " " << e2diskcount[maxcentbins]/(double)ecount[maxcentbins];
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << e2diskcount[icent]/(double)ecount[icent];
    cout << " " << minbiascorrection;
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << biascorrections[icent];
    cout << " " << endl;

    cout << "FINALPARSE_E2DISKNBD " << bestmu << " " << bestk << " " << trigeffintegral << " " << e2disknbdcount[maxcentbins]/(double)ecount[maxcentbins];
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << e2disknbdcount[icent]/(double)ecount[icent];
    cout << " " << minbiascorrection;
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << biascorrections[icent];
    cout << " " << endl;

    //====================

    cout << "FINALPARSE_OVERLAP " << bestmu << " " << bestk << " " << trigeffintegral << " " << overlapcount[maxcentbins]/(double)ecount[maxcentbins];
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << overlapcount[icent]/(double)ecount[icent];
    cout << " " << minbiascorrection;
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << biascorrections[icent];
    cout << " " << endl;

    cout << "FINALPARSE_OVERLAPB " << bestmu << " " << bestk << " " << trigeffintegral << " " << overlapgauscount[maxcentbins]/(double)ecount[maxcentbins];
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << overlapgauscount[icent]/(double)ecount[icent];
    cout << " " << minbiascorrection;
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << biascorrections[icent];
    cout << " " << endl;

    cout << "FINALPARSE_OVERLAPDISK " << bestmu << " " << bestk << " " << trigeffintegral << " " << overlapdiskcount[maxcentbins]/(double)ecount[maxcentbins];
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << overlapdiskcount[icent]/(double)ecount[icent];
    cout << " " << minbiascorrection;
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << biascorrections[icent];
    cout << " " << endl;


    cout << "FINALPARSE_OVERLAPDISKNBD " << bestmu << " " << bestk << " " << trigeffintegral << " " << overlapdisknbdcount[maxcentbins]/(double)ecount[maxcentbins];
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << overlapdisknbdcount[icent]/(double)ecount[icent];
    cout << " " << minbiascorrection;
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << biascorrections[icent];
    cout << " " << endl;

    //====================

    cout << "FINALPARSE_RBAR " << bestmu << " " << bestk << " " << trigeffintegral << " " << rbarcount[maxcentbins]/(double)ecount[maxcentbins];
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << rbarcount[icent]/(double)ecount[icent];
    cout << " " << minbiascorrection;
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << biascorrections[icent];
    cout << " " << endl;

    cout << "FINALPARSE_RBARB " << bestmu << " " << bestk << " " << trigeffintegral << " " << rbargauscount[maxcentbins]/(double)ecount[maxcentbins];
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << rbargauscount[icent]/(double)ecount[icent];
    cout << " " << minbiascorrection;
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << biascorrections[icent];
    cout << " " << endl;

    cout << "FINALPARSE_RBARDISK " << bestmu << " " << bestk << " " << trigeffintegral << " " << rbardiskcount[maxcentbins]/(double)ecount[maxcentbins];
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << rbardiskcount[icent]/(double)ecount[icent];
    cout << " " << minbiascorrection;
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << biascorrections[icent];
    cout << " " << endl;


    cout << "FINALPARSE_RBARDISKNBD " << bestmu << " " << bestk << " " << trigeffintegral << " " << rbardisknbdcount[maxcentbins]/(double)ecount[maxcentbins];
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << rbardisknbdcount[icent]/(double)ecount[icent];
    cout << " " << minbiascorrection;
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << biascorrections[icent];
    cout << " " << endl;


    //==================================

    cout << "FINALPARSE_E3 " << bestmu << " " << bestk << " " << trigeffintegral << " " << e3count[maxcentbins]/(double)ecount[maxcentbins];
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << e3count[icent]/(double)ecount[icent];
    cout << " " << minbiascorrection;
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << biascorrections[icent];
    cout << " " << endl;

    cout << "FINALPARSE_E3B " << bestmu << " " << bestk << " " << trigeffintegral << " " << e3bcount[maxcentbins]/(double)ecount[maxcentbins];
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << e3bcount[icent]/(double)ecount[icent];
    cout << " " << minbiascorrection;
    for (int icent=0;icent<maxcentbins;icent++) cout << " " << biascorrections[icent];
    cout << " " << endl;
 
  } // end extra npart text with gl loop (also for eccentricities)


  
} // end routine

//====================================================================
// Return the value of the negative binomial distribution
// with parameters (mu, k) and for given n.
double NBD_getValue(int n, double mu, double k) {

  double F;
  double f;

  if (n+k > 100.0) {
    // log method for handling large numbers
    F  = TMath::LnGamma(n + k)- TMath::LnGamma(n + 1.)- TMath::LnGamma(k);
    f  = n * TMath::Log(mu/k) - (n + k) * TMath::Log(1.0 + mu/k);
    F = F+f;
    F = TMath::Exp(F);
  } else {
    F  = TMath::Gamma(n + k) / ( TMath::Gamma(n + 1.) * TMath::Gamma(k) );
    f  = n * TMath::Log(mu/k) - (n + k) * TMath::Log(1.0 + mu/k);
    f  = TMath::Exp(f);
    F *= f;
  }

  return F;  

}


