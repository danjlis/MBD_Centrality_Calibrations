/* QA_centrality */
// Author: Daniel Lis - August 15 - 2023
// Brief : This macro class gives a quick QA analysis
//         of a run in sPHENIX
//      To be run on the output of the CentralityReco Module

//void QA_FindCentralities()

#include "QA_centrality.h"

using namespace std;
namespace fs = std::filesystem;

static TH1D *hglauber{nullptr};

double NBD_getValue(double n, double mu, double k) {

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

static double NBDGlauberConv(double *x, double *par)
{

  double ihit = x[0];
  double mu = par[0];
  double k = par[1];

  double result = 0;
  for (int ib = 2; ib <= hglauber->GetNbinsX(); ib++)
    {
      int npart = hglauber->GetBinCenter(ib);
      double event_weight = hglauber->GetBinContent(ib);
      if (event_weight <= 0) continue;
      //      if (npart > 10 && ihit < (2*npart*(int) mu)) continue;
      double nbd = NBD_getValue(ihit, mu * (double) npart, k * (double) npart); 
      result += nbd*event_weight;
    }  
  return par[2]*result;
}

class QA_centrality;

/* constructor -- for the qa_info struct to hold all QA stats */
QA_centrality::QA_centrality(int silent, int debugger)
{

  silence  = silent;
  debug = debugger;

  /* Setting path for output files and input files */

  env_p = new char[200];
  sprintf(env_p,"%s",std::getenv("MBD_CENTRALITY_PATH"));
  
  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_PATH set."<<endl;
      return;
    }

  env_out = new char[200];
  sprintf(env_out,"%s",std::getenv("MBD_CENTRALITY_OUTPUT_PATH"));
  
  if(!env_out)
    {
      std::cout << "no env MBD_CENTRALITY_OUTPUT_PATH set."<<endl;
      return;
    }

  env_calib = new char[200];
  sprintf(env_calib,"%s",std::getenv("MBD_CENTRALITY_CALIB_PATH"));
  
  if(!env_calib)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }

  env_tree = new char[200];
  sprintf(env_tree,"%s",std::getenv("MBDTREELOC"));
  
  if(!env_tree)
    {
      std::cout << "no env MBDTREELOC set."<<endl;
      return;
    }

  /* struct initialization */

  qa_info.runnumber = 0;	  
  qa_info.quality = -1;
  qa_info.nevents = 0;
  qa_info.scale_factor = 0.;
  qa_info.ZDC_percentage = 0.;
  qa_info.one_ZDC_percentage = 0.;
  qa_info.min_bias = 0.;
  qa_info.vertex_cut = 0.;
  qa_info.vertex = 0.;
  qa_info.vertex_sigma = 0.;
  qa_info.charge_sum_mean_ns = 0.;
  for (int i = 0; i < 4;i++)
    {
      qa_info.charge_sum_mean[i] = 0.;
      qa_info.centrality_check_chi2[i] = 0.;
    }
  for (int i = 0; i < 4; i++)
    {
      qa_info.glauber_mu[i] = 0.;
      qa_info.glauber_k[i] = 0.;
      qa_info.glauber_chi2[i] = 0.;
      qa_info.glauber_npart[i] = 0.;
      qa_info.glauber_ecc2[i] = 0.;
      qa_info.glauber_ecc3[i] = 0.;
      qa_info.glauber_trig_eff[i] = 0.;
      qa_info.glauber_trig_eff_err[i] = 0.;
    }
  mc_map[0] = "hijing";
  mc_map[1] = "ampt";
  mc_map[2] = "epos";
  mc_map[3] = "hijing_magoff";
  mc_map[4] = "ampt_magoff";
  mc_map[5] = "epos_magoff";
    
}

void QA_centrality::Print_QA_Info(bool use_table)
{
  std::cout << " ***************** Run " << qa_info.runnumber << " ***************** " << std::endl;
  std::cout << "   Quality: "<<qa_info.quality<<std::endl; 
  std::cout << "   ZDC percent: "<<qa_info.ZDC_percentage<<" ("<<qa_info.one_ZDC_percentage<<")"<<std::endl; 
  std::cout << "   Vertex = " << qa_info.vertex << " +/- " <<qa_info.vertex_sigma << std::endl;
  std::cout << "   Charge Sum Mean: "<<std::endl;
  std::cout << "       Regular\tBalance\tScaled\tBoth " <<std::endl;
  std::cout << "       ";
  for (int i = 0; i < 4; i++) std::cout << qa_info.charge_sum_mean[i] << "\t";
  std::cout << " "<<std::endl;
  for (int i = 0; i < 4; i++) std::cout << qa_info.centrality_check_chi2[i] << "\t";
  std::cout << " "<<std::endl;
  std::cout << "   Charge Sum Scaling factor: "<< qa_info.scale_factor << std::endl;
  std::cout << "   Glauber Fit Status (scaled): "<<std::endl;
  std::cout << "       mu           =  " << qa_info.glauber_mu[0]<<" ( "<<qa_info.glauber_mu[2]<<" )"<<std::endl;
  std::cout << "       k            =  " << qa_info.glauber_k[0]<<" ( "<<qa_info.glauber_k[2]<<" )"<<std::endl;
  std::cout << "       trigger eff. =  " << qa_info.glauber_trig_eff[0]<<" ( "<<qa_info.glauber_trig_eff[2]<<" )"<<std::endl;
  std::cout << "       chi2         =  " << qa_info.glauber_chi2[0]<<" ( "<<qa_info.glauber_chi2[2]<<" )"<<std::endl;
  std::cout << "       npart         =  " << qa_info.glauber_npart[0]<<" ( "<<qa_info.glauber_npart[2]<<" )"<<std::endl;
  std::cout << "       ecc2         =  " << qa_info.glauber_ecc2[0]<<" ( "<<qa_info.glauber_ecc2[2]<<" )"<<std::endl;
  std::cout << "       ecc3         =  " << qa_info.glauber_ecc3[0]<<" ( "<<qa_info.glauber_ecc3[2]<<" )"<<std::endl;
  std::cout << "       b         =  " << qa_info.glauber_b[0]<<" ( "<<qa_info.glauber_b[2]<<" )"<<std::endl;
  std::cout << " ********************************************* " << std::endl;
  if (use_table)
    {
      // runnumber,events,ZDC,vtx,vtxsig,mean,mean_sca,scale,chi2,mu,k,trigeff,trigeff_sca_for,chi2scafor,centchi2,centchi2sca
      std::cout << " Table Entry: "<<std::endl;
      std::cout << std::fixed << std::setprecision(5);
      std::cout << qa_info.runnumber << ",";
      std::cout << qa_info.nevents << ",";
      std::cout << qa_info.ZDC_percentage << ",";
      std::cout << qa_info.one_ZDC_percentage << ",";
      std::cout << qa_info.min_bias << ",";
      std::cout << qa_info.vertex_cut << ",";
      std::cout << qa_info.vertex << ",";
      std::cout << qa_info.vertex_sigma << ",";
      std::cout << qa_info.fillpattern << ",";
      std::cout << qa_info.charge_sum_mean[0] << ",";
      std::cout << qa_info.charge_sum_mean[2] << ",";
      std::cout << qa_info.scale_factor << ",";
      std::cout << qa_info.centrality_check_chi2[0] << ",";
      std::cout << qa_info.centrality_check_chi2[2] << ",";
      std::cout << qa_info.charge_sum_mean_ns << ",";
      std::cout << qa_info.glauber_chi2[0] << ",";
      std::cout << qa_info.glauber_mu[0] << ",";
      std::cout << qa_info.glauber_k[0] << ",";
      std::cout << qa_info.glauber_trig_eff[0] << ",";
      std::cout << qa_info.glauber_trig_eff_err[0] << ",";
      std::cout << qa_info.glauber_npart[0] << ",";
      std::cout << qa_info.glauber_ecc2[0] << ",";
      std::cout << qa_info.glauber_ecc3[0] << ",";
      std::cout << qa_info.glauber_b[0] << ",";
      std::cout << qa_info.glauber_chi2[2] << ",";
      std::cout << qa_info.glauber_mu[2] << ",";
      std::cout << qa_info.glauber_k[2] << ",";
      std::cout << qa_info.glauber_trig_eff[2] << ",";
      std::cout << qa_info.glauber_trig_eff_err[2] << ",";
      std::cout << qa_info.glauber_npart[2] << ",";
      std::cout << qa_info.glauber_ecc2[2] << ",";
      std::cout << qa_info.glauber_ecc3[2] << ",";
      std::cout << qa_info.glauber_b[2];
      std::cout << ""<<std::endl;
    }
}

/* This is that one all QA centrality start command */

void QA_centrality::QA_ReferenceRun(
				    const int runnumber
				    )
{

  qa_info.runnumber = runnumber;
  // ZDC Check for percentage that meets certain cuts.
  if (!silence)
    {
      std::cout << " Going into ZDC Check" << std::endl;
    }

  if (loadTree())
    {
      return;
    }
    QA_ZDCCheck(runnumber);
  //Go through and make distributions for MBD channels (charge and time)
  QA_MBDChannels(runnumber);
  // Charge sum plots are made here with scaled and with all cuts
  QA_MakeChargeSum(runnumber);
  // Centrality Calibrations + NBD+Glauber Fit is here

  QA_MakeDivisions(runnumber, 0);  
  QA_MakeDivisions(runnumber, 1);  
  //QA_MakeDivisions(runnumber, 2);  

  QA_MakeCentralityCalibrations(runnumber, 0);
  QA_MakeCentralityCalibrations(runnumber, 1);

  //CentralityCalibrations(runnumber, 1);
  SetReferenceRun(runnumber);

  QA_CentralityCheck(runnumber, 0);
  QA_CentralityCheck(runnumber, 1);
  QA_CentralityCheck(runnumber, 2);
  return;
}

void QA_centrality::QA_RunCentralityCheck(
					  const int runnumber,
					  const int reference)
{

  qa_info.runnumber = runnumber;
  //CentralityCalibrations(runnumber, 1);
  if (loadTree())
    {
      return;
    }
  SetReferenceRun(reference);
  QA_CentralityCheck(runnumber,0);
  QA_CentralityCheck(runnumber,1);
  QA_CentralityCheck(runnumber,2);
  return;
}

void QA_centrality::QA_MC(
			  const std::string mc_generator
			  )
{

  isSim = true;
  if (strcmp(mc_generator.c_str(), "hijing") == 0)
    qa_info.runnumber = 0;
  
  else if (strcmp(mc_generator.c_str(), "ampt") == 0)
    qa_info.runnumber = 1;
  
  else if (strcmp(mc_generator.c_str(), "epos") == 0)
    qa_info.runnumber = 2;
  if (strcmp(mc_generator.c_str(), "hijing_magoff") == 0)
    qa_info.runnumber = 3;
  
  else if (strcmp(mc_generator.c_str(), "ampt_magoff") == 0)
    qa_info.runnumber = 4;
  
  else if (strcmp(mc_generator.c_str(), "epos_magoff") == 0)
    qa_info.runnumber = 5;
  else
    return;

  int runnumber = qa_info.runnumber;
  if (loadTree())
    {
      return;
    }

  std::cout << mc_map[runnumber] <<" : " <<runnumber<<std::endl;
  QA_ZDCCheck(runnumber);
  // // Go through and make distributions for MBD channels (charge and time)
  QA_MBDChannels(runnumber);
  // Charge sum plots are made here with scaled and with all cuts
  QA_MakeChargeSum(runnumber);
  // Centrality Calibrations + NBD+Glauber Fit is here

  QA_MakeDivisions(runnumber, 0);
  QA_MakeDivisions(runnumber, 1);

  QA_MakeCentralityCalibrations(runnumber, 0);
  QA_MakeCentralityCalibrations(runnumber, 1);

  QA_CentralityCheck(runnumber, 0);
  QA_CentralityCheck(runnumber, 1);
  return;
}

void QA_centrality::Start_QA_Centrality(
					const int runnumber
					)
{

  qa_info.runnumber = runnumber;
  if (loadTree())
    {
      return;
    }

  // ZDC Check for percentage that meets certain cuts.
  QA_ZDCCheck(runnumber);
  // Go through and make distributions for MBD channels (charge and time)
  QA_MBDChannels(runnumber);
  // Charge sum plots are made here with scaled and with all cuts
  QA_MakeChargeSum(runnumber);

  // Centrality Calibrations + NBD+Glauber Fit is here
  //  QA_MakeCentralityCalibrations(runnumber, 1);
  //  QA_MakeCentralityCalibrations(runnumber, 1, true, false);
  //  SetTrigEffMUK(.94, 4.39375, 0.859025);
  //  QA_MakeCentralityCalibrations(runnumber, 1, true, false);
  //  QA_MakeCentralityCalibrations(runnumber, 1, true, false);
  // Check the centrality bins

  QA_CentralityCheck(runnumber, 0);
  QA_CentralityCheck(runnumber, 1);
  //QA_CentralityCheck(runnumber, 2);

  return;
}

void QA_centrality::QA_CentralityCheck(const int runnumber, const int scaled)
{

  double scale_factor = 1;
  // Load in charge_sum distributions
  if (runnumber != reference_run && !isSim)
    {
      TFile *fchargesum = new TFile(Form("%s/mbdana_charge_sum_%d.root", env_out, runnumber), "r");
      if (!fchargesum)
	{
	  cout << "No file exists for this charge sum plot... exiting" << endl;
	  return;
	} 
      
      TFile *fchargesumref = new TFile(Form("%s/mbdana_charge_sum_%d.root", env_out, reference_run), "r");
      if (!fchargesumref)
	{
	  cout << "No file exists for this charge sum plot... exiting" << endl;
	  return;
	} 
      
      TH1D *h_run_charge_sum = (TH1D*) fchargesum->Get("h_charge_sum_min_bias_w_vertex_cut");
      if (scaled == 1)
	h_run_charge_sum = (TH1D*) fchargesum->Get("h_charge_sum_min_bias_w_vertex_cut_balanced");
      if (!h_run_charge_sum)
	{
	  std::cout << "nah fam" << std::endl;
	}
      TH1D *h_run_charge_sum_ref = (TH1D*) fchargesumref->Get("h_charge_sum_min_bias_w_vertex_cut");
      if (scaled == 1)
	h_run_charge_sum_ref = (TH1D*) fchargesumref->Get("h_charge_sum_min_bias_w_vertex_cut_balanced");
      if (!h_run_charge_sum_ref)
	{
	  std::cout << "nah fam no ref" << std::endl;
	}
      
      scale_factor = h_run_charge_sum_ref->GetMean()/h_run_charge_sum->GetMean();

      fchargesumref->Close();
      fchargesum->Close();
    }
  std::cout << " using scale factor : " << scale_factor << std::endl;
  if (debug) std::cout << __LINE__ << std::endl;

  char *path = new char[100];
  float vhigh;
  float vlow;
  float vscale;
  const int nvertexbins = 30;
  float vertex_ranges[nvertexbins+1] = {0};
  float vertex_scale_factor[nvertexbins] = {0};
  if (isSim)
    {
      sprintf(path, "%s/mbd_vertex_scale_%s.root" ,env_calib, mc_map[runnumber].c_str());
    }
  else
    {
      sprintf(path, "%s/mbd_vertex_scale_%d.root" ,env_calib, runnumber);
    }

  TFile *fcalibs = new TFile(path, "r");
  TNtuple *tss = (TNtuple*) fcalibs->Get("tn_vertexscale");
  tss->SetBranchAddress("low_vertex", &vlow);
  tss->SetBranchAddress("high_vertex", &vhigh);
  tss->SetBranchAddress("scale", &vscale);
  for (int i = 0; i < tss->GetEntries() ; i++)
    {
      tss->GetEntry(i);      
      vertex_ranges[i] = vlow;
      vertex_ranges[i+1] = vhigh;
      vertex_scale_factor[i] = (scaled == 1? vscale : 1.0);
      if (!silence) std::cout  << " vscale : " << vscale << std::endl;
    }

  fcalibs->Close();


  // Making fresh histograms
  // Load in centrality calibrations
  char pathfile[1000];
  if (!isSim)
    {
      sprintf(pathfile, "%s/mbdana_centrality%s%d.root", env_calib, (scaled? "_bal_":"_"), reference_run);
    }
  else
    {
      sprintf(pathfile, "%s/mbdana_centrality%s%s.root", env_calib, (scaled? "_bal_":"_"), mc_map[runnumber].c_str());
    }
  char pathvtxfile[1000];
  if (!isSim)
    {
      sprintf(pathvtxfile, "%s/mbdana_centrality_vtx%s%d.root", env_calib, (scaled? "_bal_":"_"), reference_run);
    }
  else
    {
      sprintf(pathvtxfile, "%s/mbdana_centrality%s%s.root", env_calib, (scaled? "_bal_":"_"), mc_map[runnumber].c_str());
    }

  TFile *fcentralitycalib = new TFile(pathfile, "r");
  if (!fcentralitycalib)
    {
      cout << "No calib file exists for the centrality in run " << reference_run << ".. exiting"<<endl;
      return;
    }
  
  float centbin, low, high;
  TNtuple *tn_centrality = (TNtuple*) fcentralitycalib->Get("tn_centrality");
  if (!tn_centrality)
    {
      cout << "No calib tntuple exists for the centrality in run " << reference_run << ".. exiting"<<endl;
      return;
    }
  float centrality_bins[100];
  for (int i = 0 ; i < 100; i++)
    centrality_bins[i] = 0.;
  
  tn_centrality->SetBranchAddress("bin", &centbin);
  tn_centrality->SetBranchAddress("low", &low);
  tn_centrality->SetBranchAddress("high", &high);
  int cent_divs1 = tn_centrality->GetEntries();
  for (int i = 0; i < cent_divs1; i++)
    {
      tn_centrality->GetEntry(i);
      if (!silence) std::cout  <<" Calib Entry " << i <<" : "<<  low <<  " -- " << high << std::endl;
      centrality_bins[i] = low;
    }
  

  // TFile *fcentralitycalibvtx = new TFile(pathvtxfile, "r");

  // if (!fcentralitycalibvtx)
  //   {
  //     cout << "No calib file exists for the centrality in run " << reference_run << ".. exiting"<<endl;
  //     return;
  //   }
  // float centrality_bins_vtx[nvertexbins][100];
  // for (int iv =0 ; iv < nvertexbins; iv++)
  //   {
  //     TNtuple *tn_centralityvtx = (TNtuple*) fcentralitycalibvtx->Get(Form("tn_centrality_%d", iv));
  //     if (!tn_centralityvtx)
  // 	{
  // 	  cout << "No calib tntuple exists for the centrality in run " << reference_run << ".. exiting"<<endl;
  // 	  return;
  // 	}
      
  //     tn_centralityvtx->SetBranchAddress("bin", &centbin);
  //     tn_centralityvtx->SetBranchAddress("low", &low);
  //     tn_centralityvtx->SetBranchAddress("high", &high);
  

  //     for (int i = 0 ; i < 100; i++)
  // 	centrality_bins_vtx[iv][i] = 0.;

  // // filling centrality bins

  //     int cent_divs = tn_centralityvtx->GetEntries();
  //     for (int i = 0; i < cent_divs; i++)
  // 	{
  // 	  tn_centralityvtx->GetEntry(i);
  // 	  if (!silence) std::cout  <<" Calib Entry " << i <<" : "<<  low <<  " -- " << high << std::endl;
  // 	  centrality_bins_vtx[iv][i] = low;
  // 	}
  //   }
  // fcentralitycalibvtx->Close();
  fcentralitycalib->Close();

  const int nvertexbinscheck = 17;
  double vertex_ranges_check[nvertexbinscheck+1] = {-60, -50, -40, -30, -20, -15, -10, -5, -2, 2, 5, 10, 15, 20, 30, 40, 50, 60};
  TH1D *h_cent_bin = new TH1D("hcent_bins","", ndivs, -0.5, (float) ndivs - 0.5);
  TH1D *h_cent_bin_sca = new TH1D("hcent_bins_sca","", ndivs, -0.5, (float) ndivs - 0.5);
  TH1D *h_vertex = new TH1D("h_vertex", ";Vertex [cm];Fraction in CentBin", nvertexbinscheck, vertex_ranges_check);
  TH1D *h_cent_vtx[20];
  TH1D *h_cent_vertex[20];
  TH1D *h_cent_vtx_sca[20];

  TH2D *h_cent_running = new TH2D("h_cent_running", ";Event/100000; centbin", 100, 0, 100, 100, 0, 100);
  for (int i = 0; i < 20; i++)
    {
      h_cent_vertex[i] = new TH1D(Form("h_cent_vertex_%d", i), ";Vertex [cm]", 120, -60, 60);
      h_cent_vtx[i] = new TH1D(Form("h_cent_vtx_%d", i), ";Vertex [cm];Fraction in CentBin", nvertexbinscheck, vertex_ranges_check);
      h_cent_vtx_sca[i] = new TH1D(Form("h_cent_vtx_sca_%d", i), ";Vertex [cm];Fraction in CentBin", nvertexbinscheck, vertex_ranges_check);
    }

  std::vector<std::pair<float, float>> v_mbd_charge_sum_scale{};
  std::vector<std::pair<float, float>> v_mbd_charge_sum{};

  //// vertex
  int hits_n = 0;
  int hits_s = 0;      
  float nsum = 0;
  float ssum = 0;
  int hits_n2 = 0;
  int hits_s2 = 0;      
  float nsum2 = 0;
  float ssum2 = 0;
  if (!silence)   cout <<"NEvents = "<<(nevents? nevents:m_ttree->GetEntries())<<std::endl;
  int nev = (nevents > 0 && nevents < m_ttree->GetEntries()?nevents:m_ttree->GetEntries());
  for (int i = 0 ; i < (debug ? 10 :nev); i++)
    {
      //      bool minbias = false;
      hits_n2 = 0;
      hits_s2 = 0;
      ssum2 = 0;
      nsum2 = 0;
      hits_n = 0;
      hits_s = 0;
      ssum = 0;
      nsum = 0;
      m_ttree->GetEntry(i);

      if (!(((gl1_scaled >> 14) & 0x1 ) == 0x1)) continue;

      double scale = 1;
      int vbin = -1;
      if (scaled)
	{
	  for (int iv = 0; iv < nvertexbins; iv++)
	    {
	      if (mbd_vertex >= vertex_ranges[iv] &&  mbd_vertex < vertex_ranges[iv+1])
		{
		  scale = vertex_scale_factor[iv];
		  vbin = iv;
		  break;
		}
	    }
	  if (vbin == -1 ) continue;
	}
      

      for (int ich = 0 ; ich < 64; ich++)
	{
	  if (countbefore)
	    {
	      ssum+=mbd_charge[ich]*scale*scale_factor;
	      nsum += mbd_charge[ich + 64]*scale*scale_factor;
	      ssum2+=mbd_charge[ich];
	      nsum2 += mbd_charge[ich + 64];
	    }
	  if (mbd_time[ich] < 25. && mbd_charge[ich] > cthresh)
	    {
	      hits_s++;
	      if (!countbefore) ssum+=mbd_charge[ich]*scale*scale_factor;
	    }
	  if (mbd_time[ich+64] < 25. && mbd_charge[ich+64] > cthresh)
	    {
	      hits_n++;
	      if (!countbefore) nsum += mbd_charge[ich + 64]*scale*scale_factor;
	    } 
	  if (mbd_time[ich] < 25. && mbd_charge[ich] > cthresh)
	    {
	      hits_s2++;
	      if (!countbefore) ssum2+=mbd_charge[ich];
	    }
	  if (mbd_time[ich+64] < 25. && mbd_charge[ich+64] > cthresh)
	    {
	      hits_n2++;
	      if (!countbefore) nsum2 += mbd_charge[ich + 64];
	    } 
	}
      //      if (!minbias) continue;
      bool mymb = true;
      if (hits_s < 2 || hits_n < 2) mymb = false;//continue; 
      if (fabs(mbd_vertex) > z_cut) mymb = false;//continue; 
      if (!isSim && ssum > charge_sum_south_cut && nsum < charge_sum_north_cut) mymb = false;//continue; 
      if (hasZDC && !isSim && (zdc_sum[0] <= zdc_cut || zdc_sum[1] <= zdc_cut)) mymb = false;//continue; 
      if ((ssum + nsum) > m_maxsumcut) mymb = false;//continue; 
      if (!mymb) continue;

      int centbin = 0;
      double totalsum = (nsum + ssum);
      for (int ic = 0; ic < 100; ic++)
	{
	  if (scaled < 2)
	    {
	      if (totalsum < centrality_bins[ic])
		continue;
	    }
	  else
	    {
	      if (totalsum < centrality_bins[ic])
		continue;
	    }
	  centbin = ic;
	  break;
	}

      int cent5bin = centbin/5;
      //h_cent_running->Fill(i/100000, centbin);      
      h_cent_bin->Fill(centbin);
      h_cent_vertex[cent5bin]->Fill(mbd_vertex);
      h_vertex->Fill(mbd_vertex);
      h_cent_vtx[cent5bin]->Fill(mbd_vertex);
    }
  
  for (int i = 0; i < 20; i++)
    {
      for (int ib = 1; ib <= h_vertex->GetNbinsX(); ib++)
	{
	  h_cent_vtx[i]->SetBinContent(ib, h_cent_vtx[i]->GetBinContent(ib)/h_vertex->GetBinContent(ib));
	  h_cent_vtx[i]->SetBinError(ib, h_cent_vtx[i]->GetBinError(ib)/h_vertex->GetBinContent(ib));
	}
    }

  TF1 *flatline = new TF1("flatline","[0]",-0.5, (float) ndivs - .5);
  float match = 1./((float)ndivs);
  flatline->SetParameter(0, match);

  h_cent_bin->Scale(1./(static_cast<float>(h_cent_bin->Integral())));
  h_cent_bin->Fit(flatline,"NDORQ","", -0.5, 93.5);

  float chi2 = flatline->GetChisquare()/flatline->GetNDF();
  float level = flatline->GetParameter(0);;


  if (!silence)
    {
      for (int i = 0; i < h_cent_bin->GetNbinsX();i++)
	{
	  cout << i*(100/ndivs) <<"% : "<<h_cent_bin->GetBinContent(i+1) <<endl;
	}
      cout << " ************************** " <<endl;
      cout << "  Run "<< runnumber<<endl;
      cout << "    Avg / Chi2         =  "<< level <<" / "<<chi2<<endl;
      cout << " ************************** " <<endl;
    }

  qa_info.centrality_check_chi2[2*(scaled?1:0) + (scaled/2)] = chi2;
  std::string xtra = "_";
  if (scaled == 1)
    {
      xtra = "_bal_";
    }
  if(scaled == 2)
    {
      xtra = "_vtx_";
    }
  std::string filename = Form("%s/mbdana_centrality_check%s%d.root", env_out,   xtra.c_str(), runnumber);
  if (isSim)
    {
      filename = Form("%s/mbdana_centrality_check%s%s.root", env_out,   xtra.c_str(), mc_map[runnumber].c_str());
    }
  TFile *fout = new TFile(filename.c_str(),"recreate");
  h_vertex->Write();

  h_cent_bin->Write();
  h_cent_bin_sca->Write();
  h_cent_running->Write();
  for (int i = 0; i < 20; i++)
    {
      h_cent_vtx_sca[i]->Write();
      h_cent_vtx[i]->Write();
      h_cent_vertex[i]->Write();

    }
  fout->Close();
  
  
  delete h_cent_bin;
  delete h_cent_bin_sca;
  delete h_vertex;
  delete h_cent_running;
  for (int i = 0; i < 20; i++)
    {
      delete h_cent_vertex[i];
      delete h_cent_vtx[i];
      delete h_cent_vtx_sca[i];
    }

  return;
}


void QA_centrality::QA_MakeChargeSum(const int runnumber)
{


  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }
  if(!env_tree)
    {
      std::cout << "no env MBDTREELOC set."<<endl;
      return;
    }
 
 
  // Get the mbd charge sum plot to shift to.
  float mbd_charge_factor = 1.;
  float mbd_charge_factor_bal = 1.;
  
  if (reference_run && !isSim)
    {
      if (!silence) cout << "reference_run: " <<reference_run<< std::endl;
      TFile *fout = new TFile(Form("%s/mbdana_charge_sum_%d.root", env_out, reference_run), "r");
      if (!fout) return;
      TH1D *h_shift = (TH1D*) fout->Get("h_charge_sum_min_bias_w_vertex_cut");
      if (!h_shift) return;
      h_shift->SetName("h_charge_sum_shift");
      TH1D *h_balshift = (TH1D*) fout->Get("h_charge_sum_min_bias_w_vertex_cut_balanced");
      if (!h_balshift) return;
      h_balshift->SetName("h_charge_sum_balshift");

      mbd_charge_factor = h_shift->GetMean();
      mbd_charge_factor_bal = h_balshift->GetMean();

      if (!silence) std::cout << "MBD Scale Factor: "<<mbd_charge_factor << " ( " << mbd_charge_factor_bal << " ) " <<std::endl;

      fout->Close();
      
    }
  
  const int nvertexbins = 10;
  double vertex_ranges[nvertexbins + 1] = {-60, -30, -20, -10, -5, -0, 5, 10, 20, 30, 60};//, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 20, 30, 60};

  int finefactor=4;
  int nbin = 2500;
  int maxrange = 2500;
  if (isSim)
    { 
      nbin = 2500;
      maxrange = 2500;
    }
  TH2D *h_mbd_hits_ns = new TH2D("h_mbd_hits_ns", ";nhit S; nhit N", 65, -0.5, 64.5, 65, -0.5, 64.5);

  TH2D *h2_charge_sum_v_zdc = new TH2D("h_charge_sum_v_zdc","", nbin/10, -0.5, (float)maxrange - 0.5, 2000, 0, 20000);
  TH2D *h2_charge_sum_v_zdc_mb = new TH2D("h_charge_sum_v_zdc_mb","", nbin/10, -0.5, (float)maxrange - 0.5, 1200, 0, 12000);

  TProfile *hp_charge_sum_v_zdc = new TProfile("hp_charge_sum_v_zdc","", 200, -0.5,  1999.5);
  TProfile *hp_charge_sum_v_zdc_mb = new TProfile("hp_charge_sum_v_zdc_mb","", 200, -0.5, 1999.5);

  TH2D *h2_charge_sum_ns = new TH2D("h2_charge_sum_ns","", nbin, -0.5, (float)maxrange - 0.5, nbin, -0.5, (float)maxrange - 0.5);  
  TH2D *h2_charge_sum_vtx_ns = new TH2D("h2_charge_sum_vtx_ns","", nbin, -0.5, (float)maxrange - 0.5, nbin, -0.5, (float)maxrange - 0.5);  
  TH2D *h2_charge_sum_min_bias_ns = new TH2D("h2_charge_sum_min_bias_ns","", nbin, -0.5, (float)maxrange - 0.5, nbin, -0.5, (float)maxrange - 0.5);  
  TH2D *h2_charge_sum_min_bias_w_vertex_cut_ns = new TH2D("h2_charge_sum_min_bias_w_vertex_cut_ns","", nbin, -0.5, (float)maxrange - 0.5, nbin, -0.5, (float)maxrange - 0.5);  
  TH2D *h2_charge_sum_ns_cut = new TH2D("h2_charge_sum_ns_cut","", nbin, -0.5, (float)maxrange - 0.5, nbin, -0.5, (float)maxrange - 0.5);  

  TH2D *h_charge_sum_v_vtx = new TH2D("h_charge_sum_v_vtx","", nvertexbins, vertex_ranges, nbin, -0.5, (float)maxrange - 0.5);
  TH2D *h_charge_sum_v_vtx_vtx = new TH2D("h_charge_sum_v_vtx_vtx","",  nvertexbins, vertex_ranges,nbin, -0.5, (float)maxrange - 0.5);
  TH2D *h_charge_sum_v_vtx_min_bias = new TH2D("h_charge_sum_v_vtx_min_bias","",nvertexbins, vertex_ranges, nbin, -0.5, (float)maxrange - 0.5);
  TH2D *h_charge_sum_v_vtx_min_bias_w_vertex_cut = new TH2D("h_charge_sum_v_vtx_min_bias_w_vertex_cut","",nvertexbins, vertex_ranges, nbin, -0.5, (float)maxrange - 0.5);

  TH2D *h_charge_sum_balanced_v_vtx = new TH2D("h_charge_sum_balanced_v_vtx","", nvertexbins, vertex_ranges, nbin, -0.5, (float)maxrange - 0.5);
  TH2D *h_charge_sum_balanced_v_vtx_vtx = new TH2D("h_charge_sum_balanced_v_vtx_vtx","",  nvertexbins, vertex_ranges,nbin, -0.5, (float)maxrange - 0.5);
  TH2D *h_charge_sum_balanced_v_vtx_min_bias = new TH2D("h_charge_sum_balanced_v_vtx_min_bias","",nvertexbins, vertex_ranges, nbin, -0.5, (float)maxrange - 0.5);
  TH2D *h_charge_sum_balanced_v_vtx_min_bias_w_vertex_cut = new TH2D("h_charge_sum_balanced_v_vtx_min_bias_w_vertex_cut","",nvertexbins, vertex_ranges, nbin, -0.5, (float)maxrange - 0.5);

  TH2D *h_charge_sum_north_v_vtx = new TH2D("h_charge_sum_north_v_vtx","", nvertexbins, vertex_ranges, nbin, -0.5, (float)maxrange - 0.5);
  TH2D *h_charge_sum_north_v_vtx_vtx = new TH2D("h_charge_sum_north_v_vtx_vtx","",  nvertexbins, vertex_ranges,nbin, -0.5, (float)maxrange - 0.5);
  TH2D *h_charge_sum_north_v_vtx_min_bias = new TH2D("h_charge_sum_north_v_vtx_min_bias","",nvertexbins, vertex_ranges, nbin, -0.5, (float)maxrange - 0.5);
  TH2D *h_charge_sum_north_v_vtx_min_bias_w_vertex_cut = new TH2D("h_charge_sum_north_v_vtx_min_bias_w_vertex_cut","",nvertexbins, vertex_ranges, nbin, -0.5, (float)maxrange - 0.5);

  TH2D *h_charge_sum_north_balanced_v_vtx = new TH2D("h_charge_sum_north_balanced_v_vtx","", nvertexbins, vertex_ranges, nbin, -0.5, (float)maxrange - 0.5);
  TH2D *h_charge_sum_north_balanced_v_vtx_vtx = new TH2D("h_charge_sum_north_balanced_v_vtx_vtx","",  nvertexbins, vertex_ranges,nbin, -0.5, (float)maxrange - 0.5);
  TH2D *h_charge_sum_north_balanced_v_vtx_min_bias = new TH2D("h_charge_sum_north_balanced_v_vtx_min_bias","",nvertexbins, vertex_ranges, nbin, -0.5, (float)maxrange - 0.5);
  TH2D *h_charge_sum_north_balanced_v_vtx_min_bias_w_vertex_cut = new TH2D("h_charge_sum_north_balanced_v_vtx_min_bias_w_vertex_cut","",nvertexbins, vertex_ranges, nbin, -0.5, (float)maxrange - 0.5);

  TH2D *h_charge_sum_south_v_vtx = new TH2D("h_charge_sum_south_v_vtx","", nvertexbins, vertex_ranges, nbin, -0.5, (float)maxrange - 0.5);
  TH2D *h_charge_sum_south_v_vtx_vtx = new TH2D("h_charge_sum_south_v_vtx_vtx","",  nvertexbins, vertex_ranges,nbin, -0.5, (float)maxrange - 0.5);
  TH2D *h_charge_sum_south_v_vtx_min_bias = new TH2D("h_charge_sum_south_v_vtx_min_bias","",nvertexbins, vertex_ranges, nbin, -0.5, (float)maxrange - 0.5);
  TH2D *h_charge_sum_south_v_vtx_min_bias_w_vertex_cut = new TH2D("h_charge_sum_south_v_vtx_min_bias_w_vertex_cut","",nvertexbins, vertex_ranges, nbin, -0.5, (float)maxrange - 0.5);

  TH2D *h_charge_sum_south_balanced_v_vtx = new TH2D("h_charge_sum_south_balanced_v_vtx","", nvertexbins, vertex_ranges, nbin, -0.5, (float)maxrange - 0.5);
  TH2D *h_charge_sum_south_balanced_v_vtx_vtx = new TH2D("h_charge_sum_south_balanced_v_vtx_vtx","",  nvertexbins, vertex_ranges,nbin, -0.5, (float)maxrange - 0.5);
  TH2D *h_charge_sum_south_balanced_v_vtx_min_bias = new TH2D("h_charge_sum_south_balanced_v_vtx_min_bias","",nvertexbins, vertex_ranges, nbin, -0.5, (float)maxrange - 0.5);
  TH2D *h_charge_sum_south_balanced_v_vtx_min_bias_w_vertex_cut = new TH2D("h_charge_sum_south_balanced_v_vtx_min_bias_w_vertex_cut","",nvertexbins, vertex_ranges, nbin, -0.5, (float)maxrange - 0.5);

  TH1D *h_charge_sum_trigvtx[3];
  TH1D *h_charge_sum_vtx_trigvtx[3];// = new TH1D("h_charge_sum_vtx","",  nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_min_bias_trigvtx[3];// = new TH1D("h_charge_sum_min_bias","", nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_min_bias_w_vertex_cut_trigvtx[3];// = new TH1D("h_charge_sum_min_bias_w_vertex_cut","", nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_balanced_trigvtx[3];
  TH1D *h_charge_sum_balanced_vtx_trigvtx[3];// = new TH1D("h_charge_sum_balanced_vtx","",  nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_balanced_min_bias_trigvtx[3];// = new TH1D("h_charge_sum_balanced_min_bias","", nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_balanced_min_bias_w_vertex_cut_trigvtx[3];// = new TH1D("h_charge_sum_balanced_min_bias_w_vertex_cut","", nbin, -0.5, (float)maxrange - 0.5);

  for (int i = 0; i < 3; i++)
    {
      h_charge_sum_trigvtx[i] = new TH1D(Form("h_charge_sum_trigvtx_%d", i),"", nbin, -0.5, (float)maxrange - 0.5);
      h_charge_sum_vtx_trigvtx[i] = new TH1D(Form("h_charge_sum_vtx_trigvtx_%d", i),"",  nbin, -0.5, (float)maxrange - 0.5);
      h_charge_sum_min_bias_trigvtx[i] = new TH1D(Form("h_charge_sum_min_bias_trigvtx_%d", i),"", nbin, -0.5, (float)maxrange - 0.5);
      h_charge_sum_min_bias_w_vertex_cut_trigvtx[i] = new TH1D(Form("h_charge_sum_min_bias_w_vertex_cut_trigvtx_%d", i),"", nbin, -0.5, (float)maxrange - 0.5);
      h_charge_sum_balanced_trigvtx[i] = new TH1D(Form("h_charge_sum_balanced_trigvtx_%d", i),"", nbin, -0.5, (float)maxrange - 0.5);
      h_charge_sum_balanced_vtx_trigvtx[i] = new TH1D(Form("h_charge_sum_balanced_vtx_trigvtx_%d", i),"",  nbin, -0.5, (float)maxrange - 0.5);
      h_charge_sum_balanced_min_bias_trigvtx[i] = new TH1D(Form("h_charge_sum_balanced_min_bias_trigvtx_%d", i),"", nbin, -0.5, (float)maxrange - 0.5);
      h_charge_sum_balanced_min_bias_w_vertex_cut_trigvtx[i] = new TH1D(Form("h_charge_sum_balanced_min_bias_w_vertex_cut_trigvtx_%d", i),"", nbin, -0.5, (float)maxrange - 0.5);

    }
  TH1D *h_charge_sum_min_bias_w_vertex_cut_byevent[25];// = new TH1D("h_charge_sum_min_bias_w_vertex_cut","", nbin, -0.5, (float)maxrange - 0.5);
  for (int i = 0; i < 25; i++)
    {
      h_charge_sum_min_bias_w_vertex_cut_byevent[i] = new TH1D(Form("h_charge_sum_min_bias_w_vertex_cut_byevent_%d", i),"", nbin, -0.5, (float)maxrange - 0.5);
    }

  TH1D *h_charge_sum = new TH1D("h_charge_sum","", nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_vtx = new TH1D("h_charge_sum_vtx","",  nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_min_bias = new TH1D("h_charge_sum_min_bias","", nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_min_bias_w_vertex_cut = new TH1D("h_charge_sum_min_bias_w_vertex_cut","", nbin, -0.5, (float)maxrange - 0.5);

  TH1D *h_charge_sum_balanced = new TH1D("h_charge_sum_balanced","", nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_vtx_balanced = new TH1D("h_charge_sum_vtx_balanced","",  nbin, -0.5, (float)maxrange - 0.5);
  //TH1D *h_charge_sum_min_bias_balanced = new TH1D("h_charge_sum_min_bias_balanced","",  nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_min_bias_w_vertex_cut_balanced = new TH1D("h_charge_sum_min_bias_w_vertex_cut_balanced","", nbin, -0.5, (float)maxrange - 0.5);

  TH1D *h_charge_sum_scaled = new TH1D("h_charge_sum_scaled","", nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_vtx_scaled = new TH1D("h_charge_sum_vtx_scaled","",  nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_min_bias_scaled = new TH1D("h_charge_sum_min_bias_scaled","",  nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_min_bias_w_vertex_cut_scaled = new TH1D("h_charge_sum_min_bias_w_vertex_cut_scaled","", nbin, -0.5, (float)maxrange - 0.5);

  TH1D *h_charge_sum_balanced_scaled = new TH1D("h_charge_sum_balanced_scaled","", nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_vtx_balanced_scaled = new TH1D("h_charge_sum_vtx_balanced_scaled","",  nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_min_bias_balanced_scaled = new TH1D("h_charge_sum_min_bias_balanced_scaled","",  nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_min_bias_w_vertex_cut_balanced_scaled = new TH1D("h_charge_sum_min_bias_w_vertex_cut_balanced_scaled","", nbin, -0.5, (float)maxrange - 0.5);

  TH1D *h_charge_sum_fine = new TH1D("h_charge_sum_fine","", nbin*finefactor, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_fine_vtx = new TH1D("h_charge_sum_fine_vtx","",  nbin*finefactor, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_fine_min_bias = new TH1D("h_charge_sum_fine_min_bias","", nbin*finefactor, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_fine_min_bias_w_vertex_cut = new TH1D("h_charge_sum_fine_min_bias_w_vertex_cut","", nbin*finefactor, -0.5, (float)maxrange - 0.5);

  TH1D *h_charge_sum_fine_balanced = new TH1D("h_charge_sum_fine_balanced","", nbin*finefactor, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_fine_vtx_balanced = new TH1D("h_charge_sum_fine_vtx_balanced","",  nbin*finefactor, -0.5, (float)maxrange - 0.5);
  //TH1D *h_charge_sum_fine_min_bias_balanced = new TH1D("h_charge_sum_fine_min_bias_balanced","",  nbin*finefactor, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_fine_min_bias_w_vertex_cut_balanced = new TH1D("h_charge_sum_fine_min_bias_w_vertex_cut_balanced","", nbin*finefactor, -0.5, (float)maxrange - 0.5);

  TH1D *h_charge_sum_fine_scaled = new TH1D("h_charge_sum_fine_scaled","", nbin*finefactor, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_fine_vtx_scaled = new TH1D("h_charge_sum_fine_vtx_scaled","",  nbin*finefactor, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_fine_min_bias_scaled = new TH1D("h_charge_sum_fine_min_bias_scaled","",  nbin*finefactor, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_fine_min_bias_w_vertex_cut_scaled = new TH1D("h_charge_sum_fine_min_bias_w_vertex_cut_scaled","", nbin*finefactor, -0.5, (float)maxrange - 0.5);

  TH1D *h_charge_sum_fine_balanced_scaled = new TH1D("h_charge_sum_fine_balanced_scaled","", nbin*finefactor, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_fine_vtx_balanced_scaled = new TH1D("h_charge_sum_fine_vtx_balanced_scaled","",  nbin*finefactor, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_fine_min_bias_balanced_scaled = new TH1D("h_charge_sum_fine_min_bias_balanced_scaled","",  nbin*finefactor, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_fine_min_bias_w_vertex_cut_balanced_scaled = new TH1D("h_charge_sum_fine_min_bias_w_vertex_cut_balanced_scaled","", nbin*finefactor, -0.5, (float)maxrange - 0.5);

  TH1D *h_charge_sum_ns[2];
  TH1D *h_charge_sum_vtx_ns[2];
  TH1D *h_charge_sum_min_bias_w_vertex_cut_ns[2];

  for (int i = 0; i < 2; i++)
    {
      h_charge_sum_ns[i] = new TH1D(Form("h_charge_sum_%s", (i ? "n":"s")),"", nbin, -0.5, (float)maxrange - 0.5);
      h_charge_sum_vtx_ns[i] = new TH1D(Form("h_charge_sum_vtx_%s", (i ? "n":"s")),"", nbin, -0.5, (float)maxrange - 0.5);
      h_charge_sum_min_bias_w_vertex_cut_ns[i] = new TH1D(Form("h_charge_sum_min_bias_w_vertex_cut_%s", (i ? "n":"s")),"", nbin, -0.5, (float)maxrange - 0.5);
    }


  //// vertex
  int hits_n = 0;
  int hits_s = 0;      
  float nsum = 0;
  float ssum = 0;

  unsigned short mbd_bit = 10;
  if (!silence)   cout <<"NEvents = "<<(nevents? nevents:m_ttree->GetEntries())<<std::endl;
  int  nev = (nevents > 0 && nevents < m_ttree->GetEntries()?nevents:m_ttree->GetEntries());
  for (int i = 0 ; i < (debug ? 10 : nev); i++)
    {
      hits_n = 0;
      hits_s = 0;      
      nsum = 0;
      ssum = 0;
      m_ttree->GetEntry(i);
      if (!(((gl1_scaled >> 14) & 0x1 ) == 0x1)) continue;
      for (int ich = 0 ; ich < 64; ich++)
	{
	  if (countbefore)
	    {
	      nsum += mbd_charge[ich+64];
	      ssum += mbd_charge[ich];
	    }
	  if (fabs(mbd_time[ich]) < 25. && mbd_charge[ich] > cthresh)
	    {
	      hits_s++;
	      if (!countbefore) ssum += mbd_charge[ich];
	    }
	  if (fabs(mbd_time[ich+64]) < 25. && mbd_charge[ich+64] > cthresh)
	    {
	      hits_n++;
	      if (!countbefore) nsum += mbd_charge[ich+64];
	    }
	} 
      if (!silence) if (i%20000 == 0) std::cout <<" ---------   " << hits_s<< " + " <<hits_n<< "    ----"<<std::endl;
      // don't continue if this does not work.

      h_mbd_hits_ns->Fill(hits_s, hits_n);
      if (hits_s < 2 || hits_n < 2) continue; 

      h2_charge_sum_ns->Fill(ssum, nsum);
      h_charge_sum->Fill(ssum + nsum);
      h_charge_sum_fine->Fill(ssum + nsum);
      h_charge_sum_v_vtx->Fill(mbd_vertex, ssum + nsum);
      h_charge_sum_south_v_vtx->Fill(mbd_vertex, ssum);
      h_charge_sum_north_v_vtx->Fill(mbd_vertex, nsum);

      for (unsigned short ibit = 0; ibit < 3; ibit++)
	{
	  if ((gl1_scaled >> (mbd_bit + 2 + ibit)) & 0x1U)
	    h_charge_sum_trigvtx[ibit]->Fill(ssum + nsum);
	}
      
      for (int j = 0; j < 2; j++) h_charge_sum_ns[j]->Fill(mbd_charge_sum[j]);
      if (!silence) if (i%20000 == 0) std::cout <<" ---------   " << mbd_vertex<< "          ----"<<std::endl;		      


      if (fabs(mbd_vertex) < z_cut)
	{
	  h_charge_sum_vtx->Fill(ssum + nsum);
	  h_charge_sum_fine_vtx->Fill(ssum + nsum);
	  h2_charge_sum_vtx_ns->Fill(ssum,nsum);
	  h_charge_sum_v_vtx_vtx->Fill(mbd_vertex, ssum + nsum);
	  h_charge_sum_south_v_vtx_vtx->Fill(mbd_vertex, ssum);
	  h_charge_sum_north_v_vtx_vtx->Fill(mbd_vertex, nsum);
	  for (unsigned short ibit = 0; ibit < 3; ibit++)
	    {
	      if ((gl1_scaled >> (mbd_bit + 2 + ibit)) & 0x1U)
		h_charge_sum_vtx_trigvtx[ibit]->Fill(ssum + nsum);
	    }
	}

      if (!isSim && ssum > charge_sum_south_cut && nsum < charge_sum_north_cut) continue;
      if (hasZDC && !isSim && (zdc_sum[0] <= zdc_cut || zdc_sum[1] <= zdc_cut)) continue;

      h_charge_sum_min_bias->Fill(ssum + nsum);
      h_charge_sum_fine_min_bias->Fill(ssum + nsum);
      h2_charge_sum_min_bias_ns->Fill(ssum,nsum);
      h_charge_sum_v_vtx_min_bias->Fill(mbd_vertex, ssum + nsum);
      h_charge_sum_south_v_vtx_min_bias->Fill(mbd_vertex, ssum);
      h_charge_sum_north_v_vtx_min_bias->Fill(mbd_vertex, nsum);
      for (unsigned short ibit = 0; ibit < 3; ibit++)
	{
	  if ((gl1_scaled >> (mbd_bit + 2 + ibit)) & 0x1U)
	    h_charge_sum_min_bias_trigvtx[ibit]->Fill(ssum + nsum);
	}

      if (fabs(mbd_vertex) > z_cut) continue;

      if (hasZDC && !isSim) 
	{
	  h2_charge_sum_v_zdc->Fill(ssum+nsum, zdc_sum[0] + zdc_sum[1]);
	  hp_charge_sum_v_zdc->Fill(ssum+nsum, zdc_sum[0] + zdc_sum[1]);
	}

      if (hasZDC && !isSim) 
	{
	  h2_charge_sum_v_zdc_mb->Fill(ssum+nsum, zdc_sum[0] + zdc_sum[1]);
	  hp_charge_sum_v_zdc_mb->Fill(ssum+nsum, zdc_sum[0] + zdc_sum[1]);
	}
      //      if (!minbias) continue;
      h_charge_sum_v_vtx_min_bias_w_vertex_cut->Fill(mbd_vertex, ssum + nsum);
      h_charge_sum_south_v_vtx_min_bias_w_vertex_cut->Fill(mbd_vertex, ssum);
      h_charge_sum_north_v_vtx_min_bias_w_vertex_cut->Fill(mbd_vertex, nsum);
      h2_charge_sum_min_bias_w_vertex_cut_ns->Fill(ssum,nsum);
      h_charge_sum_min_bias_w_vertex_cut->Fill(ssum + nsum);
      h_charge_sum_fine_min_bias_w_vertex_cut->Fill(ssum + nsum);
      int ii = i/1000000;
      if (ii < 25)
	h_charge_sum_min_bias_w_vertex_cut_byevent[ii]->Fill(ssum + nsum);
      for (unsigned short ibit = 0; ibit < 3; ibit++)
	{
	  if ((gl1_scaled >> (mbd_bit + 2 + ibit)) & 0x1U)
	    h_charge_sum_min_bias_w_vertex_cut_trigvtx[ibit]->Fill(ssum + nsum);
	}
      

      for (int j = 0; j < 2; j++) h_charge_sum_min_bias_w_vertex_cut_ns[j]->Fill((j?nsum:ssum));      
      if (!silence) if (i%20000 == 0) std::cout <<" ---------   " << ssum + nsum << "          ----"<<std::endl;		      
    }

  double vertex_scale_factor[nvertexbins] = {0};
  double vertex_scale_factor_vtx[nvertexbins] = {0};
  double vertex_n_scale_factor[nvertexbins] = {0};
  double vertex_s_scale_factor[nvertexbins] = {0};
  for (int ib = 0; ib < nvertexbins; ib++)
    { 
      if (isSim)
	{
	  vertex_scale_factor[ib] = 1.0;
	  vertex_scale_factor_vtx[ib] = 1.0;
	  vertex_n_scale_factor[ib] = 1.0;
	  vertex_s_scale_factor[ib] = 1.0;
	  continue;
	}
      TH1D *h = (TH1D*) h_charge_sum_v_vtx_min_bias_w_vertex_cut->ProjectionY(Form("hRealMBD_%d", ib), ib+1, ib+1);
      vertex_scale_factor[ib] = h->GetMean();
      TH1D *h1 = (TH1D*) h_charge_sum_v_vtx_vtx->ProjectionY(Form("hRealMBDvtx_%d", ib), ib+1, ib+1);
      vertex_scale_factor_vtx[ib] = h1->GetMean();
      TH1D *hn = (TH1D*) h_charge_sum_north_v_vtx_min_bias_w_vertex_cut->ProjectionY(Form("hRealMBD_north_%d", ib), ib+1, ib+1);
      TH1D *hs = (TH1D*) h_charge_sum_south_v_vtx_min_bias_w_vertex_cut->ProjectionY(Form("hRealMBD_south_%d", ib), ib+1, ib+1);
      vertex_n_scale_factor[ib] = hn->GetMean();
      vertex_s_scale_factor[ib] = hs->GetMean();
      if (!silence) std::cout << ib << " --> " << vertex_n_scale_factor[ib]/vertex_s_scale_factor[ib] << std::endl;
    }
  double average_middle=vertex_scale_factor[nvertexbins/2];
  double average_middle_vtx=vertex_scale_factor_vtx[nvertexbins/2];
  double average_middle_s=vertex_s_scale_factor[nvertexbins/2];
  double average_middle_n=vertex_n_scale_factor[nvertexbins/2];
  if (!isSim)
    {
      for (int ib = 0; ib < nvertexbins; ib++)
	{
	  vertex_scale_factor[ib] =  average_middle/(vertex_scale_factor[ib]);
	  vertex_scale_factor_vtx[ib] =  average_middle_vtx/(vertex_scale_factor_vtx[ib]);
	  vertex_n_scale_factor[ib] =  average_middle_n/(vertex_n_scale_factor[ib]);
	  vertex_s_scale_factor[ib] =  average_middle_s/(vertex_s_scale_factor[ib]);
	}
    }

  for (int i = 0 ; i < (debug ? 10 : nev); i++)
    {
      //minbias = false;
      hits_n = 0;
      hits_s = 0;      
      nsum = 0;
      ssum = 0;
      m_ttree->GetEntry(i);
      if (!(((gl1_scaled >> 14) & 0x1 ) == 0x1)) continue;
      double vscale = 0;
      for (int iv = 0; iv < nvertexbins; iv++)
	{
	  if (mbd_vertex < vertex_ranges[iv+1])
	    {
	      vscale = vertex_scale_factor[iv];
	      break;
	    }
	}

      for (int ich = 0 ; ich < 64; ich++)
	{
	  if (countbefore)
	    {
	      ssum += mbd_charge[ich]*vscale;
	      nsum += mbd_charge[ich+64]*vscale;
	    }
	  if (mbd_time[ich] < 25. && mbd_charge[ich] > cthresh)
	    {
	      if (!countbefore) ssum += mbd_charge[ich]*vscale;
	      hits_s++;
	    }
	  if (mbd_time[ich+64] < 25. && mbd_charge[ich+64] > cthresh)
	    {
	      hits_n++;
	      if (!countbefore) nsum += mbd_charge[ich+64]*vscale;
	    }
	} 
      // don't continue if this does not work.
      if (hits_s < 2 || hits_n < 2) continue; 
      if (ssum > charge_sum_south_cut && nsum < charge_sum_north_cut) continue;      
      if (fabs(mbd_vertex) > z_cut) continue;

      h_charge_sum_vtx_balanced->Fill(nsum + ssum);
      h_charge_sum_fine_vtx_balanced->Fill(nsum + ssum);

      if (hasZDC && (zdc_sum[0] <= zdc_cut || zdc_sum[1] <= zdc_cut)) continue;
      //if (!minbias) continue;
      h_charge_sum_balanced_v_vtx_min_bias_w_vertex_cut->Fill(mbd_vertex, ssum + nsum);
      h_charge_sum_min_bias_w_vertex_cut_balanced->Fill(ssum + nsum);
      h_charge_sum_fine_min_bias_w_vertex_cut_balanced->Fill(ssum + nsum);
      h_charge_sum_north_balanced_v_vtx_min_bias_w_vertex_cut->Fill(mbd_vertex,nsum);
      h_charge_sum_south_balanced_v_vtx_min_bias_w_vertex_cut->Fill(mbd_vertex, ssum);

      for (unsigned short ibit = 0; ibit < 3; ibit++)
	{
	  if ((gl1_scaled >> (mbd_bit + 2 + ibit)) & 0x1U)
	    h_charge_sum_balanced_min_bias_w_vertex_cut_trigvtx[ibit]->Fill(ssum + nsum);
	}

    }


  float scale = 0;
  float scale_bal = 0;

  if (reference_run && !isSim)
    {
      if (!silence) std::cout  << "Scaling charge sum."<< std::endl;

      scale = mbd_charge_factor/(h_charge_sum_min_bias_w_vertex_cut->GetMean());
      scale_bal = mbd_charge_factor_bal/(h_charge_sum_min_bias_w_vertex_cut_balanced->GetMean());
      if (!silence) std::cout  << "Scaling by " << scale << " ( "<<scale_bal << ") " << std::endl;
      int  nev = (nevents > 0 && nevents < m_ttree->GetEntries()?nevents:m_ttree->GetEntries());
      for (int i = 0 ; i < (debug ? 10 : nev); i++)
	{
	  hits_n = 0;
	  hits_s = 0;      
	  ssum=0;
	  nsum=0;

	  m_ttree->GetEntry(i);
	  if (!(((gl1_scaled >> 14) & 0x1 ) == 0x1)) continue;
	  for (int ich = 0 ; ich < 64; ich++)
	    {
	      if (countbefore)
		{
		  ssum+=scale*mbd_charge[ich];
		  nsum += scale*mbd_charge[ich+64];
		}
	      if (mbd_time[ich] < 25. && mbd_charge[ich] > cthresh)
		{
		  hits_s++;
		  if (!countbefore) ssum+=scale*mbd_charge[ich];
		}
	      if (mbd_time[ich+64] < 25. && mbd_charge[ich+64] > cthresh)
		{
		  if (!countbefore) nsum += scale*mbd_charge[ich+64];
		  hits_n++;
		}	
	    } 
	  double total_sum = nsum+ssum;
	  // don't continue if this does not work.
	  if (hits_s < 2 || hits_n < 2) continue; 
	  h_charge_sum_scaled->Fill(total_sum);
	  h_charge_sum_fine_scaled->Fill(total_sum);
	  if (fabs(mbd_vertex) > z_cut) continue;
	  h_charge_sum_vtx_scaled->Fill(total_sum);
	  h_charge_sum_fine_vtx_scaled->Fill(total_sum);
	  if (hasZDC && (zdc_sum[0] <= zdc_cut || zdc_sum[1] <= zdc_cut)) continue;
	  h_charge_sum_fine_min_bias_scaled->Fill(total_sum);
	  h_charge_sum_min_bias_scaled->Fill(total_sum);
	  if (ssum > charge_sum_south_cut && nsum < charge_sum_north_cut) continue;
	  //if (!minbias) continue;
	  h_charge_sum_min_bias_w_vertex_cut_scaled->Fill(total_sum);
	  h_charge_sum_fine_min_bias_w_vertex_cut_scaled->Fill(total_sum);

	}

      for (int i = 0 ; i < (debug ? 10 : nev); i++)
	{
	  hits_n = 0;
	  hits_s = 0;      
	  ssum=0;
	  nsum=0;

	  m_ttree->GetEntry(i);
	  if (!(((gl1_scaled >> 14) & 0x1 ) == 0x1)) continue;
	  double vscale = 0;
	  for (int iv = 0; iv < nvertexbins; iv++)
	    {
	      if (mbd_vertex < vertex_ranges[iv+1])
		{
		  vscale = vertex_scale_factor[iv];
		  break;
		}
	    }

	  for (int ich = 0 ; ich < 64; ich++)
	    {
	      if (countbefore)
		{
		  ssum+=scale*vscale*mbd_charge[ich];
		  nsum += scale*vscale*mbd_charge[ich+64];
		}
	      if (mbd_time[ich] < 25. && mbd_charge[ich] > cthresh)
		{
		  if (!countbefore) ssum+=scale*vscale*mbd_charge[ich];
		  hits_s++;
		}
	      if (mbd_time[ich+64] < 25. && mbd_charge[ich+64] > cthresh)
		{
		  if (!countbefore) nsum += scale*vscale*mbd_charge[ich+64];
		  hits_n++;
		}	

	    } 

	  double total_sum = nsum+ssum;
	  // don't continue if this does not work.
	  if (hits_s < 2 || hits_n < 2) continue; 
	  h_charge_sum_balanced_scaled->Fill(total_sum);
	  h_charge_sum_fine_balanced_scaled->Fill(total_sum);
	  if (fabs(mbd_vertex) > z_cut) continue;
	  h_charge_sum_vtx_balanced_scaled->Fill(total_sum);
	  h_charge_sum_fine_vtx_balanced_scaled->Fill(total_sum);
	  if (hasZDC && (zdc_sum[0] <= zdc_cut || zdc_sum[1] <= zdc_cut)) continue;
	  h_charge_sum_min_bias_balanced_scaled->Fill(total_sum);
	  h_charge_sum_fine_min_bias_balanced_scaled->Fill(total_sum);
	  if (ssum > charge_sum_south_cut && nsum < charge_sum_north_cut) continue;
	  //if (!minbias) continue;
	  h_charge_sum_min_bias_w_vertex_cut_balanced_scaled->Fill(total_sum);
	  h_charge_sum_fine_min_bias_w_vertex_cut_balanced_scaled->Fill(total_sum);
	}

    }
  if (!silence)
    {
      std::cout << " *************************************** " << std::endl;

      std::cout << "  Run " << runnumber << std::endl;
      std::cout << "    Mean  : Noraml\tBalanced" << (reference_run ? "\tScaled\tScaBa":"")<< std::endl;
      std::cout << "    Raw   : "<<h_charge_sum->GetMean()<<"\t"<<h_charge_sum_balanced->GetMean();
      if (reference_run)
	{
	  std::cout << "\t\t"<<h_charge_sum_scaled->GetMean()<<"\t"<<h_charge_sum_balanced_scaled->GetMean();
	}
      std::cout << " " <<std::endl;
      std::cout << "    vtx   : "<<h_charge_sum_vtx->GetMean()<<"\t"<<h_charge_sum_vtx_balanced->GetMean();
      if (reference_run)
	{
	  std::cout << "\t\t"<<h_charge_sum_vtx_scaled->GetMean()<<"\t"<<h_charge_sum_vtx_balanced_scaled->GetMean();
	}
      std::cout << " " <<std::endl;

      std::cout << "    mb    : "<<h_charge_sum_min_bias_w_vertex_cut->GetMean()<<"\t"<<h_charge_sum_min_bias_w_vertex_cut_balanced->GetMean();  
      if (reference_run)
	{
	  std::cout << "\t\t"<<h_charge_sum_min_bias_w_vertex_cut_scaled->GetMean()<<"\t"<<h_charge_sum_min_bias_w_vertex_cut_balanced_scaled->GetMean();
	}
      std::cout << " " <<std::endl;


      std::cout << " *************************************** " << std::endl;


    }
  qa_info.charge_sum_mean_ns = h_charge_sum_min_bias_w_vertex_cut_ns[0]->GetMean()/h_charge_sum_min_bias_w_vertex_cut_ns[1]->GetMean();
  qa_info.charge_sum_mean[0] = h_charge_sum_min_bias_w_vertex_cut->GetMean();
  qa_info.charge_sum_mean[2] = h_charge_sum_min_bias_w_vertex_cut_balanced->GetMean();
  if (reference_run)
    {
      qa_info.charge_sum_mean[1] = h_charge_sum_min_bias_w_vertex_cut_scaled->GetMean();
      qa_info.charge_sum_mean[3] = h_charge_sum_min_bias_w_vertex_cut_balanced_scaled->GetMean();
      qa_info.scale_factor = scale_bal;
    }

  
  char *pathvertex = new char[100];
  if (isSim)
    {
      sprintf(pathvertex, "%s/mbd_vertex_scale_%s.root" , env_calib, mc_map[runnumber].c_str());
    }
  else
    {
      sprintf(pathvertex, "%s/mbd_vertex_scale_%d.root" ,env_calib, runnumber);
    }

  char *pathout = new char[100];
  if (isSim)
    {
      sprintf(pathout, "%s/mbdana_charge_sum_%s.root" , env_out, mc_map[runnumber].c_str());
    }
  else
    {
      sprintf(pathout, "%s/mbdana_charge_sum_%d.root" ,env_out, runnumber);
    }

  TFile *fcalib = new TFile(pathvertex, "RECREATE");
  TNtuple *tss = new TNtuple("tn_vertexscale","scalefactor by vertex","bin:low_vertex:high_vertex:scale:nscale:sscale");
  for (int iv = 0; iv < nvertexbins;iv++)
    {
      tss->Fill(iv, vertex_ranges[iv], vertex_ranges[iv+1], vertex_scale_factor[iv], vertex_n_scale_factor[iv], vertex_s_scale_factor[iv]);
    }
  tss->Write();
  fcalib->Close();

  TFile *fout = new TFile(pathout, "RECREATE");
  TNtuple *ts = new TNtuple("tn_scale","scalefactor","run:scale");
  ts->Fill(runnumber, scale_bal);
  ts->Write();
  h_mbd_hits_ns->Write();
  h2_charge_sum_ns_cut->Write();

  h_charge_sum->Write();
  h_charge_sum_vtx->Write();
  h_charge_sum_min_bias_w_vertex_cut->Write();
  h_charge_sum_min_bias->Write();
  h_charge_sum_balanced->Write();
  h_charge_sum_vtx_balanced->Write();
  h_charge_sum_min_bias_w_vertex_cut_balanced->Write();

  h_charge_sum_balanced_scaled->Write();
  h_charge_sum_vtx_balanced_scaled->Write();
  h_charge_sum_min_bias_w_vertex_cut_balanced_scaled->Write();
  h_charge_sum_scaled->Write();
  h_charge_sum_vtx_scaled->Write();
  h_charge_sum_min_bias_w_vertex_cut_scaled->Write();

  h_charge_sum_fine->Write();
  h_charge_sum_fine_vtx->Write();
  h_charge_sum_fine_min_bias_w_vertex_cut->Write();
  h_charge_sum_fine_min_bias->Write();
  h_charge_sum_fine_balanced->Write();
  h_charge_sum_fine_vtx_balanced->Write();
  h_charge_sum_fine_min_bias_w_vertex_cut_balanced->Write();

  h_charge_sum_fine_balanced_scaled->Write();
  h_charge_sum_fine_vtx_balanced_scaled->Write();
  h_charge_sum_fine_min_bias_w_vertex_cut_balanced_scaled->Write();
  h_charge_sum_fine_scaled->Write();
  h_charge_sum_fine_vtx_scaled->Write();
  h_charge_sum_fine_min_bias_w_vertex_cut_scaled->Write();

  h2_charge_sum_ns->Write();
  h2_charge_sum_vtx_ns->Write();
  h2_charge_sum_min_bias_w_vertex_cut_ns->Write();
  h2_charge_sum_min_bias_ns->Write();

  h_charge_sum_v_vtx->Write();
  h_charge_sum_v_vtx_vtx->Write();
  h_charge_sum_v_vtx_min_bias_w_vertex_cut->Write();
  h_charge_sum_v_vtx_min_bias->Write();
  h_charge_sum_north_v_vtx->Write();
  h_charge_sum_north_v_vtx_vtx->Write();
  h_charge_sum_north_v_vtx_min_bias_w_vertex_cut->Write();
  h_charge_sum_north_v_vtx_min_bias->Write();
  h_charge_sum_south_v_vtx->Write();
  h_charge_sum_south_v_vtx_vtx->Write();
  h_charge_sum_south_v_vtx_min_bias_w_vertex_cut->Write();
  h_charge_sum_south_v_vtx_min_bias->Write();

  h_charge_sum_balanced_v_vtx->Write();
  h_charge_sum_balanced_v_vtx_vtx->Write();
  h_charge_sum_balanced_v_vtx_min_bias_w_vertex_cut->Write();
  h_charge_sum_balanced_v_vtx_min_bias->Write();
  h_charge_sum_south_balanced_v_vtx->Write();
  h_charge_sum_south_balanced_v_vtx_vtx->Write();
  h_charge_sum_south_balanced_v_vtx_min_bias_w_vertex_cut->Write();
  h_charge_sum_south_balanced_v_vtx_min_bias->Write();
  h_charge_sum_north_balanced_v_vtx->Write();
  h_charge_sum_north_balanced_v_vtx_vtx->Write();
  h_charge_sum_north_balanced_v_vtx_min_bias_w_vertex_cut->Write();
  h_charge_sum_north_balanced_v_vtx_min_bias->Write();

  for (int i = 0; i < 3; i++)
    {
      h_charge_sum_trigvtx[i]->Write();
      h_charge_sum_vtx_trigvtx[i]->Write();
      h_charge_sum_min_bias_w_vertex_cut_trigvtx[i]->Write();
      h_charge_sum_min_bias_trigvtx[i]->Write();
      h_charge_sum_balanced_trigvtx[i]->Write();
      h_charge_sum_balanced_vtx_trigvtx[i]->Write();
      h_charge_sum_balanced_min_bias_w_vertex_cut_trigvtx[i]->Write();
      h_charge_sum_balanced_min_bias_trigvtx[i]->Write();

    }
  for (int i = 0; i < 25; i++)
    {
      h_charge_sum_min_bias_w_vertex_cut_byevent[i]->Write();
    }

  h2_charge_sum_v_zdc->Write();
  h2_charge_sum_v_zdc_mb->Write();
  hp_charge_sum_v_zdc->Write();
  hp_charge_sum_v_zdc_mb->Write();

  for (int j = 0; j < 2; j++) 
    {
      h_charge_sum_ns[j]->Write();
      h_charge_sum_vtx_ns[j]->Write();
      h_charge_sum_min_bias_w_vertex_cut_ns[j]->Write();
    }

  fout->Close();
  delete h_charge_sum;
  delete h_charge_sum_vtx;
  delete h_charge_sum_min_bias_w_vertex_cut;
  delete h_charge_sum_balanced;
  delete h_charge_sum_vtx_balanced;
  delete h_charge_sum_min_bias_w_vertex_cut_balanced;
  for (int j = 0; j < 2; j++) 
    {
      delete h_charge_sum_ns[j];
      delete h_charge_sum_vtx_ns[j];
      delete h_charge_sum_min_bias_w_vertex_cut_ns[j];
    }

  delete h_charge_sum_balanced_scaled;
  delete h_charge_sum_vtx_balanced_scaled;
  delete h_charge_sum_min_bias_w_vertex_cut_balanced_scaled;
  delete h_charge_sum_scaled;
  delete h_charge_sum_vtx_scaled;
  delete h_charge_sum_min_bias_w_vertex_cut_scaled;

    
}
void QA_centrality::QA_MakeDivisions(const int runnumber, const int doVertexScaled)
{

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }

  if(!env_tree)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }
 
  char *path = new char[100];
  if (isSim)
    {
      sprintf(path, "%s/mbd_vertex_scale_%s.root" , env_calib, mc_map[runnumber].c_str());
    }
  else
    {
      sprintf(path, "%s/mbd_vertex_scale_%d.root" ,env_calib, runnumber);
    }

  TFile *fcalibs = new TFile(path, "r");
  TNtuple *tss = (TNtuple*) fcalibs->Get("tn_vertexscale");
  float vhigh;
  float vlow;
  float vscale;
  float vscale_n;
  float vscale_s;
  const int nvertexbins = tss->GetEntries();;
  if (!silence) std::cout  << nvertexbins << std::endl;

  float vertex_ranges[nvertexbins+1] = {0};
  float vertex_scale_factor[nvertexbins] = {0};
  float vertex_n_scale_factor[nvertexbins] = {0};
  float vertex_s_scale_factor[nvertexbins] = {0};
  tss->SetBranchAddress("low_vertex", &vlow);
  tss->SetBranchAddress("high_vertex", &vhigh);
  tss->SetBranchAddress("scale", &vscale);
  tss->SetBranchAddress("nscale", &vscale_n);
  tss->SetBranchAddress("sscale", &vscale_s);
  for (int i = 0; i < tss->GetEntries() ; i++)
    {
      tss->GetEntry(i);      
      vertex_ranges[i] = vlow;
      vertex_ranges[i+1] = vhigh;
      if(doVertexScaled == 1)
	{
	  vertex_scale_factor[i] = vscale;
	  vertex_n_scale_factor[i] = vscale_n;
	  vertex_s_scale_factor[i] = vscale_s;
	}
      else
	{

	  vertex_scale_factor[i] = 1;
	  vertex_n_scale_factor[i] = 1;
	  vertex_s_scale_factor[i] = 1;

	}
      std::cout  << "vscale: " << vscale << std::endl;
    }

  fcalibs->Close();

  if (isSim)
    {
      sprintf(path, "%s/mbdana/mbd_tree_%s_2025.root" , mc_map[runnumber].c_str(), mc_map[runnumber].c_str());
    }
  else
    {
      sprintf(path, "Run24/run%d/mbdana/mbd_trees_%d.root" ,runnumber, runnumber);
    }


  std::vector<float> v_mbd_charge_sum{};
  std::vector<float> v_mbd_charge_sum_vertex[nvertexbins];
  for (int i = 0; i < nvertexbins;i++)
    {
      v_mbd_charge_sum_vertex[i] = {};
    }
  std::vector<float> v_mbd_charge_sum_scaled{};

  //// vertex
  int hits_n = 0;
  int hits_s = 0;      
  float nsum = 0;
  float ssum = 0;

  if (!silence)   cout <<"NEvents = "<<(nevents? nevents:m_ttree->GetEntries())<<std::endl;
  int  nev = (nevents > 0 && nevents < m_ttree->GetEntries()?nevents:m_ttree->GetEntries());
  for (int i = 0 ; i < nev; i++)
    {
      hits_n = 0;
      hits_s = 0;
      ssum = 0;
      nsum = 0;
      m_ttree->GetEntry(i);
      if (!(((gl1_scaled >> 14) & 0x1 ) == 0x1)) continue;
      double scale = 0;
      int vbin = 0;
      for (int iv = 0; iv < nvertexbins; iv++)
	{
	  if (mbd_vertex < vertex_ranges[iv+1] && mbd_vertex >= vertex_ranges[iv])
	    {
	      scale = vertex_scale_factor[iv];
	      vbin = iv;
	      break;
	    }
	}

      for (int ich = 0 ; ich < 64; ich++)
	{
	  if (countbefore)
	    {
	      ssum+=mbd_charge[ich]*scale;
	      nsum += mbd_charge[ich + 64]*scale;
	    }
	  if (mbd_time[ich] < 25. && mbd_charge[ich] > cthresh)
	    {
	      hits_s++;
	      if (!countbefore) ssum+=mbd_charge[ich]*scale;
	    }
	  if (mbd_time[ich+64] < 25. && mbd_charge[ich+64] > cthresh)
	    {
	      hits_n++;
	      if (!countbefore) nsum += mbd_charge[ich + 64]*scale;
	    } 
	}
      //if (!minbias) continue;
      if (hits_s < 2 || hits_n < 2) continue; 
      if (fabs(mbd_vertex) > z_cut) continue;
      if (!isSim && ssum > charge_sum_south_cut && nsum < charge_sum_north_cut) continue;

      if (hasZDC && !isSim && (zdc_sum[0] <= zdc_cut || zdc_sum[1] <= zdc_cut)) continue;
      if (!isSim && (ssum + nsum) > m_maxsumcut) continue; 
      v_mbd_charge_sum.push_back((ssum + nsum));
      v_mbd_charge_sum_vertex[vbin].push_back((ssum + nsum));
    }
  int size = v_mbd_charge_sum.size();
  if (!silence) std::cout  << " Starting Sort " << size << " events.";
  std::sort(v_mbd_charge_sum.begin(), v_mbd_charge_sum.end(), [] (auto a, auto b) { return a > b;} );

  for (int i =0 ; i < nvertexbins; i++)
    {
      std::sort(v_mbd_charge_sum_vertex[i].begin(), v_mbd_charge_sum_vertex[i].end(), [] (auto a, auto b) { return a > b;} );
    }
  if (!silence) std::cout  << "done."<<std::endl;

  if (!silence) std::cout  << "Divisions: "<<std::endl;
  TString   extra = "";    
  if (doVertexScaled == 1) extra = "bal_";
  if (doVertexScaled == 2) extra = "vtx_";

  TString calib_file_name = Form("%s/mbdana_centrality_%s%d.root", env_calib, extra.Data(), runnumber);
  if (isSim)
    {
      calib_file_name = Form("%s/mbdana_centrality_%s%s.root", env_calib, extra.Data(), mc_map[runnumber].c_str());
    }
  TFile *fcalib = new TFile(calib_file_name, "recreate");
  if (doVertexScaled < 2)
    {

      TNtuple *ts = new TNtuple("tn_centrality", "holds centrality divisions", "bin:low:high");
      for (int i = 0; i < divs ; i++)
	{
	  int bin = floor((float)i*(float)size/(float)divs);
	  int bin1 = floor((i+1)*(float)size/(float)divs);
	  if (i == divs - 1) bin1 = size - 1;
	  if (bin1 == size) bin1--;
	  ts->Fill(i, v_mbd_charge_sum[bin1], v_mbd_charge_sum[bin]);
	}
    }
  else 
    {
      for (int iv = 0; iv < nvertexbins; iv++)
	{
	  size = v_mbd_charge_sum_vertex[iv].size();
	  TNtuple *ts = new TNtuple(Form("tn_centrality_%d", iv), "holds centrality divisions", "bin:low:high");
	  for (int i = 0; i < divs ; i++)
	    {
	      int bin = floor((float)i*(float)size/(float)divs);
	      int bin1 = floor((i+1)*(float)size/(float)divs);
	      if (bin1 == size) bin1--;
	      if (i == divs - 1) bin1 = size - 1;
	      ts->Fill(i, v_mbd_charge_sum_vertex[iv][bin1], v_mbd_charge_sum_vertex[iv][bin]);
	    }
	}
    }
  fcalib->Write();
  fcalib->Close();

  return;
}

void QA_centrality::QA_MBDChannels(const int runnumber)
{

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }

  if(!env_tree)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }
  
  // Output from the centrality module.
  // Making fresh histograms

  //raw
  TH1D *h_charge[128];
  TH1D *h_time[128];

  //passes ZDC tube cuts
  TH1D *h_charge_min_bias[128];
  TH1D *h_time_min_bias[128];
  const int nvertexbins = 15;
  double vertex_ranges[nvertexbins + 1] = {-60, -30, -20, -15, -10, -7, -4, -1, 1, 4, 7, 10, 15, 20, 30, 60};
  TH2D *h2_charge_vertex[128];
  TH2D *h2_time_vertex[128];
  TProfile *hp_charge_vertex[128];
  TProfile *hp_time_vertex[128];
  for (int j = 0; j < 128; j++)
    {
      h_charge[j] = new TH1D(Form("h_mbd_charge_ch%d", j), "",200, 0, 10);
      h_time[j] = new TH1D(Form("h_mbd_time_ch%d", j), "",2500, -25, 25);

      h_charge_min_bias[j] = new TH1D(Form("h_mbd_charge_min_bias_ch%d", j), "",200, 0, 10);
      h_time_min_bias[j] = new TH1D(Form("h_mbd_time_min_bias_ch%d", j), "",2500, -25, 25);

      h2_charge_vertex[j] = new TH2D(Form("h2_charge_vertex_ch%d", j), "",nvertexbins, vertex_ranges, 200, 0, 10);
      h2_time_vertex[j] = new TH2D(Form("h2_time_vertex_ch%d", j), "",nvertexbins, vertex_ranges, 2500, -25, 25);

      hp_charge_vertex[j] = new TProfile(Form("hp_charge_vertex_ch%d", j), "",nvertexbins, vertex_ranges);
      hp_time_vertex[j] = new TProfile(Form("hp_time_vertex_ch%d", j), "",nvertexbins, vertex_ranges);

    }

  std::cout << "Starting,,, " <<std::endl;
  int  nev = (nevents > 0 && nevents < m_ttree->GetEntries()?nevents:m_ttree->GetEntries());
  for (int i = 0; i < nev;i++)
    {
      m_ttree->GetEntry(i);
      bool min_bias = true;

      if (!isSim) min_bias &= (zdc_sum[0] >zdc_cut && zdc_sum[1] > zdc_cut);
      min_bias &=  (fabs(mbd_vertex) < 60);
      min_bias &= (!(mbd_sum[0] > 150 && mbd_sum[1] < 10)); 

      for (int ich = 0 ; ich < 128; ich++)
	{
	  h_charge[ich]->Fill(mbd_charge[ich]);
	  h_time[ich]->Fill(mbd_time[ich]);
	  if (min_bias)
	    {
	      h_charge_min_bias[ich]->Fill(mbd_charge[ich]);
	      h_time_min_bias[ich]->Fill(mbd_time[ich]);
	      h2_charge_vertex[ich]->Fill(mbd_vertex, mbd_charge[ich]);
	      hp_charge_vertex[ich]->Fill(mbd_vertex, mbd_charge[ich]);
	      if (fabs(mbd_time[ich]) < 20)
		{
		  h2_time_vertex[ich]->Fill(mbd_vertex, mbd_time[ich]);
		  hp_time_vertex[ich]->Fill(mbd_vertex, mbd_time[ich]);
		}
	    }
	}
    }

  
  char *path_out = new char[100];
  if (isSim)
    {
      sprintf(path_out, "%s/mbdana_channels_%s.root" , env_out, mc_map[runnumber].c_str());
    }
  else
    {
      sprintf(path_out, "%s/mbdana_channels_%d.root" , env_out, runnumber);
    }
  
  TFile *fout = new TFile(path_out, "recreate");
  for (int ich = 0 ; ich < 128; ich++)
    {
      h_charge[ich]->Write();
      h_time[ich]->Write();
      h_charge_min_bias[ich]->Write();
      h_time_min_bias[ich]->Write();
      hp_charge_vertex[ich]->Write();
      hp_time_vertex[ich]->Write();
      h2_charge_vertex[ich]->Write();
      h2_time_vertex[ich]->Write();

    }

  fout->Close();			 
}


void QA_centrality::QA_ZDCCheck(const int runnumber)
{


  TH1D *h_bunch_number = new TH1D("h_bunch_number", "; Bunch Number:", 121, -0.5, 120.5);

  TH1D *h_zdc_sum_n = new TH1D("h_zdc_sum_n", ";ZDC sum [GeV]; Counts", 1000, 0, 10000);
  TH1D *h_zdc_sum_s = new TH1D("h_zdc_sum_s", ";ZDC sum [GeV]; Counts", 1000, 0, 10000);

  TH2D *h_zdc_sum_ns = new TH2D("h_zdc_sum_ns", ";ZDC sum S [GeV]; ZDC sum N [GeV]", 1000, 0, 10000, 1000, 0, 10000);
  TH2D *h_mbd_hits_ns = new TH2D("h_mbd_hits_ns", ";nhit S; nhit N", 65, -0.5, 64.5, 65, -0.5, 64.5);
  TH1D *h_zdc_energy[6];
  
  for (int i = 0; i < 6; i++)
    {
      h_zdc_energy[i] = new TH1D(Form("h_zdc_energy_%d", i), ";ZDC Energy [GeV]; Counts", 1000, 0, 5000);
    }

  //int events_per_measure = 10000;
  //int n_measurements = 0;
  TH1D *h_running_average_vertex = new TH1D("h_running_average_vertex","", 60, -60, 60);
  TGraphErrors *g_running_average_vertex = new TGraphErrors();
  g_running_average_vertex->SetName("g_running_average_vertex");
  int running_MB = 0;
  TGraphErrors *g_running_average_MB = new TGraphErrors();
  g_running_average_MB->SetName("g_running_average_MB");
  int running_ZDC = 0;
  TGraphErrors *g_running_average_ZDC = new TGraphErrors();
  g_running_average_ZDC->SetName("g_running_average_ZDC");

  TH1D *h_vertex = new TH1D("h_vertex", ";vertex [cm]; Counts", 601, -300.5, 300.5);
  TH1D *h_time_zero = new TH1D("h_time_zero", ";time zero [ns]; Counts", 501, -25.05, 25.05);


  int n_tight_vertex = 0;;
  int n_two_hits = 0;;
  int n_min_bias = 0;;
  int n_one_min_bias = 0;;
  int n_events = 0;
  int n_vertex = 0;
  int n_passed_ns_cut = 0;
  int n_zdc_in_box = 0;
  int n_nozdc_in_box = 0;
  int  nev = (nevents > 0 && nevents < m_ttree->GetEntries()?nevents:m_ttree->GetEntries());
  std::cout << "Starting,,, " <<std::endl;
  for (int i = 0; i < nev ; i++)
    {
      if (!silence)
	{
	  if (i%10000 == 0) std::cout << " Event " << i << " \r" << std::flush;
	}
      m_ttree->GetEntry(i);
      if (!(((gl1_scaled >> 14) & 0x1 ) == 0x1)) continue;
      n_events++;
      h_vertex->Fill(mbd_vertex);
      h_time_zero->Fill(mbd_time_zero);

      int hits_n = 0;
      int hits_s = 0;
      
      for (int j = 0; j < 64; j++)
	{
	  if (mbd_charge[j] > cthresh && fabs(mbd_time[j]) < 25.)
	    hits_s++;
	  if (mbd_charge[j+64] > cthresh && fabs(mbd_time[j+64]) < 25.)
	    hits_n++;
	}

      h_zdc_sum_ns->Fill(zdc_sum[0], zdc_sum[1]);
      h_mbd_hits_ns->Fill(hits_s, hits_n);
      if (!(hits_n >= 2 && hits_s >= 2) ) continue;

      n_two_hits++;
      if (fabs(mbd_vertex) > z_cut) continue;
      for (int j = 0; j < 6; j++) h_zdc_energy[j]->Fill(zdc_energy[j]);
      h_zdc_sum_n->Fill(zdc_sum[1]);
      h_zdc_sum_s->Fill(zdc_sum[0]);
      n_vertex++;

      if (mbd_sum[1] < charge_sum_north_cut && mbd_sum[0] > charge_sum_south_cut) 
	{
	  if (zdc_sum[0] > zdc_cut && zdc_sum[1] > zdc_cut)
	    {
	      n_zdc_in_box++;
	    }
	  else
	    {
	      n_nozdc_in_box++;
	    }
	  continue;
	}    

      n_passed_ns_cut++;
      running_MB++;
      if (!isSim && !(zdc_sum[0] > zdc_cut || zdc_sum[1] > zdc_cut)) continue;
      n_one_min_bias++;
      if (!isSim && !(zdc_sum[0] > zdc_cut && zdc_sum[1] > zdc_cut)) continue;
      n_min_bias++;
      h_bunch_number->Fill(bunch_number);
      running_ZDC++;

      h_running_average_vertex->Fill(mbd_vertex);
      if (fabs(mbd_vertex) > 10) continue;
      n_tight_vertex++;
    }


  if (!silence) std::cout << " done 1"<<std::endl;
  h_vertex->Fit("gaus",(silence?"Q":""),"", -60, 60);
  float mean = h_vertex->GetFunction("gaus")->GetParameter(1);

  h_vertex->Fit("gaus",(silence?"Q":""),"", mean - 20, mean+20);
  qa_info.vertex = h_vertex->GetFunction("gaus")->GetParameter(1);
  qa_info.vertex_sigma = h_vertex->GetFunction("gaus")->GetParameter(2);

  if (!silence) std::cout << "Vertex found at "<<qa_info.vertex<<std::endl;


  if ((float)n_min_bias/(float)n_passed_ns_cut < 0.1)
  {
    if (!silence) std::cout << "forcing no ZDC" << std::endl;
      hasZDC = false;
  }  
  if (forceZDC) 
    {
      hasZDC = false;
    }
  if (!silence)
    {
      std::cout << "Total events: "<< m_ttree->GetEntries()<<endl;
      std::cout << "Number of events passing Two Hits Cut             : "<<n_two_hits <<"/"<<n_events<<" = "<<static_cast<double>(n_two_hits)/static_cast<double>(n_events)<<std::endl;
      std::cout << "Number of events passing Vertex Cut               : "<<n_vertex <<"/"<<n_events<<" = "<<static_cast<double>(n_vertex)/static_cast<double>(n_events)<<std::endl;
      std::cout << "Number of events passing ChargeSum Cut            : "<<n_passed_ns_cut <<"/"<<n_events<<" = "<<static_cast<double>(n_passed_ns_cut)/static_cast<double>(n_events)<<std::endl;
      std::cout << "Number of events passing ZDC cut after all others : "<<n_min_bias <<"/"<<n_passed_ns_cut<<" = "<<static_cast<double>(n_min_bias)/static_cast<double>(n_passed_ns_cut)<<std::endl;
      std::cout << "Number of events passing 1ZDC cut after all others : "<<n_one_min_bias <<"/"<<n_passed_ns_cut<<" = "<<static_cast<double>(n_one_min_bias)/static_cast<double>(n_passed_ns_cut)<<std::endl;

      std::cout << "Number in box w/ and w.o ZDC                       :"<<n_zdc_in_box<<"("<<n_nozdc_in_box<<")/"<<n_zdc_in_box+n_nozdc_in_box<<std::endl;
      std::cout << "Number of min_bias events with vertex < 10         : "<<n_tight_vertex <<"/"<<n_min_bias<<" = "<<static_cast<double>(n_tight_vertex)/static_cast<double>(n_min_bias)<<std::endl;

    }

  qa_info.nevents = m_ttree->GetEntries();
  qa_info.vertex_cut = static_cast<double>(n_tight_vertex)/static_cast<double>(n_min_bias);
  qa_info.min_bias = static_cast<double>(n_min_bias)/static_cast<double>(n_events);
  qa_info.ZDC_percentage = static_cast<double>(n_min_bias)/static_cast<double>(n_passed_ns_cut);
  qa_info.one_ZDC_percentage = static_cast<double>(n_one_min_bias)/static_cast<double>(n_passed_ns_cut);
  qa_info.noZDC_background = static_cast<double>(n_nozdc_in_box)/static_cast<double>(n_nozdc_in_box + n_zdc_in_box);
  qa_info.ZDC_background = static_cast<double>(n_zdc_in_box)/static_cast<double>(n_nozdc_in_box + n_zdc_in_box);

  int nbunches = 0;
  double cut_events = n_min_bias/10000.;;
  for (int ibin = 1; ibin <= h_bunch_number->GetNbinsX(); ibin++)
    {
      if (cut_events < h_bunch_number->GetBinContent(ibin))
	{
	  nbunches++;
	}
    }
  qa_info.fillpattern = nbunches;
  TFile *fout = new TFile(Form("%s//mbdana_zdc_check_%d.root", env_out, runnumber),"recreate");
  for (int i = 0; i < 6; i++)
    {
      h_zdc_energy[i]->Write();
    }

  g_running_average_vertex->Write();
  g_running_average_MB->Write();
  g_running_average_ZDC->Write();

  h_bunch_number->Write();

  h_zdc_sum_n->Write();
  h_zdc_sum_s->Write();
  h_zdc_sum_ns->Write();
  h_mbd_hits_ns->Write();

  h_time_zero->Write();
  h_vertex->Write();
  fout->Close();
  
  std::cout << "Moving On."<< std::endl;
}

void QA_centrality::QA_MakeCentralityCalibrations(const int runnumber, const bool doVertexScaled, const bool use_shifted)
{

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }
  ///////////// 
  //===============================================================================================
  bool npartminusone = false;
  bool npartscaling = true;

  //double forcetrigfrac = 100*trigeff;

  if (!silence) std::cout  << "ndivs: "<<ndivs<<std::endl;

  //sphenixconst int maxcentbins  = ndivs - 1;

  //vtxflagbool sphenix = true;
  //int vtxflag = 1;
  std::string name = Form("%s/glauber/%s", env_p, histo_file.c_str());
 
  int lowfit = 100; // lowest end of standard fit range, was 30 in PHENIX
  int highfit = 1900; // lowest end of standard fit range, was 30 in PHENIX

  // the below file contains the "hNcoll" histogram to be sampled from (step #1)
  TFile *fglauber = new TFile(name.c_str()); 
  
  hglauber = (TH1D * ) fglauber->Get("hNpart"); // note that this is only npart in the Au nucleus
  // this creates a problem because there is never an Npart = 1, always Npart >= 2
  // the below file contains the "hbbcQs6" histogram of MBD data distribution (from z-vertex range)
  //=======================================================
  char *fdataname = new char[100];
  if (isSim)
    {
      sprintf(fdataname,"%s/mbdana_charge_sum_%s.root", env_out, mc_map[runnumber].c_str());  
    }
  else
    {
      sprintf(fdataname,"%s/mbdana_charge_sum_%d.root", env_out, runnumber);  
    }

  char runnum[100];
  if (isSim) sprintf(runnum, "%s", mc_map[runnumber].c_str());
  else  sprintf(runnum, "%08d", runnumber);
  //char zselect[100] = "|z_{MBD}|< 60 cm";
  //=======================================================  
  
  TFile *fdata = new TFile(fdataname);

  TH1D *hRealMBD;
  //hRealMBD = (TH1D *) fdata->Get("h_charge_sum_vtx_balanced"); // a particular z-vertex range
  hRealMBD = (TH1D *) fdata->Get("h_charge_sum_min_bias_w_vertex_cut"); // a particular z-vertex range
  if (doVertexScaled && use_shifted)  hRealMBD = (TH1D *) fdata->Get("h_charge_sum_min_bias_w_vertex_cut_balanced_scaled"); // a particular z-vertex range
  else if (doVertexScaled)  hRealMBD = (TH1D *) fdata->Get("h_charge_sum_min_bias_w_vertex_cut_balanced"); // a particular z-vertex range
  else if (use_shifted)  hRealMBD = (TH1D *) fdata->Get("h_charge_sum_min_bias_w_vertex_cut_scaled"); // a particular z-vertex range

  TH1D *hRealMBDfine;
  //hRealMBD = (TH1D *) fdata->Get("h_charge_sum_vtx_balanced"); // a particular z-vertex range
  hRealMBDfine = (TH1D *) fdata->Get("h_charge_sum_fine_min_bias_w_vertex_cut"); // a particular z-vertex range
  if (doVertexScaled && use_shifted)  hRealMBDfine = (TH1D *) fdata->Get("h_charge_sum_fine_min_bias_w_vertex_cut_balanced_scaled"); // a particular z-vertex range
  else if (doVertexScaled)  hRealMBDfine = (TH1D *) fdata->Get("h_charge_sum_fine_min_bias_w_vertex_cut_balanced"); // a particular z-vertex range
  else if (use_shifted)  hRealMBDfine = (TH1D *) fdata->Get("h_charge_sum_fine_min_bias_w_vertex_cut_scaled"); // a particular z-vertex range

  TH2D *h2RealMBD = (TH2D *) fdata->Get("h_charge_sum_v_vtx_min_bias_w_vertex_cut"); // a particular z-vertex range

  if (!hRealMBD || !h2RealMBD)
    {
      if (!silence) std::cout  << " no histograms"<<std::endl;
      return;
    }
  //hRealMBD->Scale(1./hRealMBD->Integral());
  int nhistbins = hRealMBD->GetNbinsX();
  int nhistbinsfine = hRealMBDfine->GetNbinsX();
  // j.nagle - also change 199.5 to be maxrange = ((float)nhistbins) - 0.5
  float maxrange = ((float) nhistbins) - 0.5;
  // might be good to have max Ncoll / Npart as well... (CuAu = 197 + 63 = 270)
  // Au+Au case 197 + 197 = 396
  int ncollmax = 400; // really npart here
  float maxrangencoll = ((float) ncollmax) - 0.5;

  TH1D *hSimMBDHardUnbiased = new TH1D("hSimMBDHardUnbiased","hSimMBDHardUnbiased",nhistbins,-0.5,maxrange);
  TH1D *hSimMBDHardBiased = new TH1D("hSimMBDHardBiased","hSimMBDHardBiased",nhistbins,-0.5,maxrange);

  // for these two, the trigger is not required to have fired ... (MBD vs Ncoll)

  TH2D *hSimMBDNcoll = new TH2D("hSimMBDNcoll","hSimMBDNcoll",nhistbins,-0.5,maxrange,ncollmax,-0.5,maxrangencoll);
  TH2D *hSimMBDNcollHard = new TH2D("hSimMBDNcollHard","hSimMBDNcollHard",nhistbins,-0.5,maxrange,ncollmax,-0.5,maxrangencoll);

  // this now includes the trigger efficiency turn-on and can be used to calculate <Ncoll> for each centrality bin!!!
  
  TH2D *hSimMBDNcoll_wtrig = new TH2D("hSimMBDNcoll_wtrig","hSimMBDNcoll_wtrig",nhistbins,-0.5,maxrange,ncollmax,-0.5,maxrangencoll);

  // step #2
  TF1 *flatline = new TF1("flatline","[0]",lowfit,highfit);
  flatline->SetParameter(0,1.0);

  // Starting the fit
  
  double mu = 4.1;
  double k = 0.8;
  double N = 1.0;
  double biased_mu = mu;
  double biased_k = k;
  const int nvertexbins = 15;
  double cut = hRealMBD->Integral(lowfit, highfit);
  double all =   hRealMBD->Integral()*trigeff;
  double after = cut/all;
  double scalefactor = hglauber->Integral()*after/all;
  TH1D *hRealMBDScaled = (TH1D*) hRealMBD->Clone();
  hRealMBDScaled->Scale(scalefactor);
  
  // Make the function
  TF1 *fNBD = new TF1("nbd", NBDGlauberConv, 0, 2500, 3);
  fNBD->SetParameters(mu,k, N);
  fNBD->SetParNames("mu","k","N");
  fNBD->SetParLimits(0, 2, 5.5);
  fNBD->SetParLimits(1, 0.3, 3);
  fNBD->SetParLimits(2, 0.1, 1.5);
  hRealMBDScaled->Fit("nbd", "NDOR", "",100, 2000);  
  hRealMBDScaled->Fit("nbd", "NDOR", "",100, 2000);  
  hRealMBDScaled->Fit("nbd", "NDOR", "",100, 2000);  
 
  double bestmu = fNBD->GetParameter(0);
  double bestk = fNBD->GetParameter(1);
  double sigma_mu = fNBD->GetParError(0);
  double sigma_k = fNBD->GetParError(1);
  double chi2 = fNBD->GetChisquare();
  if (!silence)
    {
      cout << "=====================================================" << endl;
      cout << "Best chi2 = " << chi2 << " with mu = " << 
	bestmu << " +/- " << sigma_mu << 
	" and k = " <<
	bestk << " +/- " << sigma_k << endl;
      cout << "=====================================================" << endl;
    }

  // use the best values in the rest of the calculation
  mu = bestmu;
  k  = bestk;
  if (!silence) std::cout << "filling " << ((use_shifted? 1 : 0) + 2*(doVertexScaled?1:0)) << std::endl;
  qa_info.glauber_mu[(use_shifted? 1 : 0) + 2*(doVertexScaled?1:0)] = mu;  
  qa_info.glauber_k[(use_shifted? 1 : 0) + 2*(doVertexScaled?1:0)] = k;  
  qa_info.glauber_chi2[(use_shifted? 1 : 0) + 2*(doVertexScaled?1:0)] = chi2;  

  double alpha = 1.00;
  double particlealpha = 1.00;

  //TF1 *trigeffcurveVertex[nvertexbins];
  TH1F *hnbd = new TH1F("hnbd","hnbd",nhistbins,-0.5,maxrange);
  for (int ibbc=0;ibbc<nhistbins;ibbc++) {
    float nbdreturn = NBD_getValue(ibbc, bestmu, bestk);
    //      cout << "hnbdcheck ibbc = " << ibbc << " bestmu,k = " << bestmu << " " << bestk << " nbd = " << nbdreturn << endl;
    hnbd->SetBinContent(ibbc+1, nbdreturn);
  }

  TH1D *hSimMBD = (TH1D *) hRealMBD->Clone();
  TH1D *hSimMBDfine = (TH1D *) hRealMBDfine->Clone();
  TH1D *hSimMBDVertex[nvertexbins];
  hSimMBD->SetName("hSimMBD");
  hSimMBD->Reset();
  hSimMBDfine->SetName("hSimMBDfine");
  hSimMBDfine->Reset();

  for (int ib = 0; ib < nvertexbins; ib++)
    {
      //trigeffcurveVertex[ib] = new TF1(Form("trigeffcurve_%d", ib),"1.0-TMath::Exp(-pow((x/[0]), [1]))",0.0,maxrange);
      hSimMBDVertex[ib] = (TH1D*) hSimMBD->Clone();
      hSimMBDVertex[ib]->SetName(Form("hSimMBD_%d", ib));
      hSimMBDVertex[ib]->Reset();
    }

  //-------------------------------------------------------------------------------------------
  // now use the best values (or default values) and re-calculate things, including bias factors...
  // loop over nNcoll values
  for (int ib=2; ib<=hglauber->GetNbinsX(); ib++) {

    int    ncoll = (int) hglauber->GetBinCenter(ib);

    double event_weightfactor = hglauber->GetBinContent(ib);
    if (event_weightfactor <= 0) continue;
  
    // scale parameters by mu' = my * ncoll and k' = k * ncoll
    for (int ihit=0; ihit< 4*(20 * ncoll * (int) mu); ihit++) {
      double dhit = (double) ihit/4.;
      double nbdvalue = NBD_getValue(dhit, mu * (double) (pow(ncoll,alpha)), k * (double) (pow(ncoll,alpha)));
      hSimMBDfine->Fill((double) dhit, nbdvalue * event_weightfactor);

    }
  }    
  for (int ib=2; ib<=hglauber->GetNbinsX(); ib++) {

    int    ncoll = (int) hglauber->GetBinCenter(ib);

    double event_weightfactor = hglauber->GetBinContent(ib);
    if (event_weightfactor <= 0) continue;
  
    double hardprocess_weightfactor = pow(ncoll,particlealpha); // assume for now that hard process scales with ncoll

    // scale parameters by mu' = my * ncoll and k' = k * ncoll
    for (int ihit=0; ihit< (20 * ncoll * (int) mu); ihit++) {

      double nbdvalue = NBD_getValue(ihit, mu * (double) (pow(ncoll,alpha)), k * (double) (pow(ncoll,alpha)));
      hSimMBD->Fill((double) ihit, nbdvalue * event_weightfactor);
      hSimMBDNcoll->Fill((double) ihit, (double) ncoll, nbdvalue * event_weightfactor);
      hSimMBDHardUnbiased->Fill((double) ihit, nbdvalue * event_weightfactor * hardprocess_weightfactor);

    }
    
    if (false)
      {
	// now try biased case

	// IHIT REPRESENTS NUMBER OF HITS IN THE MBD DETECTOR - SO NEED TO LOOP UP TO HIGH ENOUGH VALUES TO GET THE TAIL
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
	    hSimMBDHardBiased->Fill((double) ihit, 
				    nbdvalue * hardnbdvalue * event_weightfactor * hardprocess_weightfactor);
	    hSimMBDNcollHard->Fill((double) ihit, (double) ncoll, 
				   nbdvalue * hardnbdvalue * event_weightfactor * hardprocess_weightfactor);
	  }
	}
      }

  } // end loop over possible Ncoll values

  double scalefactorfine = hRealMBD->Integral(0,-1,"width")/hRealMBDfine->Integral(0,-1,"width");
  hRealMBDfine->Scale(scalefactorfine);
  hSimMBD->Scale(hRealMBD->Integral(lowfit,nhistbins)/hSimMBD->Integral(lowfit,nhistbins));
  double simscalefactorfine = hSimMBD->Integral(0,-1,"width")/hSimMBDfine->Integral(0,-1,"width");
  hSimMBDfine->Scale(simscalefactorfine);


  TH1D *hRatio = new TH1D("hRatio","hRatio",nhistbins,-0.5, maxrange);
  for (int i=1;i<=nhistbins;i++) {
    if (hSimMBD->GetBinContent(i) > 0 && hRealMBD->GetBinContent(i) > 0) {
      hRatio->SetBinContent(i, hRealMBD->GetBinContent(i)/hSimMBD->GetBinContent(i));
      hRatio->SetBinError(i, (sqrt(hRealMBD->GetBinContent(i))/hSimMBD->GetBinContent(i)));
    }
  }
  TH1D *hRatiofine = new TH1D("hRatiofine","hRatio",nhistbinsfine,-0.5, maxrange);
  for (int i=1;i<=nhistbinsfine;i++) {
    if (hSimMBDfine->GetBinContent(i) > 0 && hRealMBDfine->GetBinContent(i) > 0) {
      hRatiofine->SetBinContent(i, hRealMBDfine->GetBinContent(i)/hSimMBDfine->GetBinContent(i));
      hRatiofine->SetBinError(i, (sqrt(hRealMBDfine->GetBinContent(i))/hSimMBDfine->GetBinContent(i)));
    }
  }


  TF1 *trigeffcurve = new TF1("trigeffcurve","1.0-TMath::Exp(-pow((x/[0]), [1]))",0.0,maxrange);
  trigeffcurve->SetParameters(0.732,0.532); // just initial value guesses
  trigeffcurve->SetParLimits(0,0.2,20.0);    
  trigeffcurve->SetParLimits(1,0.2,5.0);
  hRatio->Fit(trigeffcurve,"NDOR","",1.0,80.0);
  // also fit a flat line above say 20 -> 140 ==> but just constraint it exactly at 1.0
  flatline->SetParameter(0,1.0);
  hRatio->Fit(flatline,"NDOR","",(double)lowfit,highfit);
  
  // trigger efficiency from integral comparison
  double err_real; double err_sim;
  double trigeffintegral = hRealMBD->IntegralAndError(1, nhistbins, err_real)/hSimMBD->IntegralAndError(1, nhistbins, err_sim);
  double trigeff_err = trigeffintegral*sqrt(TMath::Power(err_real/hRealMBD->Integral(1, nhistbins), 2) + TMath::Power(err_sim/hSimMBD->Integral(1, nhistbins), 2));
  if (!silence)  cout << "Trigger Efficiency from Integrals = " << trigeffintegral << endl;

  // calculate eccentricities


      
  // for (int iv = 0; iv < nvertexbins; iv++)
  //   {
  //     double mmu = h_mu->GetBinContent(iv+2);
  //     double kk = h_k->GetBinContent(iv+2);
  //     for (int ib=2; ib<=hglauber->GetNbinsX(); ib++) 
  // 	{
	  
  // 	  int    ncoll = (int) hglauber->GetBinCenter(ib);
	  
  // 	  double event_weightfactor = hglauber->GetBinContent(ib);
  // 	  if (event_weightfactor <= 0) continue;
	  
  // 	  // scale parameters by mu' = my * ncoll and k' = k * ncoll
  // 	  for (int ihit=0; ihit< (20 * ncoll * (int) mmu); ihit++) 
  // 	    {
	      
  // 	      double nbdvalue = NBD_getValue(ihit, mmu * (double) (pow(ncoll,alpha)), kk * (double) (pow(ncoll,alpha)));
  // 	      hSimMBDVertex[iv]->Fill((double) ihit, nbdvalue * event_weightfactor);
  // 	    }	  
  // 	}
  //     hSimMBDVertex[iv]->Scale(hRealMBDVertex[iv]->Integral(lowfit,nhistbins)/hSimMBDVertex[iv]->Integral(lowfit,nhistbins));
  //     for (int i=1;i<=nhistbins;i++) {
  // 	if (hSimMBDVertex[iv]->GetBinContent(i) > 0 && hRealMBDVertex[iv]->GetBinContent(i) > 0) {
  // 	  hRatioVertex[iv]->SetBinContent(i, hRealMBDVertex[iv]->GetBinContent(i)/hSimMBDVertex[iv]->GetBinContent(i));
  // 	  hRatioVertex[iv]->SetBinError(i, (sqrt(hRealMBDVertex[iv]->GetBinContent(i))/hSimMBDVertex[iv]->GetBinContent(i)));
  // 	}
  //     }
      
  //     trigeffcurveVertex[iv]->SetParameters(0.732,0.532); // just initial value guesses
  //     trigeffcurveVertex[iv]->SetParLimits(0,0.2,20.0);    
  //     trigeffcurveVertex[iv]->SetParLimits(1,0.2,5.0);
  //     hRatioVertex[iv]->Fit(trigeffcurveVertex[iv],"NDOR","",1.0,80.0);
  //     // also fit a flat line above say 20 -> 140 ==> but just constraint it exactly at 1.0
  //     flatline->SetParameter(0,1.0);
  //     hRatioVertex[iv]->Fit(flatline,"NDOR","",(double)lowfit,highfit);
      
  //     // trigger efficiency from integral comparison
  //     double err_realv; double err_simv;
  //     double trigeffintegralv = hRealMBDVertex[iv]->IntegralAndError(1, nhistbins, err_realv)/hSimMBDVertex[iv]->IntegralAndError(1, nhistbins, err_simv);
  //     double trigeff_errv = trigeffintegralv*sqrt(TMath::Power(err_realv/hRealMBDVertex[iv]->Integral(1, nhistbins), 2) + TMath::Power(err_simv/hSimMBDVertex[iv]->Integral(1, nhistbins), 2));
  //     if (!silence)  cout << "Trigger Efficiency from Integrals = " << trigeffintegralv << endl;
      
  //   }
  
  
  qa_info.glauber_trig_eff[(use_shifted? 1: 0)+2*(doVertexScaled?1:0)] = trigeffintegral;  
  qa_info.glauber_trig_eff_err[(use_shifted? 1: 0)+2*(doVertexScaled?1:0)] = trigeff_err;  
  

  TH1D *hSimMBDwTrig = new TH1D("hSimMBDwTrig","hSimMBDwTrig",nhistbins,-0.5,maxrange);
  for (int i=1;i<=nhistbins;i++) {
    if (i==1) {
      hSimMBDwTrig->SetBinContent(i,0.0); // no chance to fire the trigger if no MBD south hits
    } else {
      hSimMBDwTrig->SetBinContent(i,hSimMBD->GetBinContent(i)*trigeffcurve->Eval(hSimMBD->GetBinCenter(i)));
    }
  }

  for (int i=1;i<=nhistbins;i++) { // loop over MBD hits
    for (int j=1;j<=ncollmax;j++) { // loop over Ncoll values
      if (i==1) {
	hSimMBDNcoll_wtrig->SetBinContent(i,j,0.0);
      } else {
	hSimMBDNcoll_wtrig->SetBinContent(i,j,
					  hSimMBDNcoll->GetBinContent(i,j)*trigeffcurve->Eval(hSimMBD->GetBinCenter(i)));
      }
    }
  }

  std::string   vextra = "";    
  if (doVertexScaled) vextra = "bal_";
  TString   extra = "";    
  if (doVertexScaled || use_shifted) extra = Form("%s%s_", (use_shifted ? "sca" :""), (doVertexScaled ? "bal":""));

  TString calib_file_name = Form("%s/mbdana_centrality_%s%d.root", env_calib, vextra.c_str(), runnumber);
  if (isSim)
    {
      calib_file_name = Form("%s/mbdana_centrality_%s%s.root", env_calib, vextra.c_str(), mc_map[runnumber].c_str());
    }

  float cent_high, cent_low;
  TFile *fcalib = new TFile(calib_file_name, "r");
  TNtuple *ts = (TNtuple*) fcalib->Get("tn_centrality");
  float centrality_high[100] = {0};
  float centrality_low[100] = {0};
  ts->SetBranchAddress("low", &cent_low);
  ts->SetBranchAddress("high", &cent_high);
  for (int i = 0; i < ts->GetEntries() ; i++)
    {
      ts->GetEntry(i);      
      centrality_high[i] = cent_high;
      centrality_low[i] = cent_low;
    }

  fcalib->Close();

  // TODO: Now the same case for HARD collisions

  // NPart for each centrality bin
  TF1 *fflat = new TF1("fflat","1.0",0.0,1.0);

  std::string name_tntuple = Form("%s/SOFTX-D-15-00001/%s", env_p, ntuple_file.c_str());

  TFile *ftglauber = new TFile(name_tntuple.c_str(), "r");
  if (!ftglauber) {
    std::cerr << "no tntuple file" << std::endl;
    
  }
  TNtuple *tauau = (TNtuple*) ftglauber->Get(ntuple_name.c_str());
  if (!tauau) {
    std::cerr << "no tntuple" << std::endl;
  }

  Float_t npart;
  Float_t B;
  Float_t Ecc[5];
  Float_t Psi[5];
  std::array<double, 100> a_ecc[5];
  std::array<double, 100> a_psi[5];
  std::array<double, 100> a_b;
  std::array<double, 100> a_count;

  a_count.fill(0);
  a_b.fill(0);

  double sum_ecc[5] = {0};;
  double sum_b = {0};
  double counter = 0.0;
  tauau->SetBranchAddress("Npart", &npart);
  tauau->SetBranchAddress("B", &B);

  for (int i = 0; i < 5; i++)
    {
      a_ecc[i].fill(0);
      a_psi[i].fill(0);
      tauau->SetBranchAddress(Form("Ecc%d", i+1), &Ecc[i]);
      tauau->SetBranchAddress(Form("Psi%d", i+1), &Psi[i]);
    }

  TH1D *h_ecc2[100];
  TH1D *h_ecc3[100];
  TH1D *h_b[100];
  TH1D *h_b_all = new TH1D("h_b_all",";B;", 100 ,0,20); 
  TH1D *h_ecc2_all = new TH1D("h_ecc2_all", ";Ecc2;", 100,0,1); 
  TH1D *h_ecc3_all = new TH1D("h_ecc3_all", ";Ecc3;", 100,0,1); 
  for (int i = 0; i < 100; i++)
    {
      h_b[i] = new TH1D(Form("h_b_%d", i), ";Cent;B", 100 ,0,20); 
      h_ecc2[i] = new TH1D(Form("h_ecc2_%d", i), ";Cent;Ecc2", 100,0,1); 
      h_ecc3[i] = new TH1D(Form("h_ecc3_%d", i), ";Cent;Ecc3", 100,0,1); 
    }
  ULong64_t nentries = tauau->GetEntries();
  if (nentries > 20000) nentries = 20000;
      
  for (int i = 0; i < nentries; i++)
    {
      tauau->GetEntry(i);

      double fakembd = 0;

      for (int in = 0; in < (int) npart ; in++) fakembd += hnbd->GetRandom();

      // did the trigger fire
      bool   trigfired = false;
      if (fflat->GetRandom() < trigeffcurve->Eval(fakembd)) trigfired = true;

      if (!trigfired) continue;

      int centbin = 0;
      for (int idiv = 0; idiv < ndivs; idiv++)
	{

	  if (fakembd > centrality_low[idiv])
	    {
	      centbin = idiv;
	      break;
	    }
	}
      h_b_all->Fill(B);
      h_ecc2_all->Fill(Ecc[1]);
      h_ecc3_all->Fill(Ecc[2]);
      sum_b+=B;
      counter += 1.0;
      a_count[centbin]+=1.0;
      a_b[centbin]+=B;
      h_b[centbin]->Fill(B);
      h_ecc2[centbin]->Fill(Ecc[1]);
      h_ecc3[centbin]->Fill(Ecc[2]);
      for (int ie = 0; ie < 5; ie++)
	{
	  a_ecc[ie][centbin] += (double) Ecc[ie];
	  a_psi[ie][centbin] += (double) Psi[ie];
	  sum_ecc[ie] += (double) Ecc[ie];
	}
    }

  for (int i = 0; i < 100; i++)
    {
      a_b[i] /= a_count[i];
      for (int ie = 0; ie < 5; ie++)
	{
	  a_ecc[ie][i] /= a_count[i];
	  a_psi[ie][i] /= a_count[i];
	}
    }


  double npartstore[100] = {0};

  TH1D *h_npart_cent[100];

  for (int icent = 0; icent < 95; icent++)
    {
      h_npart_cent[icent] = new TH1D(Form("h_npart_cent_%d",icent),"", ncollmax, -0.5, maxrangencoll);

      for (int ibbc = 1+ floor(centrality_low[icent]); ibbc <= 1 + (floor(centrality_high[icent])) ; ibbc++)
	{
	  for (int inpart = 1; inpart <= ncollmax;inpart++) 
	    {
	      h_npart_cent[icent]->Fill((double)inpart - 1.0, hSimMBDNcoll_wtrig->GetBinContent(ibbc, inpart));
	    }
	}
	  
      npartstore[icent] = h_npart_cent[icent]->GetMean();
      if (!silence) std::cout << " Centbin "<<icent<<" <NPart> = " << npartstore[icent] << std::endl;
    }

  if (!silence) std::cout  << " filling npart total"<<std::endl;
  TH1D *h_npart_total = new TH1D("h_npart_total","", ncollmax, -0.5, maxrangencoll);
  for (int ibbc = 1; ibbc <= nhistbins; ibbc++)
    {
      for (int inpart = 1 ; inpart <= ncollmax;inpart++) h_npart_total->Fill((double)inpart - 1.0, hSimMBDNcoll->GetBinContent(ibbc, inpart));
    }

  double npartall = h_npart_total->GetMean();

  if (!silence) 
    {
      std::cout<<"Centrality 0-100% has <Npart> = " << npartall << std::endl;
    }
  qa_info.glauber_npart[(use_shifted? 1: 0) + (doVertexScaled?1:0)*2] = npartall;  
  qa_info.glauber_ecc2[(use_shifted? 1: 0) + (doVertexScaled?1:0)*2] = h_ecc2_all->GetMean();  
  qa_info.glauber_ecc3[(use_shifted? 1: 0) + (doVertexScaled?1:0)*2] = h_ecc3_all->GetMean();  
  qa_info.glauber_b[(use_shifted? 1: 0) + (doVertexScaled?1:0)*2] = h_b_all->GetMean();  


  TString calib_file_name2 = Form("%s/mbdana_npart_%s%d.root", env_calib, extra.Data(), runnumber);
  if (isSim)
    {
      calib_file_name2 = Form("%s/mbdana_npart_%s%s.root", env_calib, extra.Data(), mc_map[runnumber].c_str());
    }
  if (systematics)
    {
      calib_file_name2 = Form("%s/mbdana_npart_sys_%s_%s%d.root", env_calib, sysname.c_str(), extra.Data(), runnumber);
    }
  if (!silence) std::cout  << " filling npart file "<<calib_file_name2<<std::endl;
  TFile *fcalib2 = new TFile(calib_file_name2, "recreate");
  TNtuple *ts2 = new TNtuple("tn_npart", "holds npart divisions", "bin:cent_low:cent_high:npart:ecc2:ecc3:b");
  for (int i = 0; i < 95 ; i++)
    {
      ts2->Fill(i+1, centrality_low[i], centrality_high[i], npartstore[i], a_ecc[1][i], a_ecc[2][i], a_b[i]);
    }

  fcalib2->Write();
  fcalib2->Close();

  if (!silence) std::cout  << " filling file "<<std::endl;
  TString fname = Form("%s/mbdana_centrality_trigeff_%s%d.root", env_out, extra.Data(), runnumber);;
  if (isSim)
    {
      fname = Form("%s/mbdana_centrality_trigeff_%s%s.root", env_out, extra.Data(), mc_map[runnumber].c_str());;
    }
  if (!silence) std::cout  << " filling file "<<std::endl;
  TFile *fout = new TFile(fname,"recreate");


  hSimMBDwTrig->Write();
  hSimMBD->Write();
  hRealMBD->Write();
  hRealMBDfine->Write();
  hRatio->Write();
  hRatiofine->Write();
  trigeffcurve->Write();
  h_npart_total->Write();
  h_b_all->Write();
  h_ecc2_all->Write();
  h_ecc3_all->Write();
  
  for (int icent = 0; icent < 95;icent++) 
    {
      h_npart_cent[icent]->Write();
      h_b[icent]->Write();
      h_ecc2[icent]->Write();
      h_ecc3[icent]->Write();
    }
  fout->Close();


} // end routine

int QA_centrality::loadTree()
{
  
  m_file = new TFile(Form("%s/mbd_trees_%d.root", env_tree, qa_info.runnumber), "r");

  if (!m_file)
    {
      std::cout << " No Tree File found " <<std::endl;
      return 1;
    }

  m_ttree = (TTree*) m_file->Get("T");
  //REMOVE
  if (!m_ttree)
    {
      return 1;
    }
  
  m_ttree->SetBranchAddress("gl1_scaled", &gl1_scaled);
  m_ttree->SetBranchAddress("gl1_live", &gl1_live);    
  m_ttree->SetBranchAddress("minbias", &minbias);
  m_ttree->SetBranchAddress("bunch_number", &bunch_number);
  m_ttree->SetBranchAddress("zdc_sum",zdc_sum);
  m_ttree->SetBranchAddress("zdc_energy",zdc_energy);
  m_ttree->SetBranchAddress("mbd_charge_sum",mbd_charge_sum);
  m_ttree->SetBranchAddress("mbd_charge",mbd_charge);
  m_ttree->SetBranchAddress("mbd_time",mbd_time);
  m_ttree->SetBranchAddress("mbd_vertex",&mbd_vertex);
  return 0;
}
  
  // with parameters (mu, k) and for given n.

