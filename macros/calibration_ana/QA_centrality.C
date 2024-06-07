/* QA_centrality */
// Author: Daniel Lis - August 15 - 2023
// Brief : This macro class gives a quick QA analysis
//         of a run in sPHENIX
//      To be run on the output of the CentralityReco Module

//void QA_FindCentralities()

#include "QA_centrality.h"

using namespace std;
namespace fs = std::filesystem;

class QA_centrality;

/* constructor -- for the qa_info struct to hold all QA stats */
QA_centrality::QA_centrality(int silent, int debugger)
{

  silence  = silent;
  debug = debugger;

/* Setting path for output files and input files */

  env_p = new char[200];
  sprintf(env_p,"%s",std::getenv("MBD_CENTRALITY_CALIB_PATH"));
  
  if(!env_p)
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
  qa_info.min_bias = 0.;
  qa_info.vertex_cut = 0.;
  qa_info.vertex = 0.;
  qa_info.vertex_sigma = 0.;
  qa_info.charge_sum_mean_ns = 0.;
  for (int i = 0; i < 4;i++)
    {
      qa_info.charge_sum_mean[i] = 0.;
      qa_info.centrality_divs_chi2[i] = 0.;
    }
  for (int i = 0; i < 2; i++)
    {
      qa_info.glauber_mu[i] = 0.;
      qa_info.glauber_k[i] = 0.;
      qa_info.glauber_chi2[i] = 0.;
      qa_info.glauber_npart[i] = 0.;
      qa_info.glauber_trig_eff[i] = 0.;
      qa_info.glauber_trig_eff_err[i] = 0.;
      qa_info.glauber_chi2_forced[i] = 0.;
      qa_info.glauber_trig_eff_forced[i] = 0.;
      qa_info.glauber_trig_eff_forced_err[i] = 0.;
      qa_info.glauber_npart_forced[i] = 0.;
    }
  mc_map[0] = "hijing";
  mc_map[1] = "ampt";
  mc_map[2] = "epos";
    
}

void QA_centrality::Print_QA_Info(bool use_table)
{
  std::cout << " ***************** Run " << qa_info.runnumber << " ***************** " << std::endl;
  std::cout << "   Quality: "<<qa_info.quality<<std::endl; 
  std::cout << "   ZDC percent: "<<qa_info.ZDC_percentage<<std::endl; 
  std::cout << "   Vertex = " << qa_info.vertex << " +/- " <<qa_info.vertex_sigma << std::endl;
  std::cout << "   Charge Sum Mean: "<<std::endl;
  std::cout << "       Regular\tBalance\tScaled\tBoth " <<std::endl;
  std::cout << "       ";
  for (int i = 0; i < 4; i++) std::cout << qa_info.charge_sum_mean[i] << "\t";
  std::cout << " "<<std::endl;
  std::cout << "   Charge Sum Scaling factor: "<< qa_info.scale_factor << std::endl;
  std::cout << "   Glauber Fit Status (scaled): "<<std::endl;
  std::cout << "       mu           =  " << qa_info.glauber_mu[0]<<" ( "<<qa_info.glauber_mu[1]<<" )"<<std::endl;
  std::cout << "       k            =  " << qa_info.glauber_k[0]<<" ( "<<qa_info.glauber_k[1]<<" )"<<std::endl;
  std::cout << "       trigger eff. =  " << qa_info.glauber_trig_eff[0]<<" ( "<<qa_info.glauber_trig_eff[1]<<" )"<<std::endl;
  std::cout << "       chi2         =  " << qa_info.glauber_chi2[0]<<" ( "<<qa_info.glauber_chi2[1]<<" )"<<std::endl;
  std::cout << "       npart         =  " << qa_info.glauber_npart[0]<<" ( "<<qa_info.glauber_npart[1]<<" )"<<std::endl;
  std::cout << "     ********************************"<<std::endl;
  std::cout << "       Forced (scaled) :"<<std::endl;
  std::cout << "         trigger eff. =  " << qa_info.glauber_trig_eff_forced[0] << " ( " << qa_info.glauber_trig_eff_forced[1] << " )"<<std::endl;
  std::cout << "         chi2         =  " << qa_info.glauber_chi2_forced[0] << " ( " << qa_info.glauber_chi2_forced[1] << " )"<<std::endl;
  std::cout << "       npart         =  " << qa_info.glauber_npart_forced[0]<<" ( "<<qa_info.glauber_npart_forced[1]<<" )"<<std::endl;
  std::cout << "     ********************************"<<std::endl;
  std::cout << "   Centrality Div Fit:"<<std::endl;
  std::cout << "       Regular\tBalance\tScaled\tBoth " <<std::endl;
  std::cout << "       ";

  for (int i = 0; i < 4; i++) std::cout << qa_info.centrality_divs_chi2[i] << "\t";
  std::cout << " "<<std::endl;
  std::cout << " ********************************************* " << std::endl;
  if (use_table)
    {
      // runnumber,events,ZDC,vtx,vtxsig,mean,mean_sca,scale,chi2,mu,k,trigeff,trigeff_sca_for,chi2scafor,centchi2,centchi2sca
      std::cout << " Table Entry: "<<std::endl;
      std::cout << std::fixed << std::setprecision(5);
      std::cout << qa_info.runnumber << ", ";
      std::cout << qa_info.nevents << ", ";
      std::cout << qa_info.ZDC_percentage << ", ";
      std::cout << qa_info.min_bias << ", ";
      std::cout << qa_info.vertex_cut << ", ";
      std::cout << qa_info.vertex << ", ";
      std::cout << qa_info.vertex_sigma << ", ";
      std::cout << qa_info.charge_sum_mean[0] << ", ";
      std::cout << qa_info.charge_sum_mean[2] << ", ";
      std::cout << qa_info.scale_factor << ", ";
      std::cout << qa_info.charge_sum_mean_ns << ", ";
      std::cout << qa_info.glauber_chi2[1] << ", ";
      std::cout << qa_info.glauber_mu[1] << ", ";
      std::cout << qa_info.glauber_k[1] << ", ";
      std::cout << qa_info.glauber_trig_eff[1] << ", ";
      std::cout << qa_info.glauber_trig_eff_err[1] << ", ";
      std::cout << qa_info.glauber_trig_eff_forced[1] << ", ";
      std::cout << qa_info.glauber_trig_eff_forced_err[1] << ", ";
      std::cout << qa_info.glauber_npart[1] << ", ";
      std::cout << qa_info.glauber_npart_forced[1] << ", ";
      std::cout << qa_info.glauber_chi2_forced[1] <<", ";
      std::cout << qa_info.centrality_divs_chi2[2] <<", ";
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
  QA_ZDCCheck(runnumber);
// Go through and make distributions for MBD channels (charge and time)
  QA_MBDChannels(runnumber);
// Charge sum plots are made here with scaled and with all cuts
  QA_MakeChargeSum(runnumber);
// Centrality Calibrations + NBD+Glauber Fit is here
  
  QA_MakeCentralityCalibrations(runnumber, 1);

  QA_MakeDivisions(runnumber);  
  //CentralityCalibrations(runnumber, 1);
  SetReferenceRun(runnumber);
  QA_CentralityCheck(runnumber);
  return;
}
void QA_centrality::QA_MC(
			  const std::string mc_generator
		   )
{
  cthresh = 0.25;
  isSim = true;
  if (strcmp(mc_generator.c_str(), "hijing") == 0)
    qa_info.runnumber = 0;
  
  else if (strcmp(mc_generator.c_str(), "ampt") == 0)
    qa_info.runnumber = 1;
  
  else if (strcmp(mc_generator.c_str(), "epos") == 0)
    qa_info.runnumber = 2;
  else
      return;

  int runnumber = qa_info.runnumber;
  std::cout << mc_map[runnumber] <<" : " <<runnumber<<std::endl;

// Go through and make distributions for MBD channels (charge and time)
  QA_MBDChannels(runnumber);
// Charge sum plots are made here with scaled and with all cuts
  QA_MakeChargeSum(runnumber);
// Centrality Calibrations + NBD+Glauber Fit is here
  //QA_MakeCentralityCalibrations(runnumber, 1);
  QA_MakeDivisions(runnumber);

  return;
}

void QA_centrality::Start_QA_Centrality(
					const int runnumber
		   )
{

  qa_info.runnumber = runnumber;
// ZDC Check for percentage that meets certain cuts.
  QA_ZDCCheck(runnumber);
// Go through and make distributions for MBD channels (charge and time)
  QA_MBDChannels(runnumber);
// Charge sum plots are made here with scaled and with all cuts
  QA_MakeChargeSum(runnumber);

// Centrality Calibrations + NBD+Glauber Fit is here
//  QA_MakeCentralityCalibrations(runnumber, 1);
//  QA_MakeCentralityCalibrations(runnumber, 1, true, false);
  SetTrigEffMUK(.91, 4.39375, 0.859025);
  QA_MakeCentralityCalibrations(runnumber, 1, true, false);
//  QA_MakeCentralityCalibrations(runnumber, 1, true, false);
// Check the centrality bins

  QA_CentralityCheck(runnumber);

  return;
}

void QA_centrality::QA_CentralityCheck(const int runnumber)
{


  // Load in charge_sum distributions

  TFile *fchargesum = new TFile(Form("%s/output/plots/mbdana_charge_sum_%d.root", env_p, runnumber), "r");
  if (!fchargesum)
    {
      cout << "No file exists for this charge sum plot... exiting" << endl;
      return;
    } 
  TFile *fchargesumref = new TFile(Form("%s/output/plots/mbdana_charge_sum_%d.root", env_p, reference_run), "r");
  if (!fchargesumref)
    {
      cout << "No file exists for this charge sum plot... exiting" << endl;
      return;
    } 

  TH1D *h_run_charge_sum = (TH1D*) fchargesum->Get("h_charge_sum_min_bias_w_vertex_cut");
  if (!h_run_charge_sum)
    {
      std::cout << "nah fam" << std::endl;
    }
  TH1D *h_run_charge_sum_ref = (TH1D*) fchargesumref->Get("h_charge_sum_min_bias_w_vertex_cut");
  if (!h_run_charge_sum_ref)
    {
      std::cout << "nah fam no ref" << std::endl;
    }

  float scale_factor = h_run_charge_sum_ref->GetMean()/h_run_charge_sum->GetMean();
  if (debug) std::cout << __LINE__ << std::endl;

  char *path = new char[100];
  if (isSim)
    {
      sprintf(path, "%s/mbdana/mbd_trees_%s.root" , mc_map[runnumber].c_str(), mc_map[runnumber].c_str());
    }
  else
    {
      sprintf(path, "run%d/mbdana/mbd_trees_%d.root" ,runnumber, runnumber);
    }
  
  TFile *file = new TFile(Form("%s/output/%s", env_tree, path), "r");

  if (!file)
    {
      std::cout << " No Tree File found " <<std::endl;
      return;
    }

  // Making fresh histograms


  float zdc_sum[2];
  float mbd_vertex;
  float mbd_charge_sum[2];
  float mbd_charge[128];
  float mbd_time[128];


  std::vector<float> v_mbd_charge_sum{};
  TTree *t = (TTree*)file->Get("T");

  t->SetBranchAddress("zdc_sum",zdc_sum);
  t->SetBranchAddress("mbd_charge_sum",mbd_charge_sum);
  t->SetBranchAddress("mbd_charge",mbd_charge);
  t->SetBranchAddress("mbd_time",mbd_time);
  t->SetBranchAddress("mbd_vertex",&mbd_vertex);

  //// vertex
  int hits_n = 0;
  int hits_s = 0;      
  float nsum = 0;
  float ssum = 0;
  if (!silence)   cout <<"NEvents = "<<t->GetEntries()<<std::endl;

  for (int i = 0 ; i < (debug ? 10 : t->GetEntries()); i++)
    {
      bool minbias = false;
      hits_n = 0;
      hits_s = 0;
      ssum = 0;
      nsum = 0;
      t->GetEntry(i);

      if (!silence) if (i%20000 == 0) std::cout <<" ------------------------------- Event: "<<i<< " -------------------------------"<<std::endl;
      if (!silence) if (i%20000 == 0) std::cout <<" ---------   " << mbd_charge_sum[0] << " + " <<mbd_charge_sum[1]<< "    ----"<<std::endl;
      for (int ich = 0 ; ich < 64; ich++)
	{
	  if (mbd_time[ich] < 25. && mbd_charge[ich] > cthresh)
	    {
	      hits_s++;
	      ssum+=mbd_charge[ich];
	    }
	  if (mbd_time[ich+64] < 25. && mbd_charge[ich+64] > cthresh)
	    {
	      hits_n++;
	      nsum += mbd_charge[ich + 64];
	    } 
	}
      if (hits_s < 2 || hits_n < 2) continue; 
      if (fabs(mbd_vertex) > z_cut) continue;
      if (!isSim && ssum > charge_sum_south_cut && nsum < charge_sum_north_cut) continue;
      if (hasZDC && !isSim && (zdc_sum[0] <= zdc_cut || zdc_sum[1] <= zdc_cut)) continue;
      v_mbd_charge_sum.push_back(ssum + nsum);
    }
  int size = v_mbd_charge_sum.size();
  std::sort(v_mbd_charge_sum.begin(), v_mbd_charge_sum.end());
  // Load in centrality calibrations

  TFile *fcentralitycalib = new TFile(Form("%s/calib/mbdana_centrality_%d.root", env_p, reference_run), "r");
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
  
  tn_centrality->SetBranchAddress("bin", &centbin);
  tn_centrality->SetBranchAddress("low", &low);
  tn_centrality->SetBranchAddress("high", &high);
  
  float centrality_bins[100];
  for (int i = 0 ; i < 100; i++)
    centrality_bins[i] = 0.;

  // filling centrality bins

  int cent_divs = tn_centrality->GetEntries();
  for (int i = 0; i < cent_divs; i++)
    {
      tn_centrality->GetEntry(i);
      std::cout <<" Calib Entry " << i <<" : "<<  low << std::endl;
      centrality_bins[i] = low;
    }


  TH1D *h_cent_bin = new TH1D("hcent_bins","", ndivs, -0.5, (float) ndivs - 0.5);
  TH1D *h_cent_bin_sca = new TH1D("hcent_bins_sca","", ndivs, -0.5, (float) ndivs - 0.5);


  int i_cent_bin = 0;
  int count = 0;
  int i_cent_bin_sca = 0;
  int count_sca = 0;
  std::reverse(v_mbd_charge_sum.begin(), v_mbd_charge_sum.end());
  std::cout << "scale factor " << scale_factor<<std::endl;
  for (auto &sum : v_mbd_charge_sum)
    {
      if (sum*scale_factor < centrality_bins[i_cent_bin_sca])
	{
	  h_cent_bin_sca->SetBinContent(i_cent_bin_sca, count_sca);
	  h_cent_bin_sca->SetBinError(i_cent_bin_sca, sqrt(count_sca));
	  std::cout << i_cent_bin_sca << " Found "<<count_sca <<std::endl;
	  count_sca = 0;
	  i_cent_bin_sca++;
	}
      count_sca++;
    }

  for (auto &sum : v_mbd_charge_sum)
    {
      if (sum < centrality_bins[i_cent_bin])
	{
	  h_cent_bin->SetBinContent(i_cent_bin, count);
	  h_cent_bin->SetBinError(i_cent_bin, sqrt(count));
	  std::cout << i_cent_bin << " Found "<<count <<std::endl;
	  count = 0;
	  i_cent_bin++;
	}
      count++;	  
    }

  TF1 *flatline = new TF1("flatline","[0]",-0.5, (float) ndivs - .5);
  float match = 1./((float)ndivs);
  flatline->SetParameter(0, match);

  h_cent_bin->Scale(1./(static_cast<float>(h_cent_bin->Integral())));
  h_cent_bin->Fit(flatline,"NDORQ","");

  float chi2 = flatline->GetChisquare()/flatline->GetNDF();
  float level = flatline->GetParameter(0);;


  flatline->SetParameter(0,match);

  h_cent_bin_sca->Scale(1./(static_cast<float>(h_cent_bin_sca->Integral())));
  h_cent_bin_sca->Fit(flatline,"QNDOR","");

  float chi2_sca = flatline->GetChisquare()/flatline->GetNDF();
  float level_sca = flatline->GetParameter(0);;
  flatline->SetParameter(0,match);

  flatline->SetParameter(0,match);

  if (!silence)
    {
      for (int i = 0; i < h_cent_bin->GetNbinsX();i++)
	{
	  cout << i*(100/ndivs) <<"% : "<<h_cent_bin->GetBinContent(i+1) <<" -- " << h_cent_bin_sca->GetBinContent(i+1)<<endl;
	}
      cout << " ************************** " <<endl;
      cout << "  Run "<< runnumber<<endl;
      cout << "    Avg / Chi2         =  "<< level <<" / "<<chi2<<endl;
      cout << "    Avg / Chi2 sca     =  "<< level_sca <<" / "<< chi2_sca<<endl;
      cout << " ************************** " <<endl;
    }

  qa_info.centrality_divs_chi2[0] = chi2;
  qa_info.centrality_divs_chi2[2] = chi2_sca;

  TFile *fout = new TFile(Form("%s/output/plots/mbdana_centrality_check_%d.root", env_p, runnumber), "recreate");
  h_cent_bin->Write();
  h_cent_bin_sca->Write();

  fout->Close();

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
  
  if (reference_run)
    {
      if (!silence) cout << "reference_run: " <<reference_run<< std::endl;
      TFile *fout = new TFile(Form("%s/output/plots/mbdana_charge_sum_%d.root", env_p, reference_run), "r");
      if (!fout) return;
      TH1D *h_shift = (TH1D*) fout->Get("h_charge_sum_min_bias_w_vertex_cut");
      if (!h_shift) return;
      h_shift->SetName("h_charge_sum_shift");
      mbd_charge_factor = h_shift->GetMean();
      if (!silence) std::cout << "MBD Scale Factor: "<<mbd_charge_factor<<std::endl;
      TH1D *h_shift2_n = (TH1D*) fout->Get("h_charge_sum_min_bias_w_vertex_cut_n");
      TH1D *h_shift2_s = (TH1D*) fout->Get("h_charge_sum_min_bias_w_vertex_cut_s");
      mbd_charge_factor_bal = h_shift2_n->GetMean() / h_shift2_s->GetMean();
    }
  
  float low_t = -20.;
  float high_t = 20.;

  float binlengtht = 0.01;

  int nbins_t = floor((high_t - low_t)/binlengtht);

  int nbin = 2500;
  int nbin_fine = 1000000;
  int maxrange = 2500;
  if (isSim)
    { 
      nbin = 2500;
      maxrange = 2500;
    }

  bool minbias = false;

  TH2D *h2_charge_sum_v_zdc = new TH2D("h_charge_sum_v_zdc","", nbin/10, -0.5, (float)maxrange - 0.5, 1000, 0, 10000);
  TH2D *h2_charge_sum_v_zdc_mb = new TH2D("h_charge_sum_v_zdc_mb","", nbin/10, -0.5, (float)maxrange - 0.5, 1000, 0, 10000);

  TProfile *hp_charge_sum_v_zdc = new TProfile("hp_charge_sum_v_zdc","", 200, -0.5,  1999.5);
  TProfile *hp_charge_sum_v_zdc_mb = new TProfile("hp_charge_sum_v_zdc_mb","", 200, -0.5, 1999.5);

  TH2D *h2_charge_sum_ns = new TH2D("h2_charge_sum_ns","", nbin, -0.5, (float)maxrange - 0.5, nbin, -0.5, (float)maxrange - 0.5);  
  TH2D *h2_charge_sum_ns_cut = new TH2D("h2_charge_sum_ns_cut","", nbin, -0.5, (float)maxrange - 0.5, nbin, -0.5, (float)maxrange - 0.5);  

  TH1D *h_charge_sum = new TH1D("h_charge_sum","", nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_vtx = new TH1D("h_charge_sum_vtx","",  nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_min_bias = new TH1D("h_charge_sum_min_bias","", nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_min_bias_w_vertex_cut = new TH1D("h_charge_sum_min_bias_w_vertex_cut","", nbin, -0.5, (float)maxrange - 0.5);

  TH1D *h_charge_sum_fine = new TH1D("h_charge_sum_fine","", nbin_fine, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_fine_vtx = new TH1D("h_charge_sum_fine_vtx","",  nbin_fine, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_fine_min_bias_w_vertex_cut = new TH1D("h_charge_sum_fine_min_bias_w_vertex_cut","", nbin_fine, -0.5, (float)maxrange - 0.5);

  TH1D *h_charge_sum_balanced = new TH1D("h_charge_sum_balanced","", nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_vtx_balanced = new TH1D("h_charge_sum_vtx_balanced","",  nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_min_bias_w_vertex_cut_balanced = new TH1D("h_charge_sum_min_bias_w_vertex_cut_balanced","", nbin, -0.5, (float)maxrange - 0.5);

  TH1D *h_charge_sum_fine_balanced = new TH1D("h_charge_sum_fine_balanced","", nbin_fine, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_fine_vtx_balanced = new TH1D("h_charge_sum_fine_vtx_balanced","",  nbin_fine, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_fine_min_bias_w_vertex_cut_balanced = new TH1D("h_charge_sum_fine_min_bias_w_vertex_cut_balanced","", nbin_fine, -0.5, (float)maxrange - 0.5);

  TH1D *h_charge_sum_scaled = new TH1D("h_charge_sum_scaled","", nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_vtx_scaled = new TH1D("h_charge_sum_vtx_scaled","",  nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_min_bias_w_vertex_cut_scaled = new TH1D("h_charge_sum_min_bias_w_vertex_cut_scaled","", nbin, -0.5, (float)maxrange - 0.5);

  TH1D *h_charge_sum_fine_scaled = new TH1D("h_charge_sum_fine_scaled","", nbin_fine, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_fine_vtx_scaled = new TH1D("h_charge_sum_fine_vtx_scaled","",  nbin_fine, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_fine_min_bias_w_vertex_cut_scaled = new TH1D("h_charge_sum_fine_min_bias_w_vertex_cut_scaled","", nbin_fine, -0.5, (float)maxrange - 0.5);

  TH1D *h_charge_sum_balanced_scaled = new TH1D("h_charge_sum_balanced_scaled","", nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_vtx_balanced_scaled = new TH1D("h_charge_sum_vtx_balanced_scaled","",  nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_min_bias_w_vertex_cut_balanced_scaled = new TH1D("h_charge_sum_min_bias_w_vertex_cut_balanced_scaled","", nbin, -0.5, (float)maxrange - 0.5);

  TH1D *h_charge_sum_fine_balanced_scaled = new TH1D("h_charge_sum_fine_balanced_scaled","", nbin_fine, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_fine_vtx_balanced_scaled = new TH1D("h_charge_sum_fine_vtx_balanced_scaled","",  nbin_fine, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_fine_min_bias_w_vertex_cut_balanced_scaled = new TH1D("h_charge_sum_fine_min_bias_w_vertex_cut_balanced_scaled","", nbin_fine, -0.5, (float)maxrange - 0.5);

  TH1D *h_charge_sum_ns[2];
  TH1D *h_charge_sum_vtx_ns[2];
  TH1D *h_charge_sum_min_bias_w_vertex_cut_ns[2];

  for (int i = 0; i < 2; i++)
    {
      h_charge_sum_ns[i] = new TH1D(Form("h_charge_sum_%s", (i ? "n":"s")),"", nbin, -0.5, (float)maxrange - 0.5);
      h_charge_sum_vtx_ns[i] = new TH1D(Form("h_charge_sum_vtx_%s", (i ? "n":"s")),"", nbin, -0.5, (float)maxrange - 0.5);
      h_charge_sum_min_bias_w_vertex_cut_ns[i] = new TH1D(Form("h_charge_sum_min_bias_w_vertex_cut_%s", (i ? "n":"s")),"", nbin, -0.5, (float)maxrange - 0.5);
    }


  char *path = new char[100];
  if (isSim)
    {
      sprintf(path, "%s/mbdana/mbd_trees_%s.root" , mc_map[runnumber].c_str(), mc_map[runnumber].c_str());
    }
  else
    {
      sprintf(path, "run%d/mbdana/mbd_trees_%d.root" ,runnumber, runnumber);
    }
  
  TFile *file = new TFile(Form("%s/output/%s", env_tree, path), "r");

  if (!file)
    {
      std::cout << " No Tree File found " <<std::endl;
      return;
    }

  // Making fresh histograms


  float zdc_sum[2];
  float mbd_vertex;
  float mbd_charge_sum[2];
  float mbd_charge[128];
  float mbd_time[128];


  TTree *t = (TTree*)file->Get("T");

  t->SetBranchAddress("zdc_sum",zdc_sum);
  t->SetBranchAddress("mbd_charge_sum",mbd_charge_sum);
  t->SetBranchAddress("mbd_charge",mbd_charge);
  t->SetBranchAddress("mbd_time",mbd_time);
  t->SetBranchAddress("mbd_vertex",&mbd_vertex);

  //// vertex
  int hits_n = 0;
  int hits_s = 0;      
  float nsum = 0;
  float ssum = 0;
  if (!silence)   cout <<"NEvents = "<<t->GetEntries()<<std::endl;

  for (int i = 0 ; i < (debug ? 10 : t->GetEntries()); i++)
    {
      minbias = false;
      hits_n = 0;
      hits_s = 0;      
      nsum = 0;
      ssum = 0;
      t->GetEntry(i);

      if (!silence) if (i%20000 == 0) std::cout <<" ------------------------------- Event: "<<i<< " -------------------------------"<<std::endl;
      if (!silence) if (i%20000 == 0) std::cout <<" ---------   " << mbd_charge_sum[0] << " + " <<mbd_charge_sum[1]<< "    ----"<<std::endl;
      for (int ich = 0 ; ich < 64; ich++)
	{
	  if (mbd_time[ich] < 25. && mbd_charge[ich] > cthresh)
	    {
	      hits_s++;
	      ssum += mbd_charge[ich];
	    }
	  if (mbd_time[ich+64] < 25. && mbd_charge[ich+64] > cthresh)
	    {
	      hits_n++;
	      nsum += mbd_charge[ich+64];
	    }
	} 
      if (!silence) if (i%20000 == 0) std::cout <<" ---------   " << hits_s<< " + " <<hits_n<< "    ----"<<std::endl;
      // don't continue if this does not work.
      if (hits_s < 2 || hits_n < 2) continue; 

      h2_charge_sum_ns->Fill(ssum, nsum);
      h_charge_sum->Fill(ssum + nsum);
      h_charge_sum_fine->Fill(ssum + nsum);

      for (int j = 0; j < 2; j++) h_charge_sum_ns[j]->Fill(mbd_charge_sum[j]);
      if (!silence) if (i%20000 == 0) std::cout <<" ---------   " << mbd_vertex<< "          ----"<<std::endl;		      

      if (!(ssum > charge_sum_south_cut && nsum < charge_sum_north_cut))
	{
	  h_charge_sum_min_bias->Fill(ssum + nsum);
	}
      if (fabs(mbd_vertex) > z_cut) continue;

      if (hasZDC && !isSim) 
	{
	  h2_charge_sum_v_zdc->Fill(ssum+nsum, zdc_sum[0] + zdc_sum[1]);
	  hp_charge_sum_v_zdc->Fill(ssum+nsum, zdc_sum[0] + zdc_sum[1]);
	}

      h_charge_sum_vtx->Fill(ssum + nsum);
      h_charge_sum_fine_vtx->Fill(ssum + nsum);
      for (int j = 0; j < 2; j++) h_charge_sum_vtx_ns[j]->Fill(mbd_charge_sum[j]);      

      if (!isSim && ssum > charge_sum_south_cut && nsum < charge_sum_north_cut) continue;

      h2_charge_sum_ns_cut->Fill(ssum, nsum);

      if (hasZDC && !isSim && (zdc_sum[0] <= zdc_cut || zdc_sum[1] <= zdc_cut)) continue;
      if (hasZDC && !isSim) 
	{
	  h2_charge_sum_v_zdc_mb->Fill(ssum+nsum, zdc_sum[0] + zdc_sum[1]);
	  hp_charge_sum_v_zdc_mb->Fill(ssum+nsum, zdc_sum[0] + zdc_sum[1]);
	}
      h_charge_sum_min_bias_w_vertex_cut->Fill(ssum + nsum);
      h_charge_sum_fine_min_bias_w_vertex_cut->Fill(ssum + nsum);
      for (int j = 0; j < 2; j++) h_charge_sum_min_bias_w_vertex_cut_ns[j]->Fill((j?nsum:ssum));      
      if (!silence) if (i%20000 == 0) std::cout <<" ---------   " << ssum + nsum << "          ----"<<std::endl;		      
    }


  float ns_scale = h_charge_sum_ns[1]->GetMean()/h_charge_sum_ns[0]->GetMean();

  float ns_scale_vtx = h_charge_sum_vtx_ns[1]->GetMean()/h_charge_sum_vtx_ns[0]->GetMean();
  float ns_scale_min_bias_w_vertex_cut = h_charge_sum_min_bias_w_vertex_cut_ns[1]->GetMean()/h_charge_sum_min_bias_w_vertex_cut_ns[0]->GetMean();
  ns_scale_min_bias_w_vertex_cut /= mbd_charge_factor_bal; 
  for (int i = 0 ; i < (debug ? 10 : t->GetEntries()); i++)
    {
      minbias = false;
      hits_n = 0;
      hits_s = 0;      
      nsum = 0;
      ssum = 0;
      t->GetEntry(i);

      if (!silence) if (i%20000 == 0) std::cout <<" ------------------------------- Event: "<<i<< " -------------------------------"<<std::endl;
      if (!silence) if (i%20000 == 0) std::cout <<" ---------   " << mbd_charge_sum[0] << " + " <<mbd_charge_sum[1]<< "    ----"<<std::endl;
      for (int ich = 0 ; ich < 64; ich++)
	{
	  if (mbd_time[ich] < 25. && mbd_charge[ich] > cthresh)
	    {
	      hits_s++;
	      ssum += mbd_charge[ich];
	    }
	  if (mbd_time[ich+64] < 25. && mbd_charge[ich+64] > cthresh)
	    {
	      hits_n++;
	      nsum += mbd_charge[ich+64];
	    }
	} 
// don't continue if this does not work.
      if (hits_s < 2 || hits_n < 2) continue; 
      
      h_charge_sum_balanced->Fill(nsum + ssum*ns_scale);
      h_charge_sum_fine_balanced->Fill(nsum + ssum*ns_scale);
		      
      if (fabs(mbd_vertex) > z_cut) continue;

      h_charge_sum_vtx_balanced->Fill(nsum + ssum*ns_scale_vtx);
      h_charge_sum_fine_vtx_balanced->Fill(nsum + ssum*ns_scale_vtx);
      if (hasZDC && (zdc_sum[0] <= zdc_cut || zdc_sum[1] <= zdc_cut)) continue;
      if (ssum > charge_sum_south_cut && nsum < charge_sum_north_cut) continue;
      h_charge_sum_min_bias_w_vertex_cut_balanced->Fill(nsum + ssum*ns_scale_min_bias_w_vertex_cut);
      h_charge_sum_fine_min_bias_w_vertex_cut_balanced->Fill(nsum + ssum*ns_scale_min_bias_w_vertex_cut);
    }

  float scale;
  float scale_bal;
  float scale_vtx;
  float scale_vtx_bal;
  float scale_mb_vtx;
  float scale_mb_vtx_bal;

  if (reference_run && !isSim)
    {
      std::cout << "Scaling charge sum."<< std::endl;

      scale = mbd_charge_factor/(h_charge_sum->GetMean());
      scale_bal = mbd_charge_factor_bal/(h_charge_sum_balanced->GetMean());
      scale_vtx = mbd_charge_factor/(h_charge_sum_vtx->GetMean());
      scale_vtx_bal = mbd_charge_factor_bal/(h_charge_sum_vtx_balanced->GetMean());
      scale_mb_vtx = mbd_charge_factor/(h_charge_sum_min_bias_w_vertex_cut->GetMean());
      scale_mb_vtx_bal = mbd_charge_factor_bal/(h_charge_sum_min_bias_w_vertex_cut_balanced->GetMean());
      int count[4] = {0};
      std::cout << "Scaling by " << scale_mb_vtx << std::endl;
      for (int i = 0 ; i < (debug ? 10 : t->GetEntries()); i++)
	{
	  hits_n = 0;
	  hits_s = 0;      
	  ssum=0;
	  nsum=0;

	  t->GetEntry(i);

	  if (!silence) if (i%20000 == 0) std::cout <<" ------------------------------- Event: "<<i<< " -------------------------------"<<std::endl;

	  for (int ich = 0 ; ich < 64; ich++)
	    {
	      if (mbd_time[ich] < 25. && mbd_charge[ich] > cthresh)
		{
		  ssum+=mbd_charge[ich];
		  hits_s++;
		}
	      if (mbd_time[ich+64] < 25. && mbd_charge[ich+64] > cthresh)
		{
		  nsum += mbd_charge[ich+64];
		  hits_n++;
		}	

	    } 

	  // don't continue if this does not work.
	  if (hits_s < 2 || hits_n < 2) continue; 
	  count[0]++;
	  h_charge_sum_scaled->Fill(scale*(ssum + nsum));
	  h_charge_sum_fine_scaled->Fill(scale*(ssum + nsum));
	  h_charge_sum_balanced_scaled->Fill(scale_bal*(ssum + nsum*ns_scale));
	  h_charge_sum_fine_balanced_scaled->Fill(scale_bal*(ssum + nsum*ns_scale));
		      
	  if (fabs(mbd_vertex) > z_cut) continue;
	  count[1]++;
	  h_charge_sum_vtx_scaled->Fill(scale_vtx*(ssum + nsum));
	  h_charge_sum_fine_vtx_scaled->Fill(scale_vtx*(ssum + nsum));
	  h_charge_sum_vtx_balanced_scaled->Fill(scale_vtx_bal*(ssum + nsum*ns_scale_vtx));
	  h_charge_sum_fine_vtx_balanced_scaled->Fill(scale_vtx_bal*(ssum + nsum*ns_scale_vtx));

	  if (ssum > charge_sum_south_cut && nsum < charge_sum_north_cut) continue;
	  count[2]++;
	  if (hasZDC && (zdc_sum[0] <= zdc_cut || zdc_sum[1] <= zdc_cut)) continue;
	  count[3]++;
	  h_charge_sum_min_bias_w_vertex_cut_scaled->Fill(scale_mb_vtx*(ssum + nsum));
	  h_charge_sum_fine_min_bias_w_vertex_cut_scaled->Fill(scale_mb_vtx*(ssum + nsum));
	  h_charge_sum_min_bias_w_vertex_cut_balanced_scaled->Fill(scale_mb_vtx_bal*(ssum + nsum*ns_scale_min_bias_w_vertex_cut));
	}
      std::cout << " Cuts passed : " << count[0] << " -- > " <<count[1] << " -- > " <<count[2] << " -- > " <<count[3] << " "<<std::endl;

 
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
  qa_info.charge_sum_mean[1] = h_charge_sum_min_bias_w_vertex_cut_balanced->GetMean();
  if (reference_run)
    {
      qa_info.charge_sum_mean[2] = h_charge_sum_min_bias_w_vertex_cut_scaled->GetMean();
      qa_info.charge_sum_mean[3] = h_charge_sum_min_bias_w_vertex_cut_balanced_scaled->GetMean();
      qa_info.scale_factor = scale_mb_vtx;
    }
  char *pathout = new char[100];
  if (isSim)
    {
      sprintf(pathout, "%s/output/plots/mbdana_charge_sum_%s.root" , env_p, mc_map[runnumber].c_str());
    }
  else
    {
      sprintf(pathout, "%s/output/plots/mbdana_charge_sum_%d.root" ,env_p, runnumber);
    }

  TFile *fout = new TFile(pathout, "RECREATE");
  TNtuple *ts = new TNtuple("tn_scale","scalefactor","run:scale");
  ts->Fill(runnumber, scale_mb_vtx);
  ts->Write();

  h2_charge_sum_ns->Write();
  h2_charge_sum_ns_cut->Write();
  h_charge_sum->Write();
  h_charge_sum_vtx->Write();
  h_charge_sum_min_bias_w_vertex_cut->Write();
  h_charge_sum_min_bias->Write();
  h_charge_sum_balanced->Write();
  h_charge_sum_vtx_balanced->Write();
  h_charge_sum_min_bias_w_vertex_cut_balanced->Write();
  h_charge_sum_fine->Write();
  h_charge_sum_fine_vtx->Write();
  h_charge_sum_fine_min_bias_w_vertex_cut->Write();
  h_charge_sum_fine_balanced->Write();
  h_charge_sum_fine_vtx_balanced->Write();
  h_charge_sum_fine_min_bias_w_vertex_cut_balanced->Write();

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


  h_charge_sum_fine_balanced_scaled->Write();
  h_charge_sum_fine_vtx_balanced_scaled->Write();
  h_charge_sum_fine_min_bias_w_vertex_cut_balanced_scaled->Write();
  h_charge_sum_fine_scaled->Write();
  h_charge_sum_fine_vtx_scaled->Write();
  h_charge_sum_fine_min_bias_w_vertex_cut_scaled->Write();
  h_charge_sum_balanced_scaled->Write();
  h_charge_sum_vtx_balanced_scaled->Write();
  h_charge_sum_min_bias_w_vertex_cut_balanced_scaled->Write();
  h_charge_sum_scaled->Write();
  h_charge_sum_vtx_scaled->Write();
  h_charge_sum_min_bias_w_vertex_cut_scaled->Write();

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

  delete h_charge_sum_fine;
  delete h_charge_sum_fine_vtx;
  delete h_charge_sum_fine_min_bias_w_vertex_cut;
  delete h_charge_sum_fine_balanced;
  delete h_charge_sum_fine_vtx_balanced;
  delete h_charge_sum_fine_min_bias_w_vertex_cut_balanced;
  delete h_charge_sum_fine_balanced_scaled;
  delete h_charge_sum_fine_vtx_balanced_scaled;
  delete h_charge_sum_fine_min_bias_w_vertex_cut_balanced_scaled;
  delete h_charge_sum_fine_scaled;
  delete h_charge_sum_fine_vtx_scaled;
  delete h_charge_sum_fine_min_bias_w_vertex_cut_scaled;
    
}
void QA_centrality::QA_MakeDivisions(const int runnumber)
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
      sprintf(path, "%s/mbdana/mbd_trees_%s.root" , mc_map[runnumber].c_str(), mc_map[runnumber].c_str());
    }
  else
    {
      sprintf(path, "run%d/mbdana/mbd_trees_%d.root" ,runnumber, runnumber);
    }
  
  TFile *file = new TFile(Form("%s/output/%s", env_tree, path), "r");

  if (!file)
    {
      std::cout << " No Tree File found " <<std::endl;
      return;
    }

  // Making fresh histograms


  float zdc_sum[2];
  float mbd_vertex;
  float mbd_charge_sum[2];
  float mbd_charge[128];
  float mbd_time[128];


  std::vector<float> v_mbd_charge_sum{};
  TTree *t = (TTree*)file->Get("T");

  t->SetBranchAddress("zdc_sum",zdc_sum);
  t->SetBranchAddress("mbd_charge_sum",mbd_charge_sum);
  t->SetBranchAddress("mbd_charge",mbd_charge);
  t->SetBranchAddress("mbd_time",mbd_time);
  t->SetBranchAddress("mbd_vertex",&mbd_vertex);

  //// vertex
  int hits_n = 0;
  int hits_s = 0;      
  float nsum = 0;
  float ssum = 0;
  if (!silence)   cout <<"NEvents = "<<t->GetEntries()<<std::endl;

  for (int i = 0 ; i < (debug ? 10 : t->GetEntries()); i++)
    {
      bool minbias = false;
      hits_n = 0;
      hits_s = 0;
      ssum = 0;
      nsum = 0;
      t->GetEntry(i);

      if (!silence) if (i%20000 == 0) std::cout <<" ------------------------------- Event: "<<i<< " -------------------------------"<<std::endl;
      if (!silence) if (i%20000 == 0) std::cout <<" ---------   " << mbd_charge_sum[0] << " + " <<mbd_charge_sum[1]<< "    ----"<<std::endl;
      if (isSim)
	{
	  v_mbd_charge_sum.push_back(mbd_charge_sum[0] + mbd_charge_sum[1]);
	}

      for (int ich = 0 ; ich < 64; ich++)
	{
	  if (mbd_time[ich] < 25. && mbd_charge[ich] > cthresh)
	    {
	      hits_s++;
	      ssum+=mbd_charge[ich];
	    }
	  if (mbd_time[ich+64] < 25. && mbd_charge[ich+64] > cthresh)
	    {
	      hits_n++;
	      nsum += mbd_charge[ich + 64];
	    } 
	}
      if (hits_s < 2 || hits_n < 2) continue; 
      if (fabs(mbd_vertex) > z_cut) continue;
      if (!isSim && ssum > charge_sum_south_cut && nsum < charge_sum_north_cut) continue;

      if (hasZDC && !isSim && (zdc_sum[0] <= zdc_cut || zdc_sum[1] <= zdc_cut)) continue;
      v_mbd_charge_sum.push_back(ssum + nsum);
    }
  int size = v_mbd_charge_sum.size();
  std::cout << " Starting Sort " << size << " events.";
  std::sort(v_mbd_charge_sum.begin(), v_mbd_charge_sum.end());
  std::cout << "done."<<std::endl;


  std::cout << "Divisions: "<<std::endl;

  TString calib_file_name = Form("%s/calib/mbdana_centrality_%d.root", env_p, runnumber);

  TFile *fcalib = new TFile(calib_file_name, "recreate");
  TNtuple *ts = new TNtuple("tn_centrality", "holds centrality divisions", "bin:low:high");
  int divs = 100;
  if (!isSim) divs = 91;
  for (int i = 1; i < 100 ; i++)
    {
      int bin = floor(i*size/divs);
      int bin1 = floor((i-1)*size/divs);
      ts->Fill(i, v_mbd_charge_sum[size - 1 - bin], v_mbd_charge_sum[size - 1 - bin1]);
    }

  fcalib->Write();
  fcalib->Close();

  for (int i = 0; i < divs;i++)
    {
      int bin = floor(i*size/divs);
      std::cout << " Division " << i+1 <<": " << v_mbd_charge_sum[size - 1 - bin]<<std::endl;
    }
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

  char *path = new char[100];
  if (isSim)
    {
      sprintf(path, "%s/mbdana/mbd_trees_%s.root" , mc_map[runnumber].c_str(), mc_map[runnumber].c_str());
    }
  else
    {
      sprintf(path, "run%d/mbdana/mbd_trees_%d.root" ,runnumber, runnumber);
    }
  
  TFile *file = new TFile(Form("%s/output/%s", env_tree, path), "r");
  
  if (!file)
    {
      std::cout << " No Tree File found " <<std::endl;
      return;
    }

  // Making fresh histograms

  //raw
  TH1D *h_charge[128];
  TH1D *h_time[128];

  //passes ZDC tube cuts
  TH1D *h_charge_min_bias[128];
  TH1D *h_time_min_bias[128];

  for (int j = 0; j < 128; j++)
    {
      h_charge[j] = new TH1D(Form("h_mbd_charge_ch%d", j), "",200, 0, 10);
      h_time[j] = new TH1D(Form("h_mbd_time_ch%d", j), "",2500, -25, 25);

      h_charge_min_bias[j] = new TH1D(Form("h_mbd_charge_min_bias_ch%d", j), "",200, 0, 10);
      h_time_min_bias[j] = new TH1D(Form("h_mbd_time_min_bias_ch%d", j), "",2500, -25, 25);
    }

  float zdc_sum[2];
  float mbd_vertex;
  float mbd_charge[128];
  float mbd_time[128];
  float mbd_sum[2];
  
  TTree *t = (TTree*)file->Get("T");

  t->SetBranchAddress("zdc_sum", zdc_sum);
  t->SetBranchAddress("mbd_charge",mbd_charge);
  t->SetBranchAddress("mbd_vertex",&mbd_vertex);
  t->SetBranchAddress("mbd_time",mbd_time);
  t->SetBranchAddress("mbd_charge_sum",mbd_sum);

  for (int i = 0; i < t->GetEntries();i++)
    {
      t->GetEntry(i);
      bool min_bias = true;
      if (!isSim) min_bias &= (zdc_sum[0] >zdc_cut && zdc_sum[1] > zdc_cut);
      min_bias &=  (fabs(mbd_vertex - qa_info.vertex) < 60);
      min_bias &= (!(mbd_sum[0] > 150 && mbd_sum[1] < 10)); 

      for (int ich = 0 ; ich < 128; ich++)
	{
	  h_charge[ich]->Fill(mbd_charge[ich]);
	  h_time[ich]->Fill(mbd_time[ich]);
	  if (min_bias)
	    {
	      h_charge_min_bias[ich]->Fill(mbd_charge[ich]);
	      h_time_min_bias[ich]->Fill(mbd_time[ich]);
	    }
	}
    }

  char *path_out = new char[100];
  if (isSim)
    {
      sprintf(path_out, "%s/output/plots/mbdana_channels_%s.root" , env_p, mc_map[runnumber].c_str());
    }
  else
    {
      sprintf(path_out, "%s/output/plots/mbdana_channels_%d.root" , env_p, runnumber);
    }
  
  TFile *fout = new TFile(path_out, "recreate");
  for (int ich = 0 ; ich < 128; ich++)
    {
      h_charge[ich]->Write();
      h_time[ich]->Write();
      h_charge_min_bias[ich]->Write();
      h_time_min_bias[ich]->Write();
    }

  fout->Close();			 
}


void QA_centrality::QA_ZDCCheck(const int runnumber)
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

  TFile *f = new TFile(Form("%s/output/run%d/mbdana/mbd_trees_%d.root", env_tree, runnumber, runnumber), "r");
  TTree *t = (TTree*) f->Get("T");

  TH1D *h_zdc_sum_n = new TH1D("h_zdc_sum_n", ";ZDC sum [GeV]; Counts", 1000, 0, 10000);
  TH1D *h_zdc_sum_s = new TH1D("h_zdc_sum_s", ";ZDC sum [GeV]; Counts", 1000, 0, 10000);

  TH2D *h_zdc_sum_ns = new TH2D("h_zdc_sum_ns", ";ZDC sum S [GeV]; ZDC sum N [GeV]", 1000, 0, 10000, 1000, 0, 10000);
  TH2D *h_mbd_hits_ns = new TH2D("h_mbd_hits_ns", ";nhit S; nhit N", 65, -0.5, 64.5, 65, -0.5, 64.5);
  TH1D *h_zdc_energy[6];
  
  for (int i = 0; i < 6; i++)
    {
      h_zdc_energy[i] = new TH1D(Form("h_zdc_energy_%d", i), ";ZDC Energy [GeV]; Counts", 1000, 0, 5000);
    }

  TH1D *h_vertex = new TH1D("h_vertex", ";vertex [cm]; Counts", 601, -300.5, 300.5);
  TH1D *h_time_zero = new TH1D("h_time_zero", ";time zero [ns]; Counts", 501, -25.05, 25.05);

  float zdc_sum[2];
  float mbd_sum[2];
  float zdc_energy[6];
  float mbd_charge[128];
  float mbd_time[128];

  float mbd_vertex;
  float mbd_time_zero;
  t->SetBranchAddress("zdc_sum", zdc_sum);
  t->SetBranchAddress("zdc_energy", zdc_energy);
  t->SetBranchAddress("mbd_charge", mbd_charge);
  t->SetBranchAddress("mbd_time", mbd_time);
  t->SetBranchAddress("mbd_charge_sum", mbd_sum);
  t->SetBranchAddress("mbd_vertex", &mbd_vertex);
  t->SetBranchAddress("mbd_time_zero", &mbd_time_zero);
  int n_tight_vertex = 0;;
  int n_two_hits = 0;;
  int n_min_bias = 0;;
  int n_events = 0;
  int n_vertex = 0;
  int n_passed_ns_cut = 0;
  for (int i = 0; i < t->GetEntries(); i++)
    {
      t->GetEntry(i);
      h_vertex->Fill(mbd_vertex);
      h_time_zero->Fill(mbd_time_zero);
    }
  h_vertex->Fit("gaus",(silence?"Q":""),"", -60, 60);
  float mean = h_vertex->GetFunction("gaus")->GetParameter(1);

  h_vertex->Fit("gaus",(silence?"Q":""),"", mean - 20, mean+20);
  qa_info.vertex = h_vertex->GetFunction("gaus")->GetParameter(1);
  qa_info.vertex_sigma = h_vertex->GetFunction("gaus")->GetParameter(2);

  if (!silence) std:cout << "Vertex found at "<<qa_info.vertex<<std::endl;

  for (int i = 0; i < t->GetEntries(); i++)
    {
      t->GetEntry(i);
      n_events++;

      int hits_n = 0;
      int hits_s = 0;

      for (int j = 0; j < 64; j++)
	{
	  if (mbd_charge[j] > 0. && mbd_time[j] < 25.)
	    hits_s++;
	  if (mbd_charge[j+64] > 0. && mbd_time[j+64] < 25.)
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
      if (mbd_sum[1] < charge_sum_north_cut && mbd_sum[0] > charge_sum_south_cut) continue;
      n_passed_ns_cut++;
      if (!(zdc_sum[0] > zdc_cut && zdc_sum[1] > zdc_cut)) continue;
      n_min_bias++;
      if (fabs(mbd_vertex) > 10) continue;
      n_tight_vertex++;
    }

  if (n_min_bias/n_passed_ns_cut < 0.1) hasZDC = false;
  if (!silence)
    {
      std::cout << "Total events: "<< t->GetEntries()<<endl;
      std::cout << "Number of events passing Two Hits Cut             : "<<n_two_hits <<"/"<<n_events<<" = "<<static_cast<double>(n_two_hits)/static_cast<double>(n_events)<<std::endl;
      std::cout << "Number of events passing Vertex Cut               : "<<n_vertex <<"/"<<n_events<<" = "<<static_cast<double>(n_vertex)/static_cast<double>(n_events)<<std::endl;
      std::cout << "Number of events passing ChargeSum Cut            : "<<n_passed_ns_cut <<"/"<<n_events<<" = "<<static_cast<double>(n_passed_ns_cut)/static_cast<double>(n_events)<<std::endl;
      std::cout << "Number of events passing ZDC cut after all others : "<<n_min_bias <<"/"<<n_passed_ns_cut<<" = "<<static_cast<double>(n_min_bias)/static_cast<double>(n_passed_ns_cut)<<std::endl;
      std::cout << "Number of min_bias events with vertex < 10         : "<<n_tight_vertex <<"/"<<n_min_bias<<" = "<<static_cast<double>(n_tight_vertex)/static_cast<double>(n_min_bias)<<std::endl;
    }

  qa_info.nevents = t->GetEntries();
  qa_info.vertex_cut = static_cast<double>(n_tight_vertex)/static_cast<double>(n_min_bias);
  qa_info.min_bias = static_cast<double>(n_min_bias)/static_cast<double>(n_events);
  qa_info.ZDC_percentage = static_cast<double>(n_min_bias)/static_cast<double>(n_passed_ns_cut);

  TFile *fout = new TFile(Form("%s//output/plots/mbdana_zdc_check_%d.root", env_p, runnumber),"recreate");
  for (int i = 0; i < 6; i++)
    {
      h_zdc_energy[i]->Write();
    }

  h_zdc_sum_n->Write();
  h_zdc_sum_s->Write();
  h_zdc_sum_ns->Write();
  h_mbd_hits_ns->Write();

  h_time_zero->Write();
  h_vertex->Write();
  fout->Close();
}

void QA_centrality::QA_MakeCentralityCalibrations(const int runnumber, const bool doVertexSelection, const bool use_shifted , const bool use_balanced)
{

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }
  ///////////// 
  double mu = 4.4;  //defaults that are used if not determining these parameters again
  double k  = 2.14;
  if (isSim)
    {
      mu *= 130;
    }
  bool flag_determineNBDparameters = true;  
  if (set_mu && set_k) 
    {
      mu = set_mu;  //defaults that are used if not determining these parameters again
      k  = set_k;
      flag_determineNBDparameters = false;
    }
  //===============================================================================================
  bool npartminusone = false;
  bool npartscaling = true;
  bool flag_determineBiasFactors = false;
  
  double forcetrigfrac = 100*trigeff;
  double alpha = 1.00;
  double particlealpha = 1.00;
  double biasfactor = 1.55;
  double biased_mu = mu * biasfactor; // larger number of hits
  double biased_k  =  k * biasfactor;
  std::cout << "ndivs: "<<ndivs<<std::endl;
  const int maxcentbins  = ndivs - 1;

  bool sphenix = true;
  int vtxflag = 1;
  std::string name = Form("%s/glauber/lime_glauber_auau200_100k_histo.root", env_p);
  //  std::string name = Form("%s/glauber/lemon_glauber_auau200_100k_histo.root", env_p);
  int lowfit = 100; // lowest end of standard fit range, was 30 in PHENIX
  int highfit = 2000; // lowest end of standard fit range, was 30 in PHENIX
  if (isSim)
    {
      lowfit = 10000;
      highfit = 250000;
    }
  // utyy
  //===============================================================================================

  //===============================================================================================
  // the below file contains the "hNcoll" histogram to be sampled from (step #1)
  TFile *fglauber = new TFile(name.c_str()); 
  TH1D  *hNcoll;
  hNcoll = (TH1D * ) fglauber->Get("hNpart"); // note that this is only npart in the Au nucleus
    // this creates a problem because there is never an Npart = 1, always Npart >= 2
    // so bias correction procedure has a flaw; so in this case shift down by one
  
  // the below file contains the "hbbcQs6" histogram of BBC data distribution (from z-vertex range)
  //=======================================================
  char *fdataname = new char[100];
  if (isSim)
    {
      sprintf(fdataname,"%s/output/plots/mbdana_charge_sum_%s.root", env_p, mc_map[runnumber].c_str());  
    }
  else
    {
      sprintf(fdataname,"%s/output/plots/mbdana_charge_sum_%d.root", env_p, runnumber);  
    }
  //  if (sphenix) sprintf(fdataname,"../output/run%d/trees_%d.root", runnumber, runnumber);
  //if (sphenix) fdataname = "./fout_sphenix2023_auau200_23722_both_z20.root";
  //  if (sphenix) fdataname = "./fout_sphenix2023_auau200_22979_both_z20.root";

  char runnum[100];
  if (isSim) sprintf(runnum, "%s", mc_map[runnumber].c_str());
  else  sprintf(runnum, "000%d", runnumber);
  char zselect[100] = "|z_{MBD}|< 60 cm";
  //=======================================================  
  
  TFile *fdata = new TFile(fdataname);

  TH1D *hRealBBC;
  TH1D *hRealBBC_fine;
  //if (which==0) hRealBBC = (TH1D *) fdata->Get("h_mbd_charge_ns_w_zdc_cut_w_mbd_cut_and_vertex_10"); // a particular z-vertex range
  
  hRealBBC = (TH1D *) fdata->Get("h_charge_sum_min_bias_w_vertex_cut"); // a particular z-vertex range
  hRealBBC_fine = (TH1D *) fdata->Get("h_charge_sum_fine_min_bias_w_vertex_cut"); // a particular z-vertex range

  if (!isSim)
    {
      if (qa_info.ZDC_percentage < 0.5) 
	{
	  hRealBBC = (TH1D *) fdata->Get("h_charge_sum_min_bias_w_vertex_cut"); // a particular z-vertex range
	  hRealBBC_fine = (TH1D *) fdata->Get("h_charge_sum_fine_min_bias_w_vertex_cut"); // a particular z-vertex range
	}
      else if (use_shifted && use_balanced) 
	{
	  hRealBBC = (TH1D *) fdata->Get("h_charge_sum_min_bias_w_vertex_cut_balanced_scaled"); // a particular z-vertex range
	  hRealBBC_fine = (TH1D *) fdata->Get("h_charge_fine_sum_min_bias_w_vertex_cut_balanced_scaled"); // a particular z-vertex range
	}
      else if (use_shifted)
	{
	  hRealBBC = (TH1D *) fdata->Get("h_charge_sum_min_bias_w_vertex_cut_scaled"); // a particular z-vertex range
	  hRealBBC_fine = (TH1D *) fdata->Get("h_charge_fine_sum_min_bias_w_vertex_cut_scaled"); // a particular z-vertex range
	}
      else if (use_balanced) 
	{
	  hRealBBC = (TH1D *) fdata->Get("h_charge_sum_min_bias_w_vertex_cut_balanced"); // a particular z-vertex range
	  hRealBBC_fine = (TH1D *) fdata->Get("h_charge_fine_sum_fine_min_bias_w_vertex_cut_balanced_scaled"); // a particular z-vertex range
	}
    }

  if (!hRealBBC)
    {
      std::cout << " no histograms"<<std::endl;
      return;
    }
  //hRealBBC->Scale(1./hRealBBC->Integral());
  int nhistbins = hRealBBC->GetNbinsX();
  int nhistbins_fine = hRealBBC_fine->GetNbinsX();
  // j.nagle - also change 199.5 to be maxrange = ((float)nhistbins) - 0.5
  float maxrange = ((float) nhistbins) - 0.5;
  // might be good to have max Ncoll / Npart as well... (CuAu = 197 + 63 = 270)
  // Au+Au case 197 + 197 = 396
  int ncollmax = 400; // really npart here
  float maxrangencoll = ((float) ncollmax) - 0.5;

  TH1D *hSimBBCHardUnbiased = new TH1D("hSimBBCHardUnbiased","hSimBBCHardUnbiased",nhistbins,-0.5,maxrange);
  TH1D *hSimBBCHardBiased = new TH1D("hSimBBCHardBiased","hSimBBCHardBiased",nhistbins,-0.5,maxrange);

  // for these two, the trigger is not required to have fired ... (BBC vs Ncoll)

  TH2D *hSimBBCNcoll = new TH2D("hSimBBCNcoll","hSimBBCNcoll",nhistbins,-0.5,maxrange,ncollmax,-0.5,maxrangencoll);
  TH2D *hSimBBCNcoll_fine = new TH2D("hSimBBCNcoll_fine","hSimBBCNcoll",nhistbins,-0.5,maxrange,ncollmax,-0.5,maxrangencoll);
  TH2D *hSimBBCNcollHard = new TH2D("hSimBBCNcollHard","hSimBBCNcollHard",nhistbins,-0.5,maxrange,ncollmax,-0.5,maxrangencoll);

  // this now includes the trigger efficiency turn-on and can be used to calculate <Ncoll> for each centrality bin!!!
  
  TH2D *hSimBBCNcoll_wtrig = new TH2D("hSimBBCNcoll_wtrig","hSimBBCNcoll_wtrig",nhistbins,-0.5,maxrange,ncollmax,-0.5,maxrangencoll);
  TH2D *hSimBBCNcoll_fine_wtrig = new TH2D("hSimBBCNcoll_fine_wtrig","hSimBBCNcoll_wtrig",nhistbins_fine,-0.5,maxrange,ncollmax,-0.5,maxrangencoll);

  // step #2
  TH1D *hSimBBC = new TH1D("hSimBBC","hSimBBC",nhistbins,-0.5,maxrange);

  TF1 *flatline = new TF1("flatline","[0]",lowfit,highfit);
  flatline->SetParameter(0,1.0);

  //-----------------------------------------------------------
  double lowmu = 3.0;
  double highmu = 5.0;
  double lowk = 0.2;
  double highk = 2.0;
  if (isSim)
    {
      lowmu = lowmu*130;
      highmu = highmu*130;
    }
  int nmusteps = 20;
  int nksteps = 20;

  double binsizemu = (highmu-lowmu)/((double)nmusteps);
  double binsizek = (highk-lowk)/((double)nksteps);

  double bestmu = mu; double bestk = k;
  TH2D *hmuk_chi2 = new TH2D("hmuk_chi2","hmuk_chi2",nmusteps,lowmu-0.5*binsizemu,highmu+0.5*binsizemu,nksteps,lowk-0.5*binsizek,highk+0.5*binsizek);
  TH2D *hmuk_chi2_fine;
  if (flag_determineNBDparameters) {

    for (int imu=0;imu<nmusteps;imu++) {
      for (int ik=0;ik<nksteps;ik++) {

	double mutemp = hmuk_chi2->GetXaxis()->GetBinCenter(imu+1);
	double ktemp = hmuk_chi2->GetYaxis()->GetBinCenter(ik+1);

	hSimBBC->Reset();
	// loop over nNcoll values, 2 for nPart
	for (int ib=2; ib<=hNcoll->GetNbinsX(); ib++) {

	  int ncoll = (int) hNcoll->GetBinCenter(ib);
	  double event_weightfactor = hNcoll->GetBinContent(ib);
	  if (event_weightfactor <= 0) continue;

	  for (int ihit=0; ihit< (20 * ncoll * (int) mutemp); ihit++) {

	    // TRY THIS FOR AUAU 
	    if (ncoll > 10 && ihit > (2 * ncoll * (int) mu)) continue;

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
	hRatio->Fit(flatline,"NDORQ","",(double)lowfit,highfit);

	if (!silence) cout << "Parameters mu = " << mutemp << " k = " << ktemp << " chi^2 = " << flatline->GetChisquare() << " and prob = " << TMath::Prob(flatline->GetChisquare(),100) << endl;	

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

    double lowmufine = bestmu - binsizemu;
    double highmufine = bestmu + binsizemu;
    double lowkfine = bestk - binsizek;
    double highkfine = bestk + binsizek;

    double binsizemufine = (highmufine-lowmufine)/((double)nmusteps);
    double binsizekfine = (highkfine-lowkfine)/((double)nksteps);

    hmuk_chi2_fine = new TH2D("hmuk_chi2","hmuk_chi2",nmusteps,lowmufine-0.5*binsizemufine,highmufine+0.5*binsizemufine,nksteps,lowkfine-0.5*binsizekfine,highkfine+0.5*binsizekfine);
    for (int imu=0;imu<nmusteps;imu++) {
      for (int ik=0;ik<nksteps;ik++) {

	double mutemp = hmuk_chi2_fine->GetXaxis()->GetBinCenter(imu+1);
	double ktemp = hmuk_chi2_fine->GetYaxis()->GetBinCenter(ik+1);

	hSimBBC->Reset();
	// loop over nNcoll values, 2 for nPart
	for (int ib=2; ib<=hNcoll->GetNbinsX(); ib++) {

	  int ncoll = (int) hNcoll->GetBinCenter(ib);
	  double event_weightfactor = hNcoll->GetBinContent(ib);
	  if (event_weightfactor <= 0) continue;

	  for (int ihit=0; ihit< (20 * ncoll * (int) mutemp); ihit++) {

	    // TRY THIS FOR AUAU 
	    if (ncoll > 10 && ihit > (2 * ncoll * (int) mu)) continue;

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
	hRatio->Fit(flatline,"NDORQ","",(double)lowfit,highfit);

	if (!silence) cout << "Parameters mu = " << mutemp << " k = " << ktemp << " chi^2 = " << flatline->GetChisquare() << " and prob = " << TMath::Prob(flatline->GetChisquare(),100) << endl;	

	hRatio->Delete();

	hmuk_chi2_fine->SetBinContent(imu+1,ik+1,flatline->GetChisquare());

      } // end loop over k steps
    } // end loop over mu steps

    // calculate the best value and +/- 1 sigma ranges
    double lowestchi2_fine = 99999999.;
    for (int imu=0;imu<nmusteps;imu++) {
      for (int ik=0;ik<nksteps;ik++) {
	if (hmuk_chi2_fine->GetBinContent(imu+1,ik+1) < lowestchi2_fine && hmuk_chi2_fine->GetBinContent(imu+1,ik+1) != 0.0) {
	  lowestchi2_fine = hmuk_chi2_fine->GetBinContent(imu+1,ik+1);
	  bestmu = hmuk_chi2_fine->GetXaxis()->GetBinCenter(imu+1);
	  bestk  = hmuk_chi2_fine->GetYaxis()->GetBinCenter(ik+1);
	}
      }
    }
    double bestmu_minusonesigma_fine = 99999.; double bestmu_plusonesigma_fine = -99999.;
    double bestk_minusonesigma_fine = 99999.; double bestk_plusonesigma_fine = -99999.;
    for (int imu=0;imu<nmusteps;imu++) {
      for (int ik=0;ik<nksteps;ik++) {
	if (hmuk_chi2_fine->GetBinContent(imu+1,ik+1) < lowestchi2_fine+1.0 && hmuk_chi2_fine->GetBinContent(imu+1,ik+1) != 0.0) {

	  if (hmuk_chi2_fine->GetXaxis()->GetBinCenter(imu+1) < bestmu_minusonesigma_fine) 
	    bestmu_minusonesigma_fine = hmuk_chi2_fine->GetXaxis()->GetBinCenter(imu+1);
	  if (hmuk_chi2_fine->GetXaxis()->GetBinCenter(imu+1) > bestmu_plusonesigma_fine) 
	    bestmu_plusonesigma_fine = hmuk_chi2_fine->GetXaxis()->GetBinCenter(imu+1);

	  if (hmuk_chi2_fine->GetYaxis()->GetBinCenter(ik+1) < bestk_minusonesigma_fine) 
	    bestk_minusonesigma_fine = hmuk_chi2_fine->GetYaxis()->GetBinCenter(ik+1);
	  if (hmuk_chi2_fine->GetYaxis()->GetBinCenter(ik+1) > bestk_plusonesigma_fine) 
	    bestk_plusonesigma_fine = hmuk_chi2_fine->GetYaxis()->GetBinCenter(ik+1);

	}
      }
    }

    // end loop over parameters
    if (!silence)
      {
	cout << "=====================================================" << endl;
	cout << "Best chi2 = " << lowestchi2_fine << " with mu = " << 
	  bestmu << " +/- " << bestmu-bestmu_minusonesigma_fine << " / " << bestmu_plusonesigma_fine-bestmu << 
	  " and k = " <<
	  bestk << " +/- " << bestk-bestk_minusonesigma_fine << " / " << bestk_plusonesigma_fine-bestk << endl;
	cout << "=====================================================" << endl;
      }

    // use the best values in the rest of the calculation
    mu = bestmu;
    k  = bestk;
    biased_mu = mu * biasfactor;
    biased_k  =  k * biasfactor;
    if (!use_balanced)
      {
	qa_info.glauber_mu[(use_shifted? 1 : 0)] = mu;  
	qa_info.glauber_k[(use_shifted? 1 : 0)] = k;  
	qa_info.glauber_chi2[(use_shifted? 1 : 0)] = lowestchi2_fine;  
      } // end section for determining NBD values


  }


  hSimBBC->Reset();
  //-------------------------------------------------------------------------------------------
  // now use the best values (or default values) and re-calculate things, including bias factors...
  // loop over nNcoll values
  for (int ib=2; ib<=hNcoll->GetNbinsX(); ib++) {

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
      hSimBBCNcoll_fine->Fill((double) ihit, (double) ncoll, nbdvalue * event_weightfactor);
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


  hSimBBC->Scale(hRealBBC->Integral(lowfit,nhistbins)/hSimBBC->Integral(lowfit,nhistbins));

  TH1D *hRatio = new TH1D("hRatio","hRatio",nhistbins,-0.5, maxrange);
  for (int i=1;i<=nhistbins;i++) {
    if (hSimBBC->GetBinContent(i) > 0 && hRealBBC->GetBinContent(i) > 0) {
      hRatio->SetBinContent(i, hRealBBC->GetBinContent(i)/hSimBBC->GetBinContent(i));
      hRatio->SetBinError(i, (sqrt(hRealBBC->GetBinContent(i))/hSimBBC->GetBinContent(i)));
    }
  }

  TF1 *trigeffcurve = new TF1("trigeffcurve","1.0-TMath::Exp(-pow((x/[0]), [1]))",0.0,maxrange);
  trigeffcurve->SetParameters(0.732,0.532); // just initial value guesses
  trigeffcurve->SetParLimits(0,0.2,20.0);    
  trigeffcurve->SetParLimits(1,0.2,5.0);
  hRatio->Fit(trigeffcurve,"NDORQ","",1.0,80.0);
  // also fit a flat line above say 20 -> 140 ==> but just constraint it exactly at 1.0
  flatline->SetParameter(0,1.0);
  hRatio->Fit(flatline,"NDORQ","",(double)lowfit,highfit);
  
  // trigger efficiency from integral comparison
  double err_real; double err_sim;
  double trigeffintegral = hRealBBC->IntegralAndError(1, nhistbins, err_real)/hSimBBC->IntegralAndError(1, nhistbins, err_sim);
  double trigeff_err = trigeff*sqrt(TMath::Power(err_real/hRealBBC->Integral(1, nhistbins), 2) + TMath::Power(err_sim/hSimBBC->Integral(1, nhistbins), 2));
  if (!silence)  cout << "Trigger Efficiency from Integrals = " << trigeffintegral << endl;
  if (flag_determineNBDparameters) {
    qa_info.glauber_trig_eff[(use_shifted? 1: 0)] = trigeffintegral;  
    qa_info.glauber_trig_eff_err[(use_shifted? 1: 0)] = trigeff_err;  
  }
  else {
    qa_info.glauber_trig_eff_forced[(use_shifted? 1: 0)] = trigeffintegral;  
    qa_info.glauber_trig_eff_forced_err[(use_shifted? 1: 0)] = trigeff_err;  
    qa_info.glauber_chi2_forced[(use_shifted? 1: 0)] = flatline->GetChisquare();
  }

  TH1D *hSimBBCwTrig = new TH1D("hSimBBCwTrig","hSimBBCwTrig",nhistbins,-0.5,maxrange);
  for (int i=1;i<=nhistbins;i++) {
    if (i==1) {
      hSimBBCwTrig->SetBinContent(i,0.0); // no chance to fire the trigger if no BBC south hits
    } else {
      hSimBBCwTrig->SetBinContent(i,hSimBBC->GetBinContent(i)*trigeffcurve->Eval(hSimBBC->GetBinCenter(i)));
    }
  }

  for (int i=1;i<=nhistbins;i++) { // loop over BBC hits
    for (int j=1;j<=ncollmax;j++) { // loop over Ncoll values
      if (i==1) {
	hSimBBCNcoll_wtrig->SetBinContent(i,j,0.0);
      } else {
	hSimBBCNcoll_wtrig->SetBinContent(i,j,
		 hSimBBCNcoll->GetBinContent(i,j)*trigeffcurve->Eval(hSimBBC->GetBinCenter(i)));
      }
    }
  }


  TString calib_file_name = Form("%s/calib/mbdana_centrality_%d.root", env_p, runnumber);

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

  fcalib->Write();
  fcalib->Close();

  // TODO: Now the same case for HARD collisions

  // NPart for each centrality bin


  double npartstore[100] = {0};

  TH1D *h_npart_cent[100];

  for (int icent = 0; icent < 91; icent++)
    {
      h_npart_cent[icent] = new TH1D(Form("h_npart_cent_%d",icent),"", ncollmax, -0.5, maxrangencoll);

      for (int ibbc = 1+ floor(centrality_low[icent]); ibbc <= 1 + (floor(centrality_high[icent])) ; ibbc++)
	{
	  for (int inpart = 1; inpart <= ncollmax;inpart++) 
	    {
	      h_npart_cent[icent]->Fill((double)inpart - 1.0, hSimBBCNcoll_wtrig->GetBinContent(ibbc, inpart));
	    }
	}
	  
      npartstore[icent] = h_npart_cent[icent]->GetMean();
      if (!silence) std::cout << " Centbin "<<icent<<" <NPart> = " << npartstore[icent] << std::endl;
    }

  std::cout << " filling npart total"<<std::endl;
  TH1D *h_npart_total = new TH1D("h_npart_total","", ncollmax, -0.5, maxrangencoll);
  for (int ibbc = 1; ibbc <= nhistbins; ibbc++)
    {
      for (int inpart = 1 ; inpart <= ncollmax;inpart++) h_npart_total->Fill((double)inpart - 1.0, hSimBBCNcoll->GetBinContent(ibbc, inpart));
    }

  double npartall = h_npart_total->GetMean();

  if (!silence) 
    {
      std::cout<<"Centrality 0-100% has <Npart> = " << npartall << std::endl;
    }
  if (flag_determineNBDparameters) qa_info.glauber_npart[(use_shifted? 1: 0)] = npartall;  
  else {
    qa_info.glauber_npart_forced[(use_shifted? 1: 0)] = npartall;
  }


  TString calib_file_name2 = Form("%s/calib/mbdana_npart_%d.root", env_p, runnumber);
  std::cout << " filling npart file "<<calib_file_name2<<std::endl;
  TFile *fcalib2 = new TFile(calib_file_name2, "recreate");
  TNtuple *ts2 = new TNtuple("tn_npart", "holds npart divisions", "bin:cent_low:cent_high:npart");
  for (int i = 0; i < 91 ; i++)
    {
      ts2->Fill(i+1, centrality_low[i], centrality_high[i], npartstore[i]);
    }

  fcalib2->Write();
  fcalib2->Close();

  std::cout << " filling file "<<std::endl;
  TString extra = Form("%s%s%s_", (use_shifted ? "sca" :""), (use_balanced ? "bal":""), (flag_determineNBDparameters ? "":"forc"));
  if (!use_shifted && !use_balanced && flag_determineNBDparameters) extra = "";
  TString fname = Form("%s/output/plots/mbdana_centrality_trigeff_%s%d.root", env_p, extra.Data(), runnumber);;
  TFile *fout = new TFile(fname,"recreate");
  hmuk_chi2->Write();
  hmuk_chi2_fine->Write();
  hSimBBCwTrig->Write();
  hSimBBC->Write();
  hRealBBC->Write();
  hRatio->Write();
  trigeffcurve->Write();
  h_npart_total->Write();

  for (int icent = 0; icent < 91;icent++) h_npart_cent[icent]->Write();

  fout->Close();


} // end routine

// with parameters (mu, k) and for given n.
double QA_centrality::NBD_getValue(int n, double mu, double k) {

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
