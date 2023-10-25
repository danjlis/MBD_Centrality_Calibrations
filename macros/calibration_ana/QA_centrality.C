/* QA_centrality */
// Author: Daniel Lis - August 15 - 2023
// Brief : This macro class gives a quick QA analysis
//         of a run in sPHENIX
//      To be run on the output of the CentralityReco Module

const int DEBUG =  0;
const int SILENCE = 1;

//void QA_FindCentralities()
#include <string>
#include <string.h>
#include <iostream>
#include <filesystem>
using namespace std;
namespace fs = std::filesystem;

#include "/sphenix/user/samfred/commissioning/macros/calib/mbd_info.h"
#include "/sphenix/user/samfred/commissioning/macros/style.h"

// const int mbd_ring_index[64] =
//   {2, 2, 2, 1, 1, 2, 1, 0,
//    0, 2, 1, 0, 2, 1, 0, 1,
//    0, 2, 1, 0, 2, 1, 0, 2,
//    1, 0, 0, 2, 1, 1, 2, 2,
//    2, 2, 2, 1, 1, 2, 1, 0,
//    0, 2, 1, 0, 2, 1, 0, 1,
//    0, 2, 1, 0, 2, 1, 0, 2,
//    1, 0, 0, 2, 1, 1, 2, 2};
const int centrality_map[20] = {1999, 1499, 1291, 1102, 937, 790, 660, 547, 449, 363, 289, 227, 174, 130, 94, 66, 45, 0, 0, 0};


int GetCentBin(float sum);
double NBD_getValue(int n, double mu, double k);
std::pair<double, double> getFitFunction(TH1D* h, double *params, float minimumx = -25., float maximumx = 25., bool ForceSingle = false);
std::pair<int, string> tokenize_path(string path);

void QA_MBDCalibrations(const int runnumber, const int doPerRunCalibration, const int makeNewHistograms);
void QA_MakeChargeSum(const int runnumber, const int loadRunCalibration = 1, const int shift_runnumber = 0);
void QA_MBDTimeCalibrations(const int runnumber, const int useChargeTime = 0);
void QA_MakeCentralityCalibrations(const int runnumber, const bool doVertexSelection = 1, const bool use_shifted = false);
void QA_MakeCentralityComparison();
void QA_CentralityCheck(const int runnumber, const int centrality_runnumber = 21813);

void QA_centrality(
		   const int runnumber,
		   const int loadRunCalibration = 1,
		   const int doPerRunCalibration = 1,
		   const int doFindCentralities = 1,
		   const int doVertexSelection = 1
		   )
{
  //QA_MBDCalibrations(runnumber, doPerRunCalibration, 1);  
  //  QA_MBDTimeCalibrations(runnumber);  
  //QA_MakeChargeSum(runnumber, loadRunCalibration, 21813);
  QA_MakeChargeSum(runnumber, 1, 21813);
  QA_CentralityCheck(runnumber, 21813);
  return;
}

void QA_MakeCentralityComparison()
{


  const char *env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }

  char *calib_path = new char[100];

  sprintf(calib_path,"%s/calib/", env_p);


  std::vector<std::pair<int, string>> v_run_path;
  for (const auto & entry : fs::directory_iterator(calib_path))
    {
      std::pair<int, string> run_path =tokenize_path(entry.path());
      if (run_path.first)
	{
	  std::cout << run_path.first <<" : "<<run_path.second<<std::endl;
	  v_run_path.push_back(run_path);
	}
    }


  int number_of_runs = v_run_path.size();

  std::cout << " Making graphs"<<endl;

  TH1D *h_cent[19];
  TGraph *g_cent[19];

  TMultiGraph *bigone = new TMultiGraph();
  bigone->SetName("bigone");

  int nbins = 3000;
  for (int i = 0; i < 19; i++)
    {
      h_cent[i] = new TH1D(Form("h_cent_%d", i), "", nbins, -0.5, nbins - 0.5);
      g_cent[i] = new TGraph();
      g_cent[i]->Set(number_of_runs);
      //      g_cent[i]->SetName(Form("g_cent_%d", i));
      //bigone->Add(g_cent[i], Form("g_cent_%d", i));
    }

  std::cout << "Starting loop through runs"<<std::endl;
  for (int irun = 0; irun < number_of_runs; irun++)
    {
      int runnumber = v_run_path.at(irun).first;
      std::string centrality_path = v_run_path.at(irun).second;

      TFile *fin = new TFile(centrality_path.c_str());
      if (!fin)
	{
	  std::cout << centrality_path << " not found!"<<std::endl;
	  continue;
	}
      TNtuple *tin = (TNtuple*) fin->Get("tn_centrality");

      float centbin;
      float lowbin, highbin;
      tin->SetBranchAddress("bin", &centbin);
      tin->SetBranchAddress("low", &lowbin);
      tin->SetBranchAddress("high", &highbin);


      int entries = tin->GetEntries();
      std::cout << "Run " << irun << std::endl;
      for (int icent = 0; icent < entries; icent++)
	{
	  tin->GetEntry(icent);
	  if (icent >= 19) break;
	  g_cent[icent]->SetPoint(irun, irun, lowbin); 
	  h_cent[icent]->Fill(lowbin); 
	}
    }

  for (int i = 0; i < 19; i++)
    {
      bigone->Add(g_cent[i], Form("g_cent_%d", i));
      // Histogram

      double mean = h_cent[i]->GetMean();
      int bottom = -1;
      int top = -1;
      for(int j = 1; j <= nbins; j++)
	{
	  if (h_cent[i]->GetBinContent(j) && bottom == -1)
	    {
	      bottom = static_cast<int>(h_cent[i]->GetBinCenter(j));
	    }
	  if (h_cent[i]->GetBinContent(nbins + 1 - j) && top == -1)
	    {
	      top = static_cast<int>(h_cent[i]->GetBinCenter(nbins + 1 - j));
	    }

	}
      if (top <= 0 & bottom <= 0) continue;

      TH1D *h = (TH1D*) h_cent[i]->Clone();
      h->SetName("h");
      h_cent[i] = new TH1D(Form("h_cent_%d", i), "", top - bottom + 10 , bottom - 10, top + 10);

      for(int j = 1; j <= nbins; j++)
	{
	  h_cent[i]->Fill(h->GetBinCenter(j), h->GetBinContent(j));
	}
      delete h;
    }

  TFile *fout = new TFile(Form("%s/output/plots/centcomparison.root", env_p), "recreate");
  
  bigone->Write();
  for (int i = 0; i < 19; i++)
    {
      g_cent[i]->SetName(Form("g_cent_%d", i));
      g_cent[i]->Write();
      h_cent[i]->Write();
    }

  fout->Close();

  return;
}

void QA_CentralityCheck(const int runnumber, const int centrality_runnumber = 21813)
{

  const char *env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }

  // Load in charge_sum distributions

  TFile *fchargesum = new TFile(Form("%s/output/plots/mbd_charge_sum_%d.root", env_p, runnumber), "r");
  if (!fchargesum)
    {
      cout << "No file exists for this charge sum plot... exiting" << endl;
      return;
    } 
  TH1D *h_charge_sum_shifted = (TH1D*) fchargesum->Get("h_charge_sum_min_bias_w_vertex_30_shifted");
  TH1D *h_charge_sum = (TH1D*) fchargesum->Get("h_charge_sum_min_bias_w_vertex_30");
  if (!h_charge_sum)
    {
      cout << "No histogram exists for this charge sum plot... exiting" << endl;
      return;
    } 

  // Load in centrality calibrations

  TFile *fcentralitycalib = new TFile(Form("%s/calib/calib_centrality_%d.root", env_p, centrality_runnumber), "r");
  if (!fcentralitycalib)
    {
      cout << "No calib file exists for the centrality in run " << centrality_runnumber << ".. exiting"<<endl;
      return;
    }
  
  float centbin, low, high;
  TNtuple *tn_centrality = (TNtuple*) fcentralitycalib->Get("tn_centrality");
  if (!tn_centrality)
    {
      cout << "No calib tntuple exists for the centrality in run " << centrality_runnumber << ".. exiting"<<endl;
      return;
    }
  
  tn_centrality->SetBranchAddress("bin", &centbin);
  tn_centrality->SetBranchAddress("low", &low);
  tn_centrality->SetBranchAddress("high", &high);
  
  float centrality_bins[20];
  for (int i = 0 ; i < 20; i++)
    centrality_bins[i] = 0.;

  // filling centrality bins

  int cent_divs = tn_centrality->GetEntries();
  for (int i = 0; i < cent_divs; i++)
    {
      tn_centrality->GetEntry(i);

      centrality_bins[i] = low;
    }

  // now we go through the charge sum and get a percent of the run that is in each of the bins;

  TH1D *h_cent_bin = new TH1D("hcent_bins","", 20, -0.5, 19.5);
  TH1D *h_cent_bin_shifted = new TH1D("hcent_bins_shifted","", 20, -0.5, 19.5);

  int nbins = h_charge_sum->GetNbinsX();
  int i_cent_bin = 0;
  for (int i = nbins; i > 0; i--)
    {
      while (h_charge_sum->GetBinCenter(i) < centrality_bins[i_cent_bin]) i_cent_bin++;
      h_cent_bin->Fill(i_cent_bin, h_charge_sum->GetBinContent(i));
    }

  nbins = h_charge_sum_shifted->GetNbinsX();
  i_cent_bin = 0;
  for (int i = nbins; i > 0; i--)
    {
      while (h_charge_sum_shifted->GetBinCenter(i) < centrality_bins[i_cent_bin]) i_cent_bin++;
      h_cent_bin_shifted->Fill(i_cent_bin, h_charge_sum_shifted->GetBinContent(i));
    }
  for (int i = 1; i <= 20; i++)
    {
      h_cent_bin->SetBinError(i, sqrt(h_cent_bin->GetBinContent(i)));
      h_cent_bin_shifted->SetBinError(i, sqrt(h_cent_bin_shifted->GetBinContent(i)));
    }

  TF1 *flatline = new TF1("flatline","[0]",-0.5, 19.5);
  flatline->SetParameter(0,0.05);

  h_cent_bin->Scale(1./h_cent_bin->Integral());
  h_cent_bin->Fit(flatline,"NDORQ","");

  float chi2 = flatline->GetChisquare()/flatline->GetNDF();

  flatline->SetParameter(0,0.05);

  h_cent_bin_shifted->Scale(1./h_cent_bin_shifted->Integral());
  h_cent_bin_shifted->Fit(flatline,"QNDOR","");

  float chi2_shifted = flatline->GetChisquare()/flatline->GetNDF();


  cout << " ************************** " <<endl;
  cout << "  Run "<< runnumber<<endl;
  cout << "    Cent Chi2       =  "<< chi2<<endl;
  cout << "    Cent Chi2 shift =  "<< chi2_shifted<<endl;

  cout << " ************************** " <<endl;

  TFile *fout = new TFile(Form("%s/output/plots/centrality_check_%d.root", env_p, runnumber), "recreate");
  h_cent_bin->Write();
  h_cent_bin_shifted->Write();
  fout->Close();

  return;
}


void QA_MakeChargeSum(const int runnumber, const int loadRunCalibration = 1, const int shift_runnumber = 0)
{
  int before = 2;
  int start_looking = 10;
  int min_range = 8;
  TFile *calibfile;
  float gain_corr[128];
  TFile *tcalibfile;
  float time_shift_corr[128];
  float toa_shift_corr[128];
  float time_scale_corr[128];
  float g, ch, sc, sh, toa_sh;

  double tthresh = 16;
  double cthresh = 0.4;
  int central_cut = 4;
  float sigma_cut = 1.5;


  const char *env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }

  if (loadRunCalibration)
    {
      calibfile = new TFile(Form("%s/calib/calib_mbd_%d.root", env_p, runnumber), "r");
        
      if (!calibfile)
	{
	  std::cout << " No calibration File found " <<std::endl;
	  return;
	}
  
      // 64 gain corrections on each side
      //  0 -  63 -- North
      // 64 - 127 -- South
      
    
      TNtuple *s = (TNtuple*) calibfile->Get("mbd_calib");
      s->SetBranchAddress("peak", &g);
      
      for (int i = 0 ; i < 128; i ++)
	{
	  s->GetEntry(i);
	  gain_corr[i] = g;
	}
      tcalibfile = new TFile(Form("%s/calib/t0_calib_mbd_%d.root", env_p, runnumber), "r");

      if (!tcalibfile)
	{
	  std::cout << " No calibration File found " <<std::endl;
	  return;
	}
  
      // 64 gain corrections on each side
      //  0 -  63 -- North
      // 64 - 127 -- South
      
      TNtuple *ts = (TNtuple*) tcalibfile->Get("ttree");
      ts->SetBranchAddress("shift", &sh);
      ts->SetBranchAddress("scale", &sc);
      ts->SetBranchAddress("toa_shift", &toa_sh);
      ts->SetBranchAddress("channel", &ch);
      
      for (int i = 0 ; i < 128; i ++)
 	{
	  time_shift_corr[i] = -999;
	  toa_shift_corr[i] = -999;
	  time_scale_corr[i] = 1;
	}
      for (int i = 0 ; i < ts->GetEntries(); i ++)
 	{
	  ts->GetEntry(i);
	  int ich = static_cast<int>(ch);
	  
	  time_shift_corr[ich] = sh;
	  toa_shift_corr[ich] = toa_sh;
	  time_scale_corr[ich] = sc;
	}

    }

  // Get the mbd charge sum plot to shift to.
  float mbd_charge_factor = 1.;
  
  if (shift_runnumber)
    {
      TFile *fout = new TFile(Form("%s/output/plots/mbd_charge_sum_%d.root", env_p, shift_runnumber), "r");
      if (!fout) return;
      TH1D *h_shift = (TH1D*) fout->Get("h_charge_sum_min_bias_w_vertex_30");
      if (!h_shift) return;
      h_shift->SetName("h_charge_sum_shift");
      mbd_charge_factor = h_shift->GetMean();
    }

  //ku66tt66ttttyy7y
  float low_t = -20.;
  float high_t = 20.;

  float binlengtht = 0.01;

  int nbins_t = floor((high_t - low_t)/binlengtht);

  int nbin = 2500;
  int maxrange = 2500;
  bool minbias = false;

  TH1D *h_vertex = new TH1D("h_vertex", "", nbin, -250, 250);

  TH1D *h_vertex_c[20];
  for (int i = 0; i < 20; i++)  {
    h_vertex_c[i]= new TH1D(Form("h_vertex_c%d", i), "", nbin, -250, 250);
  }  
  TH1D *h_vertex_w_center = new TH1D("h_vertex_w_center", "", nbin, -250, 250);
  TH1D *h_vertex_w_mean = new TH1D("h_vertex_w_mean", "", nbin, -250, 250);
  TH1D *h_vertex_w_t = new TH1D("h_vertex_w_t", "", nbin, -250, 250);
  TH1D *h_vertex_wo_t = new TH1D("h_vertex_wo_t", "", nbin, -250, 250);
  TH1D *h_time_0 = new TH1D("h_time_0", "", nbin, -25, 25);
  TH2D *h_vtx_time_0 = new TH2D("h_vtx_time_0", "",nbin/10, -250, 250, nbin/10, -25, 25);

  TH2D *h_vtx_diff_nhit = new TH2D("h_vtx_diff_nhit", "",nbin/10, -250, 250, 129, -64.5, 64.5);

  TH2D *h_vtx_time_0_center = new TH2D("h_vtx_time_0_center", "",nbin/10, -250, 250, nbin/10, -25, 25);

  TH2D *h_vtx_diff_nhit_center = new TH2D("h_vtx_diff_nhit_center", "",nbin/10, -250, 250, 129, -64.5, 64.5);
  TH2D *h_vtx_diff_nhit_center_t = new TH2D("h_vtx_diff_nhit_center_t", "",nbin/10, -250, 250, 129, -64.5, 64.5);


  TH2D *h_rms_by_time_n = new TH2D("h_rms_by_time_n","", 60, -5, 25, 100, 0, 25);
  TH2D *h_rms_by_time_s = new TH2D("h_rms_by_time_s","", 60, -5, 25, 100, 0, 25);
  TH2D *h_mean_north_south = new TH2D("h_mean_north_south", ";<t>_{N};<t>_{S}", 60, -5, 25, 60, -5, 25);
  TH1D *h_occuppancy = new TH1D("h_occuppancy","", 128, -0.5, 127.5);
  TH1D *h_voccuppancy = new TH1D("h_voccuppancy","", 128, -0.5, 127.5);

  TEfficiency *h_ring_occuppancy[3]; 
  TEfficiency *h_ring_voccuppancy[3]; 
  for (int i = 0; i < 3; i++)
    {
      h_ring_voccuppancy[i] = new TEfficiency(Form("h_voccuppancy_ring_%d", i),"", 128, -0.5, 127.5);
      h_ring_occuppancy[i]= new TEfficiency(Form("h_occuppancy_ring_%d", i),"", 128, -0.5, 127.5);
    }
  TH1D *h_charge_sum = new TH1D("h_charge_sum","", nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_min_bias = new TH1D("h_charge_sum_min_bias","",  nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_min_bias_w_vertex_30 = new TH1D("h_charge_sum_min_bias_w_vertex_30","", nbin, -0.5, (float)maxrange - 0.5);
  TH1D *h_charge_sum_min_bias_w_vertex_30_shifted = new TH1D("h_charge_sum_min_bias_w_vertex_30_shifted","", nbin, -0.5, (float)maxrange - 0.5);

  TH1D *h_time_w_t[128];
  TH1D *h_time_wo_t[128];

  for (int ich = 0 ; ich < 128; ich++)
    {
      h_time_w_t[ich]  = new TH1D(Form("h_time_w_t_%d", ich), "", nbins_t, low_t, high_t);
      h_time_wo_t[ich]  = new TH1D(Form("h_time_wo_t_%d", ich), "", nbins_t, low_t, high_t);
    }


  TH1D *h_zdc_sum_n = new TH1D("h_zdc_sum_n", ";ZDC sum [GeV]; Counts", 1000, 0, 10000);
  TH1D *h_zdc_sum_s = new TH1D("h_zdc_sum_s", ";ZDC sum [GeV]; Counts", 1000, 0, 10000);
  TH1D *h_zdc_sum_n_prime = new TH1D("h_zdc_sum_n_prime", ";ZDC sum [GeV]; Counts", 1000, 0, 10000);
  TH1D *h_zdc_sum_s_prime = new TH1D("h_zdc_sum_s_prime", ";ZDC sum [GeV]; Counts", 1000, 0, 10000);

  TH1D *h_zdc_sum_n_w_vertex_30 = new TH1D("h_zdc_sum_n_w_vertex_30", ";ZDC sum [GeV]; Counts", 1000, 0, 10000);
  TH1D *h_zdc_sum_s_w_vertex_30 = new TH1D("h_zdc_sum_s_w_vertex_30", ";ZDC sum [GeV]; Counts", 1000, 0, 10000);
  TH1D *h_zdc_sum_n_prime_w_vertex_30 = new TH1D("h_zdc_sum_n_prime_w_vertex_30", ";ZDC sum [GeV]; Counts", 1000, 0, 10000);
  TH1D *h_zdc_sum_s_prime_w_vertex_30 = new TH1D("h_zdc_sum_s_prime_w_vertex_30", ";ZDC sum [GeV]; Counts", 1000, 0, 10000);

  // TAxis *axis_n = h_zdc_sum_n->GetXaxis();
  // TAxis *axis_s = h_zdc_sum_s->GetXaxis();
  // TAxis *axis_np = h_zdc_sum_n_prime->GetXaxis();
  // TAxis *axis_sp = h_zdc_sum_s_prime->GetXaxis();

  // const int bins = axis_n->GetNbins();

  // Axis_t from = axis_n->GetXmin();
  // Axis_t to = axis_n->GetXmax();
  // Axis_t width = (to - from) / bins;
  // Axis_t *new_bins = new Axis_t[bins + 1];

  // for (int i = 0; i <= bins; i++) {
  //   new_bins[i] = TMath::Power(10, from + i * width);
  // }

  // axis_n->Set(bins, new_bins);
  // axis_s->Set(bins, new_bins);
  // axis_np->Set(bins, new_bins);
  // axis_sp->Set(bins, new_bins);
  // h_zdc_sum_n_w_vertex_30->GetXaxis()->Set(bins, new_bins);
  // h_zdc_sum_s_w_vertex_30->GetXaxis()->Set(bins, new_bins);
  // h_zdc_sum_n_prime_w_vertex_30->GetXaxis()->Set(bins, new_bins);
  // h_zdc_sum_s_prime_w_vertex_30->GetXaxis()->Set(bins, new_bins);

  // delete[] new_bins;


  TFile *file = new TFile(Form("%s/output/run%d/mbdcalibana/mbd_calib_trees_%d.root", env_p, runnumber, runnumber), "r");
  if (!file)
    {
      std::cout << " No Tree File found " <<std::endl;
      return;
    }

  // Making fresh histograms


  float zdc_sum[2];
  float zdc_sum_prime[2];
  
  float mbd_charge[128];
  float mbd_time[128];
  float mbd_charge_raw[128];
  float mbd_time_raw[128];

  float z_vertex;
  float time_0;  

  TTree *t = (TTree*)file->Get("T");

  t->SetBranchAddress("zdc_sum",zdc_sum);
  t->SetBranchAddress("zdc_sum_prime",zdc_sum_prime);
  t->SetBranchAddress("mbd_charge_raw",mbd_charge_raw);
  t->SetBranchAddress("mbd_time_raw",mbd_time_raw);

  double charge_sum = 0;

  //// vertex
  int hits_n = 0;
  int hits_s = 0;      
  int hits_n_t = 0;
  int hits_s_t = 0;      

  std::vector<float> time_sum_n;
  std::vector<float> time_sum_s;
  float sum_n = 0.;
  float sum_s = 0.;
  float sum_n2 = 0.;
  float sum_s2 = 0.;

  cout <<"NEvents = "<<t->GetEntries()<<std::endl;

  for (int i = 0 ; i < (DEBUG ? 10 : t->GetEntries()); i++)
    {
      minbias = false;
      hits_n = 0;
      hits_s = 0;      
      hits_n_t = 0;
      hits_s_t = 0;      
      charge_sum = 0;
      time_sum_n.clear();
      time_sum_s.clear();
      sum_n = 0.;
      sum_s = 0.;
      sum_n2 = 0.;
      sum_s2 = 0.;

      t->GetEntry(i);
      if (i%20000 == 0) std::cout <<" ------------------------------- Event: "<<i<< " -------------------------------"<<std::endl;

      for (int ich = 0 ; ich < 128; ich++)
	{

	  mbd_time[ich] = ((25. - mbd_time_raw[ich] * (9.0 / 5000.) - time_shift_corr[ich]))*time_scale_corr[ich];
	  mbd_charge[ich] = mbd_charge_raw[ich]*gain_corr[ich];

	  if (mbd_charge[ich] > cthresh)
	    {
	      charge_sum += mbd_charge[ich];	      	      

	      // get hit charge channels
	      if (ich/64) hits_n++;
	      else hits_s++;

	      // if bad channel or the time is greater than 10 ns after 0-point don't count the time
	      if (ich == 56 || ich == 56+64 || fabs(mbd_time[ich]) > 10) continue;

	      float timme = mbd_time[ich];

	      if (ich/64) 
		{
		 		  
		  hits_n_t++;
		  time_sum_n.push_back(timme);
		  sum_n += timme;;
		  sum_n2 += (timme*timme);

		}
	      else
		{
		  
		  hits_s_t++;
		  time_sum_s.push_back(timme);
		  
		  sum_s += timme;
		  sum_s2 += (timme*timme);

		}

	    }
	  
	} 
      
      sort(time_sum_n.begin(), time_sum_n.end());
      sort(time_sum_s.begin(), time_sum_s.end());

      minbias = true;
      float mean_north = 999;
      float mean_south = 999;

      // does event meet cut?
      if (hits_s_t >= central_cut){

	// get the mean
	mean_south = sum_s/static_cast<float>(hits_s_t);

	// get rms
	float rms_s = sqrt(sum_s2/static_cast<float>(hits_s_t) - TMath::Power(mean_south, 2));
	int nhit_s_center = 0;
	float sum_s_center = 0.;
	
	// get rid of times outside of RMS*1.5 range
	for (unsigned int is = 0; is < time_sum_s.size(); is++)
	  {
	    if (fabs(time_sum_s.at(is) - mean_south) < sigma_cut*rms_s )
	      {
		sum_s_center += time_sum_s.at(is);
		nhit_s_center++;
	      }
	  }

	// Get the mean again
	float mean_south_center = sum_s_center/static_cast<float>(nhit_s_center);

	
	mean_south = mean_south_center;
	
      }

      else if (hits_s >=2 && (hits_s_t >= 1)){

	mean_south = sum_s/static_cast<float>(hits_s_t);

      }
      else minbias = false;

      if (hits_n_t >=central_cut){
	
	mean_north = sum_n/static_cast<float>(hits_n_t);
	
	float rms_n = sqrt(sum_n2/static_cast<float>(hits_n_t) - TMath::Power(mean_north, 2));
	int nhit_n_center = 0;
	float sum_n_center = 0.;
	
	for (unsigned int ino = 0; ino < time_sum_n.size(); ino++)
	  {
	    if (fabs(time_sum_n.at(ino) - mean_north) < sigma_cut*rms_n )
	      {
		sum_n_center += time_sum_n.at(ino);
		nhit_n_center++;
	      }
	  }
	
	float mean_north_center = sum_n_center/static_cast<float>(nhit_n_center);
       
	mean_north = mean_north_center;
      }

      else if (hits_n >=2 && hits_n_t >= 1){

	mean_north = sum_n/static_cast<float>(hits_n_t);

	h_vertex_w_mean->Fill(z_vertex);
      }

      else minbias = false;

      if (mean_north != 999 && mean_south != -999) {

	z_vertex = 15*(mean_north - mean_south);
	time_0 = (mean_north + mean_south)/2.;
      }
      else 
	{
	  z_vertex = 999;
	  time_0 = 999;
	}
      if (hits_s_t >= central_cut && hits_n_t >= central_cut)
	{
	  h_vtx_time_0_center->Fill(z_vertex, time_0);
	  h_vtx_diff_nhit_center->Fill(z_vertex, hits_s - hits_n);
	  h_vtx_diff_nhit_center_t->Fill(z_vertex, hits_s_t - hits_n_t);
	}
      h_vtx_diff_nhit->Fill(z_vertex, hits_s - hits_n);
      h_vtx_time_0->Fill(z_vertex, time_0);
      h_vertex->Fill(z_vertex);
      h_mean_north_south->Fill(mean_north, mean_south);
      if (time_0 < -1) h_vertex_wo_t->Fill(z_vertex);
      else h_vertex_w_t->Fill(z_vertex);
      //h_vertex_w_t->Fill(z_vertex);
      
      int cent_bin =GetCentBin(charge_sum);
      
      h_vertex_c[cent_bin]->Fill(z_vertex);
      h_time_0->Fill(time_0);
      h_charge_sum->Fill(charge_sum);
      h_zdc_sum_n->Fill(zdc_sum[1]);      
      h_zdc_sum_s->Fill(zdc_sum[0]);      
      h_zdc_sum_n_prime->Fill(zdc_sum_prime[1]);      
      h_zdc_sum_s_prime->Fill(zdc_sum_prime[0]);      
      if (!minbias) continue;

      h_charge_sum_min_bias->Fill(charge_sum);
  
      if (TMath::Abs(z_vertex) > 30) continue;
  
      h_charge_sum_min_bias_w_vertex_30->Fill(charge_sum);
      h_zdc_sum_n_w_vertex_30->Fill(zdc_sum[1]);      
      h_zdc_sum_s_w_vertex_30->Fill(zdc_sum[0]);      
      h_zdc_sum_n_prime_w_vertex_30->Fill(zdc_sum_prime[1]);      
      h_zdc_sum_s_prime_w_vertex_30->Fill(zdc_sum_prime[0]);      

    }


  if (shift_runnumber)
    {
      mbd_charge_factor *= 1./(h_charge_sum_min_bias_w_vertex_30->GetMean());
      for (int i = 0 ; i < (DEBUG ? 10 : t->GetEntries()); i++)
	{
	  minbias = false;
	  hits_n = 0;
	  hits_s = 0;      
	  hits_n_t = 0;
	  hits_s_t = 0;      
	  charge_sum = 0;
	  time_sum_n.clear();
	  time_sum_s.clear();
	  sum_n = 0.;
	  sum_s = 0.;
	  sum_n2 = 0.;
	  sum_s2 = 0.;

	  t->GetEntry(i);

	  for (int ich = 0 ; ich < 128; ich++)
	    {

	      mbd_time[ich] = ((25. - mbd_time_raw[ich] * (9.0 / 5000.) - time_shift_corr[ich]))*time_scale_corr[ich];
	      mbd_charge[ich] = mbd_charge_raw[ich]*gain_corr[ich];

	      if (mbd_charge[ich] > cthresh)
		{
		  charge_sum += mbd_charge[ich];	      	      

		  // get hit charge channels
		  if (ich/64) hits_n++;
		  else hits_s++;

		  // if bad channel or the time is greater than 10 ns after 0-point don't count the time
		  if (ich == 56 || ich == 56+64 || fabs(mbd_time[ich]) > 10) continue;

		  float timme = mbd_time[ich];

		  if (ich/64) 
		    {
		 		  
		      hits_n_t++;
		      time_sum_n.push_back(timme);
		      sum_n += timme;;
		      sum_n2 += (timme*timme);

		    }
		  else
		    {
		  
		      hits_s_t++;
		      time_sum_s.push_back(timme);
		  
		      sum_s += timme;
		      sum_s2 += (timme*timme);

		    }

		}
	  
	    } 
      
	  sort(time_sum_n.begin(), time_sum_n.end());
	  sort(time_sum_s.begin(), time_sum_s.end());
	  minbias = true;
	  float mean_north = 999;
	  float mean_south = 999;

	  // does event meet cut?
	  if (hits_s_t >= central_cut){

	    // get the mean
	    mean_south = sum_s/static_cast<float>(hits_s_t);

	    // get rms
	    float rms_s = sqrt(sum_s2/static_cast<float>(hits_s_t) - TMath::Power(mean_south, 2));
	    int nhit_s_center = 0;
	    float sum_s_center = 0.;
	
	    // get rid of times outside of RMS*1.5 range
	    for (unsigned int is = 0; is < time_sum_s.size(); is++)
	      {
		if (fabs(time_sum_s.at(is) - mean_south) < sigma_cut )
		  {
		    sum_s_center += time_sum_s.at(is);
		    nhit_s_center++;
		  }
	      }

	    // Get the mean again
	    float mean_south_center = sum_s_center/static_cast<float>(nhit_s_center);

	
	    mean_south = mean_south_center;
	
	  }

	  else if (hits_s >=2 && (hits_s_t >= 1)){

	    mean_south = sum_s/static_cast<float>(hits_s_t);

	  }
	  else minbias = false;

	  if (hits_n_t >=central_cut){
	
	    mean_north = sum_n/static_cast<float>(hits_n_t);
	
	    float rms_n = sqrt(sum_n2/static_cast<float>(hits_n_t) - TMath::Power(mean_north, 2));
	    int nhit_n_center = 0;
	    float sum_n_center = 0.;
	
	    for (unsigned int ino = 0; ino < time_sum_n.size(); ino++)
	      {
		if (fabs(time_sum_n.at(ino) - mean_north) < sigma_cut )
		  {
		    sum_n_center += time_sum_n.at(ino);
		    nhit_n_center++;
		  }
	      }
	
	    float mean_north_center = sum_n_center/static_cast<float>(nhit_n_center);
       
	    mean_north = mean_north_center;
	  }

	  else if (hits_n >=2 && hits_n_t >= 1){

	    mean_north = sum_n/static_cast<float>(hits_n_t);

	    h_vertex_w_mean->Fill(z_vertex);
	  }
	  else minbias = false;

	  if (mean_north != 999 && mean_south != -999) {

	    z_vertex = 15*(mean_north - mean_south);
	    time_0 = (mean_north + mean_south)/2.;
	  }
	  else 
	    {
	      z_vertex = 999;
	      time_0 = 999;
	    }
      
	  if (!minbias) continue;
	  if (TMath::Abs(z_vertex) > 30) continue;
  
	  h_charge_sum_min_bias_w_vertex_30_shifted->Fill(charge_sum*mbd_charge_factor);
	}


    }
  TFile *fout = new TFile(Form("%s/output/plots/mbd_charge_sum_%d.root", env_p, runnumber), "RECREATE");

  h_charge_sum->Write();
  h_charge_sum_min_bias->Write();
  h_charge_sum_min_bias_w_vertex_30->Write();
  if (shift_runnumber)   h_charge_sum_min_bias_w_vertex_30_shifted->Write();

  std::cout << " ******************************* " << std::endl;
  std::cout << "  Run " << runnumber <<std::endl;
  std::cout << "    Mean Charge Sum = " << h_charge_sum_min_bias_w_vertex_30->GetMean() << std::endl;
  if (shift_runnumber) std::cout << "    Mean Shifted    = " << h_charge_sum_min_bias_w_vertex_30_shifted->GetMean() << std::endl;
  h_vertex->Fit("gaus","Q","", -30, 30);
  std::cout << "    Vertex Mean     = " << h_vertex->GetFunction("gaus")->GetParameter(1) << std::endl;
  std::cout << "    Vertex std.     = " << h_vertex->GetFunction("gaus")->GetParameter(2) << std::endl;
  
  std::cout << " ******************************* " << std::endl;

  h_vertex->Write();
  h_vertex_w_center->Write();
  h_vertex_w_mean->Write();
  h_vertex_w_t->Write();
  h_vertex_wo_t->Write();
  h_vtx_diff_nhit->Write();
  h_vtx_time_0->Write();
  h_vtx_time_0_center->Write();
  h_vtx_diff_nhit_center_t->Write();
  h_vtx_diff_nhit_center->Write();
  h_zdc_sum_n->Write();
  h_zdc_sum_s->Write();
  h_zdc_sum_n_prime->Write();
  h_zdc_sum_s_prime->Write();

  h_zdc_sum_n_w_vertex_30->Write();
  h_zdc_sum_s_w_vertex_30->Write();
  h_zdc_sum_n_prime_w_vertex_30->Write();
  h_zdc_sum_s_prime_w_vertex_30->Write();
	  	  
  for (int i = 0; i < 20;  i++)
    {
      h_vertex_c[i]->Write();
    }
  h_rms_by_time_n->Write();
  h_rms_by_time_s->Write();
  h_mean_north_south->Write();
  for (int i = 0; i < 128; i++)
    {
      h_time_w_t[i]->Write();
      h_time_wo_t[i]->Write();
    }
  h_time_0->Write();
  h_occuppancy->Write();
  h_voccuppancy->Write();
  for (int i = 0; i < 3; i++)
    {
      h_ring_voccuppancy[i]->Write();
      h_ring_occuppancy[i]->Write();
    }
  fout->Close();
    
}

void QA_MBDCalibrations(const int runnumber, const int doPerRunCalibration, const int makeNewHistograms)
{
  
  // Making the calibration file for the MBD

  int before = 2;
  int start_looking =20;
  int min_range = 15;
  float gain_corr[128];
  float g;
  double chisquare[128];

  const char *env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }
  
  // Output from the centrality module.

  TFile *file = new TFile(Form("%s/output/run%d/mbdcalibana/mbd_calib_trees_%d.root", env_p, runnumber, runnumber), "r");
  if (!file)
    {
      std::cout << " No Tree File found " <<std::endl;
      return;
    }

  // Making fresh histograms

  TH1D *h_charge[128];
  TH1D *h_charge_raw[128];
  TH1D *hempty = new TH1D("hempty","", 2000, -0.5, 1999.5);
  for (int j = 0; j < 128; j++)
    {
      h_charge[j] = new TH1D(Form("h_mbd_charge_ch%d", j), "",500, 0, 10);
      h_charge_raw[j] = new TH1D(Form("h_mbd_charge_raw_ch%d", j),"", 2000, -0.5, 1999.5);
    }
  if (!makeNewHistograms)
    {
      TFile *fin2 = new TFile(Form("%s/output/plots/mbd_calib_plots_%d.root", env_p, runnumber), "r");
      for (int ich = 0 ; ich < 128; ich++)
	{
	  h_charge[ich] = (TH1D*) fin2->Get(Form("h_mbd_charge_ch%d", ich));
	  h_charge_raw[ich] = (TH1D*) fin2->Get(Form("h_mbd_charge_raw_ch%d", ich));
	}
    }
  float zdc_sum[2];
  float mbd_charge_raw[128];
  float mbd_time_raw[128];
  
  TTree *t = (TTree*)file->Get("T");

  t->SetBranchAddress("zdc_sum", zdc_sum);
  t->SetBranchAddress("mbd_charge_raw",mbd_charge_raw);
  t->SetBranchAddress("mbd_time_raw",mbd_time_raw);

  if (makeNewHistograms)
    {
      for (int i = 0; i < t->GetEntries();i++)
	{
	  if (i%10000 == 0) cout << i << "zdc: "<< zdc_sum[0] <<" " << zdc_sum[1]<<endl;
	  t->GetEntry(i);
	  
	  //	  	  if (zdc_sum[0] < 40 || zdc_sum[1] < 40) continue;
	  for (int ich = 0 ; ich < 128; ich++)
	    {
	      h_charge_raw[ich]->Fill(mbd_charge_raw[ich]);
	    }
	}
    }

  TH1D *h_peaks = new TH1D("h_peaks","", 64, -0.5, 1.5);
  TH1D *h_peaks_raw = new TH1D("h_peaks_raw","", 64, 0, 500);

  TF1 *f_exp = new TF1("f_exp","[0]*TMath::Exp(-1*(x - [1])/[2])", 0, 2000);
  TF1 *f_lan_w_exp = new TF1("lan_w_exp","[0]*TMath::Landau(x,[1],[2],3) + [3]*TMath::Exp(-1*(x - [4])/[5])", 0, 2000);
  TF1 *f_lan_w_gausexp = new TF1("lan_w_gausexp","[0]*TMath::Landau(x,[1],[2],3) + [3]*TMath::Exp(-1*(x - [4])/[5])+gaus(6)", 2000);


  TFile *fcalibout;
  TNtuple *tn;

  fcalibout = new TFile(Form("%s/calib/calib_mbd_%d.root", env_p, runnumber),"RECREATE");
  tn = new TNtuple("mbd_calib","mbd_calib","channel:peak:width");


  int rebinning_levels[5] = {2, 4, 5, 8, 10};

  for (int ich = 0 ; ich < 128; ich++){
    
    std::cout << "Raw Fitting :: Channel "<<ich;

    double found_max = -999;
    double found_min = -999;
    int tries = 0;
    while (tries < 5)
      {
    
	double local_min = -1;
	double local_max = -1;
	double max_value = 0;
	double min_value = 999999;
	int good = 0;
	int low_bin = 0;
	double upper;
	int high_bin = 0;
	TH1D *hsmooth = (TH1D*) h_charge_raw[ich]->Clone();
	hsmooth->Smooth();
	hsmooth->Smooth();
	// Finding the dip
	if (found_max < 0)
	  {
	    double threshold = 1.2 - tries*0.3;;
	    for (int ic = start_looking; ic < hsmooth->GetNbinsX()-1; ic++)
	      {
		high_bin = ic;
		max_value = hsmooth->GetBinContent(ic);
		for (int icc = ic; icc < hsmooth->GetNbinsX(); icc++)
		  {
		    if (hsmooth->GetBinContent(icc) > max_value)
		      {
			high_bin = icc;
			max_value = hsmooth->GetBinContent(icc);
		      }
		  }
		if (high_bin == ic) continue;
		if (max_value > threshold*hsmooth->GetBinContent(ic)) break;
	      }
	    
	    local_max = hsmooth->GetBinCenter(high_bin);
	    
	    for (int ic = start_looking; ic <high_bin; ic++)
	      {
		if (hsmooth->GetBinContent(ic) < min_value)
		  {
		    low_bin = ic;
		    local_min = hsmooth->GetBinCenter(ic);
		    min_value = hsmooth->GetBinContent(ic);
		  }
	      }
	  }
	else
	  {
	    hempty = (TH1D*) h_charge_raw[ich]->Clone();
 	    hempty->Reset();
	    hempty->Fill(found_max);
	    high_bin = hempty->GetMaximumBin();
	    local_max = hempty->GetBinCenter(high_bin);
	    max_value = hempty->GetBinContent(high_bin);
 	    hempty->Reset();
	    hempty->Fill(found_min);
	    low_bin = hempty->GetMaximumBin();
	    local_min = hempty->GetBinCenter(low_bin);
	    min_value = hempty->GetBinContent(low_bin);
 	  }

	std::cout <<"  Found Max/Min =  "<< local_max << " / " << local_min  <<endl;
	std::cout << " "<< endl;
	// rebinning
	if (local_max - local_min < 5 ||local_max > 1500) 
	  {
	    tries++;
	    continue;
	  }
	found_max = local_max;
	found_min = local_min;
	int peak = static_cast<double>(h_charge_raw[ich]->GetBinContent(high_bin));
	int rebin_i = -1;
	std::cout << "Rbin: 1  --> "<<peak<<endl;; 
	while (peak <= 500.)
	  {
	    rebin_i++;

	    if (rebin_i == 5) break;

	    TH1D *hn = (TH1D*) h_charge_raw[ich]->Clone();
	    hn->Rebin(rebinning_levels[rebin_i]);
	    hempty = (TH1D*) hn->Clone();
	    hempty->Reset();
	    hempty->Fill(local_max, 1);
	    peak = hn->GetBinContent(hempty->GetMaximumBin());
	    high_bin = hempty->GetMaximumBin();
	    std::cout << "Rbin: "<< rebinning_levels[rebin_i]<<"  --> "<<peak<<endl;; 
	    hempty->Reset();
	    hempty->Fill(local_min, 1);
	    low_bin = hempty->GetMaximumBin();

       	    if (peak <= 1000) break;

	  }
	if (rebin_i != -1)
	  {
	    h_charge_raw[ich]->Rebin(rebinning_levels[rebin_i]);
	    hsmooth = (TH1D*) h_charge_raw[ich]->Clone();
	    hsmooth->Smooth();
	    hsmooth->Smooth();
	  }
	// if (local_max > 600){
	//   local_max = 150;
	// }
	// else if (local_max < 40)
	//   {
	// 	local_max = 41;
	//   }
	// if (local_min > 400) 
	//   {
	// 	local_min = 30;
	// 	local_max = 150;
	//   }
	// else if (local_min < 20)
	//   {
	// 	local_min = 20;		
	//   }

	// Fit tail first

	f_exp->SetParameter(0, 200);
	f_exp->SetParameter(1, 400);
	f_exp->SetParameter(2, 350);
	
	h_charge_raw[ich]->Fit("f_exp","Q","", local_max + 300, 2000);
	if (tries == 0)
	  {
	    std::cout << " " << std::endl;
	    std::cout << "  Exponential tail fit:"<<endl;
	    std::cout << "    Norm: "<<f_exp->GetParameter(0)<<"\t Shift: "<<f_exp->GetParameter(1)<<"\t Tail: "<<f_exp->GetParameter(2)<<std::endl;
	    std::cout << "    Chi2/NDF: "<<f_exp->GetChisquare()/f_exp->GetNDF()<<std::endl;
	    std::cout << " " << std::endl;
	  }
	f_lan_w_exp->SetParLimits(0, 1, 1300000);
	f_lan_w_exp->SetParLimits(1, 0.1, 2000);
	f_lan_w_exp->SetParLimits(2, 0.1, 2000);
	f_lan_w_gausexp->SetParLimits(0, 1, 1300000);
	f_lan_w_gausexp->SetParLimits(1, 0.1, 2000);
	f_lan_w_gausexp->SetParLimits(2, 0.1, 2000);

	f_lan_w_exp->SetParameter(0, 30000);
	f_lan_w_exp->SetParameter(1, local_max);
	f_lan_w_exp->SetParameter(2, 12);
	f_lan_w_exp->SetParameter(3, f_exp->GetParameter(0));
	f_lan_w_exp->SetParameter(4, f_exp->GetParameter(1));
	f_lan_w_exp->SetParameter(5, f_exp->GetParameter(2));

	//	    f_lan_w_exp->SetParLimits(4, -100000, 0);
	    
	upper = local_min + 5 * (local_max - local_min);
	if (local_max - local_min < 10) upper = 4*local_max;
	if (tries == 0)
	  {
	    std::cout << " " << std::endl;
	    std::cout << "  Smooth:" << std::endl;
	  }
	hsmooth->Fit("lan_w_exp",(tries == 0? "":"Q"),"",local_min, upper);
	if (fabs(f_lan_w_exp->GetParameter(1) - local_max) > 20){
	  f_lan_w_exp->SetParameter(1, local_max);
	  f_lan_w_gausexp->SetParLimits(1, local_max - 20, local_max+20);
	  f_lan_w_exp->SetParLimits(1, local_max - 30, local_max+30);
	  for (int tt = 0; tt < 5; tt++) hsmooth->Fit("lan_w_exp","Q","",local_min, upper);
	}
	if (tries == 0)
	  {
	    std::cout << "Chi2/NDF: " << f_lan_w_exp->GetChisquare()/f_lan_w_exp->GetNDF()<<std::endl;
	    std::cout << " " << std::endl;	
	    std::cout << "  Real:" << std::endl;
	  }
	h_charge_raw[ich]->Fit("lan_w_exp",(tries == 0?"": "Q"),"",local_min, upper);
	if (tries == 0)std::cout << "Chi2/NDF: " << f_lan_w_exp->GetChisquare()/f_lan_w_exp->GetNDF()<<std::endl;

	f_lan_w_gausexp->SetParameter(0, f_lan_w_exp->GetParameter(0));
	f_lan_w_gausexp->SetParameter(1, f_lan_w_exp->GetParameter(1));
	f_lan_w_gausexp->SetParameter(2, f_lan_w_exp->GetParameter(2));
	f_lan_w_gausexp->SetParameter(3, f_lan_w_exp->GetParameter(3));
	f_lan_w_gausexp->SetParameter(4, f_lan_w_exp->GetParameter(4));
	f_lan_w_gausexp->SetParameter(5, f_lan_w_exp->GetParameter(5));

	float value = 0;
	float value_x = 0;

	int ibin_mark = low_bin - 1;
	for (int ilower = 1; ilower < low_bin - 1; ilower++)
	  {
	    ibin_mark = low_bin - ilower;
	    value = hsmooth->GetBinContent(ibin_mark);
	    value_x = hsmooth->GetBinCenter(ibin_mark);
	    if (hsmooth->GetBinContent(ibin_mark) > 1.2*hsmooth->GetBinContent(high_bin))
	      {
		break;
	      }
	  }

	float value_x_i = value_x;
	h_charge_raw[ich]->Fit("gaus", "","",1, value_x);

	f_lan_w_gausexp->SetParameter(6, h_charge_raw[ich]->GetFunction("gaus")->GetParameter(0));
	f_lan_w_gausexp->SetParameter(7, h_charge_raw[ich]->GetFunction("gaus")->GetParameter(1));
	f_lan_w_gausexp->SetParameter(8, h_charge_raw[ich]->GetFunction("gaus")->GetParameter(2));
	    
	int fitss = 0;
	if (tries == 0)
	  {
	    std::cout << " " << std::endl;
	    std::cout << " Real:" << std::endl;
	  }

	h_charge_raw[ich]->Fit("lan_w_gausexp",(tries == 0?"": "Q"),"",value_x_i, upper);

	if (tries == 0) std::cout << "Chi2/NDF: " << f_lan_w_gausexp->GetChisquare()/f_lan_w_gausexp->GetNDF()<<std::endl;
	while (f_lan_w_gausexp->GetChisquare()/f_lan_w_gausexp->GetNDF() > 3 && fitss < 10)
	  {

	    h_charge_raw[ich]->Fit("lan_w_gausexp","Q","",value_x_i, upper);
	    if (tries == 0) std::cout << "Fit "<<fitss<<" Chi2/NDF: " << f_lan_w_gausexp->GetChisquare()/f_lan_w_gausexp->GetNDF()<<std::endl;
	    // value_x = 0;
	    // value = 0;
	    // for (int ilower = 1; ilower < low_bin - 1; ilower++)
	    //   {
	    //     ibin_mark = low_bin - ilower;
	    //     value = hsmooth->GetBinContent(ibin_mark);
	    //     value_x = hsmooth->GetBinCenter(ibin_mark);
	    //     if (hsmooth->GetBinContent(ibin_mark) > (1.4 - 0.1*fitss)*hsmooth->GetBinContent(high_bin))
	    //       {
	    // 	break;
	    //       }
	    //   }
	    // hsmooth->Fit("gaus", "","",value_x_i, value_x);
	    // f_lan_w_gausexp->SetParameter(6, hsmooth->GetFunction("gaus")->GetParameter(0));
	    // f_lan_w_gausexp->SetParameter(7, hsmooth->GetFunction("gaus")->GetParameter(1));
	    // f_lan_w_gausexp->SetParameter(8, hsmooth->GetFunction("gaus")->GetParameter(2));
		
	    // h_charge_raw[ich]->Fit("lan_w_gausexp","Q","",value_x, upper);
	    fitss++;
	    //	    std::cout << fitss <<" - chi2/ndf = " << f_lan_w_gausexp->GetChisquare()/f_lan_w_gausexp->GetNDF()<<std::endl;
	  }
	if (tries)	std::cout << "Fit "<<fitss<<" Chi2/NDF: " << f_lan_w_gausexp->GetChisquare()/f_lan_w_gausexp->GetNDF()<<std::endl;
	if (f_lan_w_gausexp->GetChisquare()/f_lan_w_gausexp->GetNDF() < 3. + tries*.5)
	  {
	    std::cout <<" " << local_min <<"/"<< local_max<<"/"<< upper<<" -- Peak at "<<f_lan_w_gausexp->GetParameter(1)<<endl;
	    gain_corr[ich] = 1./((double) f_lan_w_gausexp->GetParameter(1));

	    break;
	  }
	std::cout <<" " << local_min <<"/"<< local_max<<"/"<< upper<<" -- Peak at "<<f_lan_w_gausexp->GetParameter(1)<<endl;
	gain_corr[ich] = 1./((double) f_lan_w_gausexp->GetParameter(1));
	tries++;
	std::cout << "try "<<tries<<std::endl;

      }

    chisquare[ich] = f_lan_w_gausexp->GetChisquare()/f_lan_w_gausexp->GetNDF();
  }

  if (makeNewHistograms)
    {
      for (int i = 0; i < t->GetEntries();i++)
	{
	  if (i%10000 == 0) cout << i << endl;
	  t->GetEntry(i);
	  //	  if (zdc_sum[0]< 40 || zdc_sum[1] < 40) continue;
	  for (int ich = 0 ; ich < 128; ich++)
	    {
	      h_charge[ich]->Fill(mbd_charge_raw[ich]*gain_corr[ich]);
	    }
	}
    }    

      
  for (int ich = 0 ; ich < 128; ich++){

    std::cout << "Recalib Fitting :: Channel "<<ich;

    
    double local_min = 0;
    double local_max = 0;
    double max_value = 0;
    int good = 0;
    int low_bin = 0;
    int high_bin = 0;
    float upper;
    TH1D *hsmooth = (TH1D*)h_charge[ich]->Clone();
    hsmooth->Smooth();
    hsmooth->Smooth();
    // Finding the dip
    int tries = 0;
    while (tries < 5)
      {
	for (int ic = start_looking; ic < hsmooth->GetNbinsX(); ic++)
	  {
	    good = 0;
	    for (int j = -1*min_range; j <= min_range; j++)
	      {
		if (j == 0) continue;

		if (hsmooth->GetBinContent(ic+1) <= hsmooth->GetBinContent(ic+1 +j))
		  good++;

	      }
	    if (good == 2*min_range)
	      {
		local_min = hsmooth->GetBinLowEdge(ic - before);
		low_bin = ic - before;
		break;
	      }
	  }

	// Find local Max after the dip

	for (int ic = low_bin; ic < hsmooth->GetNbinsX(); ic++)
	  {
	    if (hsmooth->GetBinContent(ic+1) > max_value)
	      {
		high_bin = ic+1; 
		max_value = hsmooth->GetBinContent(high_bin);
		local_max = hsmooth->GetBinLowEdge(high_bin);
	      }
	  }

	if (local_max > 2){
	  local_max = 1;
	  local_min = 0.2;
	}
	else if (local_max < 0.3){
	  local_max = 1;
	  local_min = 0.2;
	}
	h_charge[ich]->Fit("f_exp","Q","", local_max + 2, 10);
	TF1 *f = (TF1*) h_charge_raw[ich]->GetFunction("lan_w_gausexp");
	f_lan_w_exp->SetParLimits(0, 1, 1300000);
	f_lan_w_exp->SetParLimits(1, 0.1, 2);
	f_lan_w_exp->SetParLimits(2, 0.001, 2);

	f_lan_w_exp->SetParameter(0, 10000);
	f_lan_w_exp->SetParameter(1, 1);
	f_lan_w_exp->SetParameter(2, .2);
	f_lan_w_exp->SetParameter(3, f_exp->GetParameter(0));
	f_lan_w_exp->SetParameter(4, f_exp->GetParameter(1));
	f_lan_w_exp->SetParameter(5, f_exp->GetParameter(2));
	upper = 1.5 - tries*0.05;
      
	hsmooth->Fit("lan_w_exp","Q","",local_min, upper);

	for (int t = 0; t< 5; t++)	h_charge[ich]->Fit("lan_w_exp","Q","",local_min, upper);

	std::cout <<" try " <<tries <<endl;
	if (f_lan_w_exp->GetParameter(1) > 1.4 || f_lan_w_exp->GetParameter(1) < 0.7) 
	  {
	    tries++;
	    continue;
	  }
	if (f_lan_w_exp->GetChisquare()/f_lan_w_exp->GetNDF() < 5)
	  {
	    break;
	  }
	tries++;
      }
    h_peaks->Fill(f_lan_w_exp->GetParameter(1));

	    //    gain_corr[ich] = gain_corr[ich]*f_lan_w_exp->GetParameter(1);

    std::cout <<" -- low/guess/high = " << local_min <<"/"<< local_max<<"/"<< upper<<" -- Peak at "<<f_lan_w_exp->GetParameter(1)<<endl;
    tn->Fill(ich, gain_corr[ich], ((double) f_lan_w_exp->GetParameter(2)));
  }
  fcalibout->Write();
  fcalibout->Close();

  // As function of number of tubes hit, RMS of the time dists;
  TFile *fout = new TFile(Form("%s/output/plots/mbd_calib_plots_%d.root", env_p, runnumber), "RECREATE");
  for (int ich = 0 ; ich < 128; ich++)
    {
      h_charge[ich]->Write();
      h_charge_raw[ich]->Write();
    }

  h_peaks->Write();
  h_peaks_raw->Write();
  fout->Close();

  double avg_chisquare = 0.;
  for (int i = 0; i < 128; i++)
    {
      avg_chisquare += (chisquare[i]/128.);
    }

  std::cout << " ************************************ " << std::endl;
  std::cout << "  Run " << runnumber << std::endl;
  std::cout << "  Avg. Chi^2/NDF = "<<avg_chisquare << std::endl;
  std::cout << " ************************************ " << std::endl;

}


void QA_ZDCCheck(const int runnumber)
{

  const char *env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }

  TFile *f2 = new TFile(Form("%s/output/run%d/mbdcalibana/mbd_calib_trees_%d.root", env_p, runnumber, runnumber), "r");
  TTree *t2 = (TTree*) f2->Get("T");

  TH1D *h_zdc_energy[6];
  TH1D *h_zdc_energy_prime[6];
  for (int i = 0; i < 6;i++) 
    {
      h_zdc_energy[i] = new TH1D(Form("h_zdc_energy_%d", i), ";ZDC Energy [GeV]; Counts", 500, 0, 5000);
      h_zdc_energy_prime[i] = new TH1D(Form("h_zdc_energy_prime_%d", i), ";ZDC Energy [GeV]; Counts", 500, 0, 5000);
    }
  TH1D *h_zdc_sum_n = new TH1D("h_zdc_sum_n", ";ZDC sum [GeV]; Counts", 1000, 0, 10000);
  TH1D *h_zdc_sum_s = new TH1D("h_zdc_sum_s", ";ZDC sum [GeV]; Counts", 1000, 0, 10000);
  TH1D *h_zdc_sum_n_prime = new TH1D("h_zdc_sum_n_prime", ";ZDC sum [GeV]; Counts", 1000, 0, 10000);
  TH1D *h_zdc_sum_s_prime = new TH1D("h_zdc_sum_s_prime", ";ZDC sum [GeV]; Counts", 1000, 0, 10000);

  float zdc_sum[2];
  float zdc_energy[6];

  float zdc_sum_prime[2];
  float zdc_energy_prime[6];

  t2->SetBranchAddress("zdc_sum", zdc_sum);
  t2->SetBranchAddress("zdc_energy", zdc_energy);
  t2->SetBranchAddress("zdc_sum_prime", zdc_sum_prime);
  t2->SetBranchAddress("zdc_energy_prime", zdc_energy_prime);

  int n_min_bias = 0;;
  int n_min_bias_prime = 0;;
  int n_events = t2->GetEntries();
  
  for (int i = 0; i < t2->GetEntries(); i++)
    {
      t2->GetEntry(i);
      h_zdc_sum_n->Fill(zdc_sum[1]);
      h_zdc_sum_s->Fill(zdc_sum[0]);
      h_zdc_sum_n_prime->Fill(zdc_sum_prime[1]);
      h_zdc_sum_s_prime->Fill(zdc_sum_prime[0]);

      if (zdc_sum[0] >= 40 && zdc_sum[1] >= 40) n_min_bias++;
      if (zdc_sum_prime[0] >= 40 && zdc_sum_prime[1] >= 40) n_min_bias_prime++;
      for (int j= 0; j < 6; j++)
	{
	  h_zdc_energy[j]->Fill(zdc_energy[j]);
	  h_zdc_energy_prime[j]->Fill(zdc_energy_prime[j]);

	}
    }  

  std::cout << "Number of events passing ZDC coincidence                      : "<<n_min_bias <<"/"<<n_events<<" = "<<static_cast<double>(n_min_bias)/static_cast<double>(n_events)<<std::endl;
  std::cout << "Number of events passing ZDC coincidence with our calibrations: "<<n_min_bias_prime <<"/"<<n_events<<" = "<<static_cast<double>(n_min_bias_prime)/static_cast<double>(n_events)<<std::endl;

  TFile *fout21 = new TFile(Form("%s//output/plots/zdc_check_%d.root", env_p, runnumber),"recreate");
  h_zdc_sum_n->Write();
  h_zdc_sum_s->Write();
  h_zdc_sum_n_prime->Write();
  h_zdc_sum_s_prime->Write();

  for (int i = 0; i < 6; i++)
    {
      h_zdc_energy[i]->Write();
      h_zdc_energy_prime[i]->Write();
    }
  fout21->Close();
}

void QA_MBDTimeCalibrations(const int runnumber, const int useChargeTime = 0)
{

  const char *env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }

  TFile *f = nullptr;
  if (useChargeTime) f = new TFile(Form("%s/output/run%d/signals/processedsignal_%d.root", env_p, runnumber, runnumber), "r");
  TTree *t = nullptr;
  if (useChargeTime) t = (TTree*) f->Get("ttree");

  TFile *f2 = new TFile(Form("%s/output/run%d/mbdcalibana/mbd_calib_trees_%d.root", env_p, runnumber, runnumber), "r");
  TTree *t2 = (TTree*) f2->Get("T");

  TFile *fcalib = new TFile(Form("%s/calib/calib_mbd_%d.root", env_p, runnumber), "r");
  TNtuple *tcalib = (TNtuple*) fcalib->Get("mbd_calib");

  // get calib
  float gain_corr[128];
  float g;

  tcalib->SetBranchAddress("peak", &g);

  for (int i = 0 ; i < 128; i ++)
    {
      tcalib->GetEntry(i);
      gain_corr[i] = g;
    }

  
  float low_t = -20.;
  float high_t = 20.;
  float low_c = 0;
  float high_c = 1000;

  float binlengtht = 0.01;
  float binlengthc = 1;

  int nbins_t = floor((high_t - low_t)/binlengtht);
  int nbins_t2 = floor((high_t - low_t)/(binlengtht*10.));
  int nbins_c = floor((high_c - low_c)/(binlengthc*10.));




  TH1D *h_toa_channel[128];
  TH1D *h_tdc_channel[128];
  TH2D *h_tdc_toa_channel[128];
  TH1D *h_ptime_channel[128];
  TH2D *h_toa_ptime_channel[128];
  TH2D *h_tdc_ptime_channel[128];
  TH2D *h_peak_ptime_channel[128];
  TH2D *h_peak_tdc_channel[128];
  TH2D *h_peak_toa_channel[128];

  for (int i = 0; i < 128; i++)
    {
      h_toa_channel[i] = new TH1D(Form("h_toa_channel_%d", i), "", nbins_t, low_t, high_t);
      h_ptime_channel[i] = new TH1D(Form("h_ptime_channel_%d", i), "", nbins_t, low_t, high_t);
      h_tdc_channel[i] = new TH1D(Form("h_tdc_channel_%d", i), "", nbins_t, low_t, high_t);
      h_tdc_toa_channel[i] = new TH2D(Form("h_tdc_toa_channel_%d", i), "", nbins_t2, low_t, high_t, nbins_t2, low_t, high_t);
      h_tdc_ptime_channel[i] = new TH2D(Form("h_tdc_ptime_channel_%d", i), "", nbins_t2, low_t, high_t, nbins_t2, low_t, high_t);
      h_toa_ptime_channel[i] = new TH2D(Form("h_toa_ptime_channel_%d", i), "", nbins_t2, low_t, high_t, nbins_t2, low_t, high_t);
      h_peak_ptime_channel[i] = new TH2D(Form("h_peak_ptime_channel_%d", i), "",nbins_c, low_c, high_c, nbins_t2, low_t, high_t);
      h_peak_tdc_channel[i] = new TH2D(Form("h_peak_tdc_channel_%d", i), "", nbins_c, low_c, high_c, nbins_t2, low_t, high_t);
      h_peak_toa_channel[i] = new TH2D(Form("h_peak_toa_channel_%d", i), "",nbins_c, low_c, high_c, nbins_t2, low_t, high_t);
    }


  float peak[256];
  float peaktime[256];
  float pedestal[256];
  float time_of_arrival[256];

  if (useChargeTime)
    {
      t->SetBranchAddress("peak", peak);
      t->SetBranchAddress("pedestal", pedestal);
      t->SetBranchAddress("peaktime", peaktime);
      t->SetBranchAddress("time_of_arrival", time_of_arrival);
    }
  float mbd_time[128];
  float mbd_charge[128];
  t2->SetBranchAddress("mbd_time_raw", mbd_time);
  t2->SetBranchAddress("mbd_charge_raw", mbd_charge);

  double center_tdc[128];
  double center_peak_time[128];
  double center_toa[128];

  int number_of_events[128];
  for (int i = 0; i < 128; i++) 
    {
      number_of_events[i] = 0;
      center_tdc[i] = 0.;
      center_toa[i] = 0.;
      center_peak_time[i] = 0.;
    }
  if (!SILENCE) std::cout << "\rChecking the number of events per channel ..."<<std::flush;
  for (int i = 0; i < t2->GetEntries(); i++)
    {

      t2->GetEntry(i);
      if ( !SILENCE && i != 0 && i%(t2->GetEntries()/100) == 0) std::cout << "\rChecking the number of events per channel ... "<<floor(static_cast<float>(i)/static_cast<float>(t2->GetEntries())*100)<<"%"<<std::flush;

      for (int ich = 0; ich < 128; ich++)
	{
	  if (mbd_charge[ich]*gain_corr[ich] < 0.4) continue;
	  number_of_events[ich]++;
	}
    }  
  if (!SILENCE)
    {
      std::cout << "\rChecking the number of events per channel.. [Done]"<<std::endl;
      std::cout << "\rCalculating Mean.. "<<std::flush;
    }
  for (int i = 0; i < t2->GetEntries(); i++)
    {
      if (useChargeTime) t->GetEntry(i);
      t2->GetEntry(i);
      if (!SILENCE && i != 0 && i%(t2->GetEntries()/100) == 0) std::cout << "\rCalculating Mean ... "<<floor(static_cast<float>(i)/static_cast<float>(t2->GetEntries())*100)<<"%"<<std::flush;
      for (int ich = 0; ich < 128; ich++)
	{
	  double toa = -999;
	  double ptime = -999;
	  if(useChargeTime)
	  {
	    toa = static_cast<double>(time_of_arrival[8 + (ich/8)*16 + ich%8]);
	    ptime = static_cast<double>(peaktime[(ich/16 + 1)*8 + ich%8]);
	  }	  

	  double tdc = static_cast<double>(25. - (9./5000.)*mbd_time[ich]);


	  if (mbd_charge[ich]*gain_corr[ich] < 0.4) continue;

	  center_tdc[ich] += tdc/static_cast<double>(number_of_events[ich]);
	  if (useChargeTime)
	    {
	      if (toa > 0)	
		center_toa[ich] += toa/static_cast<double>(number_of_events[ich]);
	      if (ptime > 0)
		center_peak_time[ich] += ptime/static_cast<double>(number_of_events[ich]);
	      
	    }
	}
    }  
  if (!SILENCE)
    {
      std::cout << "\rCalculating Mean ... [Done]"<<std::endl;
      std::cout << "Mean TDC     : "<< center_tdc[8] << std::endl;
      if (useChargeTime)
	{
	  std::cout << "Mean PeakTime: "<< center_peak_time[8] << std::endl;
	  std::cout << "Mean TOA     : "<< center_toa[8] << std::endl;
	}
      
      std::cout << "\rFilling Histograms.. "<<std::flush;
    }
  for (int i = 0; i < t2->GetEntries(); i++)
    {
      if (useChargeTime) t->GetEntry(i);
      t2->GetEntry(i);
      if (!SILENCE &&i != 0 && i%(t2->GetEntries()/100) == 0) std::cout << "\rFilling Histograms ... "<<floor(static_cast<float>(i)/static_cast<float>(t2->GetEntries())*100)<<"%"<<std::flush;
      for (int ich = 0; ich < 128; ich++)
	{

	  double toa = -999;
	  double ptime = -999;
	  if (mbd_charge[ich]*gain_corr[ich] < 0.4) continue;
	  double tdc = (25. - static_cast<double>(mbd_time[ich])*(9./5000.)) - center_tdc[ich];
	  if (useChargeTime)
	    {	  
	      toa = static_cast<double>(time_of_arrival[8 + (ich/8)*16 + ich%8]) - center_toa[ich];
	      ptime = static_cast<double>(peaktime[(ich/16 + 1)*8 + ich%8]) - center_peak_time[ich];
	      h_toa_channel[ich]->Fill(toa);
	      h_tdc_toa_channel[ich]->Fill(tdc, toa);
	      h_toa_ptime_channel[ich]->Fill(toa, ptime);
	      h_tdc_ptime_channel[ich]->Fill(tdc, ptime);
	      h_ptime_channel[ich]->Fill(ptime);
	      h_peak_ptime_channel[ich]->Fill(peak[8 + (ich/8)*16 + ich%8] - pedestal[8 + (ich/8)*16 + ich%8], ptime);
	      h_peak_toa_channel[ich]->Fill(peak[8 + (ich/8)*16 + ich%8] - pedestal[8 + (ich/8)*16 + ich%8], toa);
	      h_peak_tdc_channel[ich]->Fill(peak[8 + (ich/8)*16 + ich%8] - pedestal[8 + (ich/8)*16 + ich%8], tdc);
	    }
	  h_tdc_channel[ich]->Fill(tdc);
	}
    }  
  if (!SILENCE) std::cout << "\rFilling Histograms ... [Done]"<<std::endl;
  
  TFile *ftcalib = new TFile(Form("%s/calib/t0_calib_mbd_%d.root", env_p, runnumber),"recreate");

  TNtuple *t_calib = new TNtuple("ttree","time calib","channel:shift:scale:toa_shift");
  if (!SILENCE)
    std::cout << "\rFitting 1D Histograms.. "<<std::flush;
  
  double middle_fit_toa[128];
  double middle_fit_tdc[128];

  double scale_fit_toa[128];
  double scale_fit_tdc[128];
  if (!SILENCE)
    std::cout << "Tube \t TDC shift \t TDC scale \t TOA Shift"<<endl;
  for (int ich = 0; ich < 128; ich++)
    {
      if (!SILENCE)
	std::cout << "channel : "<<ich <<std::endl;
      TF1 *fit_toa = 0;
      TF1 *fit_tdc = 0;
      TH1D *hclone = nullptr;
      int nparams = 9;
      double params[9];
      std::vector<std::pair<double, double>> fit_peak_heights;
      if (useChargeTime)
	{
	  if (!SILENCE)
	    {
	      std::cout << " "<<std::endl;
	      std::cout << "TOA "<<std::endl;
	      std::cout << " "<<std::endl;
	    }
	  std::pair<double, double> range = getFitFunction(h_toa_channel[ich], params);
	  fit_toa = new TF1("triple_gaussian","gaus(0) + gaus(3) + gaus(6)",range.first, range.second);
	  
	  if (params[6] < 0)
	    {
	      fit_toa = new TF1("double_gaussian","gaus(0) + gaus(3)",range.first, range.second);
	      nparams = 6;
	    }
	  for (int i = 0 ; i < nparams; i++)
	    {
	      fit_toa->SetParameter(i, params[i]);
	      if (i%3 == 2) fit_toa->SetParLimits(i, 0, 100);
	    }
	  
	  hclone = (TH1D*)h_toa_channel[ich]->Clone();
	  hclone->Rebin(20);
	  hclone->Fit(fit_toa->GetName(), "Q","", range.first, range.second);
	  
	  for (int i = 0; i < nparams/3 ; i++)
	    {
	      fit_peak_heights.push_back(std::make_pair(fit_toa->GetParameter(i*3 + 1), fit_toa->GetParameter(i*3)));
	    }
	  
	  std::sort(fit_peak_heights.begin(), fit_peak_heights.end(), [](auto &left, auto &right) {
	      return left.second > right.second;
	    });
	  if (!SILENCE) {
	    std::cout <<" Peak heights: "<<std::flush;
	    for (int i = 0; i < nparams/3 ; i++)
	      {
		std::cout << " " <<fit_peak_heights.at(i).first<<","<<fit_peak_heights.at(i).second;
	      }
	  }
	  middle_fit_toa[ich] = center_toa[ich] + fit_peak_heights.at(0).first;
	  
	  if (nparams == 6)
	    scale_fit_toa[ich] = fabs(fit_peak_heights.at(1).first - fit_peak_heights.at(0).first);
	  
	  scale_fit_toa[ich] = fabs(fit_peak_heights.at(1).first - fit_peak_heights.at(0).first);
	}
      else scale_fit_toa[ich] = 1;
      
      if (ich == 54 || ich == 56 || ich == 56+64 ){
	middle_fit_tdc[ich] = 0;
	scale_fit_tdc[ich] = -999;
      }
      else
	{	  
	  if (!SILENCE)
	    {
	      std::cout << " "<<std::endl;
	      std::cout << "TDC "<<std::endl;
	      std::cout << " "<<std::endl;
	    }
	  std::pair<double, double> range = getFitFunction(h_tdc_channel[ich], params, -15., 15., !useChargeTime);
	  
	  if (params[3] < 0)
	    {
	      fit_tdc = new TF1("single_gaussian","gaus(0)",range.first, range.second);
	      nparams = 3;
	    }
	  else if (params[6] < 0)
	    {
	      fit_tdc = new TF1("double_gaussian","gaus(0)+gaus(3)",range.first, range.second);
	      nparams = 6;
	    }
	  else 
	    {
	      fit_tdc = new TF1("triple_gaussian","gaus(0)+gaus(3)+gaus(6)",range.first, range.second);
	      nparams = 9;
	    }
	  for (int i = 0 ; i < nparams; i++)
	    {
	      fit_tdc->SetParameter(i, params[i]);
	      if (i%3 == 2) fit_tdc->SetParLimits(i, 0, 100);
	    }
	  hclone = (TH1D*) h_tdc_channel[ich]->Clone();
	  hclone->Rebin(20);
	  hclone->Smooth();
	  hclone->Fit(fit_tdc->GetName(), "Q","", range.first, range.second);
	  h_tdc_channel[ich]->Rebin(20);
	  h_tdc_channel[ich]->Fit(fit_tdc->GetName(), "Q+","", range.first, range.second);

	  std::vector<std::pair<double, double>> fit_tdc_peak_heights;
	  for (int i = 0; i < nparams/3 ; i++)
	    {
	      fit_tdc_peak_heights.push_back(std::make_pair(fit_tdc->GetParameter(i*3 + 1), fit_tdc->GetParameter(i*3)));
	    }
	      
	  std::sort(fit_tdc_peak_heights.begin(), fit_tdc_peak_heights.end(), [](auto &left, auto &right) {
	      return left.second > right.second;
	    });
	      
	      
	  middle_fit_tdc[ich] = center_tdc[ich] + fit_tdc_peak_heights.at(0).first;
      
	  scale_fit_tdc[ich] = 1;
	  if (useChargeTime) scale_fit_tdc[ich] = scale_fit_toa[ich] / (fabs(fit_tdc_peak_heights.at(1).first - fit_tdc_peak_heights.at(0).first));
	}
      t_calib->Fill(ich, middle_fit_tdc[ich],scale_fit_tdc[ich], middle_fit_toa[ich]);
    }
  if (!SILENCE)  std::cout << "\rFitting 1D Histograms.. [done]"<<std::endl;
  ftcalib->Write();
  ftcalib->Close();

  TFile *fout21 = new TFile(Form("%s//output/plots/fitsout_%d.root", env_p, runnumber),"recreate");
  for (int i = 0; i < 128; i++)
    {
      h_tdc_channel[i]->Write();
      if (useChargeTime)
	{
	  h_toa_channel[i]->Write();
	  h_ptime_channel[i]->Write();
	  h_tdc_toa_channel[i]->Write();
	  h_tdc_ptime_channel[i]->Write();
	  h_toa_ptime_channel[i]->Write();
	  h_peak_ptime_channel[i]->Write();
	  h_peak_tdc_channel[i]->Write();
	  h_peak_toa_channel[i]->Write();
	}
    }
  fout21->Close();

  for (int ich = 0; ich < 128; ich++)
    {
      h_toa_channel[ich]->Reset();
      h_tdc_channel[ich]->Reset();
      h_tdc_toa_channel[ich]->Reset();
      h_toa_ptime_channel[ich]->Reset();
      h_tdc_ptime_channel[ich]->Reset();
      h_peak_toa_channel[ich]->Reset();
      h_peak_tdc_channel[ich]->Reset();

    }

  if (!SILENCE) std::cout << "\rFilling Histograms Again.. "<<std::flush;
  for (int i = 0; i < t2->GetEntries(); i++)
    {
      if (useChargeTime) t->GetEntry(i);
      t2->GetEntry(i);
      if (!SILENCE && i != 0 && i%(t2->GetEntries()/100) == 0) std::cout << "\rFilling Histograms ... "<<floor(static_cast<float>(i)/static_cast<float>(t2->GetEntries())*100)<<"%"<<std::flush;
      for (int ich = 0; ich < 128; ich++)
	{


	  if (mbd_charge[ich]*gain_corr[ich] < 0.4) continue;
	  double tdc = scale_fit_tdc[ich]*((25. - static_cast<double>(mbd_time[ich])*(9./5000.)) - middle_fit_tdc[ich]);
	  if (useChargeTime)
	    {
	      double toa = static_cast<double>(time_of_arrival[8 + (ich/8)*16 + ich%8]) - middle_fit_toa[ich];
	      double ptime = static_cast<double>(peaktime[(ich/16 + 1)*8 + ich%8]) - center_peak_time[ich];
	      
	      h_toa_channel[ich]->Fill(toa);
	      h_tdc_toa_channel[ich]->Fill(tdc, toa);
	      h_toa_ptime_channel[ich]->Fill(toa, ptime);
	      h_tdc_ptime_channel[ich]->Fill(tdc, ptime);
	      h_peak_toa_channel[ich]->Fill(peak[8 + (ich/8)*16 + ich%8] - pedestal[8 + (ich/8)*16 + ich%8], toa);
	      h_peak_tdc_channel[ich]->Fill(peak[8 + (ich/8)*16 + ich%8] - pedestal[8 + (ich/8)*16 + ich%8], tdc);
	    }
	  h_tdc_channel[ich]->Fill(tdc);

	}
    }  
  if (!SILENCE)
    std::cout << "\rFilling Histograms ... [Done]"<<std::endl;


  TFile *fout2 = new TFile(Form("%s/output/plots/histout_%d.root", env_p, runnumber),"recreate");
  for (int i = 0; i < 128; i++)
    {
      h_tdc_channel[i]->Write();
      if (useChargeTime)
	{
	  h_toa_channel[i]->Write();
	  h_ptime_channel[i]->Write();
	  h_tdc_toa_channel[i]->Write();
	  h_tdc_ptime_channel[i]->Write();
	  h_toa_ptime_channel[i]->Write();
	  h_peak_ptime_channel[i]->Write();
	  h_peak_tdc_channel[i]->Write();
	  h_peak_toa_channel[i]->Write();
	}
    }
  fout2->Close();
}

std::pair<double, double> getFitFunction(TH1D* h, double *params, float minimumx = -25., float maximumx = 25., bool ForceSingle  = false){

  TF1 *f;
  TH1D *h_rebinned = (TH1D*) h->Clone();
  h_rebinned->Rebin(20);
  h_rebinned->Smooth();
  // finding the fit range   
  
  int bin_low = -1;
  double x_low = -999;
  int bin_high = -1;
  double x_high = -999;
  double threshold = 0.9;
  bool foundPeak = false;
  bool foundValley = false;
  std::vector<std::pair<double,double>> v_max = {};
  std::vector<double> v_minx = {};
  double maxx = -999.99;
  double maxy = -999.99;
  double miny = -999.99;
  double minx = -999.99;
  double maxvalue;
  for (int ib = 1; ib < h_rebinned->GetNbinsX()-1; ib++)
    {
      // find first bin with something in it
      if (bin_low == -1 && h_rebinned->GetBinContent(ib) > 0)
	{
	  bin_low = ib;
	  x_low = h_rebinned->GetBinLowEdge(ib);
	  v_minx.push_back(x_low);
		  
	}

      else if (bin_low == -1) continue;

      // find the peak;
      if (!foundPeak)
	{

	  // save the maximum
	  if (h_rebinned->GetBinContent(ib) > maxy)
	    {
	      maxy = h_rebinned->GetBinContent(ib);
	      maxx = h_rebinned->GetBinCenter(ib);
	    }
	  else if (h_rebinned->GetBinContent(ib) < maxy*threshold)
	    {
	      foundPeak = true;
	      if (foundValley)
		{
		  v_minx.push_back(minx);
		  foundValley = false;
		}

	      
	      std::pair<double, double> max = std::make_pair(maxx, maxy);
	      v_max.push_back(max);

	      maxx = -999.99;
	      maxy = -999.99;
	    }
	  if (foundValley)
	    {
	      if (h_rebinned->GetBinContent(ib) < miny)
		{
		  miny = h_rebinned->GetBinContent(ib - 1);
		  minx = h_rebinned->GetBinCenter(ib - 1);
		}
	    }
	  continue;
	}
      if (!foundValley)
	{
	  if (h_rebinned->GetBinContent(ib) > h_rebinned->GetBinContent(ib - 1))
	    {
	      miny = h_rebinned->GetBinContent(ib - 1);
	      minx = h_rebinned->GetBinCenter(ib - 1);
	      foundValley = true;
	      foundPeak = false;
	    }
	  if (h_rebinned->GetBinContent(ib) == 0 || h_rebinned->GetBinCenter(ib) > 18)
	    {
	      bin_high = ib;
	      x_high = h_rebinned->GetBinLowEdge(ib);
	      break;
	    }
	}
    }

  if (bin_high == -1)
    {
      bin_high = h_rebinned->GetNbinsX() - 1;
      x_high = h_rebinned->GetBinLowEdge(bin_high);
    }

  // look for most maximum maximums
  maxvalue = 2000;
  if(v_max.size() > 1 )
    {
      std::sort(v_max.begin(), v_max.end(), [](auto &left, auto &right) {
	  return left.second > right.second;
	});

      maxvalue = 0.2*v_max.at(0).second;

      std::sort(v_max.begin(), v_max.end(), [](auto &left, auto &right) {
	  return left.first < right.first;
	});
    }
  auto maxv = v_max.begin();
  while (maxv != v_max.end())
    {
      if ((*maxv).second < maxvalue)
	v_max.erase(maxv);
      else if ((*maxv).first < minimumx || (*maxv).first > maximumx)
	v_max.erase(maxv);
      else
	maxv++;
    }



  if (v_max.size() == 0)
    {
      return std::make_pair(-20, 20);
    }
  

  if (ForceSingle)
    {
      std::sort(v_max.begin(), v_max.end(), [](auto &left, auto &right) {
	  return left.second > right.second;
	});
      
      params[0] = v_max.at(0).second;
      params[1] = v_max.at(0).first;
      params[2] = 1;
      params[3] = -999.99;
      params[4] = -999.99;
      params[5] = -999.99;
      params[6] = -999.99;
      params[7] = -999.99;
      params[8] = -999.99;

      return std::make_pair( (x_low > v_max.at(0).first - 1 ? x_low :  v_max.at(0).first - 1), v_max.at(0).first + 1);
    }

  if (v_max.size() == 1)
    {
      if (v_minx.size())
	{
	  x_high = v_minx.at(0);
	  
	  while (v_minx.begin() != v_minx.end())
	    {
	      v_minx.erase(v_minx.begin());
	    }
	}
    }
  else
    {
      if(v_max.size() > 3 )
	{

	  std::sort(v_max.begin(), v_max.end(), [](auto &left, auto &right) {
	      return left.second > right.second;
	    });
	  auto maxv = v_max.begin() + 3;
	  while (maxv != v_max.end())
	    {
	      v_max.erase(maxv);
	    }	  
	  std::sort(v_max.begin(), v_max.end(), [](auto &left, auto &right) {
	      return left.first < right.first;
	    });
	}
      
      double last_x;
      bool first = false;
      bool last = false;
      auto minv = v_minx.begin();
      while (minv != v_minx.end())
	{
	  if ((*minv) > v_max.at(0).first && (*minv) < v_max.at(v_max.size() - 1).first)
	    {
	      if (!first)
		{
		  x_low = last_x;
		  first = true;
		}
	      minv++;
	      continue;
	    }
	  else if (first && !last)
	    {
	      x_high = (*minv);
	      last = true;
	    }
	  last_x = (*minv);
	  v_minx.erase(minv);

	}

    }  

  if (v_max.size() == 3)
    {
      params[0] = v_max.at(0).second;
      params[1] = v_max.at(0).first;
      params[2] = (v_minx.at(0) - x_low)/4.;
      params[3] = v_max.at(1).second;
      params[4] = v_max.at(1).first;
      params[5] = (v_minx.at(1) - v_minx.at(0))/4.;
      params[6] = v_max.at(2).second;
      params[7] = v_max.at(2).first;
      params[8] = (x_high - v_minx.at(1))/4.;
    }

  else if (v_max.size() == 2)
    {
      f = new TF1("double_gaussian", "gaus(0) + gaus(3)", x_low, x_high);
      params[0] = v_max.at(0).second;
      params[1] = v_max.at(0).first;
      params[2] = (v_minx.at(0) - x_low)/4.;
      params[3] = v_max.at(1).second;
      params[4] = v_max.at(1).first;
      params[5] = (x_high - v_minx.at(0))/4.;
      params[6] = -999.99;
      params[7] = -999.99;
      params[8] = -999.99;

    }


  
  return std::make_pair(x_low, x_high);

}

void QA_MakeCentralityCalibrations(const int runnumber, const bool doVertexSelection = 1, const bool use_shifted = false)
{

 
  const char *env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }
  ///////////// 


  //===============================================================================================
  double mu = 4.4;  //defaults that are used if not determining these parameters again
  double k  = 2.14;
  bool npartminusone = false;
  bool npartscaling = true;
  bool flag_determineBiasFactors = false;
  bool flag_determineNBDparameters = true;
  double forcetrigfrac = 100.;
  double alpha = 1.00;
  double particlealpha = 1.00;
  double biasfactor = 1.55;
  double biased_mu = mu * biasfactor; // larger number of hits
  double biased_k  =  k * biasfactor;
  const int maxcentbins  = 19;

  bool sphenix = true;
  int vtxflag = 1;
  std::string name = Form("%s/glauber/lemon_glauber_auau200_100k_histo.root", env_p);
  int lowfit = 300; // lowest end of standard fit range, was 30 in PHENIX
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

  sprintf(fdataname,"%s/output/run%d/plots/mbd_charge_sum_%d.root", env_p, runnumber, runnumber);  
  //  if (sphenix) sprintf(fdataname,"../output/run%d/trees_%d.root", runnumber, runnumber);
  //if (sphenix) fdataname = "./fout_sphenix2023_auau200_23722_both_z20.root";
  //  if (sphenix) fdataname = "./fout_sphenix2023_auau200_22979_both_z20.root";

  char runnum[100];
  sprintf(runnum, "000%d", runnumber);
  char zselect[100] = "|z_{MBD}|< 30 cm";
  //=======================================================  
  
  TFile *fdata = new TFile(fdataname);
  int which = 0;
  TH1D *hRealBBC;
  //if (which==0) hRealBBC = (TH1D *) fdata->Get("h_mbd_charge_ns_w_zdc_cut_w_mbd_cut_and_vertex_10"); // a particular z-vertex range
  if (which==0) 
    {
      if (!use_shifted) hRealBBC = (TH1D *) fdata->Get("h_charge_sum_min_bias_w_vertex_30"); // a particular z-vertex range
      else hRealBBC = (TH1D *) fdata->Get("h_charge_sum_min_bias_w_vertex_30_shifted"); // a particular z-vertex range
    }
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
  hRatio->Fit(trigeffcurve,"NDOR","",1.0,80.0);
  // also fit a flat line above say 20 -> 140 ==> but just constraint it exactly at 1.0
  flatline->SetParameter(0,1.0);
  hRatio->Fit(flatline,"NDOR","",(double)lowfit,maxrange);
  
  // trigger efficiency from integral comparison
  double trigeffintegral = hRealBBC->Integral()/hSimBBC->Integral();
  cout << "Trigger Efficiency from Integrals = " << trigeffintegral << endl;

  TH1D *hSimBBCwTrig = new TH1D("hSimBBCwTrig","hSimBBCwTrig",nhistbins,-0.5,maxrange);
  for (int i=1;i<=nhistbins;i++) {
    if (i==1) {
      hSimBBCwTrig->SetBinContent(i,0.0); // no chance to fire the trigger if no BBC south hits
    } else {
      hSimBBCwTrig->SetBinContent(i,hSimBBC->GetBinContent(i)*trigeffcurve->Eval(hSimBBC->GetBinCenter(i)));
    }
  }

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
  int centrality_got[20] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  centrality_low[maxcentbins-1] = 00;
  int ihigh;
  for (ihigh = nhistbins; ihigh > 0; ihigh--)
    {
      if (hRealBBC->GetBinContent(ihigh))
	{
	  centrality_high[0] = hRealBBC->GetBinCenter(ihigh);
	  break;
	}
    }
  int cent_bin = 0;
  for (int i=ihigh;i>=1;i--) {

    if (maxcentbins == 4) {
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > ((cent_bin + 1)*20./forcetrigfrac) && centrality_got[cent_bin] == 0) {
	centrality_low[cent_bin] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[cent_bin+1] = hSimBBCwTrig->GetBinCenter(i);
	centrality_got[cent_bin] = 1;
	cent_bin++;
      }
    } 
    if (maxcentbins==9 || maxcentbins == 8) {
      if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > ((cent_bin + 1)*10./forcetrigfrac) && centrality_got[cent_bin] == 0) {
	centrality_low[cent_bin] = hSimBBCwTrig->GetBinCenter(i);
	centrality_high[cent_bin+1] = hSimBBCwTrig->GetBinCenter(i);
	centrality_got[cent_bin] = 1;
	cent_bin++;
      }
    }
    if (maxcentbins ==19) 
      {
	if ((hSimBBCwTrig->Integral(i,nhistbins)/hSimBBCwTrig->Integral()) > ((cent_bin + 1)*5./forcetrigfrac) && centrality_got[cent_bin] == 0) {
	  centrality_low[cent_bin] = hSimBBCwTrig->GetBinCenter(i);
	  centrality_high[cent_bin+1] = hSimBBCwTrig->GetBinCenter(i);
	  centrality_got[cent_bin] = 1;
	  cent_bin++;
	}
      }
  }
  

  /////////////

  for (int j=0; j<maxcentbins; j++) {
    cout << "Centbin = " << j << " Lowcut = " << 
      (int)centrality_low[j] << 
        " Highcut = " << 
        centrality_high[j] << endl;
  }

  cout << "mean data: " << hRealBBC->GetMean()<<endl;
  cout << "mean sim : " << hSimBBCwTrig->GetMean()<<endl;
  ///////////// 
  TFile *fcalib = new TFile(Form("%s/calib/calib_centrality_%d.root", env_p, runnumber), "recreate");
  TNtuple *ts = new TNtuple("tn_centrality", "holds centrality divisions", "bin:low:high");
  for (int i = 0; i < maxcentbins ; i++)
    {
      ts->Fill(i, centrality_low[i], centrality_high[i]);
    }
  fcalib->Write();
  fcalib->Close();


  TFile *fout = new TFile(Form("%s/output/plots/mbd_centrality_trigeff_%d.root", env_p, runnumber),"recreate");

  hSimBBCwTrig->Write();
  hSimBBC->Write();
  hRealBBC->Write();
  hRatio->Write();
  trigeffcurve->Write();

  fout->Close();

			  
} // end routine

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
int GetCentBin(float sum)
{
  for (int i = 0; i < 20; i++)
    {
      if (sum > centrality_map[i])
	{
	  return i;
	}
    }
  return 19;
}

std::pair<int, string> tokenize_path(string path)
{

  if (path.find("centrality") == std::string::npos)
    {
      return make_pair(0, "tnull");
    }

  string dash = "/";
  string under = "_";
  string dot = ".";


  string name, part, r;
  int run;
  string spath = path;
  size_t pos = 0;
  while ((pos = spath.find(dash)) != std::string::npos)
    {
      spath.erase(0, pos + dash.length());
    }
  
  pos = 0;
  while ((pos = spath.find(dot)) != std::string::npos)
    {
      name = spath.substr(0, pos);
      break;
    }

  while ((pos = name.find(under)) != std::string::npos)
    {
      part = name.substr(0, pos);

      name.erase(0, pos + dash.length());
    }



  if (strcmp(part.c_str(), "centrality"))
    return std::make_pair(0, "tnull");

  run = stoi(name);

  return make_pair(run, path);
}
