/* QA_centrality */
// Author: Daniel Lis - August 15 - 2023
// Brief : This macro class gives a quick QA analysis
//         of a run in sPHENIX
//      To be run on the output of the CentralityReco Module

const int DEBUG =  0;
//void QA_FindCentralities()


;
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


void QA_MBDCalibrations(const int runnumber, const int doPerRunCalibration, const int makeNewHistograms);
void QA_MakeChargeSum(const int runnumber, const int loadRunCalibration);


std::pair<double, double> getFitFunction(TH1D* h, double *params, float minimumx = -25., float maximumx = 25.);
void QA_MBDTimeCalibrations(const int runnumber);

void QA_centrality(
		   const int runnumber,
		   const int loadRunCalibration = 1,
		   const int doPerRunCalibration = 1,
		   const int doFindCentralities = 1,
		   const int doVertexSelection = 1
		   )
{
  QA_MBDCalibrations(runnumber, doPerRunCalibration, 1);  
  QA_MBDTimeCalibrations(runnumber);  
  QA_MakeChargeSum(runnumber, loadRunCalibration);
  return;
}

void QA_MakeChargeSum(const int runnumber, const int loadRunCalibration)
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
  int central_cut = 6;
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
      
      tcalibfile->Print();
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
  TH1D *h_vertex_w_t = new TH1D("h_vertex_w_t", "", nbin, -250, 250);
  TH1D *h_vertex_wo_t = new TH1D("h_vertex_wo_t", "", nbin, -250, 250);
  TH1D *h_time_0 = new TH1D("h_time_0", "", nbin, -25, 25);

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

  TH1D *event_with_weird_time[100][2];

  TH2Poly *h_timemap_tot[100][2];
  TH2Poly *h_hitmap_tot[100][2];
  for (int i = 0 ; i < 100; i++)
    {
      event_with_weird_time[i][0] = new TH1D(Form("event_%d_with_weird_time_s",i),"", 120, -5, 25);
      event_with_weird_time[i][1] = new TH1D(Form("event_%d_with_weird_time_n", i),"", 120, -5, 25);
      for (int j = 0; j < 2; j++)
	{
	  h_hitmap_tot[i][j] = new TH2Poly(Form("h_hitmap_%s_tot_%d", (j?"n":"s"), i),";x (cm);y (cm);Charge   ",-width,width,-height,height);
	  h_timemap_tot[i][j] = new TH2Poly(Form("h_timemap_%s_tot_%d", (j?"n":"s"), i),";x (cm);y (cm);Time   ",-width,width,-height,height);
	  mbd_honeycomb(h_hitmap_tot[i][j],-width,-height,a,10,11);
	  mbd_honeycomb(h_timemap_tot[i][j],-width,-height,a,10,11);
	}
    }
  int iweird = 0;
  int weird_event[100];
  float weird_time[100][128];
  float weird_charge[100][128];

  TH1D *h_time_w_t[128];
  TH1D *h_time_wo_t[128];

  for (int ich = 0 ; ich < 128; ich++)
    {
      h_time_w_t[ich]  = new TH1D(Form("h_time_w_t_%d", ich), "", nbins_t, low_t, high_t);
      h_time_wo_t[ich]  = new TH1D(Form("h_time_wo_t_%d", ich), "", nbins_t, low_t, high_t);
    }

  TFile *file = new TFile(Form("%s/output/run%d/trees_%d.root", env_p, runnumber, runnumber), "r");
  if (!file)
    {
      std::cout << " No Tree File found " <<std::endl;
      return;
    }

  // Making fresh histograms

  float zdc_sum[2];
  
  float mbd_charge[128];
  float mbd_time[128];
  float mbd_charge_raw[128];
  float mbd_time_raw[128];

  float z_vertex;
  float time_0;  
  TTree *t = (TTree*)file->Get("T");
  t->SetBranchAddress("zdc_sum_low",zdc_sum);
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

  for (int i = 0 ; i < (DEBUG ? 10 : t->GetEntries()); i++)
    {
      if (DEBUG) std::cout <<" ------------------------------- Event: "<<i<< " -------------------------------"<<std::endl;
      if (DEBUG) std::cout << " Channel\tCharge\tTDC"<<endl; 
      for (int ich = 0 ; ich < 128; ich++)
	{
	  mbd_time[ich] = ((25. - mbd_time_raw[ich] * (9.0 / 5000.) - time_shift_corr[ich]))*time_scale_corr[ich];
	  mbd_charge[ich] = mbd_charge_raw[ich]*gain_corr[ich];

	  if (DEBUG) std::cout << " "<<ich<<"\t"<<mbd_charge[ich]<<"\t"<< mbd_time[ich]<<endl; 
	}
      for (int ich = 0 ; ich < 64; ich++)
	{


	  if (mbd_charge[ich] > cthresh)
	    {
	      hits_s++;
	      charge_sum += mbd_charge[ich];	      
	      if (!(ich == 56 || fabs(mbd_time[ich]) > 15))
		{ 
		  
		  float timme = mbd_time[ich];
		  hits_s_t++;
		  time_sum_s.push_back(timme);
		  
		  sum_s += timme;
		  sum_s2 += (timme*timme);
		}      
	    }
	  if (mbd_charge[ich + 64] > cthresh)
	    {
	      hits_n++;
	      charge_sum += mbd_charge[ich+64];	      
	      if (!(ich == 56 || fabs(mbd_time[ich+64]) > 15))
		{ 
		  
		  float timme = mbd_time[ich+64];
		  
		  
		  hits_n_t++;
		  time_sum_n.push_back(timme);
		  sum_n += timme;;
		  sum_n2 += (timme*timme);

		}
	    }
      
	} 
      
      sort(time_sum_n.begin(), time_sum_n.end());
      sort(time_sum_s.begin(), time_sum_s.end());
      float mean_north;
      float mean_south;
      if (hits_s_t >= central_cut && hits_n_t >=central_cut){
	minbias = true;
	mean_north = sum_n/static_cast<float>(hits_n_t);
	mean_south = sum_s/static_cast<float>(hits_s_t);

	float rms_n = sqrt(sum_n2/static_cast<float>(hits_n_t) - TMath::Power(mean_north, 2));
	float rms_s = sqrt(sum_s2/static_cast<float>(hits_s_t) - TMath::Power(mean_south, 2));
	int nhit_n_center = 0;
	int nhit_s_center = 0;
	float sum_n_center = 0.;
	float sum_s_center = 0.;
	for (unsigned int ino = 0; ino < time_sum_n.size(); ino++)
	  {
	    if (fabs(time_sum_n.at(ino) - mean_north) < sigma_cut )
	      {
		sum_n_center += time_sum_n.at(ino);
		nhit_n_center++;
	      }
	  }
	for (unsigned int is = 0; is < time_sum_s.size(); is++)
	  {
	    if (fabs(time_sum_s.at(is) - mean_south) < sigma_cut )
	      {
		sum_s_center += time_sum_s.at(is);
		nhit_s_center++;
	      }
	  }
	float mean_north_center = sum_n_center/static_cast<float>(nhit_n_center);
	float mean_south_center = sum_s_center/static_cast<float>(nhit_s_center);
	mean_north = mean_north_center;
	mean_south = mean_south_center;
	z_vertex = 15.*(mean_north_center - mean_south_center);
	time_0 = (mean_north_center + mean_south_center)/2.;

	h_rms_by_time_n->Fill(time_0, rms_n);
	h_rms_by_time_s->Fill(time_0, rms_s);


      }

      else if (hits_s >=2 && hits_n >=2 && (hits_s_t > 1 && hits_n_t > 1)){

	minbias = true;
	mean_north = sum_n/static_cast<float>(hits_n);
	mean_south = sum_s/static_cast<float>(hits_s);

	z_vertex = 15*(mean_north - mean_south);
	time_0 = (mean_north + mean_south)/2.;

      }
      else z_vertex = -999;


      h_vertex->Fill(z_vertex);
      h_mean_north_south->Fill(mean_north, mean_south);
      
      //h_vertex_w_t->Fill(z_vertex);
      
      int cent_bin =GetCentBin(charge_sum);
      
      h_vertex_c[cent_bin]->Fill(z_vertex);
      h_time_0->Fill(time_0);
      h_charge_sum->Fill(charge_sum);
      
      if (!minbias) continue;

      h_charge_sum_min_bias->Fill(charge_sum);
  
      if (TMath::Abs(z_vertex) > 30) continue;
  
      h_charge_sum_min_bias_w_vertex_30->Fill(charge_sum);
    }

  TFile *fout = new TFile(Form("%s/output/plots/mbd_charge_sum_%d.root", env_p, runnumber), "RECREATE");

  TTree *tweird = new TTree("tweird", "holds weird events and such");
  int wevent;
  float wcharge[128];
  float wtime[128];
  tweird->Branch("event", &wevent, "event/I");
  tweird->Branch("charge", wcharge, "charge[128]/F");
  tweird->Branch("time", wtime, "time[128]/F");
  for (int i = 0; i < iweird;i++)
    {
      wevent = weird_event[i];
      for (int j = 0; j < 128; j++)
	{
	  wcharge[j] = weird_charge[i][j];
	  wtime[j] = weird_time[i][j];
	}
      tweird->Fill();
    }
  fout->Write();
  h_charge_sum->Write();
  h_charge_sum_min_bias->Write();
  h_charge_sum_min_bias_w_vertex_30->Write();
  h_vertex->Write();
  h_vertex_w_t->Write();
  h_vertex_wo_t->Write();
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
  for (int i = 0; i < iweird;i++)
    {      
      for (int j = 0; j < 2; j++)
	{
	  h_hitmap_tot[i][j]->Write();
	  h_timemap_tot[i][j]->Write();
	  event_with_weird_time[i][j]->Write();
	}
    }
  fout->Close();
  
  
}

void QA_MBDCalibrations(const int runnumber, const int doPerRunCalibration, const int makeNewHistograms)
{
  
  // Making the calibration file for the MBD

  int before = 2;
  int start_looking = 10;
  int min_range = 8;
  float gain_corr[128];
  float g;

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
  for (int j = 0; j < 128; j++)
    {
      h_charge[j] = new TH1D(Form("h_mbd_charge_ch%d", j), "",500, 0, 10);
      h_charge_raw[j] = new TH1D(Form("h_mbd_charge_raw_ch%d", j),"", 600, -0.5, 599.5);
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

  t->SetBranchAddress("zdc_sum_low", zdc_sum);
  t->SetBranchAddress("mbd_charge_raw",mbd_charge_raw);
  t->SetBranchAddress("mbd_time_raw",mbd_time_raw);

  if (makeNewHistograms)
    {
      for (int i = 0; i < t->GetEntries();i++)
	{
	  if (i%10000 == 0) cout << i << "zdc: "<< zdc_sum[0] <<" " << zdc_sum[1]<<endl;
	  t->GetEntry(i);
	  
	  if (zdc_sum[0] < 40 || zdc_sum[1] < 40) continue;
	  for (int ich = 0 ; ich < 128; ich++)
	    {
	      h_charge_raw[ich]->Fill(mbd_charge_raw[ich]);
	    }
	}
    }

  TH1D *h_peaks = new TH1D("h_peaks","", 64, -0.5, 1.5);
  TH1D *h_peaks_raw = new TH1D("h_peaks_raw","", 64, 0, 500);

  TF1 *f_lan_w_exp = new TF1("lan_w_exp","[0]*TMath::Landau(x,[1],[2],3) + expo(3)");
  TF1 *f_lan_w_gausexp = new TF1("lan_w_gausexp","[0]*TMath::Landau(x,[1],[2],3) + expo(3) + gaus(5)");

  f_lan_w_exp->SetParLimits(0, 100, 10000000);
  f_lan_w_exp->SetParLimits(1, 30, 400);
  f_lan_w_exp->SetParLimits(2, 0.1, 4000);
  f_lan_w_gausexp->SetParLimits(0, 100, 10000000);
  f_lan_w_gausexp->SetParLimits(1, 30, 400);
  f_lan_w_gausexp->SetParLimits(2, 0.1, 4000);

  TFile *fcalibout;
  TNtuple *tn;

  fcalibout = new TFile(Form("%s/calib/calib_mbd_%d.root", env_p, runnumber),"RECREATE");
  tn = new TNtuple("mbd_calib","mbd_calib","channel:peak:width");

  for (int ich = 0 ; ich < 128; ich++){

    std::cout << "Raw Fitting :: Channel "<<ich;
    
    double local_min = 0;
    double local_max = 0;
    double max_value = 0;
    int good = 0;
    int low_bin = 0;
    double upper;
    int high_bin = 0;
    TH1D *hsmooth = (TH1D*) h_charge_raw[ich]->Clone();
    hsmooth->Smooth();
    hsmooth->Smooth();
    // Finding the dip
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
	    local_min = hsmooth->GetBinLowEdge(ic);
	    low_bin = ic;
	    break;
	  }
      }
    if (low_bin == 0)
      {
	low_bin = start_looking;
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
    
    int tries = 0;
    while (tries < 5)
      {

	

	if (tries == 0)
	  {
	    if (local_max > 600){
	      local_max = 150;
	    }
	    else if (local_max < 40)
	      {
		local_max = 41;
	      }
	    if (local_min > 400) 
	      {
		local_min = 30;
		local_max = 150;
	      }
	    else if (local_min < 20)
	      {
		local_min = 20;
		
	      }

	    f_lan_w_exp->SetParLimits(1, (local_max - 20 < 15?15:local_max - 20), 500);
	    f_lan_w_exp->SetParameter(0, 10000);

	    f_lan_w_exp->SetParameter(1, local_max);
	    f_lan_w_exp->SetParameter(2, 10);
	    f_lan_w_exp->SetParameter(3, 1);
	    f_lan_w_exp->SetParameter(4, -0.001);
	    
	    upper = local_min + (local_max - local_min < 20 ? 200 : 4*(local_max - local_min));
	    hsmooth->Fit("lan_w_exp","Q","",local_min, upper);
	    std::cout << "Chi2/NDF: " << f_lan_w_exp->GetChisquare()/f_lan_w_exp->GetNDF()<<std::endl;
	    h_charge_raw[ich]->Fit("lan_w_exp","Q","",local_min, upper);
	    std::cout << "Chi2/NDF: " << f_lan_w_exp->GetChisquare()/f_lan_w_exp->GetNDF()<<std::endl;
	    f_lan_w_gausexp->SetParameter(0, f_lan_w_exp->GetParameter(0));
	    f_lan_w_gausexp->SetParameter(1, f_lan_w_exp->GetParameter(1));
	    f_lan_w_gausexp->SetParameter(2, f_lan_w_exp->GetParameter(2));
	    f_lan_w_gausexp->SetParameter(3, f_lan_w_exp->GetParameter(3));
	    f_lan_w_gausexp->SetParameter(4, f_lan_w_exp->GetParameter(4));


	    float value = 0;
	    float value_x = 0;
	    int bin_min = hsmooth->GetBin(local_min);//f_lan_w_exp->GetParameter(1));
	    int binpeak = hsmooth->GetBin(f_lan_w_exp->GetParameter(1));
	    int ibin_mark = bin_min - 1;
	    for (int ilower = 1; ilower < bin_min - 1; ilower++)
	      {
	    	ibin_mark = bin_min - ilower;
	    	value = hsmooth->GetBinContent(ibin_mark);
	    	value_x = hsmooth->GetBinCenter(ibin_mark);
	    	if (hsmooth->GetBinContent(ibin_mark) > 1.5*hsmooth->GetBinContent(binpeak))
	    	  {
	    	    break;
	    	  }
	      }

	    
	    float value_x_i = value_x;
	    hsmooth->Fit("gaus", "","",h_charge_raw[ich]->GetMaximumBin(), value_x);
	    f_lan_w_gausexp->SetParameter(5, hsmooth->GetFunction("gaus")->GetParameter(0));
	    f_lan_w_gausexp->SetParameter(6, hsmooth->GetFunction("gaus")->GetParameter(1));
	    f_lan_w_gausexp->SetParameter(7, hsmooth->GetFunction("gaus")->GetParameter(2));
	    
	    int fitss = 0;
	    h_charge_raw[ich]->Fit("lan_w_gausexp","","",value_x, upper);
	    while (f_lan_w_gausexp->GetChisquare()/f_lan_w_gausexp->GetNDF() > 3 && fitss < 5)
	      {
		value_x = 0;
		value = 0;
		for (int ilower = 1; ilower < bin_min - 1; ilower++)
		  {
		    ibin_mark = bin_min - ilower;
		    value = hsmooth->GetBinContent(ibin_mark);
		    value_x = hsmooth->GetBinCenter(ibin_mark);
		    if (hsmooth->GetBinContent(ibin_mark) > (1.3 - 0.1*fitss)*hsmooth->GetBinContent(binpeak))
		      {
			break;
		      }
		  }
		hsmooth->Fit("gaus", "","",value_x_i, value_x);
		f_lan_w_gausexp->SetParameter(5, hsmooth->GetFunction("gaus")->GetParameter(0));
		f_lan_w_gausexp->SetParameter(6, hsmooth->GetFunction("gaus")->GetParameter(1));
		f_lan_w_gausexp->SetParameter(7, hsmooth->GetFunction("gaus")->GetParameter(2));
		
		h_charge_raw[ich]->Fit("lan_w_gausexp","Q","",value_x, upper);
		fitss++;
		std::cout << fitss <<" - chi2/ndf = " << f_lan_w_gausexp->GetChisquare()/f_lan_w_gausexp->GetNDF()<<std::endl;
	      }
	    if (f_lan_w_gausexp->GetChisquare()/f_lan_w_gausexp->GetNDF() < 3)
	      {
		break;
	      }
	    
	  }
	else
	  {

	    local_max = hsmooth->GetBinLowEdge(high_bin);
	    f_lan_w_exp->SetParameter(0, 10000);

	    f_lan_w_exp->SetParameter(1, local_max);
	    f_lan_w_exp->SetParameter(2, 10);
	    f_lan_w_exp->SetParameter(3, 1);
	    f_lan_w_exp->SetParameter(4, 0.001);

	    upper = local_min + (local_max - local_min < 20 ? 200 : 4*(local_max - local_min));

	    hsmooth->Fit("lan_w_exp","Q","",local_min, upper);
	    
	    h_charge_raw[ich]->Fit("lan_w_exp","Q","",local_min, upper);

	    if (f_lan_w_exp->GetChisquare()/f_lan_w_exp->GetNDF() < 5)
	      {
		break;
	      }

	  }
	std::cout << "try "<<tries<<std::endl;
	tries++;
      }
    
    std::cout <<" " << local_min <<"/"<< local_max<<"/"<< upper<<" -- Peak at "<<f_lan_w_exp->GetParameter(1)<<endl;
    gain_corr[ich] = 1./((double) f_lan_w_exp->GetParameter(1));
  }    

  if (makeNewHistograms)
    {
      for (int i = 0; i < t->GetEntries();i++)
	{
	  if (i%10000 == 0) cout << i << endl;
	  t->GetEntry(i);
	  if (zdc_sum[0]< 40 || zdc_sum[1] < 40) continue;
	  for (int ich = 0 ; ich < 128; ich++)
	    {
	      h_charge[ich]->Fill(mbd_charge_raw[ich]*gain_corr[ich]);
	    }
	}
    }    

      
  for (int ich = 0 ; ich < 128; ich++){

    std::cout << "Recalib Fitting :: Channel "<<ich;

    f_lan_w_exp->SetParLimits(0, 100, 50000);
    f_lan_w_exp->SetParLimits(1, 0.2, 1.5);
    f_lan_w_exp->SetParLimits(2, 0.01, 0.8);

    
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

	if (tries == 0)
	  {
	    f_lan_w_exp->SetParLimits(1, (local_max - 0.5 < 0.2?0.2:local_max - 0.5), local_max + 0.5);
	    f_lan_w_exp->SetParameter(0, 10000);
	    f_lan_w_exp->SetParameter(1, local_max);
	    f_lan_w_exp->SetParameter(2, 0.12);
	    f_lan_w_exp->SetParameter(3, 1);
	    f_lan_w_exp->SetParameter(4, 0.001);
	    
	    upper = local_min + (local_max - local_min < 20 ? 200 : 4*(local_max - local_min));
	  }
	else 
	  {
	    f_lan_w_exp->SetParLimits(1, (local_max - 0.2 < 0.2?0.2:local_max - 0.2), local_max + 0.2);
	    f_lan_w_exp->SetParameter(0, 10000);
	    f_lan_w_exp->SetParameter(1, local_max);
	    f_lan_w_exp->SetParameter(2, 0.12);
	    f_lan_w_exp->SetParameter(3, 2);
	    f_lan_w_exp->SetParameter(4, 0.001);
	    
	    upper = local_min + (local_max - local_min < 20 ? 200 : 4*(local_max - local_min));

	  }
	hsmooth->Fit("lan_w_exp","Q","",local_min, upper);
	h_charge[ich]->Fit("lan_w_exp","Q","",local_min, upper);

	std::cout <<" try " <<tries <<endl;
	if (f_lan_w_exp->GetChisquare()/f_lan_w_exp->GetNDF() < 5)
	  {
	    break;
	  }
	tries++;
      }
    h_peaks->Fill(f_lan_w_exp->GetParameter(1));

    gain_corr[ich] = gain_corr[ich]*f_lan_w_exp->GetParameter(1);

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
}

void QA_MBDTimeCalibrations(const int runnumber)
{

  const char *env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }

  TFile *f = new TFile(Form("%s/output/run%d/signals/processedsignal_%d.root", env_p, runnumber, runnumber), "r");
  TTree *t = (TTree*) f->Get("ttree");

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
      std::cout << i <<": "<< g <<" "<<gain_corr[i]<<std::endl;
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

  t->SetBranchAddress("peak", peak);
  t->SetBranchAddress("pedestal", pedestal);
  t->SetBranchAddress("peaktime", peaktime);
  t->SetBranchAddress("time_of_arrival", time_of_arrival);

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
  std::cout << "\rChecking the number of events per channel ..."<<std::flush;
  for (int i = 0; i < t->GetEntries(); i++)
    {

      t2->GetEntry(i);
      if (1 != 0 && i%(t->GetEntries()/100) == 0) std::cout << "\rChecking the number of events per channel ... "<<floor(static_cast<float>(i)/static_cast<float>(t->GetEntries())*100)<<"%"<<std::flush;

      for (int ich = 0; ich < 128; ich++)
	{
	  if (mbd_charge[ich]*gain_corr[ich] < 0.4) continue;
	  number_of_events[ich]++;
	}
    }  

  std::cout << "\rChecking the number of events per channel.. [Done]"<<std::endl;
  std::cout << "\rCalculating Mean.. "<<std::flush;
  for (int i = 0; i < t->GetEntries(); i++)
    {
      t->GetEntry(i);
      t2->GetEntry(i);
      if (i != 0 && i%(t->GetEntries()/100) == 0) std::cout << "\rCalculating Mean ... "<<floor(static_cast<float>(i)/static_cast<float>(t->GetEntries())*100)<<"%"<<std::flush;
      for (int ich = 0; ich < 128; ich++)
	{

	  double toa = static_cast<double>(time_of_arrival[8 + (ich/8)*16 + ich%8]);
	  double tdc = static_cast<double>(25. - (9./5000.)*mbd_time[ich]);

	  double ptime = static_cast<double>(peaktime[(ich/16 + 1)*8 + ich%8]);
	  if (mbd_charge[ich]*gain_corr[ich] < 0.4) continue;

	  center_tdc[ich] += tdc/static_cast<double>(number_of_events[ich]);
	  if (toa > 0)
	    center_toa[ich] += toa/static_cast<double>(number_of_events[ich]);
	  if (ptime > 0)
	    center_peak_time[ich] += ptime/static_cast<double>(number_of_events[ich]);


	}
    }  
  std::cout << "\rCalculating Mean ... [Done]"<<std::endl;
  std::cout << "Mean TOA     : "<< center_toa[8] << std::endl;
  std::cout << "Mean TDC     : "<< center_tdc[8] << std::endl;
  std::cout << "Mean PeakTime: "<< center_peak_time[8] << std::endl;

  std::cout << "\rFilling Histograms.. "<<std::flush;
  for (int i = 0; i < t->GetEntries(); i++)
    {
      t->GetEntry(i);
      t2->GetEntry(i);
      if (i != 0 && i%(t->GetEntries()/100) == 0) std::cout << "\rFilling Histograms ... "<<floor(static_cast<float>(i)/static_cast<float>(t->GetEntries())*100)<<"%"<<std::flush;
      for (int ich = 0; ich < 128; ich++)
	{


	  if (mbd_charge[ich]*gain_corr[ich] < 0.4) continue;
	  double toa = static_cast<double>(time_of_arrival[8 + (ich/8)*16 + ich%8]) - center_toa[ich];
	  double tdc = (25. - static_cast<double>(mbd_time[ich])*(9./5000.)) - center_tdc[ich];
	  double ptime = static_cast<double>(peaktime[(ich/16 + 1)*8 + ich%8]) - center_peak_time[ich];
	  h_toa_channel[ich]->Fill(toa);
	  h_tdc_channel[ich]->Fill(tdc);
	  h_tdc_toa_channel[ich]->Fill(tdc, toa);
	  h_toa_ptime_channel[ich]->Fill(toa, ptime);
	  h_tdc_ptime_channel[ich]->Fill(tdc, ptime);
	  h_ptime_channel[ich]->Fill(ptime);
	  h_peak_ptime_channel[ich]->Fill(peak[8 + (ich/8)*16 + ich%8] - pedestal[8 + (ich/8)*16 + ich%8], ptime);
	  h_peak_toa_channel[ich]->Fill(peak[8 + (ich/8)*16 + ich%8] - pedestal[8 + (ich/8)*16 + ich%8], toa);
	  h_peak_tdc_channel[ich]->Fill(peak[8 + (ich/8)*16 + ich%8] - pedestal[8 + (ich/8)*16 + ich%8], tdc);

	}
    }  
  
  std::cout << "\rFilling Histograms ... [Done]"<<std::endl;
  
  TFile *ftcalib = new TFile(Form("%s/calib/t0_calib_mbd_%d.root", env_p, runnumber),"recreate");

  TNtuple *t_calib = new TNtuple("ttree","time calib","channel:shift:scale:toa_shift");

  std::cout << "\rFitting 1D Histograms.. "<<std::flush;

  double middle_fit_toa[128];
  double middle_fit_tdc[128];

  double scale_fit_toa[128];
  double scale_fit_tdc[128];

  std::cout << "Tube \t TDC shift \t TDC scale \t TOA Shift"<<endl;
  for (int ich = 0; ich < 128; ich++)
    {
      std::cout << "channel : "<<ich <<std::endl;
      TF1 *fit_toa = 0;
      TF1 *fit_tdc = 0;
      TH1D *hclone = nullptr;
      int nparams = 9;
      double params[9];
      std::cout << " "<<std::endl;
      std::cout << "TOA "<<std::endl;
      std::cout << " "<<std::endl;
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

      std::vector<std::pair<double, double>> fit_peak_heights;
      for (int i = 0; i < nparams/3 ; i++)
	{
	  fit_peak_heights.push_back(std::make_pair(fit_toa->GetParameter(i*3 + 1), fit_toa->GetParameter(i*3)));
	}

      std::sort(fit_peak_heights.begin(), fit_peak_heights.end(), [](auto &left, auto &right) {
	  return left.second > right.second;
	});
      std::cout <<" Peak heights: "<<std::flush;
      for (int i = 0; i < nparams/3 ; i++)
	{
	  std::cout << " " <<fit_peak_heights.at(i).first<<","<<fit_peak_heights.at(i).second;
	}
      middle_fit_toa[ich] = center_toa[ich] + fit_peak_heights.at(0).first;

      if (nparams == 6)
	scale_fit_toa[ich] = fabs(fit_peak_heights.at(1).first - fit_peak_heights.at(0).first);

      scale_fit_toa[ich] = fabs(fit_peak_heights.at(1).first - fit_peak_heights.at(0).first);

      if (ich == 54 || ich == 56 || ich == 56+64 ){
	middle_fit_tdc[ich] = 0;
	scale_fit_tdc[ich] = -999;
      }
      else
	{
	  std::cout << " "<<std::endl;
	  std::cout << "TDC "<<std::endl;
	  std::cout << " "<<std::endl;

	  range = getFitFunction(h_tdc_channel[ich], params, -15., 15.);
	  fit_tdc = new TF1("triple_gaussian","gaus(0) + gaus(3) + gaus(6)",range.first, range.second);
      
	  nparams = 9;
	  if (params[6] < 0)
	    {
	      fit_tdc = new TF1("double_gaussian","gaus(0) + gaus(3)",range.first, range.second);
	      nparams = 6;
	    }
	  std::cout << "TDC Guesses: "<<endl;
	  for (int i = 0 ; i < nparams; i++)
	    {
	      std:: cout << params[i] << " ";
	      fit_tdc->SetParameter(i, params[i]);
	      if (i%3 == 2) fit_tdc->SetParLimits(i, 0, 100);
	    }
	  std::cout << " "<<endl;
	  hclone = (TH1D*) h_tdc_channel[ich]->Clone();
	  hclone->Rebin(20);
	  hclone->Smooth();
	  hclone->Fit(fit_tdc->GetName(), "Q","", range.first, range.second);
	  std::cout << "TDC Guesses: "<<endl;
	  for (int i = 0 ; i < nparams; i++)
	    {
	      std:: cout << fit_tdc->GetParameter(i) << " ";
	    }
	  std::cout << " "<<endl;
	  
	  std::vector<std::pair<double, double>> fit_tdc_peak_heights;
	  for (int i = 0; i < nparams/3 ; i++)
	    {
	      fit_tdc_peak_heights.push_back(std::make_pair(fit_tdc->GetParameter(i*3 + 1), fit_tdc->GetParameter(i*3)));
	    }
	  
	  std::sort(fit_tdc_peak_heights.begin(), fit_tdc_peak_heights.end(), [](auto &left, auto &right) {
	      return left.second > right.second;
	    });
	  
	  middle_fit_tdc[ich] = center_tdc[ich] + fit_tdc_peak_heights.at(0).first;
	  scale_fit_tdc[ich] = scale_fit_toa[ich] / (fabs(fit_tdc_peak_heights.at(1).first - fit_tdc_peak_heights.at(0).first));
	}
      std::cout << ich <<" \t "<<middle_fit_tdc[ich] <<" \t "<<scale_fit_tdc[ich]<<" \t "<<middle_fit_toa[ich]<<endl;
      t_calib->Fill(ich, middle_fit_tdc[ich],scale_fit_tdc[ich], middle_fit_toa[ich]);
    }
  std::cout << "\rFitting 1D Histograms.. [done]"<<std::endl;
  ftcalib->Write();
  ftcalib->Close();

  TFile *fout21 = new TFile(Form("%s//output/hists/fitsout_%d.root", env_p, runnumber),"recreate");
  for (int i = 0; i < 128; i++)
    {
      h_tdc_channel[i]->Write();
      h_toa_channel[i]->Write();
      h_ptime_channel[i]->Write();
      h_tdc_toa_channel[i]->Write();
      h_tdc_ptime_channel[i]->Write();
      h_toa_ptime_channel[i]->Write();
      h_peak_ptime_channel[i]->Write();
      h_peak_tdc_channel[i]->Write();
      h_peak_toa_channel[i]->Write();
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

  std::cout << "\rFilling Histograms Again.. "<<std::flush;
  for (int i = 0; i < t->GetEntries(); i++)
    {
      t->GetEntry(i);
      t2->GetEntry(i);
      if (i != 0 && i%(t->GetEntries()/100) == 0) std::cout << "\rFilling Histograms ... "<<floor(static_cast<float>(i)/static_cast<float>(t->GetEntries())*100)<<"%"<<std::flush;
      for (int ich = 0; ich < 128; ich++)
	{


	  if (mbd_charge[ich]*gain_corr[ich] < 0.4) continue;
	  double toa = static_cast<double>(time_of_arrival[8 + (ich/8)*16 + ich%8]) - middle_fit_toa[ich];
	  double tdc = scale_fit_tdc[ich]*((25. - static_cast<double>(mbd_time[ich])*(9./5000.)) - middle_fit_tdc[ich]);
	  double ptime = static_cast<double>(peaktime[(ich/16 + 1)*8 + ich%8]) - center_peak_time[ich];

	  h_toa_channel[ich]->Fill(toa);
	  h_tdc_channel[ich]->Fill(tdc);
	  h_tdc_toa_channel[ich]->Fill(tdc, toa);
	  h_toa_ptime_channel[ich]->Fill(toa, ptime);
	  h_tdc_ptime_channel[ich]->Fill(tdc, ptime);
	  h_peak_toa_channel[ich]->Fill(peak[8 + (ich/8)*16 + ich%8] - pedestal[8 + (ich/8)*16 + ich%8], toa);
	  h_peak_tdc_channel[ich]->Fill(peak[8 + (ich/8)*16 + ich%8] - pedestal[8 + (ich/8)*16 + ich%8], tdc);
	}
    }  
  std::cout << "\rFilling Histograms ... [Done]"<<std::endl;


  TFile *fout2 = new TFile(Form("%s/output/hists/histout_%d.root", env_p, runnumber),"recreate");
  for (int i = 0; i < 128; i++)
    {
      h_tdc_channel[i]->Write();
      h_toa_channel[i]->Write();
      h_ptime_channel[i]->Write();
      h_tdc_toa_channel[i]->Write();
      h_tdc_ptime_channel[i]->Write();
      h_toa_ptime_channel[i]->Write();
      h_peak_ptime_channel[i]->Write();
      h_peak_tdc_channel[i]->Write();
      h_peak_toa_channel[i]->Write();

    }
  fout2->Close();
}

std::pair<double, double> getFitFunction(TH1D* h, double *params, float minimumx = -25., float maximumx = 25.){

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
  auto maxv = v_max.begin();
  while (maxv != v_max.end())
    {
      if ((*maxv).second < 1000)
	v_max.erase(maxv);
      else if ((*maxv).first < minimumx || (*maxv).first > maximumx)
	v_max.erase(maxv);
      else
	maxv++;
    }

  if (!v_max.size())
    {
      std::cout << "nothing left" << std::endl;
      return std::make_pair(-20, 20);
    }
  else if (v_max.size() == 1)
    {
      x_high = v_minx.at(0);

      while (v_minx.end() != v_minx.begin())
	{
	  v_minx.erase(v_minx.end());
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

  for (unsigned int i = 0; i < v_max.size(); i++)
    {
      std::cout << "Max: "<<v_max.at(i).first<<"/"<<v_max.at(i).second<<std::endl;
    }
  for (unsigned int i = 0; i < v_minx.size(); i++)
    {
      std::cout << "Min: "<<v_minx.at(i)<<std::endl;
    }
  std::cout << x_low<< "-->"<<x_high<<std::endl;

  if (v_max.size() == 3)
    {
      std::cout << "Three"<<endl;
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
      std::cout << "Two"<<endl;
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
