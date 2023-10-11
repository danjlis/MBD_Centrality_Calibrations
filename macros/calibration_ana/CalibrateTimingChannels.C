#include <dlUtility.h>

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

void justFit(const int runnumber)
{


  TFile *fout2 = new TFile(Form("../output/hists/histout_%d.root", runnumber),"r");
  TH1D *h_toa_channel[128];
  TH1D *h_tdc_channel[128];

  for (int i = 0 ; i < 128; i++)
    {
      h_toa_channel[i] = (TH1D*) fout2->Get(Form("h_toa_channel_%d", i));
      h_tdc_channel[i] = (TH1D*) fout2->Get(Form("h_tdc_channel_%d", i));
    }

  for (int ich = 0; ich < 128; ich++)
    {
      std::cout << " " <<endl;
      std::cout << "Channel " << ich <<endl;
      std::cout << " TOA:" <<endl;
      TF1 *fit_toa = 0;
      double params[9];
      fit_toa = new TF1("triple_gaussian","gaus(0) + gaus(3) + gaus(6)",-20, 20);
      std::pair<double, double> range = getFitFunction(h_toa_channel[ich], params);
      fit_toa = new TF1("triple_gaussian","gaus(0) + gaus(3) + gaus(6)",range.first, range.second);
      int nparams = 9;
      if (params[6] < 0)
	{
	  fit_toa = new TF1("double_gaussian","gaus(0) + gaus(3)",range.first, range.second);
	  nparams = 6;
	}
      for (int i = 0 ; i < nparams; i++)
	{
	  fit_toa->SetParameter(i, params[i]);
	}

      h_toa_channel[ich]->Rebin(20);
      h_toa_channel[ich]->Fit(fit_toa->GetName(), "","", range.first, range.second);

      std::cout << " " <<endl;
      std::cout << " TDC:" <<endl;

      if (ich == 54 || ich == 56 || ich == 120 ) continue;
      TF1 *fit_tdc = 0;
      

      range = getFitFunction(h_tdc_channel[ich], params, -20., 17.);
      fit_tdc = new TF1("triple_gaussian","gaus(0) + gaus(3) + gaus(6)",range.first, range.second);
      
      nparams = 9;
      if (params[6] < 0)
	{
	  fit_tdc = new TF1("double_gaussian","gaus(0) + gaus(3)",range.first, range.second);
	  nparams = 6;
	}
      for (int i = 0 ; i < nparams; i++)
	{
	  fit_tdc->SetParameter(i, params[i]);
	}
      fit_tdc->SetRange(range.first, range.second);
      
      fit_tdc->SetRange(range.first, range.second);
      if (!fit_tdc)
	{
	  std::cout << "No function..."<<endl;
	  continue;
	}
      h_tdc_channel[ich]->Rebin(20);

      h_tdc_channel[ich]->Fit(fit_tdc->GetName(), "","", range.first, range.second);

    }
  TFile *fout3 = new TFile(Form("../output/hists/fitsout_%d.root", runnumber),"recreate");

  for (int i = 0 ; i < 128; i++)
    {
      h_toa_channel[i]->Write();
      h_tdc_channel[i]->Write();
    }
  fout3->Close();


  return;
}

void drawTimeOfArrival(const int runnumber)
{

  TFile *f = new TFile(Form("/sphenix/user/dlis/Projects/signal_processing/output/run%d/processedsignal_%d.root", runnumber, runnumber), "r");
  TTree *t = (TTree*) f->Get("ttree");

  TFile *f2 = new TFile(Form("/sphenix/user/dlis/Projects/centrality/output/run%d/mbd_calib_trees_%d.root", runnumber, runnumber), "r");
  TTree *t2 = (TTree*) f2->Get("T");

  TFile *fcalib = new TFile(Form("/sphenix/user/dlis/Projects/centrality/calib/calib_mbd_%d.root", runnumber), "r");
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
  
  TFile *ftcalib = new TFile(Form("/sphenix/user/dlis/Projects/centrality/calib/t0_calib_mbd_%d.root", runnumber),"recreate");

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

  TFile *fout21 = new TFile(Form("../output/hists/fitsout_%d.root", runnumber),"recreate");
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


  TFile *fout2 = new TFile(Form("../output/hists/histout_%d.root", runnumber),"recreate");
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
