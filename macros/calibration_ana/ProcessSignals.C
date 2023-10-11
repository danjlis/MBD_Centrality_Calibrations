/* ProcessSignals */
// Author: Daniel Lis - September 27 - 2023
// Brief : 

int debug = 0;

float FastHalfMax(float ped, int maxx, float maxy, int * wave , const int nsamples);
float FastSteepMax(float ped, int * wave , const int nsamples);
void Peak_Spline(const int runnumber, const int segment);
std::pair<float, float> FastMax(int maxx, float * wave , const int nsamples , float *xmax, float *ymax);

void ProcessSignals(
		    const int runnumber,
		    const int segment
		    )
{

  Peak_Spline(runnumber, segment);
  return;
}

// returns time intercept with the x axis by highest and steepest line
float FastHalfMax(int maxx, float maxy, float ped, int * wave , const int nsamples)
{

  int half_x = -1;
  float fraction_threshold = 0.5;
  float threshold = fraction_threshold*(maxy - ped);;
  if (debug)
    {
      std::cout << "maxx/maxy/ped/thresh = "<<maxx<<"/"<<maxy<<"/"<<ped<<"/"<<threshold<<std::endl;
    }
  for (int i = 0 ; i < maxx; i++)
    {
      if ((static_cast<float>(wave[i]) - ped) > threshold)
	{
	  if (i == (maxx - 1) || static_cast<float>(wave[i+1] - ped) > threshold)
	    {
	      half_x = i;
	      if (debug)
		{
		  std::cout << "half_x = "<<half_x<<std::endl;
		}
	      
	      break; 
	    }
	}
    }

  if (half_x == -1)
    {
      return -999;
    }

  float dx = 17.7623*(static_cast<float>(half_x) - static_cast<float>(half_x - 1));
  float dy = static_cast<float>(wave[half_x]) - static_cast<float>(wave[half_x - 1]);
  float dt1 = static_cast<float>(wave[half_x]) - ped - threshold;

  float dt0 = 17.7623*(static_cast<float>(half_x)) - dt1*(dx/dy);

  if (debug)
    {
      std::cout <<"dx: "<<dx<<std::endl;
      std::cout <<"dy: "<<dy<<std::endl;
      std::cout <<"dt1: "<<dt1<<std::endl;
      std::cout <<"dt0: "<<dt0<<std::endl;

    }

  return dt0;
}

// returns time intercept with the x axis by highest and steepest line
float FastSteepMax(float ped, int * wave , const int nsamples)
{

  std::pair<float, float> two_samples;
  std::pair<float, float> two_brothers;
  float slope = 0.;
  float intercept = 0.;
  float interceptx = 0.;
  float maxslope = 0;
  float maxinterceptx = -999.99;
  float maxsample = 0.;

  for (int is = 0; is < nsamples - 1; is++)
    {
      two_samples = std::make_pair(static_cast<float>(is), static_cast<float>(is+1));
      two_brothers = std::make_pair(static_cast<float>(wave[is]), static_cast<float>(wave[is+1]));
      
      slope = two_brothers.second - two_brothers.first;
      intercept = two_brothers.first - slope*two_samples.first;

      interceptx = (ped - intercept)/slope;
      if (slope > maxslope)
	{
	  maxslope = slope;
	  maxinterceptx = interceptx;
	  maxsample = static_cast<float>(is);
	}
    }

  if (maxsample == 0 || maxsample == nsamples - 2) return -999.99;
  if (maxslope == 0) return -999.99;


  if (debug)
    {
      std::cout << "slope: "<< maxslope<<std::endl;
      
      std::cout << "maxsample: "<< maxsample<<std::endl;

    }
  return (maxinterceptx - maxsample)*17.7623;
}

std::pair<float, float> FastMax(int maxx, int * wave , const int nsamples)
{
  float xmax, ymax;
  int n = 3;
  double maxxd = static_cast<double>(maxx);
  double xp[3] = {maxxd - 1, maxxd, maxxd+1};
  double yp[3] = {static_cast<double>(wave[maxx-1]), static_cast<double>(wave[maxx]), static_cast<double>(wave[maxx+1])};
  TSpline3 *sp = new TSpline3("", xp, yp, n, "b2e2", 0, 0);
  double X, Y, B, C, D;
  ymax = wave[maxx];
  xmax = maxx;

  for (int i = 0; i <= 1; i++)
    {
      sp->GetCoeff(i, X, Y, B, C, D);
      if (D == 0)
	{
	  if (C < 0)
	    {
	      // TSpline is a quadratic equation

	      float root = -B / (2 * C) + X;
	      if (root >= xp[i] && root <= xp[i + 1])
		{
		  float yvalue = sp->Eval(root);
		  if (yvalue > ymax)
		    {
		      ymax = yvalue;
		      xmax = root;
		    }
		}
	    }
	}
      else
	{
	  // find x when derivative = 0
	  float root = (-2 * C + sqrt(4 * C * C - 12 * B * D)) / (6 * D) + X;
	  if (root >= xp[i] && root <= xp[i + 1])
	    {
	      float yvalue = sp->Eval(root);
	      if (yvalue > ymax)
		{
		  ymax = yvalue;
		  xmax = root;
		}
	    }
	  root = (-2 * C - sqrt(4 * C * C - 12 * B * D)) / (6 * D) + X;
	  if (root >= xp[i] && root <= xp[i + 1])
	    {
	      float yvalue = sp->Eval(root);
	      if (yvalue > ymax)
		{
		  ymax = yvalue;
		  xmax = root;
		}
	    }
	}
    }


  delete sp;
  return std::make_pair(xmax, ymax);;  

}

void Peak_Spline(const int runnumber, const int segment)
{

  
  std::ostringstream rstr;
  rstr << std::setw(8) << std::setfill('0') << runnumber;


  std::ostringstream sstr;
  sstr << std::setw(4) << std::setfill('0') << segment;

  TH1D *h_time_of_arrival[256];
  TH1D *h_peaktime[256];

  for (int i = 0; i < 256;i++)
    {
      int j = 0;

      h_time_of_arrival[i] = new TH1D(Form("h_time_of_arrival_ch_%d_%d_sigma", i, j), "", 2000, -1, 1);
      h_peaktime[i] = new TH1D(Form("h_peaktime_ch_%d_%d_sigma", i, j), "", 310, 0, 31);
    }

  TFile *file = new TFile(Form("/sphenix/user/dlis/Projects/centrality/output/run%d/waves/waveform_tree_%s.root", runnumber, rstr.str().c_str()), "r");
  if (segment != -1) file = new TFile(Form("/sphenix/user/dlis/Projects/centrality/output/run%d/waves/waveform_tree_%s_%s.root", runnumber, rstr.str().c_str(), sstr.str().c_str()) , "r");
  
  if (!file)
    {
      std::cout << " No Tree File found " <<std::endl;
      return;
    }

  TTree *waves = (TTree*) file->Get("ttree");
  if (!waves)
    {
      std::cout << " No Tree in File found " <<std::endl;
      return;
    }

  float time_of_arrival[256];
  float peaktime[256];
  float peak[256];
  float pedestal[256];

  TFile *fouttree = new TFile(Form("/sphenix/user/dlis/Projects/centrality/output/run%d/signals/signal_%s_%s.root", runnumber, rstr.str().c_str(), sstr.str().c_str()), "RECREATE");
  if (segment == -1) fouttree = new TFile(Form("/sphenix/user/dlis/Projects/centrality/macros/waveform_tree_%s.root", rstr.str().c_str()), "RECREATE");
  if (debug) fouttree = new TFile(Form("/sphenix/user/dlis/Projects/centrality/macros/debug_waveform_tree_%s.root", rstr.str().c_str()), "RECREATE");
  TTree *tout = new TTree("ttree"," a persevering date tree");

  tout->Branch("time_of_arrival", time_of_arrival, "time_of_arrival[256]/F");
  tout->Branch("peaktime", peaktime, "peaktime[256]/F");
  tout->Branch("peak", peak, "peak[256]/F");
  tout->Branch("pedestal", pedestal, "pedestal[256]/F");

  const int nsamples = 31;
  const int nchannels = 256;

  const int starting_channel = 0;
  int waveform_mbd[256][31];

  waves->SetBranchAddress("waveform_bbc", waveform_mbd);

  int nentries = waves->GetEntries();
  if (debug) nentries = 10;
  if (debug) std::cout<< "Starting loop over " << nentries <<" entries..." << std::endl;

  for (int i = 0 ; i < nentries; i++)
    {
      waves->GetEntry(i);

      if (debug) std::cout <<" Event " << i <<": "<<std::endl;

      for (int ich = starting_channel; ich < starting_channel + nchannels; ich++)
	{
	  float time = -999.99;
	  float amp = -999.99;
	  float max = -999.99;
	  int maxx = -1;
	  float ped = (waveform_mbd[ich][0] + waveform_mbd[ich][1] + waveform_mbd[ich][2])/3.;
	  if (debug) std::cout <<" Channel " << ich <<": "<<std::endl;
	  if (debug > 3) std::cout <<" Samples: ";

	  for (int is = 0; is < nsamples; is ++)
	    {
	      if (debug > 3) std::cout << waveform_mbd[ich][is] << " ";

	      if (waveform_mbd[ich][is] > max)
		{
		  max = waveform_mbd[ich][is];
		  maxx = is;
		}
	    }

	  if (debug) std::cout << " " << std::endl;

	  if (maxx == 0 || maxx == nsamples - 1)
	    {
	      amp = max;
	      time = maxx;
	    }
	  
	  std::pair<float, float> s;

	  s = FastMax(maxx, &waveform_mbd[ich][0], nsamples);

	  float eta = FastHalfMax(maxx, max, ped, &waveform_mbd[ich][0], nsamples);

	  amp = s.second;
	  time = s.first;



	  if (debug)	  std::cout << " Amplitude at Time: " << amp << " / "<< time << std::endl;
	  if (debug)	  std::cout << " ETA: " << eta <<std::endl;
	  h_time_of_arrival[ich]->Fill(eta);
	  h_peaktime[ich]->Fill(time);
	  peaktime[ich] = time * 17.7623;
	  time_of_arrival[ich] = eta;
	  peak[ich] = amp;
	  pedestal[ich] = ped;
	}
      tout->Fill();
    }
  fouttree->Write();
  fouttree->Close();
  delete fouttree;
  return;
}
