void DrawMBDVertex(int runnumber)
{
  char *env_p = new char[200];
  sprintf(env_p,"%s",std::getenv("MBD_CENTRALITY_CALIB_PATH"));


  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }

  TFile *f = new TFile(Form("%s/output/run%d/mbdana/mbd_trees_%d.root", env_p, runnumber, runnumber), "r");
  TTree *t = (TTree*) f->Get("T");

  TH1D *h_temptime_n = new TH1D("h_temptime_n", ";t_{n} [cm]; Counts", 400, -20, 20);
  TH1D *h_temptime_s = new TH1D("h_temptime_s", ";t_{n} [cm]; Counts", 400, -20, 20);
  TH1D *h_time_n = new TH1D("h_time_n", ";t_{n} [cm]; Counts", 400, -20, 20);
  TH1D *h_time_s = new TH1D("h_time_s", ";t_{s} [cm]; Counts", 400, -20, 20);
  TH1D *h_vertex = new TH1D("h_vertex", ";vertex [cm]; Counts", 601, -300.5, 300.5);
  TH1D *h_vertex_fid = new TH1D("h_vertex_fid", ";vertex [cm]; Counts", 601, -300.5, 300.5);
  TH1D *h_vertex_other = new TH1D("h_vertex_other", ";vertex [cm]; Counts", 601, -300.5, 300.5);
  TH1D *h_vertex_zdc = new TH1D("h_vertex_zdc", ";vertex [cm]; Counts", 601, -300.5, 300.5);
  TH1D *h_vertex_hits[64];
  for (int i = 0; i < 64; i++) h_vertex_hits[i] = new TH1D(Form("h_vertex_%d", i), ";vertex [cm]; Counts", 601, -300.5, 300.5);

  TH2D *h_charge_2 = new TH2D("h_charge_2","",50, 0, 1000, 50, 0, 1000);

  TProfile2D *hp_vertex = new TProfile2D("hp_vertex","", 64, 0.5, 64.5, 64, 0.5, 64.5,"s");

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

  for (int i = 0; i < t->GetEntries(); i++)
    {
      t->GetEntry(i);

      h_vertex->Fill(mbd_vertex);
      int nhit_n = 0;
      int nhit_s = 0;
      h_temptime_n->Reset();
      h_temptime_s->Reset();
      for (int j = 0; j < 64; j++)
	{
	  if (mbd_charge[j] > 0.4) 
	    {
	      nhit_s++;
	      h_temptime_s->Fill(mbd_time[j]);
	    }
	  if (mbd_charge[j+64] > 0.4) 
	    {
	      nhit_n++;
	      h_temptime_n->Fill(mbd_time[j+64]);
	    }
	}
      if (nhit_n < 2 || nhit_s < 2) continue; 
      h_vertex_hits[min(nhit_n, nhit_s) - 2]->Fill(mbd_vertex);
      hp_vertex->Fill(nhit_s, nhit_n, mbd_vertex);
      float mean_n = h_temptime_n->GetMean();
      float mean_s = h_temptime_s->GetMean();
      if (nhit_n > 50 && nhit_s > 50)
      {
	h_time_n->Fill(mean_n);
	h_time_s->Fill(h_temptime_s->GetMean());
      }
      if (mean_n < -9 || mean_s < -9 || mean_n > 7 || mean_s > 7)
	{
	  h_vertex_other->Fill(mbd_vertex);
	}
      else
	{
	  h_vertex_fid->Fill(mbd_vertex);
	}
      
      if (fabs(mbd_vertex) < 30)
	{
	  h_charge_2->Fill(mbd_sum[0], mbd_sum[1]);
	}
      if (zdc_sum[0] > 40 && zdc_sum[1] > 40) h_vertex_zdc->Fill(mbd_vertex);
    }


  TH1D *h_mean = new TH1D("h_mean","", 65, -0.5, 64.5);
  TH1D *h_err = new TH1D("h_err","", 65, -0.5, 64.5);
  for (int i = 0; i < 63; i++)
    {
      h_vertex_hits[i]->Fit("gaus","","",-30, 30);
      double mean = h_vertex_hits[i]->GetFunction("gaus")->GetParameter(1);
      double err = h_vertex_hits[i]->GetFunction("gaus")->GetParameter(2);
      h_mean->Fill(i+2, mean); 
      h_err->Fill(i+2, err); 
    }
}
