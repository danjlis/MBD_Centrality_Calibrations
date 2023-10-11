bool DEBUG = false;

const double centrality_divs[20] = {1999, 1499, 1291, 1102, 937, 790, 660, 547, 449, 363, 289, 227, 174, 130, 94, 66, 45, 0, 0, 0};

void DrawMBDTime(int runnumber = 21813)
{
  gStyle->SetOptStat(0);
  float thresh = 0.4;
  float timethresh = 15;
  float timethreshs = 12;
  float sigma_cut = 1.5;
  float central_cut = 32;
  float vshift = 28.53;
  TFile *file = new TFile(Form("../output/run%d/trees_%d.root", runnumber, runnumber), "r");

  // 3 vertex cuts with centralities between them compared.
  
  TH1D *h_time_n = new TH1D("h_time_n","", 100, -25, 25);
  TH1D *h_time_s = new TH1D("h_time_s","", 100, -25, 25);
  TH1D *h_tubes_n = new TH1D("h_tubes_n","", 65, -0.5, 64.5);
  TH1D *h_tubes_s = new TH1D("h_tubes_s","", 65, -0.5, 64.5);
  TProfile *hp_rms_by_hits_n = new TProfile("hp_rms_by_hits_n", "", 63, .5, 64.5);
  TProfile *hp_rms_by_hits_s = new TProfile("hp_rms_by_hits_s", "", 63, .5, 64.5);

  TH1D *distance_from_mean[128];
  for (int i = 0; i < 128; i++)
    {
      distance_from_mean[i] = new TH1D(Form("distance_from_mean_%d", i), "", 100, 0, 25);
    }

  // 0 - mean
  // 1 - median
  TH1D *h_time_dist[128];
  for (int i = 0; i < 128; i++)
    {
      h_time_dist[i] = new TH1D(Form("h_time_dist_ch_%d", i), "", 1800, 30, 15);
    }
  
  TProfile *vertex_by_mean = new TProfile("vertex_by_mean","", 65, -0.5, 64.5,"s");
  TProfile *vertex_by_median = new TProfile("vertex_by_median","", 65, -0.5, 64.5,"s");
  TProfile *vertex_by_mean_cover = new TProfile("vertex_by_mean_cover","", 65, -0.5, 64.5,"s");
  TProfile *vertex_by_median_cover = new TProfile("vertex_by_median_cover","", 65, -0.5, 64.5,"s");

  TProfile2D *vertex_by_mean2 = new TProfile2D("vertex_by_mean2","", 65, -0.5, 64.5, 65, -0.5, 64.5);
  TProfile2D *vertex_by_median2 = new TProfile2D("vertex_by_median2","", 65, -0.5, 64.5, 65, -0.5, 64.5);
  TProfile2D *vertex_by_mean_cover2 = new TProfile2D("vertex_by_mean_cover2","", 65, -0.5, 64.5, 65, -0.5, 64.5);
  TProfile2D *vertex_by_median_cover2 = new TProfile2D("vertex_by_median_cover2","", 65, -0.5, 64.5, 65, -0.5, 64.5);

  TProfile *diff_vertex_methods = new TProfile("diff_vertex_methods","", 65, -0.5, 64.5);
  TProfile *diff_vertex_methods_cover = new TProfile("diff_vertex_methods_cover","", 65, -0.5, 64.5);

  TH1D *h_vertex_mean = new TH1D("h_vertex_mean", "", 501, -250, 250);
  TH1D *h_vertex_mean_center = new TH1D("h_vertex_mean_center", "", 501, -250, 250);
  TH1D *h_vertex_median = new TH1D("h_vertex_median", "", 501, -250, 250);
  TH1D *h_vertex_mean_cover = new TH1D("h_vertex_mean_cover", "", 501, -250, 250);
  TH1D *h_vertex_median_cover = new TH1D("h_vertex_median_cover", "", 501, -250, 250);

  TH1D *h_time_mean = new TH1D("h_time_mean", "", 501, -50, 50);
  TH1D *h_time_mean_center = new TH1D("h_time_mean_center", "", 501, -50, 50);
  TH1D *h_time_median = new TH1D("h_time_median", "", 501, -50, 50);
  TH1D *h_time_mean_cover = new TH1D("h_time_mean_cover", "", 501, -50, 50);
  TH1D *h_time_median_cover = new TH1D("h_time_median_cover", "", 501, -50, 50);

  TH2D *h_vtx_time_mean = new TH2D("h_vtx_time_mean", "", 51, -250, 250, 70, -10, 25);
  TH2D *h_vtx_time_mean_center = new TH2D("h_vtx_time_mean_center", "", 51, -250, 250, 70, -10, 25);
  TH2D *h_vtx_time_median = new TH2D("h_vtx_time_median", "", 51, -250, 250, 10, -10, 25);
  TH2D *h_vtx_time_mean_cover = new TH2D("h_vtx_time_mean_cover", "", 51, -250, 250, 10, -10, 25);
  TH2D *h_vtx_time_median_cover = new TH2D("h_vtx_time_median_cover", "", 51, -250, 250, 10, -10, 25);

  TH1D *h_vertex_mean_centrality[20];
  TH1D *h_vertex_median_centrality[20];
  TH1D *h_vertex_mean_cover_centrality[20];
  TH1D *h_vertex_median_cover_centrality[20];
  

  TH2D *nhit_dist = new TH2D("nhit_dist","", 65, -0.5, 64.5, 65, -0.5, 64.5);
  TH2D *nhit_dist_cover = new TH2D("nhit_dist_cover","", 65, -0.5, 64.5, 65, -0.5, 64.5);

  // charge time correlation, should be uncorrelated... maybe higher is higher??
  TH2D *h2_charge_time[128];
  for (int i = 0; i < 128; i++)
    {
      h2_charge_time[i] = new TH2D(Form("h2_charge_time_ch_%d", i), "", 1600, 0, 16000, 1600, 0, 16000);
    }

  // occupancy at different centrality.
  TEfficiency *he_occupancy_centrality[20];

  // occupancy at vertex < 30 and >30.
  TEfficiency *he_occupancy_centrality_vtx_gt_30[20];
  TEfficiency *he_occupancy_centrality_vtx_lt_30[20];
  TEfficiency *he_occupancy_centrality_vtx_cover_gt_30[20];
  TEfficiency *he_occupancy_centrality_vtx_cover_lt_30[20];
  TEfficiency *he_occupancy_centrality_vtx_center_gt_30[20];
  TEfficiency *he_occupancy_centrality_vtx_center_lt_30[20];

  for (int i = 0; i < 20; i++)
    {

      h_vertex_mean_centrality[i] = new TH1D(Form("h_vertex_mean_centrality_%i", i), "", 501, -250, 250);
      h_vertex_median_centrality[i] = new TH1D(Form("h_vertex_median_centrality_%i", i), "", 501, -250, 250);
      h_vertex_mean_cover_centrality[i] = new TH1D(Form("h_vertex_mean_cover_centrality_%i", i), "", 501, -250, 250);
      h_vertex_median_cover_centrality[i] = new TH1D(Form("h_vertex_median_cover_centrality_%i", i), "", 501, -250, 250);

  
      he_occupancy_centrality[i] = new TEfficiency(Form("he_occupancy_centrality_%d", i), "", 128, -0.5, 127.5);
      he_occupancy_centrality_vtx_lt_30[i] = new TEfficiency(Form("he_occupancy_centrality_%d_vtx_lt_30", i), "", 128, -0.5, 127.5);
      he_occupancy_centrality_vtx_gt_30[i] = new TEfficiency(Form("he_occupancy_centrality_%d_vtx_gt_30", i), "", 128, -0.5, 127.5);
      he_occupancy_centrality_vtx_cover_lt_30[i] = new TEfficiency(Form("he_occupancy_centrality_%d_vtx_cover_lt_30", i), "", 128, -0.5, 127.5);
      he_occupancy_centrality_vtx_cover_gt_30[i] = new TEfficiency(Form("he_occupancy_centrality_%d_vtx_cover_gt_30", i), "", 128, -0.5, 127.5);
      he_occupancy_centrality_vtx_center_lt_30[i] = new TEfficiency(Form("he_occupancy_centrality_%d_vtx_center_lt_30", i), "", 128, -0.5, 127.5);
      he_occupancy_centrality_vtx_center_gt_30[i] = new TEfficiency(Form("he_occupancy_centrality_%d_vtx_center_gt_30", i), "", 128, -0.5, 127.5);
    }

  // distributions below 50% central
  // distributions above 50% central
  TH1D *h_time_above_50[128];
  TH1D *h_time_below_50[128];
  TH1D *h_charge_above_50[128];
  TH1D *h_charge_below_50[128];

  for (int i = 0; i < 128; i++)
    {
      h_time_above_50[i] = new TH1D(Form("h_time_above_50_ch_%d", i), "", 1600, 0, 16000);
      h_time_below_50[i] = new TH1D(Form("h_time_below_50_ch_%d", i), "", 1600, 0, 16000);
      h_charge_above_50[i] = new TH1D(Form("h_charge_above_50_ch_%d", i), "", 1600, 0, 16000);
      h_charge_below_50[i] = new TH1D(Form("h_charge_below_50_ch_%d", i), "", 1600, 0, 16000);
    }


  float mbd_charge[128];
  float mbd_time[128];
  float mbd_charge_raw[128];
  float mbd_time_raw[128];
  int isMinBias;
  float z_vertex;
  
  TTree *t = (TTree*)file->Get("T");
  
  t->SetBranchAddress("mbd_charge",mbd_charge);
  t->SetBranchAddress("mbd_time",mbd_time);
  t->SetBranchAddress("mbd_charge_raw",mbd_charge_raw);
  t->SetBranchAddress("mbd_time_raw",mbd_time_raw);
  t->SetBranchAddress("z_vertex",&z_vertex);

  int occupancy[128];
  
  for (int i = 0; i < (DEBUG ? 10: t->GetEntries()); i++)
    {
      t->GetEntry(i);
      h_time_n->Reset();
      h_time_s->Reset();
      int hits_n = 0;
      int hits_s_cover = 0;
      int hits_s = 0;      
      std::vector<float> time_sum_n;
      std::vector<float> time_sum_s_cover;
      std::vector<float> time_sum_s;
      float sum_n = 0.;
      float sum_n2 = 0.;
      float sum_s_cover = 0.;
      float sum_s = 0.;
      float sum_s2 = 0.;

      float charge_sum = 0.;

      for (int ich = 0 ; ich < 128; ich++) occupancy[ich] = 0;

      for (int ich = 0 ; ich < 64; ich++)
	{
	  
	  if (mbd_charge[ich] > thresh)
	    {
	      float time = mbd_time[ich] - 12.5;
	      hits_s++;
	      occupancy[ich] = 1;
	      h_time_s->Fill(time);
	      time_sum_s.push_back(time);
	      sum_s += time;	      
	      sum_s2 += (time*time);	      
	      charge_sum += mbd_charge[ich];
	      h2_charge_time[ich]->Fill(mbd_charge_raw[ich], mbd_time_raw[ich]);
	      h_time_dist[ich]->Fill(time);
	      if (!(ich == 56 || ich == 54))
		{
		  time_sum_s_cover.push_back(time);
		  sum_s_cover += time;
		  hits_s_cover++;
		}
	    }      
	  if (mbd_charge[ich + 64] > thresh)
	    {
	      hits_n++;
	      occupancy[ich+64] = 1;
	      float time = mbd_time[ich + 64] - 12.5;
	      h_time_n->Fill(time);
	      time_sum_n.push_back(time);
	      sum_n += time;
	      sum_n2 += time*time;
	      h_time_dist[ich+64]->Fill(time);
	      charge_sum += mbd_charge[ich + 64];
	      h2_charge_time[ich+64]->Fill(mbd_charge_raw[ich+64], mbd_time_raw[ich+64]);
	    }

      	  
	} 

      
      int centrality_bin;
      for (centrality_bin = 0; centrality_bin < 19; centrality_bin++)
	{
	  if (charge_sum > centrality_divs[centrality_bin + 1]) break;
	}
      if (DEBUG) 
	{
	  std::cout << "Event " << i <<": "<< std::endl; 
	  std::cout << "ChargeSum: "<< charge_sum << std::endl;
	  std::cout << "nHits: "<< hits_s <<" / " << hits_n << std::endl;
	  std::cout << "Centrality: "<< centrality_bin << std::endl;
	}

      for (int ich = 0; ich < 128; ich++)
	{
	  he_occupancy_centrality[centrality_bin]->Fill(occupancy[ich], ich);
	  if (centrality_bin > 10) 
	    {
	      h_time_above_50[ich]->Fill(mbd_time_raw[ich]);
	      h_charge_above_50[ich]->Fill(mbd_charge_raw[ich]);
	    }
	  else
	    {
	      h_time_below_50[ich]->Fill(mbd_time_raw[ich]);
	      h_charge_below_50[ich]->Fill(mbd_charge_raw[ich]);
	    }
	}

      nhit_dist->Fill(hits_s, hits_n);
      h_tubes_s->Fill(hits_s);
      h_tubes_n->Fill(hits_n);


      sort(time_sum_n.begin(), time_sum_n.end());
      sort(time_sum_s_cover.begin(), time_sum_s_cover.end());
      sort(time_sum_s.begin(), time_sum_s.end());


      if (hits_s >= central_cut && hits_n >=central_cut){

	float mean_north = sum_n/static_cast<float>(hits_n);
	float mean_south = sum_s/static_cast<float>(hits_s);

	float rms_n = sqrt(sum_n2/static_cast<float>(hits_n) - TMath::Power(mean_north, 2));
	float rms_s = sqrt(sum_s2/static_cast<float>(hits_s) - TMath::Power(mean_south, 2));
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

	float vertex_mean_center = 15.*(mean_north_center - mean_south_center) + vshift;
	float time_mean_center = (mean_north_center + mean_south_center)/2.;

	h_vertex_mean_center->Fill(vertex_mean_center);
	h_time_mean_center->Fill(time_mean_center);
	h_vtx_time_mean_center->Fill(vertex_mean_center, time_mean_center);

	for (int ich = 0; ich < 128; ich++)
	  {
	    if (fabs(vertex_mean_center) < 30) he_occupancy_centrality_vtx_center_lt_30[centrality_bin]->Fill(occupancy[ich], ich);
	    else he_occupancy_centrality_vtx_center_gt_30[centrality_bin]->Fill(occupancy[ich], ich);
	  }

      }

      if (hits_s >=2 && hits_n >=2){
	float median_n;
	if (time_sum_n.size()%2)
	  {
	    median_n = time_sum_n.at(time_sum_n.size()/2);
	  }
	else
	  {
	    median_n = (time_sum_n.at(time_sum_n.size()/2) + time_sum_n.at(time_sum_n.size()/2 - 1))/2.;
	  }
	float median_s;
	if (time_sum_s.size()%2)
	  {
	    median_s = time_sum_s.at(time_sum_s.size()/2);
	  }
	else
	  {
	    median_s = (time_sum_s.at(time_sum_s.size()/2) + time_sum_s.at(time_sum_s.size()/2 - 1))/2.;
	  }


	float vertex_median = 15.*(median_n - median_s);
	float vertex_mean = 15.*(sum_n/static_cast<float>(hits_n) - (sum_s/static_cast<float>(hits_s))) + vshift;

	float time_median = (median_n + median_s)/2.;
	float time_mean = (sum_n/static_cast<float>(hits_n) + (sum_s/static_cast<float>(hits_s)))/2.;

	vertex_by_median->Fill(min(hits_s, hits_n), vertex_median);
	vertex_by_mean->Fill(min(hits_s, hits_n), vertex_mean);
 	vertex_by_median2->Fill(hits_s, hits_n, vertex_median);
	vertex_by_mean2->Fill(hits_s, hits_n, vertex_mean);

	h_vertex_mean->Fill(vertex_mean);
	h_vertex_median->Fill(vertex_median);
	h_time_mean->Fill(time_mean);
	h_time_median->Fill(time_median);
	h_vtx_time_mean->Fill(vertex_mean, time_mean);
	h_vtx_time_median->Fill(vertex_median, time_median);

	h_vertex_mean_centrality[centrality_bin]->Fill(vertex_mean);
	h_vertex_median_centrality[centrality_bin]->Fill(vertex_median);

	diff_vertex_methods->Fill(min(hits_s, hits_n), vertex_median - vertex_mean);
	for (int ich = 0; ich < 128; ich++)
	  {
	    if (fabs(vertex_mean) < 30) he_occupancy_centrality_vtx_lt_30[centrality_bin]->Fill(occupancy[ich], ich);
	    else he_occupancy_centrality_vtx_gt_30[centrality_bin]->Fill(occupancy[ich], ich);
	  }
      }

      if (hits_s >=2 && hits_s_cover >=2){
	float median_s_cover;
	if (time_sum_s_cover.size()%2)
	  {
	    median_s_cover = time_sum_s_cover.at(time_sum_s_cover.size()/2);
	  }
	else
	  {
	    median_s_cover = (time_sum_s_cover.at(time_sum_s_cover.size()/2) + time_sum_s_cover.at(time_sum_s_cover.size()/2 - 1))/2.;
	  }
	float median_s;
	if (time_sum_s.size()%2)
	  {
	    median_s = time_sum_s.at(time_sum_s.size()/2);
	  }
	else
	  {
	    median_s = (time_sum_s.at(time_sum_s.size()/2) + time_sum_s.at(time_sum_s.size()/2 - 1))/2.;
	  }

	float vertex_median_cover = 15.*(median_s_cover - median_s);
	float vertex_mean_cover = 15.*(sum_s_cover/static_cast<float>(hits_s_cover) - (sum_s/static_cast<float>(hits_s)));

	float time_median_cover = (median_s_cover + median_s)/2.;
	float time_mean_cover = (sum_s_cover/static_cast<float>(hits_s_cover) + (sum_s/static_cast<float>(hits_s)))/2.;

	vertex_by_median_cover->Fill(min(hits_s, hits_s_cover), vertex_median_cover);
	vertex_by_mean_cover->Fill(min(hits_s, hits_s_cover), vertex_mean_cover);
	vertex_by_median_cover2->Fill(hits_s, hits_s_cover, vertex_median_cover);
	vertex_by_mean_cover2->Fill(hits_s, hits_s_cover, vertex_mean_cover);
	nhit_dist_cover->Fill(hits_s, hits_s_cover);
	diff_vertex_methods_cover->Fill(min(hits_s, hits_s_cover), vertex_median_cover - vertex_mean_cover);
	h_vertex_mean_cover->Fill(vertex_mean_cover);
	h_vertex_median_cover->Fill(vertex_median_cover);
	h_time_mean_cover->Fill(time_mean_cover);
	h_time_median_cover->Fill(time_median_cover);
	h_vtx_time_mean_cover->Fill(vertex_mean_cover, time_mean_cover);
	h_vtx_time_median_cover->Fill(vertex_median_cover, time_median_cover);
	h_vertex_mean_cover_centrality[centrality_bin]->Fill(vertex_mean_cover);
	h_vertex_median_cover_centrality[centrality_bin]->Fill(vertex_median_cover);

	for (int ich = 0; ich < 128; ich++)
	  {
	    if (fabs(vertex_mean_cover) < 30) he_occupancy_centrality_vtx_cover_lt_30[centrality_bin]->Fill(occupancy[ich], ich);
	    else he_occupancy_centrality_vtx_cover_gt_30[centrality_bin]->Fill(occupancy[ich], ich);
	  }

      }

      for (int ich = 0 ; ich < 64; ich++)
	{
	  
	  if (mbd_charge[ich] > thresh && mbd_time[ich] < timethresh)
	    {
	      distance_from_mean[ich]->Fill(TMath::Abs(h_time_n->GetMean() - mbd_time[ich]));
	    }      
	  if (mbd_charge[ich + 64] > thresh && mbd_time[ich+64] < timethreshs)
	    {
	      distance_from_mean[ich+64]->Fill(TMath::Abs(h_time_n->GetMean() - mbd_time[ich+64]));
	    }
	  
	  
	} 
      
      hp_rms_by_hits_s->Fill(hits_s, h_time_s->GetRMS());
      hp_rms_by_hits_n->Fill(hits_n, h_time_n->GetRMS());
    }

  TCanvas *c1 = new TCanvas("c1", "c1");
 
  hp_rms_by_hits_n->Draw("hist");
  hp_rms_by_hits_s->SetLineColor(kRed);
  hp_rms_by_hits_s->Draw("hist same");


  TCanvas *c2 = new TCanvas("c2", "c2");
  
  h_vertex_mean->Draw("hist");

  h_vertex_median->SetLineColor(kRed);
  h_vertex_median->Draw("hist same");
  h_vertex_mean_cover->SetLineStyle(4);
  h_vertex_mean_cover->Draw("hist same");
  h_vertex_median_cover->SetLineStyle(4);
  h_vertex_median_cover->SetLineColor(kRed);
  h_vertex_median_cover->Draw("hist same");
  
  TH1D *h_half_mast = new TH1D("h_half_mast","", 128, -0.5, 127.5);

  for (int i = 0; i < 128;i++)
    {
      for (int j = 0; j < distance_from_mean[i]->GetNbinsX();j++)
	{
	  if (distance_from_mean[i]->GetBinContent(j+1) < distance_from_mean[i]->GetBinContent( distance_from_mean[i]->GetMaximumBin())/2.)
	    {
	      h_half_mast->Fill(i, distance_from_mean[i]->GetBinLowEdge(j));
	      break;
	    }
	}

    }


  TFile *out = new TFile(Form("mbd_vertex_plots_%d.root", runnumber), "RECREATE");
  hp_rms_by_hits_n->Write();
  hp_rms_by_hits_s->Write();
  vertex_by_mean->Write();
  vertex_by_median->Write();
  vertex_by_mean_cover->Write();
  vertex_by_median_cover->Write();

  vertex_by_mean2->Write();
  vertex_by_median2->Write();
  vertex_by_mean_cover2->Write();
  vertex_by_median_cover2->Write();

  h_vertex_mean->Write();
  h_vertex_mean_center->Write();
  h_vertex_median->Write();
  h_vertex_mean_cover->Write();
  h_vertex_median_cover->Write();
  h_time_mean->Write();
  h_time_mean_center->Write();
  h_time_median->Write();
  h_time_mean_cover->Write();
  h_time_median_cover->Write();
  h_vtx_time_mean->Write();
  h_vtx_time_mean_center->Write();
  h_vtx_time_median->Write();
  h_vtx_time_mean_cover->Write();
  h_vtx_time_median_cover->Write();

  h_tubes_s->Write();
  h_tubes_n->Write();
  diff_vertex_methods->Write();
  nhit_dist->Write();
  nhit_dist_cover->Write();
  h_half_mast->Write();


  for (int i = 0; i < 20; i++)
    {
      h_vertex_mean_cover_centrality[i]->Write();
      h_vertex_median_cover_centrality[i]->Write();
      h_vertex_mean_centrality[i]->Write();
      h_vertex_median_centrality[i]->Write();

      he_occupancy_centrality[i]->Write();
      he_occupancy_centrality_vtx_lt_30[i]->Write();
      he_occupancy_centrality_vtx_gt_30[i]->Write();
      he_occupancy_centrality_vtx_cover_lt_30[i]->Write();
      he_occupancy_centrality_vtx_cover_gt_30[i]->Write();
      he_occupancy_centrality_vtx_center_lt_30[i]->Write();
      he_occupancy_centrality_vtx_center_gt_30[i]->Write();
    }
  for (int i = 0; i < 128; i++)
    {
      h2_charge_time[i]->Write();
      h_time_above_50[i]->Write();
      h_charge_above_50[i]->Write();
      h_time_below_50[i]->Write();
      h_charge_below_50[i]->Write();
      distance_from_mean[i]->Write();
      h_time_dist[i]->Write();
    }
  out->Close();

}
