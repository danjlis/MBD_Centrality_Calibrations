

void DrawMBDEvents(int evt)
{
  gStyle->SetOptStat(0);

  int runnumber = 21813;

  TFile *file = new TFile(Form("../output/run%d/trees_%d.root", runnumber, runnumber), "r");

  // 3 vertex cuts with centralities between them compared.
  
  TH1D *h_time_n = new TH1D("h_time_n","", 100, -25, 25);
  TH1D *h_time_s = new TH1D("h_time_s","", 100, -25, 25);
  float mbd_charge[128];
  float mbd_time[128];
  float mbd_charge_raw[128];
  float mbd_time_raw[128];
  int isMinBias;
  float z_vertex;
  
  TTree *t = (TTree*)file->Get("T");
  
  t->SetBranchAddress("mbd_charge",mbd_charge);
  t->SetBranchAddress("mbd_time",mbd_time);
  t->SetBranchAddress("mbd_charge_raw",mbd_charge);
  t->SetBranchAddress("mbd_time_raw",mbd_time);
  t->SetBranchAddress("z_vertex",&z_vertex);


  int d = 1;
  int i = evt;
  

  t->GetEntry(i);
  h_time_n->Reset();
  h_time_s->Reset();

  for (int ich = 0 ; ich < 64; ich++)
    {
      
      if (mbd_charge[ich] > 0.1)
	h_time_n->Fill(mbd_time[ich]);
      if (mbd_charge[ich + 64] > 0.1)
	h_time_s->Fill(mbd_time[ich + 64]);
    }
  
  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->Divide(2);

  c1->cd(1);
  h_time_n->Draw();
  
  c1->cd(2);
  h_time_s->Draw();
  
  

}
