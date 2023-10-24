void get_n_events(std::string f)
{
  TFile *fd = new TFile(f.c_str(), "r");

  if (!fd) 
    {

      std::cout << "nothing..."<<endl;
      return;
    }

  TTree *t = (TTree*) fd->Get("T");
  if (!t) 
    {

      std::cout << "nothing..."<<endl;
      return;
    }

  std::cout <<"Number of Events: "<<t->GetEntries()<<endl;
  
  return;
}
