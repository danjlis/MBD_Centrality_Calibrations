#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4allraw/SingleInttPoolInput.h>
#include <fun4allraw/Fun4AllStreamingInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4allraw/Fun4AllPrdfInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>

#include <fun4all/SubsysReco.h>
#include <phool/recoConsts.h>
#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

#include <mbd/MbdReco.h>
#include <centrality/InttMbd.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>

R__LOAD_LIBRARY(libmbd_io.so) 
R__LOAD_LIBRARY(libinttmbd.so) 
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libffarawmodules.so)
R__LOAD_LIBRARY(libcalo_reco.so) 

SingleInttPoolInput *sngl[9]{};

void Fun4All_INTTMBD(int runnumber = 23697, int nEvents = 0,
			   const string &input_file00 = "intt0.list",
			   const string &input_file01 = "intt1.list",
			   const string &input_file02 = "intt2.list",
			   const string &input_file03 = "intt3.list",
			   const string &input_file04 = "intt4.list",
			   const string &input_file05 = "intt5.list",
			   const string &input_file06 = "intt6.list",
			         const string &input_file07 = "intt7.list"
			   )
{
  vector<string> infile;
  infile.push_back(input_file00);
  infile.push_back(input_file01);
  infile.push_back(input_file02);
  infile.push_back(input_file03);
  infile.push_back(input_file04);
  infile.push_back(input_file05);
  infile.push_back(input_file06);
  infile.push_back(input_file07);


  std::ostringstream rstr;
  rstr << std::setw(8) << std::setfill('0') << runnumber;
  int rollover = 0;
  int verbosity = 0;
  std::ostringstream ostr;

  ostr << std::setw(4) << std::setfill('0') << rollover;
  
  const char* env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }
  char *dir = new char[100];
  if (rollover == -1)
    {
      sprintf(dir, "macros/fun4all");
    }
  else
    {
      sprintf(dir, "output/run%d/mbdana", runnumber);
    }

  const char *tree_outfile = Form("%s/%s/mbd_ana_tree_%s_%s.root", env_p, dir, rstr.str().c_str(), ostr.str().c_str());

  std::string fname1 = Form("/sphenix/lustre01/sphnxpro/commissioning/aligned/beam-%s-%s.prdf", rstr.str().c_str(), ostr.str().c_str());

  if (FILE *file = fopen(fname1.c_str(),"r")){
    fclose(file);
  }
  else
  {
    fname1 = Form("/sphenix/lustre01/sphnxpro/commissioning/aligned_prdf/beam-%s-%s.prdf", rstr.str().c_str(), ostr.str().c_str());


    if (FILE *file = fopen(fname1.c_str(),"r")){
      fclose(file);
    }
    else
      {
	std::cout << "NOOOOO ... no "<< fname1 <<std::endl;
	return;
      }
  }

  cout <<fname1<<endl;
  Fun4AllServer *se = Fun4AllServer::instance();
  recoConsts *rc = recoConsts::instance();
  //  rc->set_IntFlag("RUNNUMBER",20445);

  //===============
  // conditions DB flags
  //===============
  // ENABLE::CDB = true;
  // global tag
  rc->set_StringFlag("CDB_GLOBALTAG","2023p002");
  // // 64 bit timestamp
  rc->set_uint64Flag("TIMESTAMP",runnumber);

  Fun4AllInputManager *in2 = new Fun4AllPrdfInputManager("in");
  in2->fileopen(fname1);
  se->registerInputManager(in2);

  Fun4AllStreamingInputManager *in = new Fun4AllStreamingInputManager("Comb");
  //  in->Verbosity(10);
  int i=0;
  for (auto iter : infile)
    {
      SingleInttPoolInput *sngl= new SingleInttPoolInput("INTT_" + to_string(i));
      //    sngl->Verbosity(3);
      sngl->AddListFile(iter);
      in->registerStreamingInput(sngl,Fun4AllStreamingInputManager::INTT);
      i++;
    }
  se->registerInputManager(in);


  MbdReco *mbdreco = new MbdReco();
  se->registerSubsystem( mbdreco );

  InttMbd *it = new InttMbd("inttmbd","inttmbd.root");
  it->Verbosity(3);
  se->registerSubsystem(it);


  se->run(nEvents);

  se->End();
  delete se;
  gSystem->Exit(0);
}
#endif
