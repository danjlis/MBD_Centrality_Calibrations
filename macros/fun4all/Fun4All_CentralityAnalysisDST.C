#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)


#include <nodedump/Dumper.h>
#include <centrality/CentralityAnalysis.h>
#include <centrality/CentralityReco.h>
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4allraw/Fun4AllPrdfInputManager.h>
#include <vector>
#include <bbc/BbcReco.h>

#include <CaloWaveFormSim.h>
#include <CaloTriggerEmulator.h>
#include <CaloPacketGetter.h>
#include <calowaveformsim/MBDEmulatorTreeMaker.h>
#include <caloreco/CaloTowerBuilder.h>
#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libbbc_io.so)
R__LOAD_LIBRARY(libcalo_reco.so) 
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libcentrality.so)
R__LOAD_LIBRARY(libcentralityanalysis.so)
#endif

void Fun4All_CentralityAnalysisDST(const int runnumber, const int rollover = 0)
{
  gSystem->Load("libbbc_io");
  gSystem->Load("libg4dst");
  gSystem->Load("libcentrality");
  gSystem->Load("libcentralityanalysis");

  std::ostringstream rstr;
  rstr << std::setw(8) << std::setfill('0') << runnumber;
  int nevents = 100000;
  int verbosity = 0;
  std::ostringstream ostr;
  if (rollover == -1)
    {
      nevents = 10;
      verbosity = 10;
      ostr << std::setw(4) << std::setfill('0') << 0;
    }
  else {
    ostr << std::setw(4) << std::setfill('0') << rollover;
  }  
  const char* env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }
  char *dir = new char[100];
  if (rollover == -1)
    {
      sprintf(dir, "macros");
    }
  else
    {
      sprintf(dir, "output/run%d/centrality", runnumber);
    }

  const char *hist_outfile = Form("/gpfs02/%s/%s/centrality_reco_hist_%s_%s.root", env_p, dir, rstr.str().c_str(), ostr.str().c_str());
  const char *tree_outfile = Form("/gpfs02/%s/%s/centrality_reco_tree_%s_%s.root", env_p, dir, rstr.str().c_str(), ostr.str().c_str());

  std::string fname1 = Form("/sphenix/user/chiu/sphenix_bbc/RUN23_TEST/DST_MBD_%s-%s.root", rstr.str().c_str(), ostr.str().c_str());

  if (FILE *file = fopen(fname1.c_str(),"r")){
    fclose(file);
  }
  else
  {
    std::cout << "NOOOOO ... no "<< fname1 <<std::endl;
    return;
  }
  

  Fun4AllServer *se = Fun4AllServer::instance();

  // recoConsts *rc = recoConsts::instance();

  // //===============
  // // conditions DB flags
  // //===============
  // // ENABLE::CDB = true;
  // // global tag
  // rc->set_StringFlag("CDB_GLOBALTAG","ProdA_2023");
  // // // 64 bit timestamp
  // rc->set_uint64Flag("TIMESTAMP",stoi(fname1.substr(fname1.length()-15,5)));

  Fun4AllInputManager *in = new Fun4AllDstInputManager("in");
  in->fileopen(fname1);
  se->registerInputManager(in);


  //  DumpCdbUrlSave *ducus = new DumpCdbUrlSave("CdbUrl");
  //se->RegisterSubsystem(ducus);

  CentralityReco *cr = new CentralityReco("CentralityReco");
  cr->Verbosity(verbosity);  
  se->registerSubsystem(cr);

  CentralityAnalysis *car = new CentralityAnalysis("CentralityAnalysis", hist_outfile, tree_outfile);
  car->Verbosity(verbosity);
  car->useZDC(false);
  se->registerSubsystem(car);

// Fun4All

  se->run(nevents);
  se->End();
}
