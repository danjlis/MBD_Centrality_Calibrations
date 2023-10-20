#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <centrality/CentralityAnalysis.h>
#include <centrality/CentralityReco.h>
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4allraw/Fun4AllPrdfInputManager.h>
#include <vector>

#include <CaloWaveFormSim.h>
#include <CaloTriggerEmulator.h>
#include <CaloPacketGetter.h>
#include <calowaveformsim/MBDEmulatorTreeMaker.h>
#include <caloreco/CaloTowerBuilder.h>

R__LOAD_LIBRARY(libcalo_reco.so) 
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libcentrality.so)
R__LOAD_LIBRARY(libcentralityanalysis.so)
#endif

void Fun4All_CentralityAnalysis(const int runnumber, const int rollover = 0)
{
  gSystem->Load("libg4dst");
  gSystem->Load("libcentrality");
  gSystem->Load("libcentralityanalysis");

  std::ostringstream rstr;
  rstr << std::setw(8) << std::setfill('0') << runnumber;

  std::ostringstream ostr;
  ostr << std::setw(4) << std::setfill('0') << rollover;

  const char* env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }

  const char *hist_outfile = Form("/gpfs02/%s/output/run%d/centrality/centrality_reco_hist_%s_%s.root", env_p, runnumber, rstr.str().c_str(), ostr.str().c_str());
  const char *tree_outfile = Form("/gpfs02/%s/output/run%d/centrality/centrality_reco_tree_%s_%s.root", env_p, runnumber, rstr.str().c_str(), ostr.str().c_str());

  const char *fname1 = Form("/sphenix/user/chiu/sphenix_bbc/RUN23_TEST/DST_MBD_%s-%s.root", rstr.str().c_str(), ostr.str().c_str());

  Fun4AllServer *se = Fun4AllServer::instance();

  Fun4AllInputManager *in = new Fun4AllDstInputManager("in");
  in->fileopen(fname1);
  se->registerInputManager(in);

  CentralityReco *cr = new CentralityReco("CentralityReco");
  cr->Verbosity(50);
  cr->useZDC(false);
  se->registerSubsystem(cr);

  CentralityAnalysis *car = new CentralityAnalysis("CentralityAnalysis", hist_outfile, tree_outfile);
  car->Verbosity(10);
  car->useZDC(false);
  se->registerSubsystem(car);

// Fun4All

  se->run(10);
  se->End();
}
