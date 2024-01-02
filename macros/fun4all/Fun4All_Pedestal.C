#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4allraw/Fun4AllPrdfInputManager.h>

#include <CaloWaveFormSim.h>
#include <CaloTriggerEmulator.h>
#include <CaloPacketGetter.h>
#include <PedTreeMaker.h>

R__LOAD_LIBRARY(libpedtreemaker.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libcalopacketgetter.so)
#endif

void Fun4All_Pedestal(const int runnumber, const string sebname, const int rollover = 1)
{
  gSystem->Load("libg4dst");
  gSystem->Load("libcalowaveformsim");
  gSystem->Load("libcalopacketgetter");

  const char* env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }

  std::ostringstream rstr;
  rstr << std::setw(8) << std::setfill('0') << runnumber;

  std::ostringstream ostr;
  ostr << std::setw(4) << std::setfill('0') << rollover;

  const char *tree_outfile = Form("/%s/output/run%d/waves/waveform_%s_%s_%s.root", env_p, runnumber, sebname, rstr.str().c_str(), ostr.str().c_str());

  std::string fname1 = Form("/sphenix/lustre01/sphnxpro/commissioning/emcal/junk/pedestal_%s-%s-%s.prdf", sebname.c_str(), rstr.str().c_str(), ostr.str().c_str());

  if (FILE *file = fopen(fname1.c_str(),"r")){
    fclose(file);
  }
  else
  {
    std::cout << "NOOOOO"<<std::endl;
    return;
  }

  Fun4AllServer *se = Fun4AllServer::instance();

  CaloPacketGetter *ca = new CaloPacketGetter("CALOPACKETGETTER_CEMC","CEMC");
  ca->set_nsamples(12);
  se->registerSubsystem(ca);

  PedTreeMaker *tt1 = new PedTreeMaker("WAVEMAKER_CEMC",tree_outfile);
  //tt1->SetVerbosity(1);
  se->registerSubsystem(tt1);
  
  Fun4AllInputManager *in = new Fun4AllPrdfInputManager("in");
  in->fileopen(fname1);
  se->registerInputManager(in);

// Fun4All

  se->run(10000);
  se->End();
}
