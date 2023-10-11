#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4allraw/Fun4AllPrdfInputManager.h>

#include <CaloWaveFormSim.h>
#include <CaloTriggerEmulator.h>
#include <CaloPacketGetter.h>
#include <WaveTreeMaker.h>
#include <calowaveformsim/MBDEmulatorTreeMaker.h>

R__LOAD_LIBRARY(libwavetreemaker.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libcalowaveformsim.so)
R__LOAD_LIBRARY(libcalotriggeremulator.so)
R__LOAD_LIBRARY(libcalopacketgetter.so)
#endif

void Fun4All_WaveForm(const int runnumber, const int rollover = 1)
{
  gSystem->Load("libg4dst");
  gSystem->Load("libcalowaveformsim");
  gSystem->Load("libcalotriggeremulator");
  gSystem->Load("libcalopacketgetter");
  gSystem->Load("libtmbdemulatortreemaker");

  std::ostringstream rstr;
  rstr << std::setw(8) << std::setfill('0') << runnumber;

  std::ostringstream ostr;
  ostr << std::setw(4) << std::setfill('0') << rollover;

  const char *tree_outfile = Form("/gpfs02/sphenix/user/dlis/Projects/centrality/output/run%d/waves/waveform_tree_%s_%s.root", runnumber, rstr.str().c_str(), ostr.str().c_str());

  std::string fname1 = Form("/sphenix/lustre01/sphnxpro/commissioning/aligned/beam-%s-%s.prdf", rstr.str().c_str(), ostr.str().c_str());

  if (FILE *file = fopen(fname1.c_str(),"r")){
    fclose(file);
  }
  else
  {
    fname1 = Form("/sphenix/lustre01/sphnxpro/commissioning/aligned_prdf/beam-%s-%s.prdf", rstr.str().c_str(), ostr.str().c_str());
  }

  if (FILE *file = fopen(fname1.c_str(),"r")){
    fclose(file);
  }
  else
  {
    std::cout << "NOOOOO"<<std::endl;
    return;
  }

  Fun4AllServer *se = Fun4AllServer::instance();

  CaloPacketGetter *ca = new CaloPacketGetter("CALOPACKETGETTER_MBD","MBD");
  ca->set_nsamples(31);
  se->registerSubsystem(ca);

  WaveTreeMaker *tt1 = new WaveTreeMaker("WAVEMAKER_MBD",tree_outfile);
  //tt1->SetVerbosity(1);
  se->registerSubsystem(tt1);
  
  Fun4AllInputManager *in = new Fun4AllPrdfInputManager("in");
  in->fileopen(fname1);
  se->registerInputManager(in);

// Fun4All

  se->run(100000);
  se->End();
}
