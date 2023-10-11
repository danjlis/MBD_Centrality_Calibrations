#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <centrality/CentralityReco.h>
#include <centrality/QACentralityReco.h>
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4allraw/Fun4AllPrdfInputManager.h>
#include <vector>
#include <bbc/BbcReco.h>

#include <CaloWaveFormSim.h>
#include <CaloTriggerEmulator.h>
#include <CaloPacketGetter.h>
#include <calowaveformsim/MBDEmulatorTreeMaker.h>
#include <caloreco/CaloTowerBuilder.h>

R__LOAD_LIBRARY(libcalo_reco.so) 
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libcentrality.so)
R__LOAD_LIBRARY(libqacentrality.so)
R__LOAD_LIBRARY(libbbc_io.so)
#endif

void Fun4All_CentralityRecoMBD(const int runnumber, const int rollover = 0)
{
  gSystem->Load("libg4dst");
  gSystem->Load("libcentrality");
  gSystem->Load("libqacentrality");
  gSystem->Load("libbbc_io");

  std::ostringstream rstr;
  rstr << std::setw(8) << std::setfill('0') << runnumber;

  std::ostringstream ostr;
  ostr << std::setw(4) << std::setfill('0') << rollover;

  std::string fname1 = Form("/sphenix/lustre01/sphnxpro/commissioning/aligned/beam-%s-%s.prdf", rstr.str().c_str(), ostr.str().c_str());

  Fun4AllServer *se = Fun4AllServer::instance();

  Fun4AllInputManager *in = new Fun4AllPrdfInputManager("in");
  in->fileopen(fname1);
  se->registerInputManager(in);

  BbcReco *ca_1 = new BbcReco();
  ca_1->Verbosity(10);
  se->registerSubsystem(ca_1);

  CaloTowerBuilder *ca_2 = new CaloTowerBuilder();
  ca_2->set_detector_type(CaloTowerBuilder::ZDC);
  ca_2->set_nsamples(31);
  ca_2->set_processing_type(CaloWaveformProcessing::FAST);
  se->registerSubsystem(ca_2);

// Fun4All
  CentralityReco *cr = new CentralityReco();
  cr->Verbosity(1);
  se->registerSubsystem(cr);

  QACentralityReco *crq = new QACentralityReco("QACentralityReco","rootfile1.root","rootfile2.root");
  crq->Verbosity(1);
  se->registerSubsystem(crq);

  se->run(100);
  se->End();
}
