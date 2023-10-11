#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <qacentrality/QACentralityReco.h>
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
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
#endif

void Fun4All_CentralityRecoQA(const int runnumber, const int rollover = 0)
{
  gSystem->Load("libg4dst");
  gSystem->Load("libcentrality");

  std::ostringstream rstr;
  rstr << std::setw(8) << std::setfill('0') << runnumber;

  std::ostringstream ostr;
  ostr << std::setw(4) << std::setfill('0') << rollover;

  const char *hist_outfile = Form("/gpfs02/sphenix/user/dlis/Projects/centrality/output/run%d/centrality_reco_hist_%s_%s.root", runnumber, rstr.str().c_str(), ostr.str().c_str());
  const char *tree_outfile = Form("/gpfs02/sphenix/user/dlis/Projects/centrality/output/run%d/centrality_reco_tree_%s_%s.root", runnumber, rstr.str().c_str(), ostr.str().c_str());

  std::string fname1 = Form("/sphenix/lustre01/sphnxpro/commissioning/aligned/beam-%s-%s.prdf", rstr.str().c_str(), ostr.str().c_str());

  Fun4AllServer *se = Fun4AllServer::instance();

  CaloTowerBuilder *ca_1 = new CaloTowerBuilder();
  ca_1->set_detector_type(CaloTowerBuilder::MBD);
  ca_1->set_nsamples(31);
  ca_1->set_processing_type(CaloWaveformProcessing::FAST);
  se->registerSubsystem(ca_1);

  CaloTowerBuilder *ca_2 = new CaloTowerBuilder();
  ca_2->set_detector_type(CaloTowerBuilder::ZDC);
  ca_2->set_nsamples(31);
  ca_2->set_processing_type(CaloWaveformProcessing::FAST);
  se->registerSubsystem(ca_2);

  QACentralityReco *cr = new QACentralityReco("QACentralityReco", hist_outfile, tree_outfile);
  cr->Verbosity(0);
  se->registerSubsystem(cr);

  Fun4AllInputManager *in = new Fun4AllPrdfInputManager("in");
  in->fileopen(fname1);
  se->registerInputManager(in);

// Fun4All

  se->run(100000);
  se->End();
}
