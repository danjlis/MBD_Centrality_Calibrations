#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>


#include <g4bbc/BbcDigitization.h>
#include <bbc/BbcReco.h>
#include <calowaveformsim/TrigTreeMaker.h>
#include <MBDTriggerEmulator.h>
#include <DLUtility.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4bbc.so)

#endif

void Fun4All_CentralityRecoMBDSim()
{
  int verbosity = 10;

  std::cout << "Starting my fun4all script." <<endl;
  gSystem->Load("libg4dst");
  gSystem->Load("libg4bbc");
  gSystem->Load("libbbc_io");

  Fun4AllServer *se = Fun4AllServer::instance();

  Fun4AllInputManager *in = new Fun4AllDstInputManager("in");
  in->AddListFile("dst_truth.list");

  Fun4AllInputManager *in2 = new Fun4AllDstInputManager("in2");
  in2->AddListFile("dst_calo_cluster.list");
  
  Fun4AllInputManager *in3 = new Fun4AllDstInputManager("in3");
  in3->AddListFile("g4hits.list");


  se->registerInputManager(in3);
  se->registerInputManager(in2);
  se->registerInputManager(in);
  
  auto bbcdigi = new BbcDigitization();
  bbcdigi->Verbosity(verbosity);
  se->registerSubsystem(bbcdigi);

  auto bbcreco = new BbcReco();
  bbcreco->Verbosity(verbosity);
  se->registerSubsystem(bbcreco);


// Fun4All
  se->run(10);
  se->End();
}
