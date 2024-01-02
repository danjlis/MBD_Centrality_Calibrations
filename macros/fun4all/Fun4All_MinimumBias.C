#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>

#include <centrality/CentralityReco.h>
#include <calotrigger/MinimumBiasClassifier.h>

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <vector>
#include <globalvertex/GlobalVertexReco.h>

#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libcentrality.so)
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libglobalvertex.so)

#endif

void Fun4All_MinimumBias(const int runnumber, const int rollover = 0)
{
  gSystem->Load("libg4dst");
  gSystem->Load("libcentrality");
  gSystem->Load("libcalotrigger");

  std::ostringstream rstr;
  rstr << std::setw(8) << std::setfill('0') << runnumber;

  std::ostringstream ostr;
  ostr << std::setw(4) << std::setfill('0') << rollover;

  std::string fname1 = Form("/sphenix/lustre01/sphnxpro/commissioning/DST_ana.386_2023p003/DST_CALOR_ana.386_2023p003-%s-%s.root", rstr.str().c_str(), ostr.str().c_str());
  
  
  if (FILE *file = fopen(fname1.c_str(),"r")){
    fclose(file);
  }
  else
    {
      std::cout << "NOOOOO ... no "<< fname1 <<std::endl;
      return;
    }

  Fun4AllServer *se = Fun4AllServer::instance();


  recoConsts *rc = recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG","dlis");
  // // 64 bit timestamp                                                                                                                                
  rc->set_uint64Flag("TIMESTAMP",runnumber);


  Fun4AllInputManager *in = new Fun4AllDstInputManager("in");
  in->fileopen(fname1);
  se->registerInputManager(in);

  // Fun4All
  GlobalVertexReco *gvr = new GlobalVertexReco();  
  se->registerSubsystem(gvr);

  MinimumBiasClassifier *mbc = new MinimumBiasClassifier();
  se->registerSubsystem(mbc);

  CentralityReco *cr = new CentralityReco();
  se->registerSubsystem(cr);


  // QACentralityReco *crq = new QACentralityReco("QACentralityReco","rootfile1.root","rootfile2.root");
  // crq->Verbosity(1;
  // se->registerSubsystem(crq);

  se->run(10);
  se->End();
}
