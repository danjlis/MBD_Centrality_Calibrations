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

#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libcentrality.so)
R__LOAD_LIBRARY(libcalotrigger.so)
#endif

void Fun4All_CentralityRecoMBD(const int runnumber, const int rollover = 0)
{
  gSystem->Load("libg4dst");
  gSystem->Load("libcentrality");
  gSystem->Load("libcalotrigger");


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

  const char* env_dst = std::getenv("DST_SOURCE_PATH");

  if(!env_dst)
    {
      std::cout << "no env DST_SOURCE_PATH set."<<endl;
      return;
    }

  const char* env_dstname = std::getenv("DST_NAME");

  if(!env_dstname)
    {
      std::cout << "no env DST_NAME set."<<endl;
      return;
    }

  const char* env_cdb = std::getenv("CDB_TAG");

  if(!env_cdb)
    {
      std::cout << "no env CDB_TAG set."<<endl;
      return;
    }

  char *dir = new char[100];
  if (rollover == -1)
    {
      sprintf(dir, "macros/fun4all");
    }
  else
    {
      sprintf(dir, "output/run%d/centana", runnumber);
    }

  std::string fname1 = Form("%s/%s-%s-%s.root", env_dst, env_dstname, rstr.str().c_str(), ostr.str().c_str());
  
  
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
  rc->set_StringFlag("CDB_GLOBALTAG","ProdA_2023");
  // // 64 bit timestamp                                                                                                                                
  rc->set_uint64Flag("TIMESTAMP",runnumber);


  Fun4AllInputManager *in = new Fun4AllDstInputManager("in");
  in->fileopen(fname1);
  se->registerInputManager(in);

  // Fun4All
  CentralityReco *cr = new CentralityReco();
  cr->Verbosity(1);
  se->registerSubsystem(cr);

  MinimumBiasClassifier *mb = new MinimumBiasClassifier();
  mb->Verbosity(1);
  se->registerSubsystem(mb);

  // QACentralityReco *crq = new QACentralityReco("QACentralityReco","rootfile1.root","rootfile2.root");
  // crq->Verbosity(1;
  // se->registerSubsystem(crq);

  se->run(10);
  se->End();
}
