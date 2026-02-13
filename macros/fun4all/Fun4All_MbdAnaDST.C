#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <vector>
#include <phool/recoConsts.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <calowaveformsim/GL1TriggerSelect.h>
#include <globalvertex/GlobalVertexReco.h>
#include <centrality/CentralityReco.h>
#include <centrality/CentralityValid.h>
#include <calotrigger/MinimumBiasClassifier.h>
#include <getepinfo/GetEPinfo.h>

// #include <calotowerbuilder/CaloTowerBuilder.h>

#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

#include <centrality/MbdAna.h>

R__LOAD_LIBRARY(libcalo_reco.so) 
R__LOAD_LIBRARY(libmbd_io.so)
R__LOAD_LIBRARY(libgetepinfo.so) 
//R__LOAD_LIBRARY(libmbd.so) 
R__LOAD_LIBRARY(libmbdana.so) 
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)
#endif

void Fun4All_MbdAnaDST(const int runnumber, const int rollover = -1)
{
  gSystem->Load("libg4dst");
  gSystem->Load("libcalo_reco");
  

  std::ostringstream rstr;
  rstr << std::setw(8) << std::setfill('0') << runnumber;

  int runnumber_low = (runnumber - runnumber%100);
  std::ostringstream rlowstr;
  rlowstr << std::setw(8) << std::setfill('0') << runnumber_low;
  int runnumber_high = ( runnumber_low + 100);
  std::ostringstream rhighstr;

  rhighstr << std::setw(8) << std::setfill('0') << runnumber_high;

  
  int nevents = 0;
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
      sprintf(dir, "output/run%d/mbdana", runnumber);
    }

  char *dstfolder = new char[100];
  sprintf(dstfolder, "run_%s_%s", rlowstr.str().c_str(), rhighstr.str().c_str());

  const char *tree_outfile = Form("%s/%s/mbd_ana_tree_%s_%s.root", env_p, dir, rstr.str().c_str(), ostr.str().c_str());

  std::string fname1 = Form("%s/%s-%s-%s.root", env_dst, env_dstname, rstr.str().c_str(), ostr.str().c_str());

  const char *ep_outfile = Form("%s/%s/epinfo_tree_%s_%s.root", env_p, dir, rstr.str().c_str(), ostr.str().c_str());


  
  if (FILE *file = fopen(fname1.c_str(),"r")){
    fclose(file);
  }
  else
    {
      std::cout << "NOOOOO ... no "<< env_dstname << " in "<<fname1 <<std::endl;
      return;
    }


  Fun4AllServer *se = Fun4AllServer::instance();
  recoConsts *rc = recoConsts::instance();
  se->Verbosity(verbosity);
  //===============
  // conditions DB flags
  //===============
  // ENABLE::CDB = true;
  // global tag
  rc->set_StringFlag("CDB_GLOBALTAG",env_cdb);
  // // 64 bit timestamp
  rc->set_uint64Flag("TIMESTAMP",runnumber);


  Fun4AllInputManager *in = new Fun4AllDstInputManager("in");
  in->fileopen(fname1);
  se->registerInputManager(in);

  GL1TriggerSelect *gts = new GL1TriggerSelect("GL1TriggerSelect");
  gts->select_trigger(10);
  se->registerSubsystem(gts);

  MinimumBiasClassifier *mb = new MinimumBiasClassifier();
  mb->Verbosity(0);
  se->registerSubsystem(mb);

  CentralityReco *cr = new CentralityReco();
  cr->Verbosity(1);
  se->registerSubsystem(cr);
  

  CentralityValid *centralityvalidation = new CentralityValid("CentralityValid","cent_valid.root");
  se->registerSubsystem(centralityvalidation);

  MbdAna *mbdana = new MbdAna("MbdAna", tree_outfile);
  mbdana->Verbosity(verbosity);
  se->registerSubsystem(mbdana);

  GetEPinfo *epana = new GetEPinfo("GetEPInfo", ep_outfile);
  se->registerSubsystem(epana);


  std::cout << "DONE :)"<<std::endl;
}
