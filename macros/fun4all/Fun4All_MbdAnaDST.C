#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <vector>
#include <phool/recoConsts.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>


// #include <calotowerbuilder/CaloTowerBuilder.h>

#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

#include <centrality/MbdAna.h>

R__LOAD_LIBRARY(libcalo_reco.so) 
R__LOAD_LIBRARY(libmbd_io.so) 
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
      sprintf(dir, "macros/fun4all");
    }
  else
    {
      sprintf(dir, "output/run%d/mbdana", runnumber);
    }

  const char *tree_outfile = Form("%s/%s/mbd_ana_tree_%s_%s.root", env_p, dir, rstr.str().c_str(), ostr.str().c_str());

  std::string fname1 = Form("/sphenix/lustre01/sphnxpro/commissioning/DSTv3/DST_CALOR-%s-%s.root", rstr.str().c_str(), ostr.str().c_str());
  
  
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

  //===============
  // conditions DB flags
  //===============
  // ENABLE::CDB = true;
  // global tag
  rc->set_StringFlag("CDB_GLOBALTAG","2023p002");
  // // 64 bit timestamp
  rc->set_uint64Flag("TIMESTAMP",runnumber);


  Fun4AllInputManager *in = new Fun4AllDstInputManager("in");
  in->fileopen(fname1);
  se->registerInputManager(in);


  MbdAna *mbdana = new MbdAna("MbdAna", tree_outfile);
  mbdana->Verbosity(verbosity);
  se->registerSubsystem(mbdana);

  se->run(nevents);
  se->End();
}
