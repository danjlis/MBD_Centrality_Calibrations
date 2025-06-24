#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
#include <mbd/MbdReco.h>
#include <g4mbd/MbdDigitization.h>
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

#include <centrality/CentralityReco.h>
#include <centrality/CentralityValid.h>
#include <calotrigger/MinimumBiasClassifier.h>
#include <centrality/MbdAna.h>
#include <DLUtility.h>
#include <fun4all/Fun4AllUtils.h>
#include <frog/FROG.h>
#include <globalvertex/GlobalVertexReco.h>

R__LOAD_LIBRARY(libglobalvertex.so)
R__LOAD_LIBRARY(libdlutility.so)
R__LOAD_LIBRARY(libFROG.so)
R__LOAD_LIBRARY(libcalo_reco.so) 
//R__LOAD_LIBRARY(libmbd_io.so) 
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libg4mbd.so) 
R__LOAD_LIBRARY(libmbdana.so) 
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libcentrality.so)
R__LOAD_LIBRARY(libcentralityvalid.so)
R__LOAD_LIBRARY(libphool.so)
#endif
  
void Fun4All_MbdAnaMC( const std::string input_3) //, const std::string input_2, const std::string input_1)
{
  gSystem->Load("libg4dst");
  gSystem->Load("libFROG");
  gSystem->Load("libcalo_reco");


  int nevents = 100000;
  int verbosity = 0;
  std::ostringstream ostr;
  
  std::pair<int, int> run_seg = Fun4AllUtils::GetRunSegment(input_3);
  int runnumber = run_seg.first;
  int segment = run_seg.second;

  ostr << std::setw(5) << std::setfill('0') << segment;

  const char* env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }
  
  //const char* env_mc = "hijing";
  const char* env_mc = std::getenv("MC_GENERATOR");

  if(!env_mc)
    {
      std::cout << "no env MC_GENERATOR set."<<endl;
      return;
    }


  const char* env_cdb = std::getenv("CDB_TAG");

  if(!env_cdb)
    {
      std::cout << "no env CDB_TAG set."<<endl;
      return;
    }

  // char *dir = new char[100];
  // sprintf(dir, "output/%s/mbdana/", env_mc);

  const char *tree_outfile = Form("mbd_ana_tree_%s_%s.root",  env_mc, ostr.str().c_str());
  const char *hist_outfile = Form("mbd_ana_hist_%s_%s.root",  env_mc, ostr.str().c_str());


  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(verbosity);
  recoConsts *rc = recoConsts::instance();

  //===============
  // conditions DB flags
  //===============
  // ENABLE::CDB = true;
  // global tag
  rc->set_StringFlag("CDB_GLOBALTAG","MDC2");
  // // 64 bit timestamp
  rc->set_uint64Flag("TIMESTAMP",runnumber);


  FROG *fr = new FROG();

  //std::string filename_1 = fr->location(input_1.c_str());
  //std::string filename_2 = fr->location(input_2.c_str());
  std::string filename_3 = fr->location(input_3.c_str());

  // Fun4AllInputManager *in1 = new Fun4AllDstInputManager("in1");
  // in1->fileopen(filename_1.c_str());
  // in1->Verbosity(0);
  // se->registerInputManager(in1);
  // Fun4AllInputManager *in2 = new Fun4AllDstInputManager("in2");
  // in2->fileopen(filename_2.c_str());
  // in2->Verbosity(0);
  // se->registerInputManager(in2);
  Fun4AllInputManager *in3 = new Fun4AllDstInputManager("in3");
  in3->fileopen(filename_3.c_str());
  in3->Verbosity(0);
  se->registerInputManager(in3);

  // auto mbddigi = new MbdDigitization();
  // mbddigi->Verbosity(verbosity);
  // se->registerSubsystem(mbddigi);
  
  // auto mbdreco = new MbdReco();
  // mbdreco->Verbosity(verbosity);
  // se->registerSubsystem(mbdreco);

  GlobalVertexReco* gblvertex = new GlobalVertexReco();
  gblvertex->Verbosity(verbosity);
  se->registerSubsystem(gblvertex);


  MinimumBiasClassifier *mb = new MinimumBiasClassifier();
  mb->Verbosity(0);
  mb->setIsSim(true);
  mb->setOverwriteScale("/sphenix/user/dlis/Projects/centrality/cdb/calibrations/scales/cdb_centrality_scale_1.root");
  mb->setOverwriteVtx("/sphenix/user/dlis/Projects/centrality/cdb/calibrations/vertexscales/cdb_centrality_vertex_scale_1.root");
  se->registerSubsystem(mb);

  CentralityReco *cr = new CentralityReco();
  cr->Verbosity(0);
  cr->setOverwriteScale("/sphenix/user/dlis/Projects/centrality/cdb/calibrations/scales/cdb_centrality_scale_1.root");
  cr->setOverwriteVtx("/sphenix/user/dlis/Projects/centrality/cdb/calibrations/vertexscales/cdb_centrality_vertex_scale_1.root");
  cr->setOverwriteDivs("/sphenix/user/dlis/Projects/centrality/cdb/calibrations/divs/cdb_centrality_1.root");
  se->registerSubsystem(cr);

  CentralityValid *centralityvalidation = new CentralityValid("CentralityValid",hist_outfile);
  centralityvalidation->setIsSim(true);
  se->registerSubsystem(centralityvalidation);

  MbdAna *mbdana = new MbdAna("MbdAna", tree_outfile);
  mbdana->Verbosity(verbosity);
  mbdana->setIsSim(true);
  se->registerSubsystem(mbdana);

  se->run(nevents);
  se->End();
}
