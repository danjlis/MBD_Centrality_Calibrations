#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <frog/FROG.h>
#include <calotrigger/TriggerRunInfoReco.h>
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <vector>
#include <phool/recoConsts.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>


#include <caloreco/CaloTowerBuilder.h>
#include <fun4all/Fun4AllUtils.h>
#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>
#include <mbd/MbdReco.h>
#include <zdcinfo/ZdcReco.h>
//#include <getepinfo/GetEPinfo.h>

#include <mbdana/MbdAna.h>
#include <caloreco/CaloTowerCalib.h>

//#include <calowaveformsim/GL1TriggerSelect.h>
#include <globalvertex/GlobalVertexReco.h>
#include <centrality/CentralityReco.h>
//#include <centrality/CentralityValid.h>
#include <calotrigger/MinimumBiasClassifier.h>
#include <eventplaneinfo/EventPlaneReco.h>
#include <epd/EpdReco.h>

//R__LOAD_LIBRARY(libcentralityvalid.so)
//R__LOAD_LIBRARY(libGetEPinfo.so)
R__LOAD_LIBRARY(libepd.so)
R__LOAD_LIBRARY(libeventplaneinfo.so)
//R__LOAD_LIBRARY(libemulatortreemaker.so)
R__LOAD_LIBRARY(libcalo_reco.so) 
R__LOAD_LIBRARY(libmbd.so) 
R__LOAD_LIBRARY(libzdcinfo.so) 
//R__LOAD_LIBRARY(libmbd.so) 
R__LOAD_LIBRARY(libmbdana.so) 
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libcentrality.so)
R__LOAD_LIBRARY(libglobalvertex.so)
#endif


void Fun4All_MbdAnaDST24(const std::string filename, const int phase = 0)
{
  gSystem->Load("libg4dst");
  gSystem->Load("libcalo_reco");
  gSystem->Load("libFROG");
  gSystem->Load("libcentrality");
  //  gSystem->Load("libcentralityvalid");


  std::pair<int, int> run_seg = Fun4AllUtils::GetRunSegment(filename);
  int runnumber = run_seg.first;
  int segment = run_seg.second;

  int nevents = 100000;
  int verbosity = 0;

  CaloTowerDefs::BuilderType buildertype = CaloTowerDefs::kPRDFTowerv4;

  const char* env_p = std::getenv("MBD_CENTRALITY_CALIB_PATH");

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }


  const char* env_cdb = std::getenv("CDB_TAG");

  if(!env_cdb)
    {
      std::cout << "no env CDB_TAG set."<<endl;
      return;
    }

  char *dir = new char[100];

  sprintf(dir, "output/run%d/mbdana", runnumber);

  const char *ep_outfile = Form("epinfo_tree_%08d_%05d.root", runnumber, segment);
  
  char *calib_outfile = Form("epcalib_%08d_%05d.root", runnumber, segment);
  if (phase)
    {
      calib_outfile = Form("%s/%s/event_plane_calib_phase%d_%d.root", env_p, dir, phase, runnumber);
    }
  const char *tree_outfile = Form("mbd_ana_tree_%08d_%05d.root", runnumber, segment);

  const char *hist_outfile = Form("mbd_ana_hist_%08d_%05d.root",runnumber, segment);

  // const char *ep_outfile = Form("%s/%s/epinfo_tree_%08d_%05d.root", env_p, dir, runnumber, segment);
  
  // const char *calib_outfile = Form("%s/%s/epcalib_%08d_%05d.root", env_p, dir, runnumber, segment);

  // const char *tree_outfile = Form("%s/%s/mbd_ana_tree_%08d_%05d.root", env_p, dir, runnumber, segment);

  // const char *hist_outfile = Form("%s/%s/mbd_ana_hist_%08d_%05d.root", env_p, dir, runnumber, segment);

  
  Fun4AllServer *se = Fun4AllServer::instance();
  recoConsts *rc = recoConsts::instance();
  se->Verbosity(verbosity);
  //===============
  //===============
  // ENABLE::CDB = true;
  // global tag
  rc->set_StringFlag("CDB_GLOBALTAG","ProdA_2024");
  // // 64 bit timestamp
  rc->set_uint64Flag("TIMESTAMP",runnumber);
  CDBInterface::instance()->Verbosity(1);
  FROG *fr = new FROG();

  std::string input = fr->location(filename.c_str());
  Fun4AllInputManager *in = new Fun4AllDstInputManager("in");
  in->fileopen(filename.c_str());
  in->Verbosity(0);
  se->registerInputManager(in);

  //  TriggerRunInfoReco* trig = new TriggerRunInfoReco();
  //se->registerSubsystem(trig);

  // MBD/BBC Reconstruction
  MbdReco *mbdreco = new MbdReco();
  mbdreco->Verbosity();
  se->registerSubsystem(mbdreco);

  EpdReco *epdreco = new EpdReco();
  se->registerSubsystem(epdreco);

  CaloTowerBuilder *caZDC = new CaloTowerBuilder("ZDCBUILDER");
  caZDC->set_detector_type(CaloTowerDefs::ZDC);
  caZDC->set_builder_type(buildertype);
  caZDC->set_processing_type(CaloWaveformProcessing::FAST);
  caZDC->set_nsamples(16);
  caZDC->set_offlineflag();
  se->registerSubsystem(caZDC);

  //ZDC Reconstruction--Calib Info
  ZdcReco *zdcreco = new ZdcReco();
  zdcreco->set_zdc1_cut(0.0);
  zdcreco->set_zdc2_cut(0.0);
  se->registerSubsystem(zdcreco);

  GlobalVertexReco *gvr = new GlobalVertexReco("GlobalVertexReco");
  se->registerSubsystem(gvr);

  // GL1TriggerSelect *gts = new GL1TriggerSelect("GL1TriggerSelect");
  // //  gts->Verbosity(2);
  // gts->select_trigger(10);
  // se->registerSubsystem(gts);

  // Fun4All

  // MinimumBiasClassifier *mb = new MinimumBiasClassifier();
  // mb->Verbosity(0);
  // se->registerSubsystem(mb);

  // CentralityReco *cr = new CentralityReco();
  // cr->Verbosity(0);
  // se->registerSubsystem(cr);

  // EventPlaneReco *epreco = new EventPlaneReco();
  // epreco->set_mbd_epreco(true);
  // epreco->set_sepd_epreco(false);
  // se->registerSubsystem(epreco);

  // CentralityValid *centralityvalidation = new CentralityValid("CentralityValid",hist_outfile);
  // se->registerSubsystem(centralityvalidation);

  MbdAna *mbdana = new MbdAna("MbdAna", tree_outfile);
  mbdana->Verbosity(verbosity);
  se->registerSubsystem(mbdana);

  // GetEPinfo *epana = new GetEPinfo("GetEPInfo", ep_outfile, calib_outfile);
  // epana->SetPhase(phase);
  // se->registerSubsystem(epana);

  se->run(0);
  se->End();

  std::cout << "DONE :)"<<std::endl;
}
