#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4allraw/Fun4AllPrdfInputManager.h>
#include <vector>
#include <phool/recoConsts.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>


// #include <calotowerbuilder/CaloTowerBuilder.h>
#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloWaveformProcessing.h>
#include <caloreco/CaloTowerCalib.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawClusterPositionCorrection.h>

#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <caloreco/DeadHotMapLoader.h>

#include <caloreco/TowerInfoDeadHotMask.h>

#include <caloreco/RawClusterDeadHotMask.h>

#include <centrality/MbdCalibrationAnalysis.h>
R__LOAD_LIBRARY(libcalo_reco.so) 
R__LOAD_LIBRARY(libmbdcalibrationanalysis.so) 
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
#endif

void Fun4All_MbdCalibrationAnalysis(const int runnumber, const int rollover = 0)
{
  gSystem->Load("libg4dst");
  gSystem->Load("libcalo_reco");
  gSystem->Load("libmbdcalibrationanalysis");

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
      sprintf(dir, "macros");
    }
  else
    {
      sprintf(dir, "output/run%d/mbdcalibana", runnumber);
    }
  const char *tree_outfile = Form("%s/%s/mbd_calib_ana_tree_%s_%s.root", env_p, dir, rstr.str().c_str(), ostr.str().c_str());

  std::string fname1 = Form("/sphenix/lustre01/sphnxpro/commissioning/aligned/beam-%s-%s.prdf", rstr.str().c_str(), ostr.str().c_str());

  if (FILE *file = fopen(fname1.c_str(),"r")){
    fclose(file);
  }
  else
  {
    fname1 = Form("/sphenix/lustre01/sphnxpro/commissioning/aligned_prdf/beam-%s-%s.prdf", rstr.str().c_str(), ostr.str().c_str());


    if (FILE *file = fopen(fname1.c_str(),"r")){
      fclose(file);
    }
    else
      {
	std::cout << "NOOOOO ... no "<< fname1 <<std::endl;
	return;
      }
  }
  Fun4AllServer *se = Fun4AllServer::instance();
  recoConsts *rc = recoConsts::instance();

  //===============
  // conditions DB flags
  //===============
  // ENABLE::CDB = true;
  // global tag
  rc->set_StringFlag("CDB_GLOBALTAG","ProdA_2023");
  // // 64 bit timestamp
  rc->set_uint64Flag("TIMESTAMP",stoi(fname1.substr(fname1.length()-15,5)));


  Fun4AllInputManager *in = new Fun4AllPrdfInputManager("in");
  in->fileopen(fname1);
  se->registerInputManager(in);

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

  CaloTowerCalib *calibZDC = new CaloTowerCalib("ZDC");
  calibZDC -> set_detector_type(CaloTowerCalib::ZDC);
  calibZDC->Verbosity(verbosity);
  se -> registerSubsystem(calibZDC);

  MbdCalibrationAnalysis *cr = new MbdCalibrationAnalysis("MbdCalibrationsAnalysis",tree_outfile);
  cr->Verbosity(verbosity);
  se->registerSubsystem(cr);

// Fun4All

  se->run(nevents);
  se->End();
}
