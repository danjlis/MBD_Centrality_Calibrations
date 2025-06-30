#ifndef CENTRALITYSCALECDB
#define CENTRALITYSCALECDB

#include <cdbobjects/CDBTTree.h>
#include <sphenixnpc/CDBUtils.h>
#include <ffamodules/CDBInterface.h>

#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libffamodules.so)

R__LOAD_LIBRARY(libsphenixnpc.so)
R__LOAD_LIBRARY(libphool.so)

R__LOAD_LIBRARY(libcdbobjects.so)

int CDB_MakeCentralityScaleCDBTree(const int runnumber)
{

  
  char *env_p = new char[200];
  sprintf(env_p,"%s",std::getenv("MBD_CENTRALITY_CALIB_PATH"));

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return 1;
    }
  TString cdb_name = Form("%s/cdb/calibrations24/scales/cdb_centrality_scale_%d.root", env_p, runnumber);

  if (runnumber == 0)
    {
      cdb_name = Form("%s/cdb/calibrations24/scales/cdb_centrality_scale_hijing.root", env_p);
    }
  else if (runnumber == 1)
    {
      cdb_name = Form("%s/cdb/calibrations24/scales/cdb_centrality_scale_ampt.root", env_p);
    }
  else if (runnumber == 2)
    {
      cdb_name = Form("%s/cdb/calibrations24/scales/cdb_centrality_scale_epos.root", env_p);
    }
  else if (runnumber == 3)
    {
      cdb_name = Form("%s/cdb/calibrations24/scales/cdb_centrality_scale_hijing_magoff.root", env_p);
    }
  else if (runnumber == 4)
    {
      cdb_name = Form("%s/cdb/calibrations24/scales/cdb_centrality_scale_ampt_magoff.root", env_p);
    }
  else if (runnumber == 5)
    {
      cdb_name = Form("%s/cdb/calibrations24/scales/cdb_centrality_scale_epos_magoff.root", env_p);
    }

  CDBTTree *cdbttree = nullptr;
  

  TString calib_file_name = Form("%s/output/plots/mbdana_charge_sum_%d.root", env_p, runnumber);
  if (runnumber == 0)
    {
      calib_file_name = Form("%s/output/plots/mbdana_charge_sum_hijing.root", env_p);
    }
  else if (runnumber == 1)
    {
      calib_file_name = Form("%s/output/plots/mbdana_charge_sum_ampt.root", env_p);
    }
  else if (runnumber == 2)
    {
      calib_file_name = Form("%s/output/plots/mbdana_charge_sum_epos.root", env_p);
    }
  else if (runnumber == 3)
    {
      calib_file_name = Form("%s/output/plots/mbdana_charge_sum_hijing_magoff.root", env_p);
    }
  else if (runnumber == 4)
    {
      calib_file_name = Form("%s/output/plots/mbdana_charge_sum_ampt_magoff.root", env_p);
    }
  else if (runnumber == 5)
    {
      calib_file_name = Form("%s/output/plots/mbdana_charge_sum_epos_magoff.root", env_p);
    }


  TFile *fcalib = new TFile(calib_file_name.Data(), "r");
  if (!fcalib) 
    {
      cout << " No file " <<endl;
      return 1;
    }
  if (runnumber != 54912 && runnumber != 54280 && runnumber >  100)
    {
      TNtuple *ts = (TNtuple*) fcalib->Get("tn_scale");
      if (!ts) 
	{
	  cout << " No TNtuple " << endl;
	  return 1;
	}
    
      float scale;
      ts->SetBranchAddress("scale",&scale);
      
      cdbttree = new CDBTTree(cdb_name.Data());
      for (int i=0; i<ts->GetEntries(); i++)
	{
	  ts->GetEntry(i);
	  cdbttree->SetDoubleValue(i,"centralityscale",scale);
	}
    }
  else
    {
      
      cdbttree = new CDBTTree(cdb_name.Data());
      
      cdbttree->SetDoubleValue(0,"centralityscale",1.0);
    }
  fcalib->Close();

  cdbttree->Commit();
  cdbttree->Print();
  cdbttree->WriteCDBTTree();
  delete cdbttree;
  return 0;
}

void CDB_ReadTree(const int runnumber)
{
  char *env_p = new char[200];
  sprintf(env_p,"%s",std::getenv("MBD_CENTRALITY_CALIB_PATH"));

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    } 
  TString calib_file_name = Form("%s/cdb/calibrations24/scales/cdb_centrality_scale_%d.root", env_p, runnumber);

  CDBTTree *cdbttree = new CDBTTree(calib_file_name.Data());

  cdbttree->LoadCalibrations();

  cout << cdbttree->GetDoubleValue(0, "centralityscale") << endl;
  

  cdbttree->Print();
  delete cdbttree;
  gSystem->Exit(0);
}

// if you get errors like
// Error inserting payload test.root, msg: "DataBaseException: Global Tag is locked."
// use CDBUtils to unlock your global tag (unlockGlobalTag("<global tag>")

void CDB_CentralityScaleCDBInsert(const int runnumber)
{
  char *env_p = new char[200];
  sprintf(env_p,"%s",std::getenv("MBD_CENTRALITY_CALIB_PATH"));

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }

  TString calib_file_name = Form("%s/cdb/calibrations24/scales/cdb_centrality_scale_%d.root", env_p, runnumber);

  recoConsts *rc = recoConsts::instance();
// please choose a unique name, if it is your username it's easier to see who created it
  rc->set_StringFlag("CDB_GLOBALTAG","2024p007");
  CDBUtils *cdb = new CDBUtils(rc->get_StringFlag("CDB_GLOBALTAG"));
  cdb->createPayloadType("CentralityScale");

  cdb->insertPayload("CentralityScale",calib_file_name.Data(),runnumber);

  cout << cdb->getUrl("CentralityScale", runnumber)<<endl;;

  return;
}

void CDB_CentralityScaleCDBRead(const int runnumber)
{
  recoConsts *rc = recoConsts::instance();
// please choose a unique name, if it is your username it's easier to see who created it
  rc->set_StringFlag("CDB_GLOBALTAG","2024p007"); 
  rc->set_uint64Flag("TIMESTAMP",runnumber);
// 1000000 is the insert timestamp. Higher timestamps work, lower time stamps do not


  CDBInterface *cdb = CDBInterface::instance();
  cdb->Verbosity(1);
  //CDBInterface *cdb = new CDBInterface(rc->get_StringFlag("CDB_GLOBALTAG"));
  cout << "using insert timestamp to retrieve no end time payload" << endl;
  cout << cdb->getUrl("CentralityScale");
  CDBTTree *cdbttree = new CDBTTree(cdb->getUrl("CentralityScale"));

  cdbttree->LoadCalibrations();
  cdbttree->Print();

  gSystem->Exit(0);
  return;
}


#endif
