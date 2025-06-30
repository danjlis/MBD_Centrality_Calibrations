#ifndef CENTRALITYCDB
#define CENTRALITYCDB

#include <cdbobjects/CDBTTree.h>
#include <sphenixnpc/CDBUtils.h>
#include <ffamodules/CDBInterface.h>

#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libffamodules.so)

R__LOAD_LIBRARY(libsphenixnpc.so)
R__LOAD_LIBRARY(libphool.so)

R__LOAD_LIBRARY(libcdbobjects.so)

int CDB_MakeCentralityCDBTree(const int runnumber)
{

  
  char *env_p = new char[200];
  sprintf(env_p,"%s",std::getenv("MBD_CENTRALITY_CALIB_PATH"));

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return 1;
    }

  TString calib_file_name = Form("%s/calib/mbdana_centrality_bal_%d.root", env_p, runnumber);
  if (runnumber == 0)
    {
      calib_file_name = Form("%s/calib/mbdana_centrality_hijing.root", env_p);
    }
  else if (runnumber == 1)
    {
      calib_file_name = Form("%s/calib/mbdana_centrality_ampt.root", env_p);
    }
  else if (runnumber == 2)
    {
      calib_file_name = Form("%s/calib/mbdana_centrality_epos.root", env_p);
    }
  else if (runnumber == 3)
    {
      calib_file_name = Form("%s/calib/mbdana_centrality_hijing_magoff.root", env_p);
    }
  else if (runnumber == 4)
    {
      calib_file_name = Form("%s/calib/mbdana_centrality_ampt_magoff.root", env_p);
    }
  else if (runnumber == 5)
    {
      calib_file_name = Form("%s/calib/mbdana_centrality_epos_magoff.root", env_p);
    }

  TFile *fcalib = new TFile(calib_file_name.Data(), "r");
  if (!fcalib) 
    {
      cout << " No file " <<endl;
      return 1;
    }
  TNtuple *ts = (TNtuple*) fcalib->Get("tn_centrality");
  if (!ts) 
    {
      cout << " No TNtuple " << endl;
      return 1;
    }

  float lowbin;
  float npart;
  ts->SetBranchAddress("low",&lowbin);
  TString cdb_name = Form("%s/cdb/calibrations24/divs/cdb_centrality_%d.root", env_p, runnumber);
  if (runnumber == 0)
    {
      cdb_name = Form("%s/cdb/calibrations24/divs/cdb_centrality_hijing.root", env_p);
    }
  else if (runnumber == 1)
    {
      cdb_name = Form("%s/cdb/calibrations24/divs/cdb_centrality_ampt.root", env_p);
    }
  else if (runnumber == 2)
    {
      cdb_name = Form("%s/cdb/calibrations24/divs/cdb_centrality_epos.root", env_p);
    }
  else if (runnumber == 3)
    {
      cdb_name = Form("%s/cdb/calibrations24/divs/cdb_centrality_hijing_magoff.root", env_p);
    }
  else if (runnumber == 4)
    {
      cdb_name = Form("%s/cdb/calibrations24/divs/cdb_centrality_ampt_magoff.root", env_p);
    }
  else if (runnumber == 5)
    {
      cdb_name = Form("%s/cdb/calibrations24/divs/cdb_centrality_epos_magoff.root", env_p);
    }

  CDBTTree *cdbttree = new CDBTTree(cdb_name.Data());
  for (int i=0; i<ts->GetEntries(); i++)
  {
    ts->GetEntry(i);

    cdbttree->SetFloatValue(i,"centralitydiv",lowbin);
  }
  cdbttree->Commit();
  cdbttree->Print();
  cdbttree->WriteCDBTTree();
  delete cdbttree;
  return 0;
}
int CDB_MakeCentralityCDBTree(const int runnumber, const int divnumber)
{

  
  char *env_p = new char[200];
  sprintf(env_p,"%s",std::getenv("MBD_CENTRALITY_CALIB_PATH"));

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return 1;
    }

  TString calib_file_name = Form("%s/calib/mbdana_centrality_bal_%d.root", env_p, divnumber);
  if (runnumber == 0)
    {
      calib_file_name = Form("%s/calib/mbdana_centrality_hijing.root", env_p);
    }
  else if (runnumber == 1)
    {
      calib_file_name = Form("%s/calib/mbdana_centrality_ampt.root", env_p);
    }
  else if (runnumber == 2)
    {
      calib_file_name = Form("%s/calib/mbdana_centrality_epos.root", env_p);
    }
  else if (runnumber == 3)
    {
      calib_file_name = Form("%s/calib/mbdana_centrality_hijing_magoff.root", env_p);
    }
  else if (runnumber == 4)
    {
      calib_file_name = Form("%s/calib/mbdana_centrality_ampt_magoff.root", env_p);
    }
  else if (runnumber == 5)
    {
      calib_file_name = Form("%s/calib/mbdana_centrality_epos_magoff.root", env_p);
    }

  TFile *fcalib = new TFile(calib_file_name.Data(), "r");
  if (!fcalib) 
    {
      cout << " No file " <<endl;
      return 1;
    }
  TNtuple *ts = (TNtuple*) fcalib->Get("tn_centrality");
  if (!ts) 
    {
      cout << " No TNtuple " << endl;
      return 1;
    }

  float lowbin;
  float npart;
  ts->SetBranchAddress("low",&lowbin);
  TString cdb_name = Form("%s/cdb/calibrations24/divs/cdb_centrality_%d.root", env_p, runnumber);
  if (runnumber == 0)
    {
      cdb_name = Form("%s/cdb/calibrations24/divs/cdb_centrality_hijing.root", env_p);
    }
  else if (runnumber == 1)
    {
      cdb_name = Form("%s/cdb/calibrations24/divs/cdb_centrality_ampt.root", env_p);
    }
  else if (runnumber == 2)
    {
      cdb_name = Form("%s/cdb/calibrations24/divs/cdb_centrality_epos.root", env_p);
    }
  else if (runnumber == 3)
    {
      cdb_name = Form("%s/cdb/calibrations24/divs/cdb_centrality_hijing_magoff.root", env_p);
    }
  else if (runnumber == 4)
    {
      cdb_name = Form("%s/cdb/calibrations24/divs/cdb_centrality_ampt_magoff.root", env_p);
    }
  else if (runnumber == 5)
    {
      cdb_name = Form("%s/cdb/calibrations24/divs/cdb_centrality_epos_magoff.root", env_p);
    }

  CDBTTree *cdbttree = new CDBTTree(cdb_name.Data());
  for (int i=0; i<ts->GetEntries(); i++)
  {
    ts->GetEntry(i);

    cdbttree->SetFloatValue(i,"centralitydiv",lowbin);
  }
  cdbttree->Commit();
  cdbttree->Print();
  cdbttree->WriteCDBTTree();
  delete cdbttree;
  return 0;
}

void CDB_ReadCentralityTree(const int runnumber)
{
  char *env_p = new char[200];
  sprintf(env_p,"%s",std::getenv("MBD_CENTRALITY_CALIB_PATH"));

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }

  TString calib_file_name = Form("%s/cdb/calibrations24/divs/cdb_centrality_%d.root", env_p, runnumber);

  CDBTTree *cdbttree = new CDBTTree(calib_file_name.Data());

  cdbttree->LoadCalibrations();
  for (int i=0; i<20; i++)
    {
      cout << cdbttree->GetFloatValue(i, "centralitydiv") << endl;
    }

  cdbttree->Print();
  delete cdbttree;
  gSystem->Exit(0);
}

// if you get errors like
// Error inserting payload test.root, msg: "DataBaseException: Global Tag is locked."
// use CDBUtils to unlock your global tag (unlockGlobalTag("<global tag>")

void CDB_CentralityCDBInsert(const int runnumber)
{
  char *env_p = new char[200];
  sprintf(env_p,"%s",std::getenv("MBD_CENTRALITY_CALIB_PATH"));

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }

  TString calib_file_name = Form("%s/cdb/calibrations24/divs/cdb_centrality_%d.root", env_p, runnumber);

  recoConsts *rc = recoConsts::instance();
// please choose a unique name, if it is your username it's easier to see who created it
  rc->set_StringFlag("CDB_GLOBALTAG","2024p007");
  CDBUtils *cdb = new CDBUtils(rc->get_StringFlag("CDB_GLOBALTAG"));
  cdb->createPayloadType("Centrality");

  cdb->insertPayload("Centrality",calib_file_name.Data(),10);

  cout << cdb->getUrl("Centrality", 10)<<endl;;

  return;
}

void CDB_CentralityCDBRead(const int runnumber)
{
  recoConsts *rc = recoConsts::instance();
// please choose a unique name, if it is your username it's easier to see who created it
  rc->set_StringFlag("CDB_GLOBALTAG","2024p007"); 
  rc->set_uint64Flag("TIMESTAMP",6);
// 1000000 is the insert timestamp. Higher timestamps work, lower time stamps do not


  CDBInterface *cdb = CDBInterface::instance();
  //CDBInterface *cdb = new CDBInterface(rc->get_StringFlag("CDB_GLOBALTAG"));
  cout << "using insert timestamp to retrieve no end time payload" << endl;
  rc->set_uint64Flag("TIMESTAMP",runnumber);
  cout << cdb->getUrl("Centrality");
  CDBTTree *cdbttree = new CDBTTree(cdb->getUrl("Centrality"));

  cdbttree->LoadCalibrations();
  cdbttree->Print();

  gSystem->Exit(0);
  return;
}


#endif
