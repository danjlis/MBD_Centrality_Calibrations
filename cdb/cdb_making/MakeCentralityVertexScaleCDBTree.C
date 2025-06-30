#ifndef CENTRALITYVERTEXSCALECDB
#define CENTRALITYVERTEXSCALECDB

#include <cdbobjects/CDBTTree.h>
#include <sphenixnpc/CDBUtils.h>
#include <ffamodules/CDBInterface.h>

#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libffamodules.so)

R__LOAD_LIBRARY(libsphenixnpc.so)
R__LOAD_LIBRARY(libphool.so)

R__LOAD_LIBRARY(libcdbobjects.so)

int CDB_MakeCentralityVertexScaleCDBTree(const int runnumber)
{

  
  char *env_p = new char[200];
  sprintf(env_p,"%s",std::getenv("MBD_CENTRALITY_CALIB_PATH"));

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return 1;
    }

  TString cdb_name = Form("%s/cdb/calibrations24/vertexscales/cdb_centrality_vertex_scale_%d.root", env_p, runnumber);
  if (runnumber == 0)
    {
      cdb_name = Form("%s/cdb/calibrations24/vertexscales/cdb_centrality_vertex_scale_hijing.root", env_p);
    }
  else if (runnumber == 1)
    {
      cdb_name = Form("%s/cdb/calibrations24/vertexscales/cdb_centrality_vertex_scale_ampt.root", env_p);
    }
  else if (runnumber == 2)
    {
      cdb_name = Form("%s/cdb/calibrations24/vertexscales/cdb_centrality_vertex_scale_epos.root", env_p);
    }
  else if (runnumber == 3)
    {
      cdb_name = Form("%s/cdb/calibrations24/vertexscales/cdb_centrality_vertex_scale_hijing_magoff.root", env_p);
    }
  else if (runnumber == 4)
    {
      cdb_name = Form("%s/cdb/calibrations24/vertexscales/cdb_centrality_vertex_scale_ampt_magoff.root", env_p);
    }
  else if (runnumber == 5)
    {
      cdb_name = Form("%s/cdb/calibrations24/vertexscales/cdb_centrality_vertex_scale_epos_magoff.root", env_p);
    }

  CDBTTree *cdbttree = nullptr;
  
  TString calib_file_name = Form("%s/calib/mbd_vertex_scale_%d.root", env_p, runnumber);
  if (runnumber == 0)
    {
      calib_file_name = Form("%s/calib/mbd_vertex_scale_hijing.root", env_p);
    }
  else if (runnumber == 1)
    {
      calib_file_name = Form("%s/calib/mbd_vertex_scale_ampt.root", env_p);
    }
  else if (runnumber == 2)
    {
      calib_file_name = Form("%s/calib/mbd_vertex_scale_epos.root", env_p);
    }
  else if (runnumber == 3)
    {
      calib_file_name = Form("%s/calib/mbd_vertex_scale_hijing_magoff.root", env_p);
    }
  else if (runnumber == 4)
    {
      calib_file_name = Form("%s/calib/mbd_vertex_scale_ampt_magoff.root", env_p);
    }
  else if (runnumber == 5)
    {
      calib_file_name = Form("%s/calib/mbd_vertex_scale_epos_magoff.root", env_p);
    }

  TFile *fcalib = new TFile(calib_file_name.Data(), "r");
  if (!fcalib) 
    {
      cout << " No file " <<endl;
      return 1;
    }
  TNtuple *ts = (TNtuple*) fcalib->Get("tn_vertexscale");
  if (!ts) 
    {
      cout << " No TNtuple " << endl;
      return 1;
    }
  
  float scale;
  float lowvertex;
  float highvertex;
  ts->SetBranchAddress("scale",&scale);
  ts->SetBranchAddress("low_vertex",&lowvertex);
  ts->SetBranchAddress("high_vertex",&highvertex);

  int vertexbins = ts->GetEntries();  
  cdbttree = new CDBTTree(cdb_name.Data());

  cdbttree->SetIntValue(0,"nvertexbins",vertexbins);


  for (int i=0; i<vertexbins; i++)
    {
      ts->GetEntry(i);
      cdbttree->SetDoubleValue(i,"high_vertex",highvertex);
      cdbttree->SetDoubleValue(i,"low_vertex",lowvertex);
      cdbttree->SetDoubleValue(i,"scale",scale);
    }
  fcalib->Close();

  cdbttree->Commit();
  cdbttree->Print();
  cdbttree->WriteCDBTTree();
  delete cdbttree;
  return 0;
}

void CDB_VertexReadTree(const int runnumber)
{
  char *env_p = new char[200];
  sprintf(env_p,"%s",std::getenv("MBD_CENTRALITY_CALIB_PATH"));

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    } 
  TString calib_file_name = Form("%s/cdb/calibrations24/vertexscales/cdb_centrality_vertex_scale_%d.root", env_p, runnumber);

  CDBTTree *cdbttree = new CDBTTree(calib_file_name.Data());

  cdbttree->LoadCalibrations();
  int vertexbins =cdbttree->GetSingleIntValue("nvertexbins");
  cout << cdbttree->GetSingleIntValue("nvertexbins") << endl;
  for (int i=0; i<vertexbins; i++)
    {

      cout << cdbttree->GetDoubleValue(i,"low_vertex") << " " <<
	cdbttree->GetDoubleValue(i,"vertex_vertex") << " " << 
	cdbttree->GetDoubleValue(i,"scale") << std::endl;
    }
  

  cdbttree->Print();
  delete cdbttree;
  gSystem->Exit(0);
}

// if you get errors like
// Error inserting payload test.root, msg: "DataBaseException: Global Tag is locked."
// use CDBUtils to unlock your global tag (unlockGlobalTag("<global tag>")

void CDB_CentralityVertexScaleCDBInsert(const int runnumber)
{
  char *env_p = new char[200];
  sprintf(env_p,"%s",std::getenv("MBD_CENTRALITY_CALIB_PATH"));

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }

  TString calib_file_name = Form("%s/cdb/calibrations24/vertexscales/cdb_centrality_vertex_scale_%d.root", env_p, runnumber);

  recoConsts *rc = recoConsts::instance();
// please choose a unique name, if it is your username it's easier to see who created it
  rc->set_StringFlag("CDB_GLOBALTAG","dlis");
  CDBUtils *cdb = new CDBUtils(rc->get_StringFlag("CDB_GLOBALTAG"));

  cdb->createPayloadType("CentralityVertexScale");

  cdb->insertPayload("CentralityVertexScale",calib_file_name.Data(),runnumber);

  cout << cdb->getUrl("CentralityVertexScale", runnumber)<<endl;;

  return;
}

void CDB_CentralityVertexScaleCDBRead(const int runnumber)
{
  recoConsts *rc = recoConsts::instance();
// please choose a unique name, if it is your username it's easier to see who created it
  rc->set_StringFlag("CDB_GLOBALTAG","dlis"); 
  rc->set_uint64Flag("TIMESTAMP",runnumber);
// 1000000 is the insert timestamp. Higher timestamps work, lower time stamps do not


  CDBInterface *cdb = CDBInterface::instance();
  cdb->Verbosity(1);
  //CDBInterface *cdb = new CDBInterface(rc->get_StringFlag("CDB_GLOBALTAG"));
  cout << "using insert timestamp to retrieve no end time payload" << endl;
  cout << cdb->getUrl("CentralityVertexScale");
  CDBTTree *cdbttree = new CDBTTree(cdb->getUrl("CentralityVertexScale"));

  cdbttree->LoadCalibrations();
  cdbttree->Print();

  gSystem->Exit(0);
  return;
}


#endif
