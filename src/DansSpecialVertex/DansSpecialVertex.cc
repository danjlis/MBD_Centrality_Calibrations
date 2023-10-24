#include "DansSpecialVertex.h"
#include <globalvertex/GlobalVertexv1.h>
#include <globalvertex/GlobalVertexMapv1.h>

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <TF1.h>
#include <TFile.h>
#include <TH2.h>
#include <TNtuple.h>
#include <TSystem.h>

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

DansSpecialVertex::DansSpecialVertex(const std::string &name, const std::string &hist_name)
  : SubsysReco(name)
{
  _hist_filename = hist_name;
}

DansSpecialVertex::~DansSpecialVertex()
{

}
int DansSpecialVertex::Init(PHCompositeNode *topNode)
{
  // Histograms
  hm = new Fun4AllHistoManager("DansVertexHistos");

  h_vertex = new TH1D("h_vertex",";z[cm];counts", 3000, -300, 300);
  h_time_0 = new TH1D("h_time_0",";t0[ns];counts", 500, -25, 25);
  hm->registerHisto(h_vertex);
  hm->registerHisto(h_time_0);

  CreateNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;

}

int DansSpecialVertex::InitRun(PHCompositeNode * )
{

  if (!_runnumber)
    {
      std::cout << "There is no runnumber... exiting"<<std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }


  if (DownloadCalibs())
    {
      std::cout << "Download failed" <<std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  return Fun4AllReturnCodes::EVENT_OK;  ;
}

void DansSpecialVertex::ResetVars()
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
  for (int i = 0; i < 128; i++)
  {
    m_mbd_charge[i] = 0.;
    m_mbd_time[i] = 0.;
    m_mbd_side[i] = -999;
    m_mbd_ipmt[i] = -999;
  }
  return;
}

int DansSpecialVertex::FillVars()
{

  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
  
  for ( int ich = 0; ich < 256; ich++)
    {
      _tmp_tower = _pmts_mbd->get_tower_at_channel(ich);
      short side = ich%128;

      // 0-7 -> Time
      // 8-15 -> Charge
      short ipmt = (ich%16)%8 + (ich/16)*8;
      // 0 -> Time
      // 1 -> Charge
      short type = (ich%16)/8;

      _energy = _tmp_tower->get_energy();
      
      if (type)
	{

	  m_mbd_charge[ipmt] = _energy*charge_calibs[ipmt];
	}
      else
	{
	  m_mbd_time[ipmt] = (25. - (_energy * (9.0 / 5000.))) - time_calibs[ipmt];
	  m_mbd_side[ipmt] = side;
	  m_mbd_ipmt[ipmt] = ipmt;
	}
    }  

  return Fun4AllReturnCodes::EVENT_OK;
}

void DansSpecialVertex::CalculateVertexAndTime()
{

    //// vertex
  hits_n = 0;
  hits_s = 0;      
  hits_n_t = 0;
  hits_s_t = 0;      

  time_sum_n.clear();
  time_sum_s.clear();
  sum_n = 0.;
  sum_s = 0.;
  sum_n2 = 0.;
  sum_s2 = 0.;

  z_vertex = 999;
  time_0 = 999;

  for (int ich = 0; ich < 128; ich++)
    {
      if (m_mbd_charge[ich] > cthresh)
	{
	  // get hit charge channels
	  if (ich/64) hits_n++;
	  else hits_s++;

	  // if bad channel or the time is greater than 10 ns after 0-point don't count the time
	  if (ich == 56 || ich == 56+64 || fabs(m_mbd_time[ich]) > 10) continue;
	  
	  float timme = m_mbd_time[ich];

	  if (ich/64) 
	    {
	      
	      hits_n_t++;
	      time_sum_n.push_back(timme);
	      sum_n += timme;;
	      sum_n2 += (timme*timme);
	      
	    }
	  else
	    {
	      
	      hits_s_t++;
	      time_sum_s.push_back(timme);
	      
	      sum_s += timme;
	      sum_s2 += (timme*timme);
	      
	    }
	  
	}
      
    } 
      
  float mean_north = 999;
  float mean_south = 999;

  // does event meet cut?
  if (hits_s_t >= central_cut){
    
    // get the mean
    mean_south = sum_s/static_cast<float>(hits_s_t);
    
    // get rms
    //    float rms_s = sqrt(sum_s2/static_cast<float>(hits_s_t) - TMath::Power(mean_south, 2));
    int nhit_s_center = 0;
    float sum_s_center = 0.;
    
    // get rid of times outside of RMS*1.5 range
    for (unsigned int is = 0; is < time_sum_s.size(); is++)
      {
	if (fabs(time_sum_s.at(is) - mean_south) < sigma_cut )
	  {
	    sum_s_center += time_sum_s.at(is);
	    nhit_s_center++;
	  }
      }
    
    // Get the mean again
    float mean_south_center = sum_s_center/static_cast<float>(nhit_s_center);
    
    
    mean_south = mean_south_center;
    
  }
  
  else if (hits_s >=2 && (hits_s_t >= 1)){
    
    mean_south = sum_s/static_cast<float>(hits_s_t);
    
  }
  
  
  if (hits_n_t >=central_cut){
    
    mean_north = sum_n/static_cast<float>(hits_n_t);
    
    //    float rms_n = sqrt(sum_n2/static_cast<float>(hits_n_t) - TMath::Power(mean_north, 2));
    int nhit_n_center = 0;
    float sum_n_center = 0.;
    
    for (unsigned int ino = 0; ino < time_sum_n.size(); ino++)
      {
	if (fabs(time_sum_n.at(ino) - mean_north) < sigma_cut )
	  {
	    sum_n_center += time_sum_n.at(ino);
	    nhit_n_center++;
	  }
      }
    
    float mean_north_center = sum_n_center/static_cast<float>(nhit_n_center);
    
    mean_north = mean_north_center;
  }
  
  else if (hits_n >=2 && hits_n_t >= 1){
    
    mean_north = sum_n/static_cast<float>(hits_n_t);
    
  }
  
  
  
  if (mean_north != 999 && mean_south != 9999) {

    z_vertex = 15*(mean_north - mean_south);
    time_0 = (mean_north + mean_south)/2.;
  }
  else 
    {
      z_vertex = 999;
      time_0 = 999;
    }

  h_vertex->Fill(z_vertex);
  h_time_0->Fill(time_0);

  _my_vertex = new GlobalVertexv1();

  _my_vertex->set_x(0.0);
  _my_vertex->set_y(0.0);
  _my_vertex->set_z(z_vertex);
  _my_vertex->set_error(0,0,1.0);
  _my_vertex->set_error(1,1,1.0);
  _my_vertex->set_error(2,2,1.0);

  _my_vertex->set_t(time_0);
  _my_vertex->set_t_err(0.01);


  _my_vertex_map->insert(_my_vertex);

  if (Verbosity()) _my_vertex_map->identify();
  return;
}

int DansSpecialVertex::process_event(PHCompositeNode *topNode)
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << " :: " << __LINE__ << std::endl;
  }

  // Get Nodes from the Tree
  if (GetNodes(topNode))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Reset Arrays
  ResetVars();

  // Fill Arrays
  if (FillVars())
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Fill Arrays
  CalculateVertexAndTime();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int DansSpecialVertex::GetNodes(PHCompositeNode *topNode)
{
  if (Verbosity())
    {
      std::cout << __FILE__ << " :: " << __FUNCTION__ << " :: " << __LINE__ << std::endl;
    }
  
  _pmts_mbd = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_MBD");
  
  if (!_pmts_mbd)
    {
      std::cout << "no mbd pmts node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  _my_vertex_map = findNode::getClass<GlobalVertexMapv1>(topNode, "DansSpecialVertexMap");
  if (!_my_vertex_map)
    {
      std::cout << "no dansspecialvertexmap node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  return Fun4AllReturnCodes::EVENT_OK;  
}


void DansSpecialVertex::CreateNodes(PHCompositeNode *topNode)
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << " :: " << __LINE__ << std::endl;
  }

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing doing nothing" << std::endl;
  }

  PHNodeIterator dstIter(dstNode);

  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "GLOBAL"));
  if (!detNode)
  {
    std::cout << PHWHERE << "Detector Node missing, making one" << std::endl;
    detNode = new PHCompositeNode("GLOBAL");
    dstNode->addNode(detNode);
  }

  GlobalVertexMapv1 *vmap = new GlobalVertexMapv1();
  
  PHIODataNode<PHObject> *dansspecialvertexNode = new PHIODataNode<PHObject>(vmap, "DansSpecialVertexMap", "PHObject");
  detNode->addNode(dansspecialvertexNode);

  return;
}
 

int DansSpecialVertex::DownloadCalibs()
{
  calibfile_charge = new TFile(Form("/sphenix/user/dlis/Projects/centrality/calib/calib_mbd_%d.root", _runnumber), "r");
        
  if (!calibfile_charge)
    {
      std::cout << " No calibration File found " <<std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  
  // 64 gain corrections on each side
  //  0 -  63 -- North
  // 64 - 127 -- South
      
  float g;
  tn_charge = (TNtuple*) calibfile_charge->Get("mbd_calib");
  tn_charge->SetBranchAddress("peak", &g);
      
  for (int i = 0 ; i < 128; i ++)
    {
      tn_charge->GetEntry(i);
      charge_calibs[i] = g;
      if (Verbosity()) std::cout << i << "\t" << g << std::endl;
    }

  calibfile_time = new TFile(Form("/sphenix/user/dlis/Projects/centrality/calib/t0_calib_mbd_%d.root", _runnumber), "r");

  if (!calibfile_time)
    {
      std::cout << " No calibration File found " <<std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  
  // 64 gain corrections on each side
  //  0 -  63 -- North
  // 64 - 127 -- South
 
  float shift;
  float ichannel;
  tn_time = (TNtuple*) calibfile_time->Get("ttree");
  tn_time->SetBranchAddress("shift", &shift);
  tn_time->SetBranchAddress("channel", &ichannel);
        
  for (int i = 0 ; i < 128; i ++)
    {
      time_calibs[i] = -999;
    }
  for (int i = 0 ; i < tn_time->GetEntries(); i ++)
    {
      tn_time->GetEntry(i);
      int ich = static_cast<int>(ichannel);
      time_calibs[ich] = shift;
      if (Verbosity()) std::cout << i << "\t" << shift << std::endl;
    }

  std::cout << "hi" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;

  
}

int DansSpecialVertex::End(PHCompositeNode * /* topNode*/)
{
  hm->dumpHistos(_hist_filename.c_str());  

  return 0;
}
