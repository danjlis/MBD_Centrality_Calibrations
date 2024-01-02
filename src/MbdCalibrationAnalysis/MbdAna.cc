#include "MbdAna.h"

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

#include <mbd/MbdOut.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

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

MbdAna::MbdAna(const std::string &name, const std::string &tree_name)
  : SubsysReco(name)
{
  useZDC = true;
  _tree_filename = tree_name;
}

MbdAna::~MbdAna()
{

}
int MbdAna::Init(PHCompositeNode * /*unused*/)
{
  // Histograms

  outfile = new TFile(_tree_filename.c_str(), "RECREATE");
  // ttree
  ttree = new TTree("T", "a perservering date tree");

  ttree->Branch("mbd_charge", m_mbd_charge, "mbd_charge[128]/F");
  ttree->Branch("mbd_time", m_mbd_time, "mbd_time[128]/F");
  ttree->Branch("mbd_time_t", m_mbd_time_t, "mbd_time_t[128]/F");
  ttree->Branch("mbd_time_q", m_mbd_time_q, "mbd_time_q[128]/F");
  ttree->Branch("mbd_side", m_mbd_side, "mbd_side[128]/I");
  ttree->Branch("mbd_ipmt", m_mbd_ipmt, "mbd_ipmt[128]/I");

  ttree->Branch("mbd_vertex", &m_mbd_vertex, "mbd_vertex/F");
  ttree->Branch("mbd_vertex_err", &m_mbd_vertex_err, "mbd_vertex_err/F");
  ttree->Branch("mbd_time_zero", &m_mbd_time_zero, "mbd_time_zero/F");
  ttree->Branch("mbd_time_zero_err", &m_mbd_time_zero_err, "mbd_time_zero_err/F");

  ttree->Branch("mbd_charge_sum", m_mbd_charge_sum, "mbd_charge[2]/F");

  ttree->Branch("zdc_energy", m_zdc_energy, "zdc_energy[6]/F");
  ttree->Branch("zdc_sum", m_zdc_sum, "zdc_sum[2]/F");

  return Fun4AllReturnCodes::EVENT_OK;

}

int MbdAna::InitRun(PHCompositeNode * /*unused*/)
{

  return 0;
}

void MbdAna::ResetVars()
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
  m_mbd_vertex = 999;
  m_mbd_vertex_err = 999;
  m_mbd_time_zero = 999;
  m_mbd_time_zero_err = 999;

  for (int i = 0; i < 2; i ++)
    m_mbd_charge_sum[i] = 0.;
  for (int i = 0; i < 128; i++)
    {
      m_mbd_charge[i] = 0.;
      m_mbd_time[i] = 0.;
      m_mbd_time_t[i] = 0.;
      m_mbd_time_q[i] = 0.;
      m_mbd_side[i] = -999;
      m_mbd_ipmt[i] = -999;
    }

  for (int i = 0; i < 6; i++)
    {
      m_zdc_energy[i] = 0;
    }
  for (int i = 0; i < 2; i++)
    {
      m_zdc_sum[i]  = 0;
    }
  return;
}

int MbdAna::FillVars()
{

  if (Verbosity())
    {
      std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
    }
  
  int npmt = _pmts_mbd->get_npmt();
  for ( int ipmt = 0; ipmt < npmt; ipmt++)
    {
      _tmp_pmt = _pmts_mbd->get_pmt(ipmt);

      
      m_mbd_charge[ipmt] = _tmp_pmt->get_q();
      m_mbd_time[ipmt] = _tmp_pmt->get_time();
      m_mbd_time_t[ipmt] = _tmp_pmt->get_tt();
      m_mbd_time_q[ipmt] = _tmp_pmt->get_tq();
      m_mbd_side[ipmt] = ipmt/64;
      m_mbd_ipmt[ipmt] = ipmt;
    }


  for (int i = 0 ; i < 2 ; i++) m_mbd_charge_sum[i] = _mbd_out->get_q(i);
  
  std::vector<unsigned int> vtx_id;
  std::vector<float> vtx;
  std::vector<float> vtx_err;
  std::vector<float> time_0;
  std::vector<float> time_0_err;

  
  m_mbd_vertex_id = 1;
  m_mbd_vertex = _mbd_out->get_zvtx();
  m_mbd_vertex_err = _mbd_out->get_zvtxerr();
  m_mbd_time_zero = _mbd_out->get_t0();
  m_mbd_time_zero_err = _mbd_out->get_t0err();

  if (Verbosity())
    {
      std::cout << "MBD stuff:"<<m_mbd_vertex <<" +/- "<<m_mbd_vertex_err<<std::endl;
    }

  if (useZDC){
    int j  = 0;
    for (unsigned int i = 0; i < _towers_zdc->size(); i++)
      {
	_tmp_tower = _towers_zdc->get_tower_at_channel(i);
	_energy = _tmp_tower->get_energy();
	if (_energy!=0) {
	  m_zdc_energy[j] = _energy;
	  m_zdc_sum[j/3] += _energy;
	  j++;
	}
	if (Verbosity())
	  {
	    std::cout << i <<"\t"<<j<<"t"<<_energy<<std::endl; 
	  }
      }
  }
  

  
  return Fun4AllReturnCodes::EVENT_OK;
}

int MbdAna::process_event(PHCompositeNode *topNode)
{
  if (Verbosity())
    {
      std::cout << __FILE__ << " :: " << __FUNCTION__ << " :: " << __LINE__ << std::endl;
    }

  // Get Nodes from the Tree
  if (GetNodes(topNode))
    {
      std::cout << " get nodes didn't wrk " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  // Reset Arrays
  ResetVars();

  // Fill Arrays
  if (FillVars())
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  ttree->Fill();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int MbdAna::GetNodes(PHCompositeNode *topNode)
{
  if (Verbosity())
    {
      std::cout << __FILE__ << " :: " << __FUNCTION__ << " :: " << __LINE__ << std::endl;
    }
  
  _pmts_mbd = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
  
  if (!_pmts_mbd)
    {
      std::cout << "no mbd pmts node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  _mbd_out = findNode::getClass<MbdOut>(topNode, "MbdOut");
  
  if (!_mbd_out)
    {
      std::cout << "no mbd out node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }  

  _towers_zdc = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_ZDC");
  
  if (!_towers_zdc)
    {
      std::cout << "no zdc towers node " << std::endl;

      useZDC = false;
      //      return Fun4AllReturnCodes::ABORTRUN;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MbdAna::End(PHCompositeNode * /* topNode*/)
{
  outfile->Write();
  outfile->Close();
  

  return 0;
}
