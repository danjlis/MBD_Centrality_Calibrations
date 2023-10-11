#include "MbdCalibrationAnalysis.h"

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

MbdCalibrationAnalysis::MbdCalibrationAnalysis(const std::string &name, const std::string &tree_name)
  : SubsysReco(name)
{
  _zdc_gain_factors[0] = 1.37;
  _zdc_gain_factors[1] = 0.64;
  _zdc_gain_factors[2] = 0.44;
  _zdc_gain_factors[3] = 1.39;
  _zdc_gain_factors[4] = 0.78;
  _zdc_gain_factors[5] = 0.29;

  _tree_filename = tree_name;
}

MbdCalibrationAnalysis::~MbdCalibrationAnalysis()
{

}
int MbdCalibrationAnalysis::Init(PHCompositeNode * /*unused*/)
{
  // Histograms

  outfile = new TFile(_tree_filename.c_str(), "RECREATE");
  // ttree
  ttree = new TTree("T", "a perservering date tree");

  ttree->Branch("mbd_charge_raw", m_mbd_charge_raw, "mbd_charge_raw[128]/F");
  ttree->Branch("mbd_time_raw", m_mbd_time_raw, "mbd_time_raw[128]/F");
  ttree->Branch("mbd_side", m_mbd_side, "mbd_side[128]/I");
  ttree->Branch("mbd_ipmt", m_mbd_ipmt, "mbd_ipmt[128]/I");

  ttree->Branch("zdc_energy_low", m_zdc_energy_low, "zdc_energy_low[6]/F");
  ttree->Branch("zdc_energy_high", m_zdc_energy_high, "zdc_energy_high[6]/F");
  ttree->Branch("zdc_sum_low", m_zdc_sum_low, "zdc_sum_low[2]/F");
  ttree->Branch("zdc_sum_high", m_zdc_sum_high, "zdc_sum_high[2]/F");

  return Fun4AllReturnCodes::EVENT_OK;

}

int MbdCalibrationAnalysis::InitRun(PHCompositeNode * /*unused*/)
{

  return 0;
}

void MbdCalibrationAnalysis::ResetVars()
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
  for (int i = 0; i < 128; i++)
  {
    m_mbd_charge_raw[i] = 0.;
    m_mbd_time_raw[i] = 0.;
    m_mbd_side[i] = -999;
    m_mbd_ipmt[i] = -999;
  }

  for (int i = 0; i < 6; i++)
  {
    m_zdc_energy_low[i] = 0;
    m_zdc_energy_high[i] = 0;
  }
  for (int i = 0; i < 2; i++)
  {
    m_zdc_sum_low[i] = 0;
    m_zdc_sum_high[i] = 0;
  }
  return;
}

int MbdCalibrationAnalysis::FillVars()
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
	  m_mbd_charge_raw[ipmt] = _energy;
	}
      else
	{

	  m_mbd_time_raw[ipmt] = _energy;
	  m_mbd_side[ipmt] = side;
	  m_mbd_ipmt[ipmt] = ipmt;
	}
    }

  for (unsigned int i = 0; i < _towers_zdc->size(); i++)
  {
    _tmp_tower = _towers_zdc->get_tower_at_channel(i);
    _energy = _tmp_tower->get_energy();

    if (!(i % 2))
    {
      if ((i / 2) % 4 == 3)
      {
        m_zdc_sum_low[i / 8] = _energy;
      }
      else
      {
        m_zdc_energy_low[i / 2] = _zdc_gain_factors[i / 2] * _energy;
      }
    }
    else
    {
      if ((i / 2) % 4 == 3)
      {
        m_zdc_sum_high[i / 8] = _energy;
      }
      else
      {
        m_zdc_energy_high[i / 2] = _zdc_gain_factors[i / 2] * _energy;
      }
    }
  }


  return Fun4AllReturnCodes::EVENT_OK;
}

int MbdCalibrationAnalysis::process_event(PHCompositeNode *topNode)
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

  ttree->Fill();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int MbdCalibrationAnalysis::GetNodes(PHCompositeNode *topNode)
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
  
  _towers_zdc = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_ZDC");
  
  if (!_towers_zdc)
    {
      std::cout << "no zdc towers node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  return Fun4AllReturnCodes::EVENT_OK;  
}

int MbdCalibrationAnalysis::End(PHCompositeNode * /* topNode*/)
{
  outfile->Write();
  outfile->Close();
  

  return 0;
}
