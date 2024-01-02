#include "InttMbd.h"

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

#include <ffarawobjects/InttRawHit.h>
#include <ffarawobjects/InttRawHitContainer.h>

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

InttMbd::InttMbd(const std::string &name, const std::string &tree_name)
  : SubsysReco(name)
{
  _tree_filename = tree_name;
}

InttMbd::~InttMbd()
{

}
int InttMbd::Init(PHCompositeNode * /*unused*/)
{
  // Histograms

  outfile = new TFile(_tree_filename.c_str(), "RECREATE");
  // ttree
  ttree = new TTree("T", "a perservering date tree");

  ttree->Branch("mbd_fem_clock", &m_mbd_fem_clock, "mbd_fem_clock/s");
  ttree->Branch("mbd_xmit_clock", &m_mbd_fem_clock, "mbd_fem_clock/s");
  ttree->Branch("mbd_event", &m_mbd_event, "mbd_event/i");
  ttree->Branch("int_bco", &m_intt_bco, "intt_bco/l");

  return Fun4AllReturnCodes::EVENT_OK;

}

int InttMbd::InitRun(PHCompositeNode * /*unused*/)
{

  return 0;
}

void InttMbd::ResetVars()
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }

  m_mbd_fem_clock = std::numeric_limits<unsigned short>::signaling_NaN();
  m_mbd_xmit_clock = std::numeric_limits<unsigned short>::signaling_NaN();
  m_mbd_event = std::numeric_limits<int>::signaling_NaN();
  
  m_intt_bco = std::numeric_limits<unsigned long int>::signaling_NaN();

  return;
}

int InttMbd::FillVars()
{

  
  m_mbd_fem_clock = _mbd_out->get_femclock();;
  m_mbd_xmit_clock = _mbd_out->get_clock();
  m_mbd_event = _mbd_out->get_evt();
  

  for (unsigned int i =0 ; i< _intt_hit_container->get_nhits();i++)
    {

      _intt_hit= _intt_hit_container->get_hit(i);
      if (i == 0) m_intt_bco = _intt_hit->get_bco();
      else if (m_intt_bco != _intt_hit->get_bco())
	{
	  m_intt_bco = std::numeric_limits<unsigned long int>::signaling_NaN();
	  break;
	}
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttMbd::process_event(PHCompositeNode *topNode)
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

int InttMbd::GetNodes(PHCompositeNode *topNode)
{
  if (Verbosity())
    {
      std::cout << __FILE__ << " :: " << __FUNCTION__ << " :: " << __LINE__ << std::endl;
    }
  

  _mbd_out = findNode::getClass<MbdOut>(topNode, "MbdOut");
  
  if (!_mbd_out)
    {
      std::cout << "no mbd pmts node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  _intt_hit_container = findNode::getClass<InttRawHitContainer>(topNode, "InttRawHitContainer");
  
  if (!_intt_hit_container)
    {
      std::cout << "no intt node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }


  return Fun4AllReturnCodes::EVENT_OK;
}

int InttMbd::End(PHCompositeNode * /* topNode*/)
{
  outfile->Write();
  outfile->Close();
  

  return 0;
}
