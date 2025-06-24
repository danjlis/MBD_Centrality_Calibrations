#include "MbdAna.h"
#include <calotrigger/MinimumBiasInfo.h>
#include <calobase/TowerInfoDefs.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <zdcinfo/Zdcinfo.h>
#include <calotrigger/TriggerAnalyzer.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <ffaobjects/EventHeaderv1.h>

#include <mbd/MbdOut.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <eventplaneinfo/EventplaneinfoMap.h>
#include <eventplaneinfo/Eventplaneinfo.h>
#include <centrality/CentralityInfo.h>

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
  //  triggeranalyzer = new TriggerAnalyzer();
  outfile = new TFile(_tree_filename.c_str(), "RECREATE");
  // ttree
  ttree = new TTree("T", "a perservering date tree");
  ttree->Branch("minbias", &m_minbias, "minbias/I");

  ttree->Branch("centrality", &m_centrality, "centrality/F");
  // ttree->Branch("psi1n", &m_psi1n, "psi1n/F");
  // ttree->Branch("psi1s", &m_psi1s, "psi1s/F");
  // ttree->Branch("psi1", &m_psi1, "psi1/F");
  // ttree->Branch("psi2n", &m_psi2n, "psi2n/F");
  // ttree->Branch("psi2s", &m_psi2s, "psi2s/F");
  // ttree->Branch("psi2", &m_psi2, "psi2/F");
  // ttree->Branch("psi3n", &m_psi3n, "psi3n/F");
  // ttree->Branch("psi3s", &m_psi3s, "psi3s/F");
  // ttree->Branch("psi3", &m_psi3, "psi3/F");

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
  if (m_issim)
    {
      ttree->Branch("bimp",&m_bimp);
    }
  ttree->Branch("gl1_scaled", &m_gl1_scaled, "gl1_scaled/l");
  ttree->Branch("gl1_live", &m_gl1_live, "gl1_live/l");
  ttree->Branch("bunch_number", &m_bunch_number, "bunch_number/l");

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
  m_centrality = 999;

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
  m_bimp = 999;
  return;
}

int MbdAna::FillVars()
{

  if (Verbosity())
    {
      std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
    }
  
  if (m_issim)
    {

      m_bimp = m_eventheader->get_ImpactParameter();
    }


  int npmt = 128;//_pmts_mbd->get_npmt();
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

  if (useZDC)
    {
      for (unsigned int i = 0; i < 2; i++)
	{
	  m_zdc_sum[i] = _zdcinfo->get_zdc_energy(i);
	}
    }
  

  // if (_minimumbiasinfo)
  //   m_minbias = (_minimumbiasinfo->isAuAuMinimumBias()? 1 : 0);


  // auto pMBDS = eventplaneinfomap->get(EventplaneinfoMap::MBDS);
  // auto pMBDN = eventplaneinfomap->get(EventplaneinfoMap::MBDN);

  // //first order stuff
  // std::pair<double, double> mbdsouthQ1;
  // std::pair<double, double> mbdnorthQ1;
  // std::pair<double, double> mbdQ1;
  // mbdsouthQ1 = pMBDS->get_qvector(1);
  // mbdnorthQ1 = pMBDN->get_qvector(1);

  // m_psi1s = GetPsi(mbdsouthQ1.first, mbdsouthQ1.second, 1);
  // m_psi1n = GetPsi(mbdnorthQ1.first, mbdnorthQ1.second, 1);
  // mbdQ1 = std::make_pair(mbdsouthQ1.first + mbdnorthQ1.first, mbdsouthQ1.second + mbdnorthQ1.second);
  // m_psi1 = GetPsi(mbdQ1.first, mbdQ1.second, 1);
   
  // //second order stuff
  // std::pair<double, double> mbdsouthQ2;
  // std::pair<double, double> mbdnorthQ2;
  // std::pair<double, double> mbdQ2;
  // mbdsouthQ2 = pMBDS->get_qvector(2);
  // mbdnorthQ2 = pMBDN->get_qvector(2);

  // m_psi2s = GetPsi(mbdsouthQ2.first, mbdsouthQ2.second, 2);
  // m_psi2n = GetPsi(mbdnorthQ2.first, mbdnorthQ2.second, 2);
  // mbdQ2 = std::make_pair(mbdsouthQ2.first + mbdnorthQ2.first, mbdsouthQ2.second + mbdnorthQ2.second);
  // m_psi2 = GetPsi(mbdQ2.first, mbdQ2.second, 2);

  // //third order stuff
  // std::pair<double, double> mbdsouthQ3;
  // std::pair<double, double> mbdnorthQ3;
  // std::pair<double, double> mbdQ3;
  // mbdsouthQ3 = pMBDS->get_qvector(3);
  // mbdnorthQ3 = pMBDN->get_qvector(3);

  // m_psi3s = GetPsi(mbdsouthQ3.first, mbdsouthQ3.second, 3);
  // m_psi3n = GetPsi(mbdnorthQ3.first, mbdnorthQ3.second, 3);
  // mbdQ3 = std::make_pair(mbdsouthQ3.first + mbdnorthQ3.first, mbdsouthQ3.second + mbdnorthQ3.second);
  // m_psi3 = GetPsi(mbdQ3.first, mbdQ3.second, 3);

  
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

  if (!m_issim)
    {
      m_gl1_scaled = _gl1_packet->getTriggerVector();
      m_gl1_live = _gl1_packet->getLiveVector();
      m_bunch_number = _gl1_packet->getBunchNumber();
      // Fill Arrays
    }
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
  if (!m_issim)
    {
      _gl1_packet = findNode::getClass<Gl1Packet>(topNode, 14001);
  
      if (!_gl1_packet)
	{
	  std::cout << "no gl1packet node " << std::endl;
	  return Fun4AllReturnCodes::ABORTRUN;
	}

    }
  else
    {
      m_eventheader = findNode::getClass<EventHeader>(topNode, "EventHeader");
  
      if (!m_eventheader)
	{
	  std::cout << "no event header " << std::endl;
	  return Fun4AllReturnCodes::ABORTRUN;
	}
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

  _zdcinfo = findNode::getClass<Zdcinfo>(topNode, "Zdcinfo");
  
  if (!_zdcinfo)
    {
      //      std::cout << "no zdc towers node " << std::endl;

      useZDC = false;
      //      return Fun4AllReturnCodes::ABORTRUN;
    }


  // centralityinfo = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
  // eventplaneinfomap = findNode::getClass<EventplaneinfoMap>(topNode, "EventplaneinfoMap");
  // //  _min_bias_info = findNode::getClass<MinimumBiasInfo>(topNode, "MinimumBiasInfo");
  // _minimumbiasinfo = findNode::getClass<MinimumBiasInfo>(topNode, "MinimumBiasInfo");
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int MbdAna::End(PHCompositeNode * /* topNode*/)
{
  outfile->Write();
  outfile->Close();
  

  return 0;
}
double MbdAna::GetPsi(double Qx, double Qy, unsigned int order) 
{
  double temp;
  if ((Qx == 0.0) && (Qy == 0.0))
  {
    temp = NAN;
  }
  else
  {
    temp = atan2(Qy, Qx) / ((double) order);
  }
  return temp;
}
