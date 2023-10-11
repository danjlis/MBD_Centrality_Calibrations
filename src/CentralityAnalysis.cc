#include "CentralityAnalysis.h"

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>
#include <bbc/BbcDefs.h>
#include <bbc/BbcPmtInfoV1.h>
#include <bbc/BbcPmtInfoContainerV1.h>
#include <bbc/BbcOutV1.h>
#include <centrality/CentralityInfov2.h>
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

CentralityAnalysis::CentralityAnalysis(const std::string &name, const std::string &hist_name, const std::string &tree_name)
  : SubsysReco(name)
{
  _zdc_gain_factors[0] = 1.37;
  _zdc_gain_factors[1] = 0.64;
  _zdc_gain_factors[2] = 0.44;
  _zdc_gain_factors[3] = 1.39;
  _zdc_gain_factors[4] = 0.78;
  _zdc_gain_factors[5] = 0.29;

  const int centrality_map[20] = {1999, 1499, 1291, 1102, 937, 790, 660, 547, 449, 363, 289, 227, 174, 130, 94, 66, 45, 0, 0, 0};
  for (int i = 0; i < 20; i++)
  {
    _centrality_map[i] = centrality_map[i];
  }

  _offset = 28.52;

  const float mbd_vertex_cuts[5] = {50, 30, 20, 10, 5};
  for (int i = 0; i < 5; i++)
  {
    _mbd_vertex_cuts[i] = mbd_vertex_cuts[i];
  }

  _hist_filename = hist_name;
  _tree_filename = tree_name;
}

CentralityAnalysis::~CentralityAnalysis()
{
  delete hm;
}
int CentralityAnalysis::Init(PHCompositeNode * /*unused*/)
{
  // Histograms

  outfile = new TFile(_tree_filename.c_str(), "RECREATE");
  // ttree
  ttree = new TTree("T", "a perservering date tree");

  ttree->Branch("z_vertex", &_z_vertex);
  ttree->Branch("isMinBias", &_isMinBias);
  ttree->Branch("centile", &_centile);
  ttree->Branch("mbd_charge_sum", &_mbd_charge_sum);
  ttree->Branch("mbd_charge_sum_n", &_mbd_charge_sum_n);
  ttree->Branch("mbd_charge_sum_s", &_mbd_charge_sum_s);

  ttree->Branch("zdc_energy_sum", &_zdc_energy_sum);
  ttree->Branch("zdc_energy_sum_n", &_zdc_energy_sum_n);
  ttree->Branch("zdc_energy_sum_s", &_zdc_energy_sum_s);

  ttree->Branch("mbd_charge", m_mbd_charge, "mbd_charge[128]/F");
  ttree->Branch("mbd_time_t", m_mbd_time_t, "mbd_time_t[128]/F");
  ttree->Branch("mbd_time_q", m_mbd_time_q, "mbd_time_q[128]/F");
  ttree->Branch("mbd_side", m_mbd_side, "mbd_side[128]/I");

  ttree->Branch("zdc_energy_low", m_zdc_energy_low, "zdc_energy_low[6]/F");
  ttree->Branch("zdc_energy_high", m_zdc_energy_high, "zdc_energy_high[6]/F");
  ttree->Branch("zdc_sum_low", m_zdc_sum_low, "zdc_sum_low[2]/F");
  ttree->Branch("zdc_sum_high", m_zdc_sum_high, "zdc_sum_high[2]/F");

  // histograms

  hm = new Fun4AllHistoManager("CENTRALITY_HIST");

  h_mbd_vertex = new TH1D("h_mbd_vertex", "", 2000, -500, 500);
  hm->registerHisto(h_mbd_vertex);

  h_mbd_vertex_w_zdc_cut = new TH1D("h_mbd_vertex_zdc_cut", "", 2000, -500, 500);
  hm->registerHisto(h_mbd_vertex_w_zdc_cut);

  h_mbd_charge_ns = new TH1D("h_mbd_charge_ns", "", 2500, 0, 2500);
  h_mbd_charge_n = new TH1D("h_mbd_charge_n", "", 2500, 0, 2500);
  h_mbd_charge_s = new TH1D("h_mbd_charge_s", "", 2500, 0, 2500);
  hm->registerHisto(h_mbd_charge_ns);
  hm->registerHisto(h_mbd_charge_n);
  hm->registerHisto(h_mbd_charge_s);

  h_mbd_charge_ns_w_zdc_cut = new TH1D("h_mbd_charge_ns_w_zdc_cut", "", 2500, 0, 2500);
  h_mbd_charge_n_w_zdc_cut = new TH1D("h_mbd_charge_n_w_zdc_cut", "", 2500, 0, 2500);
  h_mbd_charge_s_w_zdc_cut = new TH1D("h_mbd_charge_s_w_zdc_cut", "", 2500, 0, 2500);
  hm->registerHisto(h_mbd_charge_ns_w_zdc_cut);
  hm->registerHisto(h_mbd_charge_n_w_zdc_cut);
  hm->registerHisto(h_mbd_charge_s_w_zdc_cut);

  h_mbd_charge_ns_w_zdc_cut_w_mbd_cut = new TH1D("h_mbd_charge_ns_w_zdc_cut_w_mbd_cut", "", 2500, 0, 2500);
  h_mbd_charge_n_w_zdc_cut_w_mbd_cut = new TH1D("h_mbd_charge_n_w_zdc_cut_w_mbd_cut", "", 2500, 0, 2500);
  h_mbd_charge_s_w_zdc_cut_w_mbd_cut = new TH1D("h_mbd_charge_s_w_zdc_cut_w_mbd_cut", "", 2500, 0, 2500);
  hm->registerHisto(h_mbd_charge_ns_w_zdc_cut_w_mbd_cut);
  hm->registerHisto(h_mbd_charge_n_w_zdc_cut_w_mbd_cut);
  hm->registerHisto(h_mbd_charge_s_w_zdc_cut_w_mbd_cut);
  for (int i = 0; i < 5; i++)
    {
      h_mbd_charge_ns_w_zdc_cut_w_mbd_cut_and_vertex[i] = new TH1D(Form("h_mbd_charge_ns_w_zdc_cut_w_mbd_cut_and_vertex_%d", static_cast<int>(_mbd_vertex_cuts[i])), "", 2500, 0, 2500);
      h_mbd_charge_n_w_zdc_cut_w_mbd_cut_and_vertex[i] = new TH1D(Form("h_mbd_charge_n_w_zdc_cut_w_mbd_cut_and_vertex_%d", static_cast<int>(_mbd_vertex_cuts[i])), "", 2500, 0, 2500);
      h_mbd_charge_s_w_zdc_cut_w_mbd_cut_and_vertex[i] = new TH1D(Form("h_mbd_charge_s_w_zdc_cut_w_mbd_cut_and_vertex_%d", static_cast<int>(_mbd_vertex_cuts[i])), "", 2500, 0, 2500);
      hm->registerHisto(h_mbd_charge_ns_w_zdc_cut_w_mbd_cut_and_vertex[i]);
      hm->registerHisto(h_mbd_charge_n_w_zdc_cut_w_mbd_cut_and_vertex[i]);
      hm->registerHisto(h_mbd_charge_s_w_zdc_cut_w_mbd_cut_and_vertex[i]);
    }
  for (int i = 0; i < 3; i++)
    {
      h_mbd_ring_charge_sum_n[i] = new TH1D(Form("h_mbd_ring_charge_sum_n_%d", i), "", 2500, 0, 2500);
      h_mbd_ring_charge_sum_s[i] = new TH1D(Form("h_mbd_ring_charge_sum_s_%d", i), "", 2500, 0, 2500);
      hm->registerHisto(h_mbd_ring_charge_sum_n[i]);
      hm->registerHisto(h_mbd_ring_charge_sum_s[i]);

      h_mbd_ring_charge_sum_n_w_zdc_cut[i] = new TH1D(Form("h_mbd_ring_charge_sum_n_w_zdc_cut_%d", i), "", 2500, 0, 2500);
      h_mbd_ring_charge_sum_s_w_zdc_cut[i] = new TH1D(Form("h_mbd_ring_charge_sum_s_w_zdc_cut_%d", i), "", 2500, 0, 2500);
      hm->registerHisto(h_mbd_ring_charge_sum_n_w_zdc_cut[i]);
      hm->registerHisto(h_mbd_ring_charge_sum_s_w_zdc_cut[i]);

      h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut[i] = new TH1D(Form("h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut_%d", i), "", 2500, 0, 2500);
      h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut[i] = new TH1D(Form("h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut_%d", i), "", 2500, 0, 2500);
      hm->registerHisto(h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut[i]);
      hm->registerHisto(h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut[i]);
      for (int j = 0; j < 5; j++)
	{
	  h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut_and_vertex[i][j] = new TH1D(Form("h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut_and_vertex_%d_%d", static_cast<int>(_mbd_vertex_cuts[j]), i), "", 2500, 0, 2500);
	  h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut_and_vertex[i][j] = new TH1D(Form("h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut_and_vertex_%d_%d", static_cast<int>(static_cast<int>(_mbd_vertex_cuts[j])), i), "", 2500, 0, 2500);
	  hm->registerHisto(h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut_and_vertex[i][j]);
	  hm->registerHisto(h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut_and_vertex[i][j]);
	}
    }

  // zdc
  h_zdc_energy_ns = new TH1D("h_zdc_energy_ns", "", 2000, 0, 7000);
  h_zdc_energy_n = new TH1D("h_zdc_energy_n", "", 2000, 0, 7000);
  h_zdc_energy_s = new TH1D("h_zdc_energy_s", "", 2000, 0, 7000);
  hm->registerHisto(h_zdc_energy_ns);
  hm->registerHisto(h_zdc_energy_n);
  hm->registerHisto(h_zdc_energy_s);

  // correlations
  h_zdc_mbd_corr_ns = new TH2D("h_zdc_mbd_corr_ns", "", 100, 0, 2500, 100, 0, 10000);
  h_zdc_mbd_corr_n = new TH2D("h_zdc_mbd_corr_n", "", 100, 0, 2500, 100, 0, 10000);
  h_zdc_mbd_corr_s = new TH2D("h_zdc_mbd_corr_s", "", 100, 0, 2500, 100, 0, 10000);
  hm->registerHisto(h_zdc_mbd_corr_ns);
  hm->registerHisto(h_zdc_mbd_corr_n);
  hm->registerHisto(h_zdc_mbd_corr_s);

  h_zdc_mbd_corr_ns_w_zdc_cut = new TH2D("h_zdc_mbd_corr_ns_w_zdc_cut", "", 100, 0, 2500, 100, 0, 10000);
  h_zdc_mbd_corr_n_w_zdc_cut = new TH2D("h_zdc_mbd_corr_n_w_zdc_cut", "", 100, 0, 2500, 100, 0, 10000);
  h_zdc_mbd_corr_s_w_zdc_cut = new TH2D("h_zdc_mbd_corr_s_w_zdc_cut", "", 100, 0, 2500, 100, 0, 10000);
  hm->registerHisto(h_zdc_mbd_corr_ns_w_zdc_cut);
  hm->registerHisto(h_zdc_mbd_corr_n_w_zdc_cut);
  hm->registerHisto(h_zdc_mbd_corr_s_w_zdc_cut);

  h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut = new TH2D("h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut", "", 100, 0, 2500, 100, 0, 10000);
  h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut = new TH2D("h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut", "", 100, 0, 2500, 100, 0, 10000);
  h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut = new TH2D("h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut", "", 100, 0, 2500, 100, 0, 10000);
  hm->registerHisto(h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut);
  hm->registerHisto(h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut);
  hm->registerHisto(h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut);

  for (int i = 0; i < 5; i++)
    {
      h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut_and_vertex[i] = new TH2D(Form("h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut_and_vertex_%d", static_cast<int>(_mbd_vertex_cuts[i])), "", 100, 0, 2500, 100, 0, 10000);
      h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut_and_vertex[i] = new TH2D(Form("h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut_and_vertex_%d", static_cast<int>(_mbd_vertex_cuts[i])), "", 100, 0, 2500, 100, 0, 10000);
      h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut_and_vertex[i] = new TH2D(Form("h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut_and_vertex_%d", static_cast<int>(_mbd_vertex_cuts[i])), "", 100, 0, 2500, 100, 0, 10000);
      hm->registerHisto(h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut_and_vertex[i]);
      hm->registerHisto(h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut_and_vertex[i]);
      hm->registerHisto(h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut_and_vertex[i]);
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityAnalysis::InitRun(PHCompositeNode * /*unused*/)
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }

  return 0;
}

void CentralityAnalysis::ResetVars()
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
  _z_vertex = -999;
  _isMinBias = false;
  _centile = -999;
  _mbd_charge_sum = 0.;
  _mbd_charge_sum_n = 0.;
  _mbd_charge_sum_s = 0.;
  _zdc_energy_sum = 0.;
  _zdc_energy_sum_n = 0.;
  _zdc_energy_sum_s = 0.;
  for (int i = 0; i < 3; i++)
  {
    _mbd_ring_charge_sum_n[i] = 0.;
    _mbd_ring_charge_sum_s[i] = 0.;
  }
  for (int i = 0; i < 128; i++)
  {
    m_mbd_charge[i] = 0.;
    m_mbd_time_t[i] = 0.;
    m_mbd_time_q[i] = 0.;
    m_mbd_ipmt[i] = 0;
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

int CentralityAnalysis::FillVars()
{

  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
  
  for ( int ipmt = 0; ipmt < BbcDefs::BBC_N_PMT; ipmt++)
    {
      _tmp_pmt = _pmts_mbd->get_tower_at_channel(ipmt);
      _side = ipmt%64;
      
      _charge = _tmp_pmt->get_q();
      _time_t = _tmp_pmt->get_tt();
      _time_q = _tmp_pmt->get_tq();

      m_mbd_charge[ipmt] = _charge;
      m_mbd_time_q[ipmt] = _time_q;
      m_mbd_time_t[ipmt] = _time_t;
      m_mbd_side[ipmt] = _side;
      m_mbd_ipmt[ipmt] = ipmt;
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

  // centrality info

  _isMinBias = _central->isMinBias();

  _centile = (_central->has_centile(CentralityInfo::PROP::mbd_NS)?_central->get_centile(CentralityInfo::PROP::mbd_NS) : -999.99);

  // bbc vertex

  _z_vertex = _bbc_out->get_zvtx();
  _tubes_hit_s = _bbc_out->get_npmt(Bbc::South);
  _tubes_hit_n = _bbc_out->get_npmt(Bbc::North);

  _mbd_charge_sum_s = _bbc_out->get_q(Bbc::South);
  _mbd_charge_sum_n = _bbc_out->get_q(Bbc::North);

  _mbd_charge_sum = _mbd_charge_sum_n + _mbd_charge_sum_s;

  _mbd_time_s = _bbc_out->get_time(Bbc::South);
  _mbd_time_n = _bbc_out->get_time(Bbc::North);

  _mbd_time_0 = _bbc_out->get_t0();


  if (Verbosity() >= 10)
  {
    std::cout << "--------- ZDC data: ----------" << std::endl;
    std::cout << "South:" << std::endl;
    for (int i = 0; i < 3; i++)
    {
      std::cout << i << " : " << m_zdc_energy_low[i] << " (" << m_zdc_energy_high[i] << ") " << std::endl;
    }
    std::cout << "Sum : " << m_zdc_sum_low[0] << " (" << m_zdc_sum_high[0] << ") " << std::endl;
    std::cout << "North:" << std::endl;
    for (int i = 0; i < 3; i++)
    {
      std::cout << i << " : " << m_zdc_energy_low[i + 3] << " (" << m_zdc_energy_high[i + 3] << ") " << std::endl;
    }
    std::cout << "Sum : " << m_zdc_sum_low[1] << " (" << m_zdc_sum_high[1] << ") " << std::endl;

    std::cout << "--------- MBD data: ----------" << std::endl;
    std::cout << "South:" << std::endl;
    for (int i = 0; i < 64; i++)
    {
      std::cout << m_mbd_ipmt[i] << " : " << m_mbd_charge[i] << "  " << m_mbd_time_t[i] <<" "<<m_mbd_time_q[i]<< ") " << std::endl;
    }
    std::cout << "North:" << std::endl;
    for (int i = 64; i < 128; i++)
    {
      std::cout << m_mbd_ipmt[i] << " : " << m_mbd_charge[i] << "  " << m_mbd_time_t[i] <<" "<<m_mbd_time_q[i]<< ") " << std::endl;
    }
    std::cout << "--------- CentralityInfo: ----------" << std::endl;
    std::cout << "   IsMinBias : "<< _isMinBias << std::endl;
    std::cout << "   Centile  : "<< _centile << std::endl;

    std::cout << "--------- BbcOut: ----------" << std::endl;
    std::cout << "   z-vertex   : "<< _z_vertex << std::endl;
    std::cout << "  nPmt S -- N : "<< _tubes_hit_s <<" -- "<< _tubes_hit_n << std::endl;
    std::cout << "  Charge Sum  : "<< _mbd_charge_sum_s <<"(S) + "<< _mbd_charge_sum_n <<" (N) = "<< _mbd_charge_sum <<  std::endl;
    std::cout << "  Time 0      : "<< _mbd_time_s <<"(S) + "<< _mbd_time_n <<" (N) = "<< _mbd_time_0 <<  std::endl;

  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void CentralityAnalysis::FillHistograms()
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
  h_mbd_vertex->Fill(_z_vertex);

  h_mbd_charge_ns->Fill(_mbd_charge_sum);
  h_mbd_charge_n->Fill(_mbd_charge_sum_n);
  h_mbd_charge_s->Fill(_mbd_charge_sum_s);

  for (int i = 0; i < 3; i++)
  {
    h_mbd_ring_charge_sum_n[i]->Fill(_mbd_ring_charge_sum_n[i]);
    h_mbd_ring_charge_sum_s[i]->Fill(_mbd_ring_charge_sum_s[i]);
  }
  // zdc
  h_zdc_energy_ns->Fill(_zdc_energy_sum);
  h_zdc_energy_n->Fill(_zdc_energy_sum_n);
  h_zdc_energy_s->Fill(_zdc_energy_sum_s);

  // correlations
  h_zdc_mbd_corr_ns->Fill(_mbd_charge_sum, _zdc_energy_sum);
  h_zdc_mbd_corr_n->Fill(_mbd_charge_sum_n, _zdc_energy_sum_n);
  h_zdc_mbd_corr_s->Fill(_mbd_charge_sum_s, _zdc_energy_sum_s);

  if (!_zdc_check)
  {
    return;
  }

  h_mbd_vertex_w_zdc_cut->Fill(_z_vertex);

  h_mbd_charge_ns_w_zdc_cut->Fill(_mbd_charge_sum);
  h_mbd_charge_n_w_zdc_cut->Fill(_mbd_charge_sum_n);
  h_mbd_charge_s_w_zdc_cut->Fill(_mbd_charge_sum_s);

  for (int i = 0; i < 3; i++)
  {
    h_mbd_ring_charge_sum_n_w_zdc_cut[i]->Fill(_mbd_ring_charge_sum_n[i]);
    h_mbd_ring_charge_sum_s_w_zdc_cut[i]->Fill(_mbd_ring_charge_sum_s[i]);
  }

  // correlations
  h_zdc_mbd_corr_ns_w_zdc_cut->Fill(_mbd_charge_sum, _zdc_energy_sum);
  h_zdc_mbd_corr_n_w_zdc_cut->Fill(_mbd_charge_sum_n, _zdc_energy_sum_n);
  h_zdc_mbd_corr_s_w_zdc_cut->Fill(_mbd_charge_sum_s, _zdc_energy_sum_s);
  if (Verbosity())
  {
    std::cout << "tubes hit in N/S : " << _tubes_hit_n << " / " << _tubes_hit_s << std::endl;
  }
  if (_tubes_hit_s < 2 || _tubes_hit_n < 2)
  {
    return;
  }

  h_mbd_charge_ns_w_zdc_cut_w_mbd_cut->Fill(_mbd_charge_sum);
  h_mbd_charge_n_w_zdc_cut_w_mbd_cut->Fill(_mbd_charge_sum_n);
  h_mbd_charge_s_w_zdc_cut_w_mbd_cut->Fill(_mbd_charge_sum_s);

  for (int i = 0; i < 3; i++)
  {
    h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut[i]->Fill(_mbd_ring_charge_sum_n[i]);
    h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut[i]->Fill(_mbd_ring_charge_sum_s[i]);
  }

  // correlations
  h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut->Fill(_mbd_charge_sum, _zdc_energy_sum);
  h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut->Fill(_mbd_charge_sum_n, _zdc_energy_sum_n);
  h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut->Fill(_mbd_charge_sum_s, _zdc_energy_sum_s);

  for (int j = 0; j < 5; j++)
  {
    if (std::abs(_z_vertex) > _mbd_vertex_cuts[j])
    {
      continue;
    }
    h_mbd_charge_ns_w_zdc_cut_w_mbd_cut_and_vertex[j]->Fill(_mbd_charge_sum);
    h_mbd_charge_n_w_zdc_cut_w_mbd_cut_and_vertex[j]->Fill(_mbd_charge_sum_n);
    h_mbd_charge_s_w_zdc_cut_w_mbd_cut_and_vertex[j]->Fill(_mbd_charge_sum_s);

    for (int i = 0; i < 3; i++)
    {
      h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut_and_vertex[i][j]->Fill(_mbd_ring_charge_sum_n[i]);
      h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut_and_vertex[i][j]->Fill(_mbd_ring_charge_sum_s[i]);
    }

    // correlations
    h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut_and_vertex[j]->Fill(_mbd_charge_sum, _zdc_energy_sum);
    h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut_and_vertex[j]->Fill(_mbd_charge_sum_n, _zdc_energy_sum_n);
    h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut_and_vertex[j]->Fill(_mbd_charge_sum_s, _zdc_energy_sum_s);
  }

  return;
}

int CentralityAnalysis::process_event(PHCompositeNode *topNode)
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

  FillHistograms();
  ttree->Fill();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityAnalysis::GetNodes(PHCompositeNode *topNode)
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << " :: " << __LINE__ << std::endl;
  }
  
  _pmts_mbd = findNode::getClass<BbcPmtInfoContainerV1>(topNode, "BbcPmtInfoContainer");
  
  if (!_pmts_mbd)
    {
      std::cout << "no mbd pmts node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  _bbc_out = findNode::getClass<BbcOutV1>(topNode, "BbcOut");
  
  if (!_bbc_out)
    {
      std::cout << "no mbd out node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  _towers_zdc = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_ZDC");

  if (!_towers_zdc)
  {
    std::cout << "no zdc towers node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _central = findNode::getClass<CentralityInfov2>(topNode, "CentralityInfo");

  if (!_central)
  {
    std::cout << "no centrality node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityAnalysis::End(PHCompositeNode * /* topNode*/)
{
  outfile->Write();
  outfile->Close();
  
  hm->dumpHistos(_hist_filename.c_str(), "RECREATE");

  return 0;
}
