#ifndef CENTRALITY_QACENTRALITYRECO_H
#define CENTRALITY_QACENTRALITYRECO_H

#include <fun4all/SubsysReco.h>

#include <limits>

// Forward declarations
class PHCompositeNode;
class TProfile;
class TFile;
class TNtuple;
class TTree;
class TowerInfo;
class TowerInfoContainer;
class TH1;
class TH2;
class MbdCalibrationAnalysis : public SubsysReco
{
 public:
  //! constructor
  explicit MbdCalibrationAnalysis(const std::string &name = "MbdCalibrationsAnalysis", const std::string &tree_name = "mbd_calib.root");

  //! destructor
  virtual ~MbdCalibrationAnalysis();

  //! full initialization
  int Init(PHCompositeNode *);
  int InitRun(PHCompositeNode *);
  int FillVars();
  int GetNodes(PHCompositeNode *);
  //! event processing method
  int process_event(PHCompositeNode *);
  
  //! end of run method
  int End(PHCompositeNode *);

  void ResetVars();

 protected:

  enum Bbc
  {
    South = 0,
    North = 1
  };


  float zdc_factors[6] = {};

  TFile *outfile = nullptr;
  TTree *ttree = nullptr;
  TowerInfo *_tmp_tower = nullptr;  
  TowerInfoContainer *_pmts_mbd = nullptr;
  TowerInfoContainer *_towers_zdc = nullptr;
  TowerInfoContainer *_towers_zdc_raw = nullptr;
  std::string _tree_filename;

  int _tubes_hit_s = std::numeric_limits<int>::signaling_NaN();;
  int _tubes_hit_n = std::numeric_limits<int>::signaling_NaN();;

  float _energy = std::numeric_limits<float>::signaling_NaN();
  float m_mbd_charge_raw[128]{};
  float m_mbd_time_raw[128]{};
  
  float m_mbd_side[128]{};
  float m_mbd_ipmt[128]{};

  float m_zdc_energy[6]{};
  float m_zdc_sum[2]{};

  float m_zdc_energy_prime[6]{};
  float m_zdc_sum_prime[2]{};

};

#endif
