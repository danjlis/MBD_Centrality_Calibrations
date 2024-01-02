#ifndef MBDANA_H
#define MBDANA_H

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
class MbdOut;
class MbdPmtHit;
class MbdPmtContainer;

class TH1;
class TH2;
class MbdAna : public SubsysReco
{
 public:
  //! constructor
  explicit MbdAna(const std::string &name = "MbdAna", const std::string &tree_name = "mbd_calib.root");

  //! destructor
  virtual ~MbdAna();

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

  bool useZDC;
  TFile *outfile = nullptr;
  TTree *ttree = nullptr;
  TowerInfo *_tmp_tower = nullptr;  
  MbdPmtHit *_tmp_pmt = nullptr;  
  MbdPmtContainer *_pmts_mbd = nullptr;
  TowerInfoContainer *_towers_zdc = nullptr;
  MbdOut *_mbd_out = nullptr;

  std::string _tree_filename;

  int _tubes_hit_s = std::numeric_limits<int>::signaling_NaN();;
  int _tubes_hit_n = std::numeric_limits<int>::signaling_NaN();;

  float _energy = std::numeric_limits<float>::signaling_NaN();
  float m_mbd_charge[128]{};
  float m_mbd_time[128]{};
  float m_mbd_time_t[128]{};
  float m_mbd_time_q[128]{};

  float m_mbd_side[128]{};
  float m_mbd_ipmt[128]{};

  unsigned int m_mbd_vertex_id = std::numeric_limits<unsigned int>::signaling_NaN();
  float m_mbd_vertex = std::numeric_limits<float>::signaling_NaN();
  float m_mbd_vertex_err = std::numeric_limits<float>::signaling_NaN();
  float m_mbd_time_zero = std::numeric_limits<float>::signaling_NaN();
  float m_mbd_time_zero_err = std::numeric_limits<float>::signaling_NaN();

  float m_mbd_charge_sum[2]{};

  float m_zdc_energy[6]{};
  float m_zdc_sum[2]{};

};

#endif
