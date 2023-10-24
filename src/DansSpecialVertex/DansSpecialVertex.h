#ifndef CENTRALITY_QACENTRALITYRECO_H
#define CENTRALITY_QACENTRALITYRECO_H

#include <fun4all/SubsysReco.h>
#include <globalvertex/GlobalVertexv1.h>
#include <globalvertex/GlobalVertexMapv1.h>

#include <limits>
#include <vector>
// Forward declarations
class PHCompositeNode;
class TFile;
class TNtuple;
class TowerInfo;
class TowerInfoContainer;
class GlobalVertexv1;
class GlobalVertexMapv1;
class Fun4AllHistoManager;
class TH1D;

class DansSpecialVertex : public SubsysReco

{
 public:
  //! constructor
  explicit DansSpecialVertex(const std::string &name = "DansSpecialVertex", const std::string &hist_name = "mbd_vertex.root");

  //! destructor
  virtual ~DansSpecialVertex();

  //! full initialization
  int Init(PHCompositeNode *);
  int InitRun(PHCompositeNode *);
  int FillVars();
  void CreateNodes(PHCompositeNode *topNode);
  int GetNodes(PHCompositeNode *);
  //! event processing method
  int process_event(PHCompositeNode *);
  void CalculateVertexAndTime();

  void SetRunNumber(int runnumber) { _runnumber = runnumber; }

  int DownloadCalibs();

  //! end of run method
  int End(PHCompositeNode *);

  void ResetVars();

 protected:

  enum Bbc
  {
    South = 0,
    North = 1
  };

  Fun4AllHistoManager *hm = nullptr;

  GlobalVertexMapv1 *_my_vertex_map = nullptr;
  GlobalVertexv1 *_my_vertex = nullptr;

  TH1D *h_vertex = nullptr;
  TH1D *h_time_0 = nullptr;

  TFile *calibfile_time = nullptr;
  TFile *calibfile_charge = nullptr;

  TNtuple *tn_time = nullptr;
  TNtuple *tn_charge = nullptr;

  TowerInfo *_tmp_tower = nullptr;  
  TowerInfoContainer *_pmts_mbd = nullptr;
  TowerInfoContainer *_towers_zdc = nullptr;

  float m_mbd_charge[128]{};
  float m_mbd_time[128]{};
  
  float m_mbd_side[128]{};
  float m_mbd_ipmt[128]{};

  float charge_calibs[128]{};
  float time_calibs[128]{};

  std::string _hist_filename;

  float _energy = std::numeric_limits<float>::signaling_NaN();

  int _runnumber = std::numeric_limits<int>::signaling_NaN();

  double tthresh = 16;
  double cthresh = 0.4;
  int central_cut = 4;
  float sigma_cut = 1.5;
  
  float z_vertex = std::numeric_limits<float>::signaling_NaN();
  float time_0 = std::numeric_limits<float>::signaling_NaN();

  int hits_n = std::numeric_limits<int>::signaling_NaN();
  int hits_s = std::numeric_limits<int>::signaling_NaN();
  int hits_n_t = std::numeric_limits<int>::signaling_NaN();
  int hits_s_t = std::numeric_limits<int>::signaling_NaN();

  std::vector<float> time_sum_n;
  std::vector<float> time_sum_s;

  float sum_n = std::numeric_limits<float>::signaling_NaN();
  float sum_s = std::numeric_limits<float>::signaling_NaN();
  float sum_n2 = std::numeric_limits<float>::signaling_NaN();
  float sum_s2 = std::numeric_limits<float>::signaling_NaN();

};

#endif
