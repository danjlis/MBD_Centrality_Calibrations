#ifndef MBDANA_H
#define MBDANA_H

#include <fun4all/SubsysReco.h>
#include <ffarawobjects/Gl1Packet.h>
#include <limits>

// Forward declarations
class PHCompositeNode;
class TriggerAnalyzer;
class EventHeader;
class TProfile;
class TFile;
class TNtuple;
class TTree;
class Zdcinfo;
class Gl1Packet;
class MbdOut;
class MbdPmtHit;
class MbdPmtContainer;
class MinimumBiasInfo;
class EventplaneinfoMap;
class Eventplaneinfo;
class CentralityInfo;
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
  void setIsSim(bool useSim){ m_issim = useSim;}  
  //! end of run method
  int End(PHCompositeNode *);
  double GetPsi(double Qx, double Qy, unsigned int order);

  void ResetVars();

 protected:

  bool m_issim{false};
  enum Bbc
  {
    South = 0,
    North = 1
  };

  bool useZDC;
  TFile *outfile = nullptr;
  TTree *ttree = nullptr;
  MbdPmtHit *_tmp_pmt = nullptr;  
  MbdPmtContainer *_pmts_mbd = nullptr;
  Zdcinfo *_zdcinfo = nullptr;
  Gl1Packet *_gl1_packet = nullptr;
  MbdOut *_mbd_out = nullptr;
  TriggerAnalyzer *triggeranalyzer{nullptr};
  MinimumBiasInfo *_minimumbiasinfo{nullptr};
  CentralityInfo *centralityinfo{nullptr};
  EventplaneinfoMap *eventplaneinfomap{nullptr};
  
  EventHeader *m_eventheader{nullptr};
  std::string _tree_filename;

  int _tubes_hit_s = std::numeric_limits<int>::signaling_NaN();;
  int _tubes_hit_n = std::numeric_limits<int>::signaling_NaN();;

  float _energy = std::numeric_limits<float>::signaling_NaN();
  float m_bimp{999};
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


  float m_psi1{0};
  float m_psi1s{0};
  float m_psi1n{0};
  
  float m_psi2{0};
  float m_psi2s{0};
  float m_psi2n{0};
  
  float m_psi3{0};
  float m_psi3s{0};
  float m_psi3n{0};

  float m_centrality{0};
  uint64_t m_bunch_number{0};
  uint64_t m_gl1_live{0};
  uint64_t m_gl1_scaled{0};
  int m_minbias{0};

};

#endif
