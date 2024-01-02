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
class InttRawHit;
class InttRawHitContainer;
class TH1;
class TH2;
class InttMbd : public SubsysReco
{
 public:
  //! constructor
  explicit InttMbd(const std::string &name = "InttMbd", const std::string &tree_name = "mbd_intt_aligned.root");

  //! destructor
  virtual ~InttMbd();

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

  TFile *outfile = nullptr;
  TTree *ttree = nullptr;

  MbdOut *_mbd_out = nullptr;  
  InttRawHit *_intt_hit = nullptr;
  InttRawHitContainer *_intt_hit_container = nullptr;
  
  std::string _tree_filename;

  unsigned short m_mbd_fem_clock = std::numeric_limits<unsigned short>::signaling_NaN();
  unsigned short m_mbd_xmit_clock = std::numeric_limits<unsigned short>::signaling_NaN();
  int m_mbd_event = std::numeric_limits<int>::signaling_NaN();
  unsigned long int  m_intt_bco = std::numeric_limits<unsigned long int>::signaling_NaN();

};

#endif
