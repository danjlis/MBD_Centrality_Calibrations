// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef MBDDECODER_H
#define MBDDECODER_H

#include <caloreco/CaloWaveformProcessing.h>
#include <caloreco/CaloTowerDefs.h>

#include <fun4all/SubsysReco.h>

#include <limits>
#include <string>

class CaloWaveformProcessing;
class PHCompositeNode;
class TowerInfoContainer;
class TowerInfoContainerv3;


class MBDdecoder : public SubsysReco
{
 public:
  explicit MBDdecoder(const std::string &name = "MBDdecoder");
  ~MBDdecoder() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  void CreateNodeTree(PHCompositeNode *topNode);

  void set_nsamples(int _nsamples)
  {
    m_nsamples = _nsamples;
    return;
  }
  void set_dataflag(bool flag)
  {
    m_isdata = flag;
    return;
  }

  void set_processing_type(CaloWaveformProcessing::process processingtype)
  {
    _processingtype = processingtype;
  }


  void set_outputNodePrefix(const std::string &name)
  {
    m_outputNodePrefix = name;
    return;
  }

 private:
  CaloWaveformProcessing *WaveformProcessing {nullptr};
  TowerInfoContainer *m_CaloInfoContainer {nullptr};  //! Calo info
  bool m_isdata {true};
  bool _bdosoftwarezerosuppression {false};
  int m_packet_low {std::numeric_limits<int>::min()};
  int m_packet_high {std::numeric_limits<int>::min()};
  int m_nsamples {16};
  int m_nchannels {192};
  int m_nzerosuppsamples {2};
  int _nsoftwarezerosuppression {40};
  CaloWaveformProcessing::process _processingtype {CaloWaveformProcessing::NONE};

  std::string m_detector {"CEMC"};
  std::string m_outputNodePrefix {"TOWERS_"};
  std::string TowerNodeName;
};

#endif  // MBDDECODER_H
