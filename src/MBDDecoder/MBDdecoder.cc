#include "MBDdecoder.h"
#include <caloreco/CaloTowerDefs.h>
#include <caloreco/CaloWaveformProcessing.h>

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoContainerv3.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <TSystem.h>

#include <climits>
#include <iostream>  // for operator<<, endl, basic...
#include <memory>    // for allocator_traits<>::val...
#include <vector>    // for vector

//____________________________________________________________________________..
MBDdecoder::MBDdecoder(const std::string &name)
: SubsysReco(name)
{
  WaveformProcessing = new CaloWaveformProcessing();
}

//____________________________________________________________________________..
MBDdecoder::~MBDdecoder()
{
  delete WaveformProcessing;
}

//____________________________________________________________________________..
int MBDdecoder::InitRun(PHCompositeNode *topNode)
{
  WaveformProcessing->set_processing_type(_processingtype);
  WaveformProcessing->set_softwarezerosuppression(_bdosoftwarezerosuppression, _nsoftwarezerosuppression);
  m_detector = "MBD";
  m_packet_low = 1001;
  m_packet_high = 1003;
  m_nchannels = 128;
  WaveformProcessing->set_template_name("CEMC_TEMPLATE");
  WaveformProcessing->set_processing_type(CaloWaveformProcessing::FAST);

  WaveformProcessing->initialize_processing();
  CreateNodeTree(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MBDdecoder::process_event(PHCompositeNode *topNode)
{
  if (!m_isdata)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  std::vector<std::vector<float>> waveforms;
  std::vector<short> clocks;
  // if we are going from prdf
  Event *_event = findNode::getClass<Event>(topNode, "PRDF");
  if (_event == nullptr)
    {
      std::cout << PHWHERE << " Event not found" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  if (_event->getEvtType() != DATAEVENT)  /// special events where we do not read out the calorimeters
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  for (int pid = m_packet_low; pid <= m_packet_high; pid++)
    {
      Packet *packet = _event->getPacket(pid);
      if (packet)
	{
	  int nchannels = packet->iValue(0, "CHANNELS");
	  short clkcounter = (0xffff & packet->iValue(0, "FEMCLOCK"));
	  clocks.push_back(clkcounter);
	  if (nchannels != m_nchannels)  // packet is corrupted and reports too many channels
	    {
	      return Fun4AllReturnCodes::ABORTEVENT;
	    }
	  // int sector = 0;

	  for (int channel = 0; channel < nchannels; channel++)
	    {
	      // mask empty channels

	      std::vector<float> waveform;
	      waveform.reserve(m_nsamples);
	      for (int samp = 0; samp < m_nsamples; samp++)
		{
		  waveform.push_back(packet->iValue(samp, channel));
		}
	      waveforms.push_back(waveform);

	      waveform.clear();
	    }

	  delete packet;
	}
    }
  
  // waveform vector is filled here, now fill our output. methods from the base class make sure
  // we only fill what the chosen container version supports
  std::vector<std::vector<float>> processed_waveforms = WaveformProcessing->process_waveform(waveforms);
  int n_channels = processed_waveforms.size();
  for (int i = 0; i < n_channels; i++)
    {
      TowerInfo *towerinfo = m_CaloInfoContainer->get_tower_at_channel(i);
      towerinfo->set_time(clocks.at(0));
      towerinfo->set_energy(processed_waveforms.at(i).at(0));
      //towerinfo->set_time_float(processed_waveforms.at(i).at(1));
      //towerinfo->set_pedestal(processed_waveforms.at(i).at(2));
      //towerinfo->set_chi2(processed_waveforms.at(i).at(3));

    }
  clocks.clear();
  waveforms.clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

void MBDdecoder::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator topNodeItr(topNode);
  // DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(topNodeItr.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cout << "PHComposite node created: DST" << std::endl;
      dstNode = new PHCompositeNode("DST");
      topNode->addNode(dstNode);
    }
  // towers
  PHNodeIterator nodeItr(dstNode);
  PHCompositeNode *DetNode;
  // enum CaloTowerDefs::DetectorSystem and TowerInfoContainer::DETECTOR are different!!!!
  TowerInfoContainer::DETECTOR DetectorEnum = TowerInfoContainer::DETECTOR::DETECTOR_INVALID;
  std::string DetectorNodeName;

  DetectorEnum = TowerInfoContainer::DETECTOR::MBD;
  DetectorNodeName = "MBD";

  DetNode = dynamic_cast<PHCompositeNode *>(nodeItr.findFirst("PHCompositeNode", DetectorNodeName));
  if (!DetNode)
    {
      DetNode = new PHCompositeNode(DetectorNodeName);
      dstNode->addNode(DetNode);
    }
  
  m_CaloInfoContainer = new TowerInfoContainerv1(DetectorEnum);

  TowerNodeName = m_outputNodePrefix + m_detector;
  PHIODataNode<PHObject> *newTowerNode = new PHIODataNode<PHObject>(m_CaloInfoContainer, TowerNodeName, "PHObject");
  DetNode->addNode(newTowerNode);
}
