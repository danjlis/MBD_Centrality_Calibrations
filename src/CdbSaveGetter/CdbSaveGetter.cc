#include "CdbSaveGetter.h"

#include <ffaobjects/CdbUrlSavev1.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllBase.h>

#include <phool/PHCompositeNode.h>
#include <phool/phool.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

CdbSaveGetter::CdbSaveGetter(const std::string &name)
  : SubsysReco(name)
{

}

CdbSaveGetter::~CdbSaveGetter()
{

}
int CdbSaveGetter::Init(PHCompositeNode * /*unused*/)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

int CdbSaveGetter::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }

  if (!GetNodes(topNode)){
    return Fun4AllReturnCodes::ABORTRUN;
  }
  std::cout << "**********************************************************************"<<std::endl;
  std::cout << "*    The Holy CDB Url Save getter bestows upon you the location of   *" << std::endl;
  std::cout << "* the magical cantations containing the power knowledge itself. Be   *" << std::endl;
  std::cout << "* careful of those that hide their mathematical formulations, and    *" << std::endl;
  std::cout << "* heed to noone the belief that their numbers are as good as yours.  *" << std::endl;
  std::cout << "*                                                                    *" << std::endl;
  std::cout << "*                                                                    *" << std::endl;
  std::cout << "*                                                                    *" << std::endl;
  std::cout << "* The file in question:                                              *" << std::endl;
  std::cout << "*                                                                    *" << std::endl;
  _cdb_save->identify();
  std::cout << "*                                                                    *" << std::endl;
  std::cout << "**********************************************************************"<<std::endl;
  std::cout << "*                                                                    *" << std::endl;
  std::cout << "*  Thank you  :)                                                     *" << std::endl;
  std::cout << "*                                                                    *" << std::endl;
  std::cout << "**********************************************************************"<<std::endl;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int CdbSaveGetter::process_event(PHCompositeNode * /*unused*/)
{
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int CdbSaveGetter::GetNodes(PHCompositeNode *topNode)
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << " :: " << __LINE__ << std::endl;
  }
  
  _cdb_save = findNode::getClass<CdbUrlSavev1>(topNode, "CdbUrl");
  
  if (!_cdb_save)
    {
      std::cout << "no cdb save node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CdbSaveGetter::End(PHCompositeNode * /* topNode*/)
{

  return 0;
}
