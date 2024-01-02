#ifndef CDBSAVEGETTER_H
#define CDBSAVEGETTER_H

#include <fun4all/SubsysReco.h>

#include <limits>

// Forward declarations
class CdbUrlSavev1;
class PHCompositeNode;
class CdbSaveGetter : public SubsysReco
{
 public:
  //! constructor
  explicit CdbSaveGetter(const std::string &name = "CdbSaveGetter");

  //! destructor
  virtual ~CdbSaveGetter();

  //! full initialization
  int Init(PHCompositeNode *);
  int InitRun(PHCompositeNode *);

  int GetNodes(PHCompositeNode *);
  //! event processing method
  int process_event(PHCompositeNode *);
  
  //! end of run method
  int End(PHCompositeNode *);

 protected:
  CdbUrlSavev1 *_cdb_save = nullptr;
};

#endif
