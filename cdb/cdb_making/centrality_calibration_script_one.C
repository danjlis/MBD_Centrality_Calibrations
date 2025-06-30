#include <fstream>
#include "MakeCentralityScaleCDBTree.C"
#include "MakeCentralityVertexScaleCDBTree.C"
#include "MakeCentralityCDBTree.C"

void centrality_calibration_script_one(const int x)
{

  std::cout <<x<<std::endl;

  CDB_MakeCentralityScaleCDBTree(x);
  CDB_MakeCentralityVertexScaleCDBTree(x);
  
  //CDB_CentralityScaleCDBInsert(x);
  //CDB_CentralityScaleCDBInsert(x);


  CDB_MakeCentralityCDBTree(x);

  //  CDB_CentralityCDBInsert(23696);
}
