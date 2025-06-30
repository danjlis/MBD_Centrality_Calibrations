#include <fstream>
#include "MakeCentralityScaleCDBTree.C"
#include "MakeCentralityVertexScaleCDBTree.C"
#include "MakeCentralityCDBTree.C"

void centrality_calibration_script(const string input="grl.list")
{

  std::ifstream runnumbers(input);

  int x, y;
  char *env_p = new char[200];
  sprintf(env_p,"%s",std::getenv("MBD_CENTRALITY_CALIB_PATH"));

  if(!env_p)
    {
      std::cout << "no env MBD_CENTRALITY_CALIB_PATH set."<<endl;
      return;
    }

  std::ofstream calib_file("calibration_file.list");
  calib_file << "runnumber,divs,scales,vertexscales"<<std::endl;
  while (runnumbers >> x >> y)
    {
      std::cout <<x<<"," << y << std::endl;

      TString cdb_name1 = Form("%s/cdb/calibrations24/divs/cdb_centrality_%d.root", env_p, x);
      TString cdb_name2 = Form("%s/cdb/calibrations24/scales/cdb_centrality_scale_%d.root", env_p, x);
      TString cdb_name3 = Form("%s/cdb/calibrations24/vertexscales/cdb_centrality_vertex_scale_%d.root", env_p, x);

      if (CDB_MakeCentralityScaleCDBTree(x)) continue;
      if (CDB_MakeCentralityVertexScaleCDBTree(x)) continue;;
      if (CDB_MakeCentralityCDBTree(x, y)) continue;

      calib_file << x << "," << cdb_name1 << "," << cdb_name2 << "," << cdb_name3  <<std::endl;
    }
}
