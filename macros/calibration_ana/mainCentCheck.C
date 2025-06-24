#include <qa_centrality/QA_centrality.h>

R__LOAD_LIBRARY(libQA_centrality.so);

int mainCentCheck(const int runnumber, const int ref = 0)
{
  gSystem->Load("libQA_centrality");
  std::cout << __LINE__ <<std::endl;
  QA_centrality *c = new QA_centrality(0);
  c->SetNDivs(100);
  int referencerun = ref;
  if (ref == 0)
    referencerun = runnumber;

  //  c->SetTrigEffMUK(.91, 3.84, 0.47);
  c->QA_RunCentralityCheck(runnumber, referencerun);
  std::cout << __LINE__ <<std::endl;
  c->Print_QA_Info(true);
  return 0;
}
