#include <qa_centrality/QA_centrality.h>

R__LOAD_LIBRARY(libQA_centrality.so);

int mainQAMC(const std::string runnumber)
{
  gSystem->Load("libQA_centrality");

  QA_centrality *c = new QA_centrality(0);
  c->SetNDivs(100);
  //  c->SetTrigEffMUK(.91, 3.84, 0.47);
  c->QA_MC(runnumber);
  c->Print_QA_Info(true);
  return 0;
}
