#include <qa_centrality/QA_centrality.h>

R__LOAD_LIBRARY(libQA_centrality.so);
int mainQA(const int runnumber)
{
  gSystem->Load("libQA_centrality");

  QA_centrality *c = new QA_centrality(1);
  c->SetReferenceRun(23696);
  c->Start_QA_Centrality(runnumber);
  c->Print_QA_Info(true);
  return 0;
}
