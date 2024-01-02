#include <qa_centrality/QA_centrality.h>

R__LOAD_LIBRARY(libQA_centrality.so);

int mainreferenceQA(const int runnumber)
{
  gSystem->Load("libQA_centrality");

  QA_centrality *c = new QA_centrality(0);
  c->QA_ReferenceRun(runnumber);
  c->Print_QA_Info(true);
  return 0;
}
