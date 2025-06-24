#include <qa_centrality/QA_centrality.h>

R__LOAD_LIBRARY(libQA_centrality.so);
int mainQA(const int runnumber, const int reference_run)
{
  gSystem->Load("libQA_centrality");

  QA_centrality *c = new QA_centrality(0);
  c->SetReferenceRun(reference_run);
  c->SetNEvents(1000000);

  c->SetCountBefore(false);
  c->SetNDivs(100);
  c->SetDivs(93);
  c->SetChargeCut(0.5);

  c->Start_QA_Centrality(runnumber);
  c->Print_QA_Info(true);
  return 0;
}
