#include <qa_centrality/QA_centrality.h>

R__LOAD_LIBRARY(libQA_centrality.so);
int mainQA_2(const int runnumber)
{
  gSystem->Load("libQA_centrality");

  QA_centrality *c = new QA_centrality(0);
  c->SetReferenceRun(23696);
  c->QA_MakeCentralityCalibrations(runnumber,1);
  c->SetTrigEffMUK(0.91, 4.05, 0.65);
  c->QA_MakeCentralityCalibrations(runnumber,1, true);
  c->Print_QA_Info(true);
  return 0;
}
