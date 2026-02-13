#include <qa_centrality/QA_centrality.h>

R__LOAD_LIBRARY(libQA_centrality.so);
int mainQA_OO(const int runnumber, const int reference_run)
{
  gSystem->Load("libQA_centrality");

  QA_centrality *c = new QA_centrality(0);
  c->SetReferenceRun(reference_run);

  c->setForceZDC(true);
  c->SetNEvents(1000000);
  c->SetCountBefore(true);
  c->SetTriggerBit(12);
  c->SetOO(true);
  c->SetMBDHitCut(1);
  c->SetNDivs(100);
  c->SetDivs(88);
  c->SetChargeCut(0.4);
  c->setNTupleFile("glau_oo_ntuple.root");
  c->setNTupleName("nt_O_O");
  c->setHistoFile("lemon_oo_hists.root");

  c->Start_QA_Centrality(runnumber);
  c->Print_QA_Info(true);
  return 0;
}
