#include <qa_centrality/QA_centrality.h>

R__LOAD_LIBRARY(libQA_centrality.so);

int mainreferenceQA_OO(const int runnumber)
{
  gSystem->Load("libQA_centrality");
  QA_centrality *c = new QA_centrality(0);
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
  //  c->SetTrigEffMUK(.91, 3.84, 0.47);
  c->QA_ReferenceRun(runnumber);
  std::cout << __LINE__ <<std::endl;
  c->Print_QA_Info(true);
  return 0;
}
