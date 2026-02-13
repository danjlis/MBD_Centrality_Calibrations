#include <qa_centrality/QA_centrality.h>

R__LOAD_LIBRARY(libQA_centrality.so);

int mainreferenceQA(const int runnumber)
{
  gSystem->Load("libQA_centrality");
  QA_centrality *c = new QA_centrality(0);
  //c->setForceZDC(true);
  c->SetNEvents(1000000);
  c->SetCountBefore(false);
  c->SetNDivs(100);
  c->SetDivs(93);
  c->SetChargeCut(0.5);
  c->SetMBDHitCut(2);
  c->SetMBDTimeCut(20);
  c->SetTriggerBit(10);
  //c->setNTupleFile("glau_auau_ntuple.root");
  //c->setHistoFile("lemon_glauber_auau200_100k_histo.root");
  //  c->SetTrigEffMUK(.91, 3.84, 0.47);
  c->QA_ReferenceRun(runnumber);
  std::cout << __LINE__ <<std::endl;
  c->Print_QA_Info(true);
  return 0;
}
