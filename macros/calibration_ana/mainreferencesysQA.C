#include <qa_centrality/QA_centrality.h>

R__LOAD_LIBRARY(libQA_centrality.so);

int mainreferencesysQA(const int runnumber)
{
  gSystem->Load("libQA_centrality");
  QA_centrality *c = new QA_centrality(1);
  c->SetCountBefore(false);
  c->SetNDivs(100);
  c->SetDivs(92);
  c->setSysName("nominal");
  c->SetChargeCut(0.5);
  c->QA_ReferenceRun(runnumber);
  std::cout << " ************ Nominal ************* " << std::endl;
  c->Print_QA_Info(true);
  c->SetDivs(89);
  c->SetChargeCut(0.5);
  c->setSysName("90");
  c->QA_ReferenceRun(runnumber);
  std::cout << " ************ 91 divs ************* " << std::endl;
  c->Print_QA_Info(true);
  c->SetDivs(95);
  c->setSysName("95");
  c->SetChargeCut(0.5);
  c->QA_ReferenceRun(runnumber);
  std::cout << " ************ 95 divs ************* " << std::endl;
  c->Print_QA_Info(true);


  c->SetDivs(92);

  c->setNTupleFile("glau_auau_ntuple_39mb.root");
  c->setHistoFile("lime_auau_hists_39mb.root");
  c->setSysName("39mb");
  c->QA_ReferenceRun(runnumber);
  std::cout << " ************ 39mb ************* " << std::endl;
  c->Print_QA_Info(true);


  c->setNTupleFile("glau_auau_ntuple_45mb.root");
  c->setHistoFile("lime_auau_hists_45mb.root");
  c->setSysName("45mb");
  c->QA_ReferenceRun(runnumber);
  std::cout << " ************ 45mb ************* " << std::endl;
  c->Print_QA_Info(true);


  c->setNTupleFile("glau_auau_ntuple_A1.root");
  c->setNTupleName("nt_AuA1_AuA1");
  c->setHistoFile("lime_auau_hists_A1.root");
  c->setSysName("A1");
  c->QA_ReferenceRun(runnumber);
  std::cout << " ************ A1 ************* " << std::endl;
  c->Print_QA_Info(true);


  c->setNTupleFile("glau_auau_ntuple_A2.root");
  c->setNTupleName("nt_AuA2_AuA2");
  c->setHistoFile("lime_auau_hists_A2.root");
  c->setSysName("A2");
  c->QA_ReferenceRun(runnumber);
  std::cout << " ************ A2 ************* " << std::endl;
  c->Print_QA_Info(true);

  return 0;
}
