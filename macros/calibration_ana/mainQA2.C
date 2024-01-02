#include "QA_centrality.h"
#include <stdio.h>

int mainQA2(void)
{
  gSystem->Load("libQA_centrality.so");

  std::cout << " this is a shared library test..."<<std::endl;;
  QA_centrality *c = new QA_centrality(0,0);
  c->QA_MakeCentralityCalibrations(23746, 1);

  return 0;
}
