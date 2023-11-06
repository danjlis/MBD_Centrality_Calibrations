#ifndef MACRO_RUNDUMP_C
#define MACRO_RUNDUMP_C


#include <nodedump/Dumper.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstInputManager.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libfun4all.so)
// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libphnodedump.so)

void BigDumper(const std::string &infile = "/sphenix/lustre01/sphnxpro/commissioning/DSTv3/DST_CALOR-00021199-0000.root", const int evts=5)
{
  gSystem->Load("libg4dst.so");
  Fun4AllServer* se = Fun4AllServer::instance();

  Dumper *dmp = new Dumper();
  gSystem->Exec("mkdir dump");
  dmp->SetOutDir("./dump");

  se->registerSubsystem(dmp);

  Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTin");
  se->registerInputManager(in);
  se->fileopen("DSTin",infile);
  se->run(evts);
  se->End();
  delete se;
}

#endif
