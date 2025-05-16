#include "g4libs.C"
#include "Simulation.h"
#include "TGeant4.h"
void Run(std::string configMacro = "Config.C")
{
  g4libs();

  gSystem->Load("libG4clhep");

  // MC application
  Simulation* appl = new Simulation();
  appl->SetRTheta(200., 20.);

  appl->InitMC(configMacro);

  TStopwatch timer;
  timer.Start();
  appl->RunMC(5);
  timer.Stop();
  timer.Print();

  delete appl;
}
