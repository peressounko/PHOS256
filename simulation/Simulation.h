#ifndef SIMULATION_H
#define SIMULATION_H

/// \file Simulation.h
/// \brief Definition of the Simulation class
///
#include "TFile.h"
#include "TTree.h"
#include "TMCVerbose.h"
#include <TGeoUniformMagField.h>
#include <TVirtualMCApplication.h>

#include "TMCRootManager.h"

// This project
#include "Hall.h"
#include "Magnet.h"
#include "Phos.h"
#include "Stack.h"
#include "MagField.h"
#include "GenBox.h"
#include "Digitizer.h"
#include "Clusterizer.h"

#include "Hit.h"
#include "Digit.h"

class Simulation : public TVirtualMCApplication
{
 public:
  Simulation();
  Simulation(const Simulation* s);
  virtual ~Simulation() = default;
  static Simulation* Instance() { return fSimulation; }

  // methods
  void InitMC(std::string configName);
  void RunMC(int nofEvents);
  void FinishRun();

  void SetMagnet(Magnet* det) { fMagnet = det; }
  void SetPHOS(Phos* det) { fPHOS = det; }
  void SetGenerator(GenBox* gen) { fGenerator = gen; }

  void SetRTheta(double r = 100, double theta = 20)
  {
    fRad = r;
    fTheta = theta;
  }

  virtual TVirtualMCApplication* CloneForWorker() const;
  virtual void InitOnWorker();
  virtual void FinishRunOnWorker();

  virtual void ConstructGeometry();
  virtual void InitGeometry() { printf("InitGeometry\n"); }
  virtual void AddParticles() {}
  virtual void AddIons() {}
  virtual void GeneratePrimaries();
  virtual void BeginEvent();
  virtual void BeginPrimary() {}
  virtual void PreTrack() {}
  virtual void Stepping();
  virtual void PostTrack() {}
  virtual void FinishPrimary();
  virtual void FinishEvent();

 private:
  static Simulation* fSimulation; //! pointer to this instance

  mutable TMCRootManager* fRootManager; //!< Root manager
  bool fIsMaster = true;                ///< If is on master thread

  double fRad = 100.;
  double fTheta = 20.;

  Stack* fStack = nullptr;       ///< VMC stack
  GenBox* fGenerator = nullptr;  ///< Primary generator
  MagField* fMagField = nullptr; ///< Uniform magnetic field

  // Calorimeter description
  Hall* fHall = nullptr;
  Magnet * fMagnet = nullptr;
  Phos* fPHOS = nullptr;

  // processors
  Digitizer* fDigitizer = nullptr;
  Clusterizer* fClusterizer = nullptr;
  // Data collections
  std::vector<Hit> fHits;
  TClonesArray* fDigits = nullptr; //
  TObjArray* fClusters = nullptr;

  TTree* fTree = nullptr;
  TFile* fOutFile = nullptr;

  ClassDef(Simulation, 1) // Interface to MonteCarlo application
};

#endif //
