#include <TGeoManager.h>
#include <TGeoUniformMagField.h>
#include <TInterpreter.h>
#include <TPDGCode.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TVector3.h>
#include <TVirtualGeoTrack.h>
#include <TVirtualMC.h>
#include "TGeant4.h"

#include <TMCRootManager.h>

#include "Phos.h"
#include "Simulation.h"
#include "GenBox.h"
#include "Stack.h"

using namespace std;

Simulation* Simulation::fSimulation = nullptr;

Simulation::Simulation() : TVirtualMCApplication(),
                           fIsMaster(true),
                           fRad(100.),
                           fTheta(20.),
                           fStack(nullptr),
                           fGenerator(nullptr),
                           fMagField(nullptr),
                           fHall(nullptr),
                           fPHOS(nullptr),
                           fDigitizer(nullptr),
                           fClusterizer(nullptr),
                           fDigits(nullptr),
                           fClusters(nullptr),
                           fTree(nullptr),
                           fOutFile(nullptr)
{
  fSimulation = this;
}

Simulation::Simulation(const Simulation* s)
{
}
//_____________________________________________________________________________
void Simulation::InitMC(std::string configName)
{
  // Initialize containers

  if (configName.size()) {
    gROOT->LoadMacro(configName.data());
    gInterpreter->ProcessLine("Config()");
    if (!gMC) {
      Fatal(
        "InitMC", "Processing Config() has failed. (No MC is instantiated.)");
    }
  }

  // Create Root manager
  if (!gMC->IsMT()) {
    fRootManager = new TMCRootManager(GetName(), TMCRootManager::kWrite);
    fRootManager->SetDebug(true);
  }

  // Create a user stack
  fStack = new Stack(1000);
  if (fRootManager) {
    fRootManager->Register("stack", "Stack", &fStack);
  }

  // Set data to MC
  gMC->SetStack(fStack);

  // Constant magnetic field (in kiloGauss)
  // field value: 0.2*tesla (= 2.0 kiloGauss) in x
  fMagField = new MagField(2., 0.0, 0);

  gMC->SetMagField(fMagField);

  // Create a primary generator
  if (!fGenerator) {
    fGenerator = new GenBox(fStack);
  } else {
    fGenerator->SetStack(fStack);
  }

  gMC->Init();
  gMC->BuildPhysics();

  fOutFile = TFile::Open("PHOSReco.root", "recreate");
  fTree = new TTree("PHOS256", "Reconstruction tree");
  fTree->Branch("MCParticles", "TClonesArray", fStack->GetParticles(), 32000, 99);
}
//_____________________________________________________________________________
void Simulation::ConstructGeometry()
{
  // Construct GEANT geometry for all detectros

  if (!gGeoManager)
    new TGeoManager("PHOS256", "PHOS256 geometry");
  if (!fHall)
    fHall = new Hall();
  fHall->CreateGeometry();

  // construct GEANT geometry
  if (!fPHOS) {
    std::cout << "Creating Default PHOS!!!" << std::endl;
    fPHOS = new Phos(fRad, fTheta); // Make configurable
  }
  fPHOS->SetHitContainer(&fHits);
  fPHOS->CreateMaterials();
  fPHOS->CreateGeometry(); // creates the geometry for GEANT
}

//_____________________________________________________________________________
void Simulation::RunMC(Int_t nofEvents)
{
  /// Run MC.
  /// \param nofEvents Number of events to be processed
  gMC->ProcessRun(nofEvents);
  FinishRun();
}

//_____________________________________________________________________________
void Simulation::FinishRun()
{
  /// Finish MC run.

  if (fRootManager) {
    fRootManager->WriteAll();
    fRootManager->Close();
  }

  fOutFile->cd();
  fTree->Write();
  fOutFile->Close();

  gGeoManager->Export("geometry.root");
}

//_____________________________________________________________________________
TVirtualMCApplication* Simulation::CloneForWorker() const
{
  // cout << "Simulation::CloneForWorker " << this << endl;
  return new Simulation(*this);
}

//_____________________________________________________________________________
void Simulation::InitOnWorker()
{
  // cout << "Simulation::InitForWorker " << this << endl;

  // Create Root manager
  fRootManager = new TMCRootManager("PHSO256", TMCRootManager::kWrite);
  // fRootManager->SetDebug(true);

  // Set data to MC
  gMC->SetStack(fStack);
  gMC->SetMagField(fMagField);

  if (fRootManager) {
    fRootManager->Register("stack", "PHOS256Stack", &fStack);
  }
}

//_____________________________________________________________________________
void Simulation::FinishRunOnWorker()
{
  // cout << "Simulation::FinishWorkerRun: " << endl;
  if (fRootManager) {
    fRootManager->WriteAll();
    fRootManager->Close();
  }
}

//_____________________________________________________________________________
void Simulation::GeneratePrimaries()
{
  /// Fill the user stack (derived from TVirtualMCStack) with primary particles.

  fGenerator->Generate();
}

//_____________________________________________________________________________
void Simulation::BeginEvent()
{
  // gMC->BeginEvent();
  /// User actions at beginning of event
  // Clear TGeo tracks (if filled)
  if (TString(gMC->GetName()) == "TGeant3TGeo" &&
      gGeoManager->GetListOfTracks() && gGeoManager->GetTrack(0) &&
      ((TVirtualGeoTrack*)gGeoManager->GetTrack(0))->HasPoints()) {
    gGeoManager->ClearTracks();
  }
  fPHOS->Reset();
  fHits.clear();
  if (!fDigitizer) {
    fDigitizer = new Digitizer();
    fDigitizer->SetHits(&fHits);
    fDigits = new TClonesArray("Digit", 10);
    fDigitizer->SetDigits(fDigits);
    fTree->Branch("Digits", &fDigits, 32000, 99);
  }
  if (!fClusterizer) {
    fClusterizer = new Clusterizer();
    fClusterizer->Init();
    fClusters = new TObjArray();
    fClusterizer->SetDigits(fDigits);
    fClusterizer->SetClusters(fClusters);
    fTree->Branch("Clusters", "TObjArray", &fClusters, 32000, 0);
  }
}

//_____________________________________________________________________________
void Simulation::Stepping()
{
  /// User actions at each step

  fPHOS->ProcessHits();
}

//_____________________________________________________________________________
void Simulation::FinishPrimary()
{
  /// User actions after finishing each primary track
  fPHOS->FinishPrimary();
}

//_____________________________________________________________________________
void Simulation::FinishEvent()
{
  /// User actions after finishing an event

  // Print info about primary particle
  if (fStack->GetNtrack()) {
    TParticle* primary = fStack->GetParticle(0);
    cout << endl
         << ">>> Event " << gMC->CurrentEvent()
         << " >>> Simulation truth : " << primary->GetPDG()->GetName() << " ("
         << primary->Px() * 1e03 << ", " << primary->Py() * 1e03 << ", "
         << primary->Pz() * 1e03 << ") MeV" << endl;
  }
  fStack->Purge(); 


  // Call detectors
  fPHOS->FinishEvent();
  if (fDigitizer) {
    fDigitizer->ProcessEvent();
  }
  if (fClusterizer) {
    fClusterizer->ProcessEvent();
  }

  // fRootManager->Fill();
  fTree->Fill();

  fStack->Reset();
}
