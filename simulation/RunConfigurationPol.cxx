/// \file RunConfigurationPol.cxx
/// \brief Implementation of the RunConfigurationPol class
///
/// Geant4 Example adapted to Virtual Monte Carlo \n

#include "RunConfigurationPol.h"
#include "G4EmLivermorePolarizedPhysics.hh"
#include "TG4ComposedPhysicsList.h"
#include "TG4Globals.h"
#include "TG4EmPhysicsList.h"
#include "TG4HadronPhysicsList.h"
#include "TG4SpecialPhysicsList.h"
#include "TG4ExtraPhysicsList.h"
#include "G4EmLivermorePolarizedPhysics.hh"
#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics.hh"
#include "G4DecayPhysics.hh"
#include <FTFP_BERT.hh>
//_____________________________________________________________________________
RunConfigurationPol::RunConfigurationPol(
  const TString& physicsList, const TString& specialProcess)
  : TG4RunConfiguration("geomRoot", physicsList, specialProcess)
{
  /// Standard constructor
  /// \param physicsList     Selection of physics
  /// \param specialProcess  Selection of the special processes
  ///
  /// The option for geometry selection has to be set here to "geomGeant4",
  /// as geometry will be defined directly via Geant4.
  /// \see More on the available option in class TG4RunConfiguration:
  /// http://ivana.home.cern.ch/ivana/g4vmc_html/classTG4RunConfiguration.html
}

G4VUserPhysicsList* RunConfigurationPol::CreatePhysicsList()
{

  TG4ComposedPhysicsList* builder = new TG4ComposedPhysicsList();

  G4VModularPhysicsList* physicsList = new FTFP_BERT();
  // Standard electromagnetic physics
  // physicsList->RegisterPhysics(new G4EmStandardPhysics(1)); //already in FTFP_BERT
  physicsList->ReplacePhysics(new G4EmLivermorePolarizedPhysics());

  G4cout << "Added FTFP_BERT+G4EmLivermorePolarizedPhysics " << G4endl;

  // decay physics
  // physicsList->RegisterPhysics(new G4DecayPhysics()); //already in FTFP_BERT

  builder->AddPhysicsList(physicsList);

  // if (hadronSelection != "") {
  //   G4cout << "Adding HadronPhysicsList " << hadronSelection << G4endl;
  //   builder->AddPhysicsList(new TG4HadronPhysicsList(hadronSelection));
  // }

  // if (extraSelection != "") {
  //   G4cout << "Adding ExtraPhysicsList " << extraSelection << G4endl;
  //   builder->AddPhysicsList(
  //     new TG4ExtraPhysicsList(extraSelection, fParameters));
  // }

  // // add option here
  // G4cout << "Adding SpecialPhysicsList " << fSpecialProcessSelection.Data()
  //        << G4endl;
  // builder->AddPhysicsList(
  //   new TG4SpecialPhysicsList(fSpecialProcessSelection.Data()));

  G4cout << "Created default PhysicsList " << G4endl;

  return builder;

  // new G4EmLivermorePolarizedPhysics();
}

// //_____________________________________________________________________________
// G4VUserDetectorConstruction* RunConfiguration::CreateDetectorConstruction()
// {
//   /// The Geant4 VMC detector construction is overridden with the detector
//   /// construction class from the Geant4 novice example N03 library.

//   return new Phos();
// }
