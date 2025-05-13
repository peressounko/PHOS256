/// \file RunConfiguration.cxx
/// \brief Implementation of the RunConfiguration class
///
/// Geant4 Example adapted to Virtual Monte Carlo \n

#include "RunConfiguration.h"

//_____________________________________________________________________________
RunConfiguration::RunConfiguration(
  const TString& physicsList, const TString& specialProcess)
  : TG4RunConfiguration("geomGeant4", physicsList, specialProcess),
    fUseLocalMagField(false)
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

// //_____________________________________________________________________________
// G4VUserDetectorConstruction* RunConfiguration::CreateDetectorConstruction()
// {
//   /// The Geant4 VMC detector construction is overridden with the detector
//   /// construction class from the Geant4 novice example N03 library.

//   return new Phos();
// }
