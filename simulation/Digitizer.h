//--------------------------------------------------------------------
//
// Description:
//       Digitizer - takes Hits and makes Digits
//
//
//--------------------------------------------------------------------

#ifndef DIGITIZER_H
#define DIGITIZER_H 1

#include <vector>
#include <TClonesArray.h>

#include "Digit.h"
#include "Geometry.h"
#include "SimParams.h"

class Digitizer
{
 public:
  Digitizer() = default;
  ~Digitizer() = default;

  void SetHits(std::vector<Hit>* hits) { fHits = hits; }
  void SetDigits(TClonesArray* digits) { fDigits = digits; }

  void ProcessEvent();

 private:
  double SimulateNoiseEnergy();                             // Simulation of noise of electronics
  double NonLinearity(const double e);                      // simulate non-lineraity
  double DigitizeEnergy(const double e);                    // Account final width of ADC
  double TimeResolution(const double time, const double e); // Apply final time resolution
  double SimulateNoiseTime();                               // calculate time in noise digit
  double SimulateLightCollection(const double lostenergy);  // Simulate Poissonian light production and collection

 private:
  std::vector<Hit>* fHits = nullptr; //! input hits

  TClonesArray* fDigits = nullptr; //! output digits

  SimParams* fSimParams; //! Class with all simulation parameters

  const Geometry* fGeom; //! Geometry parameters and methods
  ClassDefNV(Digitizer, 1);
};

#endif
