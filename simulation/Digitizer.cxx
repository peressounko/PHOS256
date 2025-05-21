//--------------------------------------------------------------------
//
// Description:
//      Digitizer - takes Hits and makes Digits
//
//
// Author List:
//                D.Peresunko, KI, 2025
//--------------------------------------------------------------------

#include <iostream>
#include "Digitizer.h"

#include <TRandom.h>

using namespace std;
using namespace TMath;

//__________________________________________________________________________
void Digitizer::ProcessEvent()
{
  // Add all Hits with energy deposition in same tower
  // Simulate electronic noise
  // Add non-linearity, digitization of energy
  // Remove digits in bad map and below threshold

  if (!fHits) {
    return;
  }
  // Reset output Array
  if (fDigits) {
    fDigits->Clear();
  } else {
    // Reset output Array
    std::cout << "ERROR: fDigits not set!!!" << std::endl;
  }
  int nDigits = 0;

  // Class with list of parameters
  if (!fSimParams) {
    fSimParams = SimParams::Instance();
  }
  if (!fGeom) {
    fGeom = Geometry::Instance();
  }

  // Find the first cell with signal
  int nHits = fHits->size();
  int nextSigId = 999999; // in case no ginal in
  auto h = fHits->begin();
  if (h != fHits->end()) {
    nextSigId = h->GetCellID();
  }

  int nTotCells = fGeom->GetNCristalsInModule();
  // go through all cells, either add noisy digit, or signal+noise
  for (Int_t cellId = 0; cellId < nTotCells; cellId++) {
    // If signal exist in this cell, add noise to it, otherwise just create noise digit
    if (cellId == nextSigId) {
      // Create new Digit
      Digit* digit = new ((*fDigits)[nDigits++]) Digit(*h);
      ++h;
      while (h != fHits->end()) {
        if (digit->CanAdd(*h)) { // Point from the same Tower
          digit->AddHit(*h);
          ++h;
        } else { // Points are sorted according to cellID. If no more points left, finish
          nextSigId = h->GetCellID();
          break;
        }
      }
      // Add Electroinc noise, apply non-linearity, digitize, de-calibrate, time resolution
      double energy = digit->GetE();
      // Emulate Poissonian light collection
      if (fSimParams->fSmearLightCollection) {
        energy = SimulateLightCollection(energy);
      }

      if (fSimParams->fSimulateNoise) {
        // Simulate electronic noise
        energy += SimulateNoiseEnergy();
        if (energy < 0) {
          energy = 0.;
        }
      }

      if (fSimParams->fApplyNonLinearity) {
        energy = NonLinearity(energy);
      }
      if (fSimParams->fApplyDigitization) {
        energy = DigitizeEnergy(energy);
      }

      // Convert energy from GeV to ADC units - as in real data
      energy /= fSimParams->fADCwidth; // V
      digit->SetE(energy);

      if (fSimParams->fApplyTimeResolution) {
        digit->SetTime(TimeResolution(digit->GetTime(), energy));
      }
    } else { // No signal in this cell,
      // Simulate noise
      if (fSimParams->fSimulateNoise) {
        double energy = SimulateNoiseEnergy();
        double time = SimulateNoiseTime();
        if (energy > fSimParams->fZSthreshold) {
          energy /= fSimParams->fADCwidth;
          new ((*fDigits)[nDigits++])
            Digit(cellId, energy, time, -1); // current sellId, energy, random time, no primary
        }
      }
    }
  }
}

//_______________________________________________________________________
double Digitizer::SimulateNoiseEnergy()
{
  // Simulation of noise of electronics
  // Here we assume, that noise is independent on signal
  // and just Gaus with fixed width
  return gRandom->Gaus(0., fSimParams->fAPDNoise);
}
//_______________________________________________________________________
double Digitizer::NonLinearity(const double e)
{

  return e * fSimParams->fCellNonLineaityC * (1. + fSimParams->fCellNonLineaityA * exp(-e / fSimParams->fCellNonLineaityB));
}
//_______________________________________________________________________
double Digitizer::DigitizeEnergy(const double e)
{
  // distretize energy if necessary
  double w = fSimParams->fADCwidth;
  return w * ceil(e / w);
}
//_______________________________________________________________________
double Digitizer::TimeResolution(const double time, const double e)
{
  // apply time resolution
  if (e <= 0)
    return 0.;
  double timeResolution = fSimParams->fTimeResolutionA + fSimParams->fTimeResolutionB / e;
  return gRandom->Gaus(time, timeResolution);
}
//_______________________________________________________________________
double Digitizer::SimulateNoiseTime()
{
  // Evaluate (random) time of noise digits as uniform in some ranges
  return gRandom->Uniform(fSimParams->fMinNoiseTime, fSimParams->fMaxNoiseTime);
}
//_______________________________________________________________________
double Digitizer::SimulateLightCollection(const double lostenergy)
{
  // Emulate production of scintillator light and its collection
  double coef = fSimParams->fLightYieldPerGeV;
  double lightYield = gRandom->Poisson(coef * lostenergy);
  return lightYield / coef;
}
