#ifndef CLUSTER_H
#define CLUSTER_H 1

// Cluster is a set of cells with common side or corner
// Some parameters like full energy, position, time, shape are calculated

#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "Digit.h"

class Cluster : public TObject
{
 public:
  Cluster() = default;

  /** Constructor with hit parameters **/
  Cluster(const Digit& digit);

  /** Destructor **/
  ~Cluster() = default;

  // Returm momentum of photon assuming it came from the provided vertex
  void GetMomentum(TLorentzVector& p, const TVector3& vertex) const;

  void AddDigit(const Digit& digit, Double_t edep = 0);
  void EvalAll(); // Evaluate cluster parameters

  void Purify(); // Remove digits below threshold

  void CorrectVertex(const TVector3& vertex);

  int GetNumberOfLocalMax(int* maxAt, float* maxAtEnergy) const; // Finds local maxima

  int GetDigitCellId(int i) const { return fDigitIDEnergy.at(i).first; } // detectorId of i-th digit

  // Get tower ID and energy of i-th cotributing digit
  void GetDigitParams(int i, int& detId, float& eDigit) const
  {
    if (i >= 0 && i < fNDigits) {
      detId = fDitisId[i];
      eDigit = fDigitsE[i];
    } else {
      detId = -1;
      eDigit = 0;
    }
  }
  // Method used in unfolding and should not be used by analysers
  void GetTransientDigitParams(int i, int& detId, float& eDigit) const
  {
    std::pair<int, float> p = fDigitIDEnergy.at(i);
    detId = p.first;
    eDigit = p.second;
  }

  // Number of MC tracks
  int GetNumberOfTracks() const { return fNPrimaries; }
  void GetMCTrack(int i, int& trackId, float& trackEdep) const
  {
    if (i >= 0 && i < fNPrimaries) {
      trackId = fPrimId[i];
      trackEdep = fPrimE[i];
    } else {
      trackId = -1;
      trackEdep = 0;
    }
  }
  float GetE() const { return fE; };

  float GetEcore() const { return fEcore; };

  float GetChi2() const { return fChi2; };

  float GetTime() const { return fTime; };

  float GetX() const { return fX; };

  float GetY() const { return fY; };

  float GetZ() const { return fZ; };

  float GetRad() const { return TMath::Sqrt(fX * fX + fY * fY + fZ * fZ); };

  int GetMultiplicity() const
  {
    return TMath::Max((int)fDigitIDEnergy.size(),
                      fNDigits); // either from std::vector or in final cluster - from array
  };

  int GetNLM() const { return fNExLM; }

  void SetNLM(int n) { fNExLM = n; }

  void SetTime(float time) { fTime = time; };

  // Dispersion parameters
  void GetLambdas(float& l1, float& l2)
  {
    l1 = fLambda1;
    l2 = fLambda2;
  }

 protected:
  void FillArrays(); // Fill arrays to store to disk

 protected:
  // Parameters used for re-calibration if necessary
  std::vector<std::pair<int, float>> fDigitIDEnergy; //! transient list of contributed digits with energies
  std::map<int, float> fMCTracks;                    //! transient trackID and energy deposit

  int fNDigits;    // Digit multiplicity
  int* fDitisId;   //[fNDigits] cellId
  float* fDigitsE; //[fNDigits] deposited energy

  int fNPrimaries; // Number of primaries
  int* fPrimId;    //[fNPrimaries] cellId
  float* fPrimE;   //[fNPrimaries] deposited energy per cell

  float fE;       // cluster energy
  float fEcore;   // cluster energy core
  float fEcore1p; // cluster energy core
  float fEcore2p; // cluster energy core

  float fTime; // cluster time

  // cluster coordinates in global system
  float fX; // x-coordinate of cluster
  float fY; // y-coordinate of cluster
  float fZ; // z-coordinate of cluster

  // Dispersion and shower shape parameters
  float fDisp;    // Dispersion
  float fChi2;    // Chi2 of a fit with EM shape
  float fLambda1; // smaller Disp axis
  float fLambda2; // larger Disp axis

  int fNExLM; // Number of local maxima or NLM in parent cluster before unfolding

  ClassDefNV(Cluster, 1);
};

#endif
