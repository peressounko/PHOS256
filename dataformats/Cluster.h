#ifndef CLUSTER_H
#define CLUSTER_H 1

#include <map>
#include "TMath.h"
#include "TObject.h"
class Digit;
class TLorentzVector;
class TVector3;

class Cluster : public TObject
{
 public:
  Cluster() = default;

  Cluster(const Digit* digit);

  virtual ~Cluster() = default;

  // Returm momentum of photon assuming it came from the provided vertex
  void GetMomentum(TLorentzVector& p) const;

  void AddDigit(const Digit* digit, Double_t edep = 0);
  void EvalAll(); // Evaluate cluster parameters

  void Purify(); // Remove digits below threshold

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
  int GetNumberOfLabels() const { return fNPrimaries; }
  void GetLabel(int i, int& label, float& labelEdep) const
  {
    if (i >= 0 && i < fNPrimaries) {
      label = fPrimId[i];
      labelEdep = fPrimE[i];
    } else {
      label = -1;
      labelEdep = 0;
    }
  }
  void UpdateLabel(int i, int newLabel) { fPrimId[i] = newLabel; }

  float GetE() const { return fE; };

  float GetEcore() const { return fEcore; };

  float GetChi2() const { return fChi2; };

  float GetTime() const { return fTime; };

  float GetLocX() const { return fLocX; }
  float GetLocZ() const { return fLocZ; }

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

  static double ShowerShape(double dx, double dz); // Parameterization of EM shower

 protected:
  void FillArrays(); // Fill arrays to store to disk

 protected:
  static constexpr int kNLMMax = 5; // Maximal number of local maxima
  // Parameters used for re-calibration if necessary
  std::vector<std::pair<int, float>> fDigitIDEnergy; //! transient list of contributed digits with energies
  std::map<int, float> fMCTracks;                    //! transient trackID and energy deposit

  int fNDigits = 0;          // Digit multiplicity
  int* fDitisId = nullptr;   //[fNDigits] cellId
  float* fDigitsE = nullptr; //[fNDigits] deposited energy

  int fNPrimaries = 0;     // Number of primaries
  int* fPrimId = nullptr;  //[fNPrimaries] cellId
  float* fPrimE = nullptr; //[fNPrimaries] deposited energy per cell

  float fE = 0;     // cluster energy
  float fEcore = 0; // cluster energy core

  float fTime = 0; // cluster time

  // cluster coordinatex in local system
  float fLocX;
  float fLocZ;

  // cluster coordinates in global system
  float fX = 0; // x-coordinate of cluster
  float fY = 0; // y-coordinate of cluster
  float fZ = 0; // z-coordinate of cluster

  // Dispersion and shower shape parameters
  float fDisp = 0;    // Dispersion
  float fChi2 = 0;    // Chi2 of a fit with EM shape
  float fLambda1 = 0; // smaller Disp axis
  float fLambda2 = 0; // larger Disp axis

  int fNExLM = 0; // Number of local maxima or NLM in parent cluster before unfolding

  ClassDefNV(Cluster, 1);
};

#endif
