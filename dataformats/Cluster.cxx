////////////////////////////////////////////////////////////////
//                                                            //
//  MpdEmcHitClusterKI                                          //
//  Cluster production for EMC                                //
//  Author List : D.Peresunko., RRCKI, 2019                   //
//                                                            //
////////////////////////////////////////////////////////////////

#include "MpdEmcClusterKI.h"
#include <iostream>
#include <numeric>
#include "FairLogger.h"
#include "MpdEmcClusterizerKI.h"
#include "MpdEmcDigitKI.h"
#include "MpdEmcGeoUtils.h"
#include "MpdEmcSimParams.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TVector3.h"

using namespace std;
using namespace TMath;

// Function used for sorting primaries
bool compE(const std::pair<int, float>& a, const std::pair<int, float>& b)
{
  return a.second > b.second;
}

// -----   Default constructor   -------------------------------------------

MpdEmcClusterKI::MpdEmcClusterKI()
  : fNDigits(0), fDitisId(nullptr), fDigitsE(nullptr), fNPrimaries(0), fPrimId(nullptr), fPrimE(nullptr), fE(0), fEcore(0), fEcore1p(0), fEcore2p(0), fTime(0), fX(0), fY(0), fZ(0), fdPhi(9999), fdZ(9999), fTrackId(-1), fDisp(0), fChi2(0), fLambda1(0), fLambda2(0), fNExLM(0)
{
}

// -----   Standard constructor   ------------------------------------------
MpdEmcClusterKI::MpdEmcClusterKI(const MpdEmcDigitKI* digit)
  : fNDigits(0), fDitisId(nullptr), fDigitsE(nullptr), fNPrimaries(0), fPrimId(nullptr), fPrimE(nullptr), fE(0), fEcore(0), fEcore1p(0), fEcore2p(0), fTime(0), fX(0), fY(0), fZ(0), fdPhi(9999), fdZ(9999), fTrackId(-1), fDisp(0), fChi2(0), fLambda1(0), fLambda2(0), fNExLM(0)
{
  AddDigit(digit);
}
// -----   Destructor   ----------------------------------------------------

MpdEmcClusterKI::~MpdEmcClusterKI()
{
  if (fDitisId) {
    delete[] fDitisId;
  }
  if (fDigitsE) {
    delete[] fDigitsE;
  }
  if (fPrimId) {
    delete[] fPrimId;
  }
  if (fPrimE) {
    delete[] fPrimE;
  }
}

void MpdEmcClusterKI::GetMomentum(TLorentzVector& p, const TVector3* vertex) const
{
  // Returm momentum of photon assuming it came from the provided vertex
  if (fX == 0. && fY == 0.) { // position not defined
    p.SetPxPyPzE(0., 0., 0., 0.);
    return;
  }
  double dx = fX, dy = fY, dz = fZ;
  if (vertex) { // calculate direction from vertex
    dx -= vertex->X();
    dy -= vertex->Y();
    dz -= vertex->Z();
  }

  Double_t r = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
  if (r > 0) {
    p.SetPxPyPzE(fE * dx / r, fE * dy / r, fE * dz / r, fE);
  } else {
    p.SetPxPyPzE(0., 0., 0., 0.);
  }
}
// -----  Print  -----------------------------------------------------------

void MpdEmcClusterKI::Print(const Option_t* opt) const
{
  cout << "MpdEmcClusterKI: " << endl;
  cout << "\tDeposited energy: " << fE << "\tMean time: " << fTime << "\tDigits: " << fNDigits
       << "   Rho cluster: " << GetRho() << "   Phi cluster: " << GetPhi() << "   Z cluster: " << fZ
       << "   dPhi cluster: " << fdPhi << "   dZ cluster: " << fdZ << "   index cluster: " << fTrackId << endl;

  cout << "\tNumber of tracks in module: " << fNPrimaries << endl;
  for (int i = 0; i < fNPrimaries; i++)
    cout << " " << fPrimId[i] << ", ";
  cout << endl;
}
//____________________________________________________________________________
void MpdEmcClusterKI::AddDigit(const MpdEmcDigitKI* digit, double energy)
{
  // Add another digit to cluster:
  // increase full energy, re-calculate time if necessary, add MC info, add digit to list of digits

  if (energy == 0) { // use full energy of digit, otherwise method is called from Unfolding.
    energy = digit->GetE();
  }
  fE += energy;
  // check if this is contribution with highest energy
  std::vector<pair<int, float>>::iterator elist = fDigitIDEnergy.begin();
  bool found = false;
  while (elist != fDigitIDEnergy.end()) {
    if ((*elist).second > energy) {
      found = true;
      break;
    }
    elist++;
  }
  if (!found) {
    fTime = digit->GetTime();
  }

  fDigitIDEnergy.emplace_back(digit->GetDetId(), energy);

  // now add new primaries
  int nPrim = digit->GetNPrimaries();
  for (int iprim = 0; iprim < nPrim; iprim++) {
    Int_t primId;
    Float_t primEdep;
    digit->GetPrimary(iprim, primId, primEdep);
    // Scale primary energy accosging to proportion
    if (digit->GetE() > 0) {
      primEdep *= energy / digit->GetE();
      std::map<Int_t, Float_t>::iterator it =
        fMCTracks.lower_bound(primId);                          // pointer, where new pair should appear, not necessarilly existing
      if ((it == fMCTracks.end()) || ((*it).first != primId)) { // do not exist yet
        fMCTracks.insert(it, {primId, primEdep});
      } else { // Add energy deposition
        (*it).second += primEdep;
      }
    }
  }
}
//____________________________________________________________________________
void MpdEmcClusterKI::EvalAll()
{
  // Calculate
  // Cluster coordinates
  // the shower dispersion, moments and Chi2
  // Fills final arrays with accociated digits and energy deposition

  // Remove cells below threshold and possibly single ceparated cells
  Purify();

  if (fE <= 0) { // Energy can not be zero in good cluster, do nothing
    return;
  }

  // Eval position and dispersion
  double wtot = 0.;
  double x = 0.;
  double y = 0.;
  double z = 0.;
  double phi = 0.;

  double dphiphi = 0.;
  double dphiz = 0.;
  double dzz = 0.;

  fDisp = 0.;
  fLambda1 = 0;
  fLambda2 = 0;

  MpdEmcGeoUtils* geom = MpdEmcGeoUtils::GetInstance();
  MpdEmcSimParams* simParams = MpdEmcSimParams::GetInstance();
  double logWeight = simParams->LogWeight();

  // 1) Find covariance matrix elements:
  //    || dxx dxz ||
  //    || dxz dzz ||

  std::vector<std::pair<int, float>>::iterator digIterator;

  digIterator = fDigitIDEnergy.begin();
  while (digIterator != fDigitIDEnergy.end()) {
    int detID = (*digIterator).first;
    float energy = (*digIterator).second;
    digIterator++;

    double xi, yi, zi, w, phii;
    geom->DetIdToGlobalPosition(detID, xi, yi, zi);
    if (energy > 0) {
      w = TMath::Max(0., logWeight + TMath::Log(energy / fE));
      x += w * xi;
      y += w * yi;
      z += w * zi;
      phii = ATan2(yi, xi);
      // in case of cluster around 2pi: avoid averaging of delta and 2pi+delta
      if (wtot > 0) { // not first digit
        double phiCurrent = phi / wtot;
        if (phii > phiCurrent + TMath::Pi()) {
          phii -= TMath::TwoPi();
        }
        if (phii < phiCurrent - TMath::Pi()) {
          phii += TMath::TwoPi();
        }
      }
      phi += w * phii;
      dphiphi += w * phii * phii;
      dphiz += w * phii * zi;
      dzz += w * zi * zi;
      wtot += w;
    }
  }
  if (wtot > 0) {
    x /= wtot;
    y /= wtot;
    z /= wtot;
    double r = TMath::Sqrt(x * x + y * y);
    phi = phi * r / wtot;
    dphiphi = dphiphi * r * r / wtot;
    dphiz = dphiz * r / wtot;
    dphiphi -= phi * phi;
    dphiz -= phi * z;
    dzz /= wtot;
    dzz -= z * z;

    fLambda2 = 0.5 * (dphiphi + dzz) + TMath::Sqrt(0.25 * (dphiphi - dzz) * (dphiphi - dzz) + dphiz * dphiz);
    fLambda1 = 0.5 * (dphiphi + dzz) - TMath::Sqrt(0.25 * (dphiphi - dzz) * (dphiphi - dzz) + dphiz * dphiz);
    if (fLambda2 > 0)
      fLambda2 = sqrt(fLambda2);
    if (fLambda1 > 0)
      fLambda1 = sqrt(fLambda1);
  } else {
    fLambda2 = fLambda1 = 0.;
  }
  fX = x;
  fY = y;
  fZ = z;

  // calculates agreement of the cluster shape with parameterization in MpdEmcClusterizer::ShowerShape
  // and core energies
  fChi2 = 0;
  int ndf = 0;
  fEcore = 0.;
  fEcore1p = 0.;
  fEcore2p = 0.;
  digIterator = fDigitIDEnergy.begin();
  while (digIterator != fDigitIDEnergy.end()) {
    int detID = (*digIterator).first;
    float energy = (*digIterator).second;
    digIterator++;

    double xi, yi, zi;
    geom->DetIdToGlobalPosition(detID, xi, yi, zi);
    double dz = zi - fZ;
    double dphi = TMath::Sqrt((xi - fX) * (xi - fX) + (yi - fY) * (yi - fY));
    double ss = MpdEmcClusterizerKI::ShowerShape(dphi, dz);

    if (sqrt(dz * dz + dphi * dphi) < simParams->EcoreR()) {
      fEcore += energy;
    }

    if (ss > simParams->EcoreCut1()) {
      fEcore1p += energy;
    }
    if (ss > simParams->EcoreCut2()) {
      fEcore2p += energy;
    }

    if (ss > simParams->Chi2radiusCut()) {
      double sigma = 0.008 * ss * fE * (1 - ss); // TODO: What is the meaning of 0.008???
      if (sigma != 0) {
        fChi2 += pow(ss * fE - energy, 2) / sigma;
        ndf++;
      }
    }
  }
  if (ndf > 0) {
    fChi2 /= ndf;
  } else {
    fChi2 = -1;
  }

  // Correct full energy
  fE = simParams->ENonLinCorrection(0) + simParams->ENonLinCorrection(1) * fE;

  // Prepare for writing
  FillArrays();
}

//____________________________________________________________________________
int MpdEmcClusterKI::GetNumberOfLocalMax(int* maxAt, float* maxAtEnergy) const
{
  // Calculates the number of local maxima in the cluster using LocalMaxCut as the minimum
  // energy difference between maximum and surrounding digits

  int n = GetMultiplicity();
  MpdEmcGeoUtils* geom = MpdEmcGeoUtils::GetInstance();
  MpdEmcSimParams* simParams = MpdEmcSimParams::GetInstance();
  float locMaxCut = simParams->LocalMaximumCut();

  bool* isLocalMax = new bool[n];
  for (int i = 0; i < n; i++)
  // isLocalMax[i] = true;
  {
    isLocalMax[i] = false;
    float en1 = fDigitIDEnergy[i].second;
    if (en1 > simParams->ClusteringThreshold())
      isLocalMax[i] = true;
  }

  for (int i = 0; i < n; i++) {
    int detId1 = fDigitIDEnergy[i].first;
    float en1 = fDigitIDEnergy[i].second;

    for (int j = i + 1; j < n; j++) {
      int detId2 = fDigitIDEnergy[j].first;
      float en2 = fDigitIDEnergy[j].second;

      if (geom->AreNeighboursVertex(detId1, detId2) == 1) {
        if (en1 > en2) {
          isLocalMax[j] = false;
          // but may be digit too is not local max ?
          if (en2 > en1 - locMaxCut) {
            isLocalMax[i] = false;
          }
        } else {
          isLocalMax[i] = false;
          // but may be digitN is not local max too?
          if (en1 > en2 - locMaxCut) {
            isLocalMax[j] = false;
          }
        }
      } // if Areneighbours
    } // digit j
  } // digit i

  int iDigitN = 0;
  for (int i = 0; i < n; i++) {
    if (isLocalMax[i]) {
      maxAt[iDigitN] = i;
      maxAtEnergy[iDigitN] = fDigitIDEnergy[i].second;
      iDigitN++;
      if (iDigitN >= simParams->NLMMax()) { // Note that size of output arrays is limited:
        LOG(error) << "Too many local maxima, cluster multiplicity " << n;
        return 0;
      }
    }
  }
  delete[] isLocalMax;
  return iDigitN;
}
//____________________________________________________________________________
void MpdEmcClusterKI::Purify()
{
  // Removes digits below threshold
  // If after purifying isolated digits remain, remove them too

  MpdEmcGeoUtils* geom = MpdEmcGeoUtils::GetInstance();
  MpdEmcSimParams* simParams = MpdEmcSimParams::GetInstance();
  auto digIterator = fDigitIDEnergy.begin();
  while (digIterator != fDigitIDEnergy.end()) {
    if (digIterator->second < simParams->DigitMinEnergy()) { // too soft, remove
      digIterator = fDigitIDEnergy.erase(digIterator);
    } else {
      digIterator++;
    }
  }

  int mult = fDigitIDEnergy.size();
  if (mult == 0) {
    return;
  }

  // Remove non-connected cells
  int index[mult];
  bool used[mult] = {false};
  int inClu = 0;
  double eMax = 0.;
  // find maximum
  for (int i = 0; i < mult; i++) {
    if (eMax < fDigitIDEnergy[i].second) {
      eMax = fDigitIDEnergy[i].second;
      index[0] = i;
      inClu = 1;
    }
  }
  used[index[0]] = true; // mark as used
  for (Int_t i = 0; i < inClu; i++) {
    int index1 = fDigitIDEnergy[i].first;
    for (int j = 0; j < mult; j++) {
      if (used[j]) // already used
        continue;
      int index2 = fDigitIDEnergy[j].first;
      if (geom->AreNeighboursVertex(index1, index2) == 1) {
        index[inClu] = j;
        inClu++;
        used[j] = true;
      }
    }
  }

  // cleanup not-connected cells
  for (int i = fDigitIDEnergy.size() - 1; i >= 0; --i) {
    if (!used[i]) {
      fDigitIDEnergy.erase(fDigitIDEnergy.begin() + i);
    }
  }
}
void MpdEmcClusterKI::FillArrays()
{
  // As root is not able to handle std::map and std::vector, Fill arrays to store to disk

  fNDigits = fDigitIDEnergy.size();
  fDitisId = new Int_t[fNDigits];
  fDigitsE = new Float_t[fNDigits];
  std::vector<std::pair<int, float>>::iterator digIterator = fDigitIDEnergy.begin();
  int i = 0;
  while (digIterator != fDigitIDEnergy.end()) {
    fDitisId[i] = (*digIterator).first;
    fDigitsE[i] = (*digIterator).second;
    digIterator++;
    i++;
  }

  // Sort map accordig to deposited energy
  std::vector<std::pair<int, float>> vec(fMCTracks.begin(), fMCTracks.end());
  std::sort(vec.begin(), vec.end(), compE);

  fNPrimaries = fMCTracks.size();
  fNPrimaries = TMath::Min(fNPrimaries, MpdEmcSimParams::GetInstance()->NPrimMax());
  fPrimId = new Int_t[fNPrimaries];
  fPrimE = new Float_t[fNPrimaries];
  i = 0;
  std::vector<std::pair<int, float>>::iterator primIterator = vec.begin();
  while (primIterator != vec.end() && i < fNPrimaries) {
    fPrimId[i] = (*primIterator).first;
    fPrimE[i] = (*primIterator).second;
    i++;
    primIterator++;
  }
}
void MpdEmcClusterKI::CorrectVertex(double zVtx)
{
  // correct cluster position for Z ccordinate of vertex
  // function is sA*sin(z/sW)+a+b*z, where sA,sW,a,b may depend on Zvtx and E

  MpdEmcSimParams* simParams = MpdEmcSimParams::GetInstance();
  double logE = TMath::Log(fE); // NB: here we assume that energy already corrected for 1/3

  double sA = simParams->ZcorrSinA(0) + simParams->ZcorrSinA(1) * logE; //-5.59366e-01 -2.67599e-02*logE
  double sW = simParams->ZcorrSinW(0) + simParams->ZcorrSinW(1) * logE;
  double a = (simParams->ZcorrA(0) + simParams->ZcorrA(1) * logE) * zVtx;
  double b = simParams->ZcorrB(0) + simParams->ZcorrB(1) * logE;
  fZ += (sA * TMath::Sin(fZ / sW) + a + b * fZ);
}
ClassImp(MpdEmcClusterKI);
