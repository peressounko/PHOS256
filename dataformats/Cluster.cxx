////////////////////////////////////////////////////////////////
//                                                            //
//  Cluster                                                   //
//  Cluster production for EMC                                //
//  Author List : D.Peresunko., RRCKI, 2025                   //
//                                                            //
////////////////////////////////////////////////////////////////

#include "Geometry.h"
#include "SimParams.h"
#include "Cluster.h"
#include <iostream>
#include <numeric>
#include "Digit.h"
#include "TLorentzVector.h"
#include "TVector3.h"

// Function used for sorting primaries
bool compE(const std::pair<int, float>& a, const std::pair<int, float>& b)
{
  return a.second > b.second;
}

// -----   Standard constructor   ------------------------------------------
Cluster::Cluster(const Digit* digit)
  : fNDigits(0), fDitisId(nullptr), fDigitsE(nullptr), fNPrimaries(0), fPrimId(nullptr), fPrimE(nullptr), fE(0), fEcore(0), fTime(0), fLocX(0), fLocZ(0), fX(0), fY(0), fZ(0), fDisp(0), fChi2(0), fLambda1(0), fLambda2(0), fNExLM(0)
{
  AddDigit(digit);
}

void Cluster::GetMomentum(TLorentzVector& p) const
{
  // Returm momentum of photon assuming it came from the provided vertex
  if (fX == 0. && fY == 0.) { // position not defined
    p.SetPxPyPzE(0., 0., 0., 0.);
    return;
  }
  double dx = fX, dy = fY, dz = fZ;

  Double_t r = std::sqrt(dx * dx + dy * dy + dz * dz);
  if (r > 0) {
    p.SetPxPyPzE(fE * dx / r, fE * dy / r, fE * dz / r, fE);
  } else {
    p.SetPxPyPzE(0., 0., 0., 0.);
  }
}
//____________________________________________________________________________
void Cluster::AddDigit(const Digit* digit, double energy)
{
  // Add another digit to cluster:
  // increase full energy, re-calculate time if necessary, add MC info, add digit to list of digits

  if (energy == 0) { // use full energy of digit, otherwise method is called from Unfolding.
    energy = digit->GetE();
  }
  fE += energy;
  // check if this is contribution with highest energy
  std::vector<std::pair<int, float>>::iterator elist = fDigitIDEnergy.begin();
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

  fDigitIDEnergy.emplace_back(digit->GetCellId(), energy);

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
void Cluster::EvalAll()
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
  double z = 0.;

  double dxx = 0.;
  double dxz = 0.;
  double dzz = 0.;

  fDisp = 0.;
  fLambda1 = 0;
  fLambda2 = 0;

  Geometry* geom = Geometry::Instance();
  const SimParams* simParams = SimParams::Instance();
  double logWeight = simParams->fLogWeight;

  // 1) Find covariance matrix elements:
  //    || dxx dxz ||
  //    || dxz dzz ||

  std::vector<std::pair<int, float>>::iterator digIterator;

  digIterator = fDigitIDEnergy.begin();
  while (digIterator != fDigitIDEnergy.end()) {
    int detID = (*digIterator).first;
    float energy = (*digIterator).second;
    digIterator++;

    double xi, zi, w, phii;
    geom->DetIdToLocalPosition(detID, xi, zi);
    if (energy > 0) {
      w = TMath::Max(0., logWeight + TMath::Log(energy / fE));
      x += w * xi;
      z += w * zi;
      dxx += w * xi * xi;
      dxz += w * xi * zi;
      dzz += w * zi * zi;
      wtot += w;
    }
  }
  if (wtot > 0) {
    x /= wtot;
    z /= wtot;
    dxx = dxx / wtot;
    dxz = dxz / wtot;
    dzz /= wtot;
    dxx -= x * x;
    dxz -= x * z;
    dzz -= z * z;

    fLambda2 = 0.5 * (dxx + dzz) + std::sqrt(0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz);
    fLambda1 = 0.5 * (dxx + dzz) - std::sqrt(0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz);
    if (fLambda2 > 0)
      fLambda2 = std::sqrt(fLambda2);
    if (fLambda1 > 0)
      fLambda1 = std::sqrt(fLambda1);
  } else {
    fLambda2 = fLambda1 = 0.;
  }
  fLocX = x;
  fLocZ = z;

  // calculates agreement of the cluster shape with parameterization in MpdEmcClusterizer::ShowerShape
  // and core energies
  fChi2 = 0;
  int ndf = 0;
  fEcore = 0.;
  digIterator = fDigitIDEnergy.begin();
  while (digIterator != fDigitIDEnergy.end()) {
    int detID = (*digIterator).first;
    float energy = (*digIterator).second;
    digIterator++;

    double xi, zi;
    geom->DetIdToLocalPosition(detID, xi, zi);
    double dz = zi - fLocZ;
    double dx = xi - fLocX;
    double ss = ShowerShape(dx, dz);

    if (sqrt(dz * dz + dx * dx) < simParams->fCoreR) {
      fEcore += energy;
    }

    if (ss > simParams->fChi2radiusCut) {
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
  fE = simParams->fCluNonLineaityA * fE + simParams->fCluNonLineaityB * exp(-fE / simParams->fCluNonLineaityC);

  // Fill Global coordinates
  //  convert local position in module to global position in World
  TVector3 globalPos;
  geom->Local2Global(fLocX, fLocZ, globalPos);

  fX = globalPos.X();
  fY = globalPos.Y();
  fZ = globalPos.Z();

  // Prepare for writing
  FillArrays();
}

//____________________________________________________________________________
int Cluster::GetNumberOfLocalMax(int* maxAt, float* maxAtEnergy) const
{
  // Calculates the number of local maxima in the cluster using LocalMaxCut as the minimum
  // energy difference between maximum and surrounding digits

  int n = GetMultiplicity();
  const Geometry* geom = Geometry::Instance();
  const SimParams* simParams = SimParams::Instance();
  float locMaxCut = simParams->fLocalMaximumCut;

  bool* isLocalMax = new bool[n];
  for (int i = 0; i < n; i++)
  // isLocalMax[i] = true;
  {
    isLocalMax[i] = false;
    float en1 = fDigitIDEnergy[i].second;
    if (en1 > simParams->fClusteringThreshold)
      isLocalMax[i] = true;
  }

  for (int i = 0; i < n; i++) {
    int detId1 = fDigitIDEnergy[i].first;
    float en1 = fDigitIDEnergy[i].second;

    for (int j = i + 1; j < n; j++) {
      int detId2 = fDigitIDEnergy[j].first;
      float en2 = fDigitIDEnergy[j].second;

      if (geom->AreNeighbours(detId1, detId2) == 1) {
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
      if (iDigitN >= kNLMMax) { // Note that size of output arrays is limited:
        std::cout << "ERROR: Too many local maxima, cluster multiplicity " << n;
        return 0;
      }
    }
  }
  delete[] isLocalMax;
  return iDigitN;
}
//____________________________________________________________________________
void Cluster::Purify()
{
  // Removes digits below threshold
  // If after purifying isolated digits remain, remove them too

  const Geometry* geom = Geometry::Instance();
  const SimParams* simParams = SimParams::Instance();
  auto digIterator = fDigitIDEnergy.begin();
  while (digIterator != fDigitIDEnergy.end()) {
    if (digIterator->second < simParams->fDigitMinEnergy) { // too soft, remove
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
      if (geom->AreNeighbours(index1, index2) == 1) {
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
void Cluster::FillArrays()
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
  fNPrimaries = std::min(fNPrimaries, SimParams::Instance()->fNPrimMax);
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
//__________________________________________________________________________
double Cluster::ShowerShape(double dx, double dz)
{
  // Shower shape of EM cluster: proportion of energy in cell as a distancs dx,dz (phi,z directions)
  // from the center of gravity of cluster.
  // Due to projective geometry may also depend on Z coordinate. TODO: explore Z-dependence

  // Parameterization from V.Riabov. TODO: verify with beam-test

  double frac = 0;
  double x = std::sqrt(dx * dx + dz * dz);

  if (x < 0.25)
    return 0.73;

  if (x < 4.5) {
    frac = 7.55666e+000 * exp(-2.35773e+000 + 1.23342e-001 * x - 2.53958e-001 * x * x + 1.41214e-002 * x * x * x);
  }
  if (x >= 4.5 && x < 8.33) {
    frac = 4.26832e-002 * exp(1.47109e+000 - 2.06712e-001 * x - 6.60220e-002 * x * x + 4.88207e-003 * x * x * x);
  }
  if (x >= 8.33 && x < 2.65 * 4) {
    frac = 4.46227e-002 * exp(1.56817e+000 - 7.41051e-001 * x / 4 - 4.80563e-001 * x / 4 * x / 4);
  }
  if (x >= 2.65 * 4) {
    frac = 1.24899e-002 * exp(3.60075e-001 - 8.15748e-001 * x / 4 - 3.74305e-002 * x / 4 * x / 4 * x / 4);
  }

  if (frac < 1e-24)
    frac = 1e-24;

  return frac;
}
