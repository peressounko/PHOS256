//-----------------------------------------------------------
//
// Author List:
//      D.Peresunko, KI, 2025
//
//-----------------------------------------------------------

// This Class' Header ------------------
#include "Clusterizer.h"
#include "Cluster.h"
#include "Digit.h"
#include "Geometry.h"
#include "SimParams.h"
// #include "CalibParams.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TMath.h"
#include <iostream>

//__________________________________________________________________________
void Clusterizer::Init()
{

  // // Temporary solution: get calibration from file TODO!!!!
  // if (!fCalibData) {

  //   TString CalFile = "$VMCWORKDIR/input/MpdEmcCalib.root";

  //   TFile* f = new TFile(CalFile);

  //   fCalibData = (MpdEmcCalibParams*)f->Get("CalibrationDef");
  //   cout << "ECAL: Read out EMC calib data" << endl;
  // }

  // Class with list of parameters
  if (!fSimParams) {
    fSimParams = SimParams::Instance();
  }
  if (!fGeom) {
    fGeom = Geometry::Instance();
  }
}

//__________________________________________________________________________
void Clusterizer::ProcessEvent()
{
  fClustersArray->Delete(); // No possibility to Clear for TObjArray
  fNumberOfClusters = 0;

  // Prepare Digits: Calibration, bad map, cleanup
  PrepareDigits();

  // Collect list of clusters
  MakeClusters();

  // Split clusters with several local maxima if necessary
  MakeUnfoldings();

  // Evaluate cluster position, dispersion etc.
  EvalClusters();
}
//__________________________________________________________________________
void Clusterizer::PrepareDigits()
{
  // Apply (de)calibration
  // TODO: Implement Energy and time calibration for data and de-calibrations for MC
  // Apply BadMap
  // TODO: Implement bad map

  // Remove digits below clustering threshold
  int n = fDigitsArray->GetEntriesFast();
  double calib = fSimParams->fADCwidth;
  for (int i = 0; i < n; i++) {
    Digit* digit = static_cast<Digit*>(fDigitsArray->UncheckedAt(i));
    if (!digit) { // already removed e.g. by bad map selection
      continue;
    }
    // Calibrate energy
    // digit->SetE(digit->GetE() * fCalibData->GetGain(digit->GetCellId()));
    digit->SetE(digit->GetE() * calib);

    if (digit->GetE() < fSimParams->fDigitMinEnergy) {
      fDigitsArray->RemoveAt(i);
    }
  }
  fDigitsArray->Compress();
  // Alredy sorted by construction, but ROOT requires this to allow search
  fDigitsArray->Sort();
}
//__________________________________________________________________________
void Clusterizer::MakeClusters()
{
  // Combine digits into cluster according to definition of neighbours in
  // MpdEmcGeoUtils::AreNighbours()
  // Note, that only digits in same sector contribute to the cluster

  // Mark all digits as unused yet
  int nDigits = fDigitsArray->GetEntriesFast();
  if (nDigits < 1) { // nothing to do
    return;
  }
  bool* digitsUsed = new bool[nDigits]{false};

  int iFirst = 0; // first index of digit which potentially can be a part of cluster
                  // e.g. first digit in this sector

  for (int i = 0; i < nDigits; i++) {
    if (digitsUsed[i])
      continue;

    const Digit* digitSeed = static_cast<Digit*>(fDigitsArray->UncheckedAt(i));

    // is this digit so energetic that start cluster?
    Cluster* clu = nullptr;
    int iDigitInCluster = 0;
    if (digitSeed->GetE() > fSimParams->fClusteringThreshold) {
      // start a new EMC RecPoint
      clu = new Cluster(digitSeed);
      fClustersArray->AddAtAndExpand(clu, fNumberOfClusters);
      fNumberOfClusters++;
      digitsUsed[i] = true;
      iDigitInCluster = 1;
    } else {
      continue;
    }

    // Now scan remaining digits in list to find neigbours of our seed
    int index = 0;
    while (index < iDigitInCluster) { // scan over digits already in cluster
      int digitInCluTowerId = clu->GetDigitCellId(index);
      index++;
      for (int j = iFirst; j < nDigits;
           j++) { // upper limit not really matters, AreNeighbour stops this loop for too far digits
        if (digitsUsed[j]) {
          continue; // look through remaining digits
        }
        const Digit* digitN = static_cast<Digit*>(fDigitsArray->UncheckedAt(j));

        // call (digit,digitN) in THAT oder !!!!!
        int ineb = fGeom->AreNeighbours(digitInCluTowerId, digitN->GetCellId());
        switch (ineb) {
          case -1: // too early (e.g. previous sector), do not look before j at subsequent passes
            iFirst = j + 1;
            continue;
          case 0: // not a neighbour
            continue;
          case 1: // are neighbours
            clu->AddDigit(digitN);
            iDigitInCluster++;
            digitsUsed[j] = true;
            continue;
          case 2: // too far from each other, stop loop
          default:
            goto nextDigit;

        } // switch
      }
    nextDigit:;
    } // loop over cluster
  } // energy theshold
  delete[] digitsUsed;
}
//__________________________________________________________________________
void Clusterizer::MakeUnfoldings()
{
  // Split cluster if several local maxima are found

  if (!fSimParams->fUnfoldClusters) {
    return;
  }
  int* maxAt = new int[fSimParams->fNLMMax]; // NLMMax:Maximal number of local maxima
  float* maxAtEnergy = new float[fSimParams->fNLMMax];

  float localMaxCut = fSimParams->fLocalMaximumCut;
  int numberOfNotUnfolded = fNumberOfClusters;
  for (int index = 0; index < numberOfNotUnfolded; index++) {
    Cluster* clu = static_cast<Cluster*>(fClustersArray->At(index));

    int nMultipl = clu->GetMultiplicity();
    int nMax = clu->GetNumberOfLocalMax(maxAt, maxAtEnergy);
    if (nMax > 1) {
      UnfoldOneCluster(clu, nMax, maxAt, maxAtEnergy);

      fClustersArray->Remove(clu);
      fClustersArray->Compress();
      index--;
      fNumberOfClusters--;
      numberOfNotUnfolded--;
    } else {
      clu->SetNLM(1); // Only one local maximum
    }
  }
  delete[] maxAt;
  delete[] maxAtEnergy;
}
//____________________________________________________________________________
void Clusterizer::UnfoldOneCluster(Cluster* iniClu, int nMax, int* digitId, float* maxAtEnergy)
{
  // Performs the unfolding of a cluster with nMax overlapping showers
  // Parameters: iniClu cluster to be unfolded
  //             nMax number of local maxima found (this is the number of new clusters)
  //             digitId: index of digits, corresponding to local maxima
  //             maxAtEnergy: energies of digits, corresponding to local maxima

  // Take initial cluster and calculate local coordinates of digits
  // To avoid multiple re-calculation of same parameters
  int mult = iniClu->GetMultiplicity();
  std::vector<double> x(mult);
  std::vector<double> z(mult);
  std::vector<double> e(mult);
  std::vector<std::vector<double>> eInClusters(mult, std::vector<double>(nMax));

  for (int idig = 0; idig < mult; idig++) {
    int detID;
    float eDigit;
    iniClu->GetTransientDigitParams(idig, detID, eDigit);
    e[idig] = eDigit;
    double lx, lz;
    fGeom->DetIdToLocalPosition(detID, lx, lz);
    x[idig] = lx;
    z[idig] = lz;
  }

  // Coordinates of centers of clusters
  std::vector<double> xMax(nMax);
  std::vector<double> zMax(nMax);
  std::vector<double> eMax(nMax);

  for (int iclu = 0; iclu < nMax; iclu++) {
    xMax[iclu] = x[digitId[iclu]];
    zMax[iclu] = z[digitId[iclu]];
    eMax[iclu] = e[digitId[iclu]];
  }

  std::vector<double> prop(nMax); // proportion of clusters in the current digit

  // Try to decompose cluster to contributions
  int nIterations = 0;
  bool insuficientAccuracy = true;
  while (insuficientAccuracy && nIterations < fSimParams->fNMaxIterations) {
    // Loop over all digits of parent cluster and split their energies between daughter clusters
    // according to shower shape
    for (int idig = 0; idig < mult; idig++) {
      double eEstimated = 0;
      for (int iclu = 0; iclu < nMax; iclu++) {
        prop[iclu] = eMax[iclu] * Cluster::ShowerShape(x[idig] - xMax[iclu],
                                                       z[idig] - zMax[iclu]);
        eEstimated += prop[iclu];
      }
      if (eEstimated == 0.) { // numerical accuracy
        continue;
      }
      // Split energy of digit according to contributions
      for (int iclu = 0; iclu < nMax; iclu++) {
        eInClusters[idig][iclu] = e[idig] * prop[iclu] / eEstimated;
      }
    }
    // Recalculate parameters of clusters and check relative variation of energy and absolute of position
    insuficientAccuracy = false; // will be true if at least one parameter changed too much
    for (int iclu = 0; iclu < nMax; iclu++) {
      double oldX = xMax[iclu];
      double oldZ = zMax[iclu];
      double oldE = eMax[iclu];
      // new energy, need for weight
      eMax[iclu] = 0;
      for (int idig = 0; idig < mult; idig++) {
        eMax[iclu] += eInClusters[idig][iclu];
      }
      xMax[iclu] = 0;
      zMax[iclu] = 0.;
      double wtot = 0.;
      for (int idig = 0; idig < mult; idig++) {
        double w = std::max(std::log(eInClusters[idig][iclu] / eMax[iclu]) + fSimParams->fLogWeight, 0.);
        xMax[iclu] += x[idig] * w;
        zMax[iclu] += z[idig] * w;
        wtot += w;
      }
      if (wtot > 0.) {
        xMax[iclu] /= wtot;
        zMax[iclu] /= wtot;
      }
      // Compare variation of parameters
      insuficientAccuracy += (std::abs(eMax[iclu] - oldE) > fSimParams->fUnfogingEAccuracy);
      insuficientAccuracy += (std::abs(xMax[iclu] - oldX) > fSimParams->fUnfogingXZAccuracy);
      insuficientAccuracy += (std::abs(zMax[iclu] - oldZ) > fSimParams->fUnfogingXZAccuracy);
    }
    nIterations++;
  }

  // Iterations finished, add new clusters
  for (int iclu = 0; iclu < nMax; iclu++) {
    Cluster* clu = new Cluster();
    fClustersArray->AddAtAndExpand(clu, fNumberOfClusters + iclu);
    clu->SetNLM(nMax);
  }
  for (int idig = 0; idig < mult; idig++) {
    int detID;
    float eDigit;
    ;
    iniClu->GetTransientDigitParams(idig, detID, eDigit);
    Digit testdigit(detID, 0, 0, 0);                     // test digit
    int jdigit = fDigitsArray->BinarySearch(&testdigit); // Look for digit with same detID
    if (jdigit == -1) {
      std::cout << "Clusterizer::UnfoldOneCluster: Can not find Digit with detID=" << detID << std::endl;
      continue;
    }
    Digit* digit = static_cast<Digit*>(fDigitsArray->At(jdigit));
    for (int iclu = 0; iclu < nMax; iclu++) {
      Cluster* clu = static_cast<Cluster*>(fClustersArray->UncheckedAt(fNumberOfClusters + iclu));
      clu->AddDigit(digit, eInClusters[idig][iclu]); // Fills geometry and MC infor from Digit,+ correct energy
    }
  }
  fNumberOfClusters += nMax;
}
//__________________________________________________________________________
void Clusterizer::EvalClusters()
{
  // Calculate cluster properties
  int n = fClustersArray->GetEntriesFast();

  for (int i = 0; i < n; i++) {
    Cluster* clu = static_cast<Cluster*>(fClustersArray->UncheckedAt(i));
    // Eval all variables: Energy, CoreEnergy, position, Dispersion,...
    clu->EvalAll();
  }
}
