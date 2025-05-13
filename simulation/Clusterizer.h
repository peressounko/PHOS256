//-----------------------------------------------------------
//
// Author List:
//          D.Peresunko, KI, 2025
//-----------------------------------------------------------

#ifndef CLUSTERIZER_H
#define CLUSTERIZER_H

// Base Class Headers ----------------
#include <map>

class TClonesArray;
class TObjArray;
class MpdEmcSimParams;
class MpdEmcGeoUtils;
class MpdEmcClusterKI;
class MpdEmcCalibParams;

class Clusterizer
{
 public:
  // Constructors/Destructors ---------
  Clusterizer() = default;
  ~Clusterizer() = default;

  void Init();

  void ProcessEvent();

  // functions to be used in special case of re-reconstruction
  void SetDigits(TClonesArray* digits) { fDigitsArray = digits; }
  TObjArray* Clusters() { return fClustersArray; }
  void SetCalibParams(CalibParams* calib) { fCalibData = calib; }

 protected:
  void PrepareDigits();  // Calibrate, Allpy BadMap, clean...
  void MakeClusters();   // Do the job
  void MakeUnfoldings(); // Find and unfold clusters with few local maxima
  void UnfoldOneCluster(MpdEmcClusterKI* iniClu, Int_t nMax, int* digitId, float* maxAtEnergy);
  // Performs the unfolding of a cluster with nMax overlapping showers
  // Parameters: iniClu cluster to be unfolded
  //             nMax number of local maxima found (this is the number of new clusters)
  //             digitId: index of digits, corresponding to local maxima
  //             maxAtEnergy: energies of digits, corresponding to local maxima
  void EvalClusters();

 private:
  TClonesArray* fDigitsArray = nullptr; //! Input digits array
  TObjArray* fClustersArray = nullptr;  //! output clusters array
  SimParams* fSimParams = nullptr;      //! Configuration parameters
  Geometry* fGeom = nullptr;            //! Geometry class
  CalibParams* fCalibData = nullptr;    //! Calibration parameters
  ClassDef(Clusterizer, 1);
};

#endif
