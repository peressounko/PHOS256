/// Class to perform reconstruction of raw data
///

#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include <vector>
#include <string>
#include "TClonesArray.h"
#include "TTree.h"
#include "TFile.h"

#include "Cell.h"
#include "Clusterizer.h"
#include "BadChannelsMap.h"
#include "Fitter.h"
#include "RawReader.h"

class Reconstruction
{

 public:
  /// \brief Constructor
  Reconstruction() = default;

  /// \brief Destructor
  virtual ~Reconstruction() = default;

  // Main function, performs reconstruction
  // parameter defines when reconstruction should be stopped:
  // d: after Digits production
  // c: after Cluster production
  void Run(std::string option = "c");

  void SetRawFileName(std::string rawfilename) { fRawfilename = rawfilename; }

 protected:
  void Init();
  void Reset();
  void WriteOutput();

  bool MakeDigits();

  float Calibrate(int absId, bool isHG, float amp);
  float CalibrateT(int absId, bool isHG, float time);

 private:
  static constexpr int NCELLS = 256;

  int fEventSelection = 7;

  bool fMakeClu = false;

  std::string fRawfilename = "";

  // processors
  RawReader* fRawReader = nullptr;
  Clusterizer* fClusterizer = nullptr;

  BadChannelsMap* fBadMap = nullptr;

  Fitter* fFitter = nullptr;

  // Data collections
  std::vector<Cell> fCells; //! transient list of cells

  TClonesArray* fDigits = nullptr; //
  TObjArray* fClusters = nullptr;

  TTree* fTree = nullptr;
  TFile* fOutFile = nullptr;

  ClassDef(Reconstruction, 1);
}; // End of Reconstruction

#endif