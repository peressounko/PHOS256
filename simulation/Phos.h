#ifndef PHOS_H
#define PHOS_H

//_________________________________________________________________________
// Implementation PHOS calorimeter
#include <vector>
#include <map>

// --- ROOT system ---

// --- PHOS256 header files ---
#include "Hit.h"

// Materials for all volumes, including Hall, PHOS256 and Magnet
enum materials { kAir = 1,
                 kPWO,
                 kAl,
                 kTyvek,
                 kFoam,
                 kSi,
                 kPCB,
                 kSteel,
                 kFiberglas,
                 kCables };

class Phos
{

 public:
  Phos() = default;

  // Will create calorimeter at distance r to crystal surface and at polar ange theta
  // with normal pointing to interaction point
  Phos(float r, float theta);
  ~Phos() = default;

  void CreateConverter(float d = 10., float h = 2.)
  {
    fCreateConterter = true;
    fConvDist = d;
    fConvThickness = h;
  }

  // prepare materials, mixtures and media for calorimeter
  void CreateMaterials();

  // construct GEANT geometry
  void CreateGeometry(); // creates the geometry for GEANT

  void SetHitContainer(std::vector<Hit>* hits) { fHits = hits; }

  // Tracking
  void Reset();
  bool ProcessHits();
  void FinishEvent();
  void FinishPrimary();

 protected:
  void CreateGeometryforEMC(void);    // creates the PHOS geometry for GEANT
  void CreateConverterGeometry(void); // creates converter in fronnt of PHOS
  void AddHit(int primary, int id, float* hits);

 private:
  bool fCreateConterter = false; // do not create converted by default
  float fConvDist = 10.;         // distance from crystall surface to converter
  float fConvThickness = 2.;
  // Simulation
  std::map<int, int> fSuperParents;  //! map of current tracks to SuperParents: entered PHOS active volumes particles
  std::vector<Hit>* fHits = nullptr; //! Collection of PHOS hits
  int fCurrentTrackID = 0;           //! current track Id
  int fCurrentCellID = 0;            //! current cell Id
  int fCurentSuperParent = 0;        //! current SuperParent ID: particle entered PHOS
  Hit* fCurrentHit = nullptr;        //! current Hit

  ClassDefNV(Phos, 1) // Implementation of PHOS manager class for layout EMC+PPSD
};

#endif // Phos_H
