#ifndef PHOS_H
#define PHOS_H

//_________________________________________________________________________
// Implementation PHOS calorimeter
#include <vector>
#include <map>

// --- ROOT system ---

// --- PHOS256 header files ---
#include "Hit.h"

class Phos
{

 public:
  Phos() = default;

  // Will create calorimeter at distance r to crystal surface and at polar ange theta
  // with normal pointing to interaction point
  Phos(float r, float theta);
  ~Phos() = default;

  // prepare materials, mixtures and media for calorimeter
  void CreateMaterials();

  // construct GEANT geometry
  void CreateGeometry(); // creates the geometry for GEANT

  // Tracking
  void Reset();
  bool ProcessHits();
  void FinishEvent();
  void FinishPrimary();

 protected:
  void CreateGeometryforEMC(void);     // creates the PHOS geometry for GEANT
  void CreateGeometryforSupport(void); // creates the Support geometry for GEANT
  void AddHit(int primary, int id, float* hits);

 private:
  // Simulation
  std::map<int, int> fSuperParents;  //! map of current tracks to SuperParents: entered PHOS active volumes particles
  std::vector<Hit>* fHits = nullptr; //! Collection of PHOS hits
  int fCurrentTrackID = 0;           //! current track Id
  int fCurrentCellID = 0;            //! current cell Id
  int fCurentSuperParent = 0;        //! current SuperParent ID: particle entered PHOS
  Hit* fCurrentHit = nullptr;        //! current Hit

  // Geometry and GEANT
  int fIdtmed[21]; //! GEANT media
  int fIdmate[21]; //! GEANT materials

  ClassDefNV(Phos, 1) // Implementation of PHOS manager class for layout EMC+PPSD
};

#endif // Phos_H
