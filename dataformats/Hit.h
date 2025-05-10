#ifndef HIT_H
#define HIT_H 1
// Hit - energy deposited by one particle entering calorimeter in one cell
// that is products of EM shower are assigned to one parent particle

#include "TObject.h"

class Hit : public TObject
{
 public:
  /** Default constructor **/
  Hit() = default;

  // cellID: calorimeter cell number
  // e: deposited energy (GeV)
  // t: time in seconds
  // label: index of parent particle
  Hit(int cellID, float e, float t, int label) : fE(e), fTime(t), fCellID(cellID), fLabel(label) {}

  ~Hit() = default;

  Hit& operator+=(Hit& h2)
  {
    if (fCellID == h2.fCellID) {
      fE += h2.fE;
      fTime = std::min(fTime, h2.fTime);
    }
    return *this;
  }

  bool operator==(const Hit& rhs) const { return ((fCellID == rhs.fCellID) && (fLabel == rhs.fLabel)); }

  bool operator<(const Hit& rhs) const
  {
    if (fCellID == rhs.fCellID) {
      return fLabel < rhs.fLabel;
    }
    return fCellID < rhs.fCellID;
  }

  // getters in alphabetinc order
  void AddEnergy(float e) { fE += e; }
  void SetEnergy(float e) { fE = e; }

  int GetCellID() const { return fCellID; }

  float GetE() const { return fE; };

  int GetLabel() const { return fLabel; };

  float GetTime() const { return fTime; };

 protected:
  float fE;    // energy
  float fTime; // hit mean time
  int fCellID; // detector id of each hit
  int fLabel;  // number of tracks, which have contribution in module
  ClassDefNV(Hit, 1);
};

#endif
