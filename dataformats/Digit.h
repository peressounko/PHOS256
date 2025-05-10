#ifndef DIGIT_H
#define DIGIT_H 1
// Digit is a sum of energy deposited by all particles in event in one cell of calorimeter (plus noise)

#include "TObject.h"
#include "Hit.h"

class Digit : public TObject
{
 public:
  /** Default constructor **/
  Digit() = default;

  Digit(Hit& h);

  Digit(int cellId, float energy, float time, int label);

  ~Digit() = default;

  Digit& operator+=(Digit& d);

  // Check, if point can be added (hits from the same Tower)
  bool CanAdd(Hit& h) const { return h.GetCellID() == fCellID; }

  // Adds point (add energy, change time if necessary, add primary)
  void AddHit(Hit& h);

  int GetCellId() const { return fCellID; }
  float GetE() const { return fE; }
  float GetTime() const { return fTime; };

  void SetE(float e) { fE = e; }
  void SetTime(float t) { fTime = t; }

  /// To allow sorting
  bool IsSortable() const { return true; }

  /// \brief Method ised for sorting Digits
  //  \param Another Digit
  //  \return
  int Compare(const TObject* obj) const;

  int GetNPrimaries() const { return fNprimary; }

  bool GetPrimary(int i, int& label, float& edep) const
  {
    if (i >= fNprimary) {
      label = -1;
      edep = 0;
      return false;
    } else {
      label = fPrimary[i];
      edep = fPrimEdep[i];
      return true;
    }
  }

 protected:
  static constexpr int kMaxLabels = 3; // number of primaries to store
  int fCellID = 0;                     // Cell index
  float fE = 0.;                       // Full energy
  float fTime = 0.;                    // minimal time
  int fNprimary = 0.;                  // Number of primaries
  int fPrimary[kMaxLabels] = {0};      // Array of primaries
  float fPrimEdep[kMaxLabels] = {0};   // Array of deposited energies

  ClassDefNV(Digit, 1);
};

#endif
