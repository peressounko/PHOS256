////////////////////////////////////////////////////////////////
//                                                            //
//  Authors : D.Peresunko, KI                                 //
//                                                            //
////////////////////////////////////////////////////////////////

#include "Digit.h"

Digit::Digit(Hit& h)
  : fCellID(h.GetCellID()), fE(h.GetE()), fTime(h.GetTime())
{
  if (fE > 0) { // only primaries with non-zero energy deposition count
    fNprimary = 1;
    fPrimary[0] = h.GetLabel();
    fPrimEdep[0] = h.GetE();
  }
}

Digit::Digit(int cellID, float energy, float time, int label)
  : fCellID(cellID), fE(energy), fTime(time)
{
  if (label >= 0 && energy > 0) {
    fNprimary = 1;
    fPrimary[0] = label;
    fPrimEdep[0] = energy;
  }
}

Digit& Digit::operator+=(Digit& d)
{
  // Adds digits (add energy, change time if necessary, add primary)
  if (fCellID != d.fCellID) {
    return *this;
  }
  if (d.fE == 0.) { // do nothing
    return *this;
  }
  if (d.fE > fE) {
    fTime = d.fTime;
  }
  fE += d.fE;
  // Primaries Check if track already exist
  for (int ii = 0; ii < d.fNprimary; ii++) {
    int iprim = d.fPrimary[ii];
    float edep = d.fPrimEdep[ii];
    bool found = false;
    for (int i = 0; i < fNprimary; ++i) {
      if (fPrimary[i] == iprim) {
        fPrimEdep[i] += edep;
        found = true;
        break;
      }
    }
    if (!found) { // remove softest if necessary
      if (fNprimary < kMaxLabels - 1) {
        fPrimary[fNprimary] = iprim;
        fPrimEdep[fNprimary] = edep;
        fNprimary++;
      } else { // check last and replace
        if (edep > fPrimEdep[kMaxLabels - 1]) {
          fPrimary[kMaxLabels - 1] = iprim;
          fPrimEdep[kMaxLabels - 1] = edep;
        }
      }
    }
    // sort if necessary
    for (int i = kMaxLabels - 2; i--;) {
      if (fPrimEdep[i + 1] > fPrimEdep[i]) {
        float tmpE = fPrimEdep[i];
        fPrimEdep[i] = fPrimEdep[i + 1];
        fPrimEdep[i + 1] = tmpE;
        int l = fPrimary[i];
        fPrimary[i] = fPrimary[i + 1];
        fPrimary[i + 1] = l;
      }
    }
  }
  return *this;
}

void Digit::AddHit(Hit& h)
{
  // Adds point (add energy, change time if necessary, add primary)
  if (fCellID != h.GetCellID()) {
    return;
  }
  if (h.GetE() == 0.) { // do nothing
    return;
  }
  if (h.GetE() > fE) {
    fTime = h.GetTime();
  }
  fE += h.GetE();
  // Check if track already exist
  int iprim = h.GetLabel();
  bool found = false;
  for (int i = 0; i < fNprimary; ++i) {
    if (fPrimary[i] == iprim) {
      fPrimEdep[i] += h.GetE();
      found = true;
      break;
    }
  }
  if (!found) { // remove softest if necessary
    if (fNprimary < kMaxLabels - 1) {
      fPrimary[fNprimary] = h.GetLabel();
      fPrimEdep[fNprimary] = h.GetE();
      fNprimary++;
    } else { // check last and replace
      if (h.GetE() > fPrimEdep[kMaxLabels - 1]) {
        fPrimary[kMaxLabels - 1] = h.GetLabel();
        fPrimEdep[kMaxLabels - 1] = h.GetE();
      }
    }
  }
  // sort if necessary
  for (int i = kMaxLabels - 2; i--;) {
    if (fPrimEdep[i + 1] > fPrimEdep[i]) {
      float tmpE = fPrimEdep[i];
      fPrimEdep[i] = fPrimEdep[i + 1];
      fPrimEdep[i + 1] = tmpE;
      int l = fPrimary[i];
      fPrimary[i] = fPrimary[i + 1];
      fPrimary[i + 1] = l;
    }
  }
}

int Digit::Compare(const TObject* obj) const
{
  const Digit* rhs = dynamic_cast<const Digit*>(obj);
  if (!rhs) {
    return 1;
  }
  if (fCellID < rhs->fCellID) {
    return -1;
  } else {
    if (fCellID == rhs->fCellID) {
      return 0;
    } else {
      return 1;
    }
  }
}

ClassImp(Digit);