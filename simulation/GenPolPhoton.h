#ifndef GENPOLPHOTON_H
#define GENPOLPHOTON_H

#include "GenBox.h"

class GenPolPhoton : public GenBox
{
 public:
  GenPolPhoton(TVirtualMCStack* stack);
  GenPolPhoton(const GenPolPhoton& origin, TVirtualMCStack* stack);
  GenPolPhoton() = default;
  virtual ~GenPolPhoton() = default;

  // methods
  virtual void Generate();

  // set methods
  virtual void SetNofPrimaries(int nofPrimaries) { fNofPrimaries = nofPrimaries; }
  void SetMomentumRange(double ptMin, double ptMax)
  {
    fPtMin = ptMin;
    fPtMax = ptMax;
  }
  void SetThetaRange(double thetamin, double thetamax)
  {
    fThetaMin = thetamin;
    fThetaMax = thetamax;
  }
  void SetPhiRange(double phimin, double phimax)
  {
    fPhiMin = phimin;
    fPhiMax = phimax;
  }
  void SetPolarization(double phi) { fPolAngle = phi; }
  void SetUniformPolarization() { fUniformPol = true; }

  void SetStack(TVirtualMCStack* stack) { fStack = stack; }

 protected:
  // data members
  bool fUniformPol = false; ///< generate uniform polarization if true
  double fPolAngle = 90.;   ///< degree

  ClassDef(GenPolPhoton, 1) // GenPolPhoton
};
#endif // GENBOX_H
