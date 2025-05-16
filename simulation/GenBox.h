#ifndef GENBOX_H
#define GENBOX_H

#include <TVirtualMCApplication.h>

class TVirtualMCStack;
class TVector3;

class A01DetectorConstruction;

class GenBox : public TObject
{
 public:
  GenBox(TVirtualMCStack* stack);
  GenBox(const GenBox& origin, TVirtualMCStack* stack);
  GenBox() = default;
  ~GenBox() = default;

  // methods
  void Generate();

  // set methods
  void SetNofPrimaries(int nofPrimaries) { fNofPrimaries = nofPrimaries; }
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
  virtual void SetPart(int part) { fIpart = part; }

  void SetStack(TVirtualMCStack* stack) { fStack = stack; }

 private:
  // data members
  TVirtualMCStack* fStack = nullptr; ///< VMC stack
  int fNofPrimaries = 0;             ///< Number of primary particles
  int fIpart = 22;                   ///< Default particle PDG
  double fPtMin = 0.1;               ///< Default particle momentum
  double fPtMax = 1.;                ///< Default particle momentum
  double fThetaMin = 10.;            ///< degree
  double fThetaMax = 20.;            ///< degree
  double fPhiMin = 80.;              ///< degree
  double fPhiMax = 100.;             ///< degree

  ClassDef(GenBox, 1) // GenBox
};
#endif // GENBOX_H
