#ifndef GENMUON_H
#define GENMUON_H

#include <TVirtualMCApplication.h>
#include <GenBox.h>

class TVirtualMCStack;
class TVector3;

class GenMuon : public GenBox
{
 public:
  GenMuon() = default;
  GenMuon(TVirtualMCStack* stack);
  //  GenMuon(const GenMuon& origin, TVirtualMCStack* stack);
  virtual ~GenMuon() = default;

  // methods
  virtual void Generate();

  // set methods
  void SetFirstTrigger(double x0 = 0., double y0 = 0, double dx = 16., double dy = 16)
  {
    fTr1X0 = x0;
    fTr1Y0 = y0;
    fTr1Dx = dx;
    fTr1Dy = dy;
  }
  void SetSecondTrigger(double x0 = 0., double y0 = 0, double dx = 16., double dy = 16)
  {
    fUseSecondTr = true;
    fTr2X0 = x0;
    fTr2Y0 = y0;
    fTr2Dx = dx;
    fTr2Dy = dy;
  }

 protected:
  static constexpr double kEmin = 0.04;
  static constexpr double kE0 = 4.29; // see M.Martemianov presentation
  static constexpr double knEn = 2.01;

  // Trigger parameters
  double fIntMax = 0.;
  double fTr1X0 = 0.;
  double fTr1Y0 = 0.;
  double fTr1Dx = 18.1;
  double fTr1Dy = 18.1;
  bool fUseSecondTr = false;
  double fDistTrigZ = 115.; // Distance between triggers
  double fTr2X0 = 0.;
  double fTr2Y0 = 0.;
  double fTr2Dx = 18.1;
  double fTr2Dy = 18.1;

  ClassDef(GenMuon, 1) // GenMuon
};
#endif // GENBOX_H
