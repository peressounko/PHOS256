#include <TMath.h>
#include "MagField.h"

//______________________________________________________________________________
MagField::MagField(Double_t Bx, Double_t By, Double_t Bz)
  : TVirtualMagField("Uniform magnetic field")
{
  /// Standard constructor
  /// \param Bx   The x component of the field value (in kiloGauss)
  /// \param By   The y component of the field value (in kiloGauss)
  /// \param Bz   The z component of the field value (in kiloGauss)

  fB[0] = Bx;
  fB[1] = By;
  fB[2] = Bz;
}

//______________________________________________________________________________
void MagField::Field(const Double_t* x, Double_t* B)
{
  /// Fill in the field value B in the given position at x.
  /// \param x   The position
  /// \param B   The field value (in kiloGauss)

  // G4VSolid* magneticSolid = new
  // G4Tubs("magneticTubs",0.,1.*m,1.*m,0.,360.*deg); G4VSolid* magneticSolid =
  // new G4Tubs("magneticTubs",0.,1.*m,1.*m,0.,360.*deg); G4RotationMatrix*
  // fieldRot = new G4RotationMatrix(); fieldRot->rotateX(90.*deg); new
  // G4PVPlacement(fieldRot,G4ThreeVector(),magneticLogical,
  //                  "magneticPhysical",worldLogical,0,0);

  // Magnetic field in the cube
  if (std::abs(x[0]) < 50. && std::abs(x[1]) < 50. && std::abs(x[2]) < 50.) {
    B[0] = fB[0];
    B[1] = fB[1];
    B[2] = fB[2];
  } else {
    B[0] = 0.;
    B[1] = 0.;
    B[2] = 0.;
  }
}
