
// #include <TGeoGlobalMagField.h>
// #include <TVirtualMC.h>
// #include "G4Box.hh"
// #include "G4LogicalVolume.hh"
// #include "G4Material.hh"
// #include "G4PVPlacement.hh"
// #include "G4NistManager.hh"

#include <TGeoElement.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>

#include "Hall.h"
#include "Phos.h" //for material numbering

//_____________________________________________________________________________
void Hall::CreateGeometry()
{
  //
  // Create the geometry of the exprimental hall
  //

  if (!gGeoManager->GetMedium("AIR")) {

    Double_t a;       // Mass of a mole in g/mole
    Double_t z;       // Atomic number
    Double_t density; // Material density in g/cm3
    TGeoElement* elN = new TGeoElement("Nitrogen", "N", z = 7., a = 14.01);
    TGeoElement* elO = new TGeoElement("Oxygen", "O", z = 8., a = 16.00);
    TGeoMixture* matAir = new TGeoMixture("AIR", kAir, density = 1.29e-03);
    matAir->AddElement(elN, 0.7);
    matAir->AddElement(elO, 0.3);

    // Paremeter for tracking media
    Double_t param[20];
    param[0] = 0;     // isvol  - Not used
    param[1] = 0;     // ifield - User defined magnetic field
    param[2] = 10.;   // fieldm - Maximum field value (in kiloGauss)
    param[3] = -20.;  // tmaxfd - Maximum angle due to field deflection
    param[4] = -0.01; // stemax - Maximum displacement for multiple scat
    param[5] = -.3;   // deemax - Maximum fractional energy loss, DLS
    param[6] = .001;  // epsil - Tracking precision
    param[7] = -.8;   // stmin
    for (Int_t i = 8; i < 20; ++i)
      param[i] = 0.;

    new TGeoMedium("AIR", kAir, matAir, param);
  }

  // Cube 5*5*5 m
  float dHall[3] = {250., 250., 250.};
  TGeoVolume* top = gGeoManager->Volume("World", "BOX", 1, dHall, 3);
  gGeoManager->SetTopVolume(top);
}
