
#include <TGeoElement.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TVirtualMC.h>

#include "Magnet.h"

//_____________________________________________________________________________
void Magnet::CreateGeometry()
{
  //
  // Create the geometry of the magnet
  //
  // std::cout << "Magnet creating geometry" << std::endl;

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

  // --- Stainless steel (let it be pure iron) ---
  if(!gGeoManager->GetMedium("STEEL")){
    double a, z, density;
    TGeoMixture* mixSteel = new TGeoMixture("STEEL", 1, density = 7.87);
    TGeoElement* elFe = new TGeoElement("Ferrum", "Fe", z = 26., a = 55.845);
    mixSteel->AddElement(elFe, 1);
    new TGeoMedium("STEEL", 2, mixSteel, param);
  }
  int imed = gGeoManager->GetMedium("STEEL")->GetId();


  // Create top coil.
  float magX = 30.;  //vertical half size
  float magY = 30.;
  float magZ = 40.;
  float paramsBox[3] = {magX,magY,magZ};
  gGeoManager->Volume("MAGTOP", "BOX", imed, paramsBox, 3); // 20-> idtmed[20]

  float h = 30.; //half of the height of the gap

  float paramsTube[5] = {0.,   //innerRadius
                        4.,   //outerRadius 
                        h+2.*magX,//height
                        0.,   //starting angle
                        360.};// spanning angle degree

  gGeoManager->Volume("MAGCOL", "TUBE", imed, paramsTube, 5); // 20-> idtmed[20]


  // --- Position  Magnet in Hall ---
// Mother frame:
//       z axis along the beam
//       x axis up
//       y axis horizontal toward calorimeter

  int idrotm = 0;
  TVirtualMC::GetMC()->Matrix(idrotm, 90.,0.,    //Double_t thetaX, Double_t phiX,
                                      90.,90.,   //Double_t thetaY, Double_t phiY
                                      0.,0.);    //Double_t thetaZ, Double_t phiZ
  float dVtx = 10.;  //distance from vtx to front surface of magnet

  gGeoManager->Node("MAGTOP", 1, "World", h + magX, 0., dVtx + magZ, idrotm, true, static_cast<double*>(nullptr));
  gGeoManager->Node("MAGTOP", 2, "World",-h - magX, 0., dVtx + magZ, idrotm, true, static_cast<double*>(nullptr));

  //position columns
  TVirtualMC::GetMC()->Matrix(idrotm, 0.,0.,    //Double_t thetaX, Double_t phiX,
                                      90.,90.,   //Double_t thetaY, Double_t phiY
                                      90.,0.);    //Double_t thetaZ, Double_t phiZ  
  float dzCol = 11.;  //distance from back magnet to column
  float dyCol=5.;
  gGeoManager->Node("MAGCOL", 1, "World", 0., magY + dyCol,  dVtx + 2*magZ - dzCol, idrotm, true, static_cast<double*>(nullptr));
  gGeoManager->Node("MAGCOL", 2, "World", 0., -magY - dyCol, dVtx + 2*magZ - dzCol, idrotm, true, static_cast<double*>(nullptr));
  gGeoManager->Node("MAGCOL", 3, "World", 0., magY + dyCol,  dVtx + dzCol, idrotm, true, static_cast<double*>(nullptr));
  gGeoManager->Node("MAGCOL", 4, "World", 0., -magY - dyCol, dVtx + dzCol, idrotm, true, static_cast<double*>(nullptr));

  std::cout << "Created Magnet geometry" << std::endl;


}
