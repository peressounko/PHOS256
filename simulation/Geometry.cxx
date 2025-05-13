//_________________________________________________________________________
// Geometry class  for PHOS : singleton
// PHOS consists of the electromagnetic calorimeter (EMCA)
// Mother frame:
//       z axis along the beam
//       x axis up
//       y axis horizontal toward calorimeter

// --- ROOT system ---
// #include "TVector3.h"
// #include "TRotation.h"
// #include "TParticle.h"
// #include <TGeoManager.h>
// #include <TGeoMatrix.h>

// --- Standard library ---
#include <iostream> // std::cout
#include <math.h>
#include <cassert>

// --- PHOS256 header files ---
#include "Geometry.h"

Geometry* Geometry::fgGeom = nullptr;

void Geometry::GetModuleAngles(float angle[3][2]) const
{
  // angle of the module positioned so than normal vector is in horizontal plane and polar angle theta to z direction
  // A rotation matrix is described to GEANT by giving the polar and azimuthal angles of the axes of the DRS
  //  (x′, y′, z′) in the MRS
  angle[0][0] = fModuleAngle[0][0];
  angle[0][1] = fModuleAngle[0][1];
  angle[1][0] = fModuleAngle[1][0];
  angle[1][1] = fModuleAngle[1][1];
  angle[2][0] = fModuleAngle[2][0];
  angle[2][1] = fModuleAngle[2][1];
}
void Geometry::GetModuleCenter(float pos[3]) const
{
  // Distance from IP to modeult center
  pos[0] = fModuleCenter[0];
  pos[0] = fModuleCenter[1];
  pos[0] = fModuleCenter[2];
}

//____________________________________________________________________________
Geometry::Geometry(float r, float theta) : fModR(r), fModTheta(theta)
{
  Init();

  fgGeom = this;
}

//____________________________________________________________________________
const Geometry* Geometry::Instance(float r, float theta)
{
  // Returns the pointer of the unique instance
  if (fgGeom && (fgGeom->fModR != r || fgGeom->fModTheta != theta)) {
    std::cout << "Geometry already initialized with another geometry parameters:" << std::endl;
    std::cout << " Old r=" << fgGeom->fModR << ", old theta=" << fgGeom->fModTheta << std::endl;
    std::cout << " New r=" << r << ", new theta=" << theta << std::endl;
    assert(false);
  } else {
    if (fgGeom) {
      return fgGeom;
    } else {
      fgGeom = new Geometry(r, theta);
    }
  }
  return fgGeom;
}

int Geometry::RelToAbsId(int moduleNumber, int strip, int cell)
{
  // calculates absolute cell Id from moduleNumber, strip (number) and cell (number)
  // PHOS layout parameters:
  const int nStrpZ = 8;                  // Number of strips along z-axis
  const int nCrystalsInModule = 16 * 16; // Total number of crystals in module
  const int nCellsXInStrip = 8;          // Number of crystals in strip unit along x-axis
  const int nZ = 16;                     // nStripZ * nCellsZInStrip

  int row = nStrpZ - (strip - 1) % nStrpZ;
  int col = (int)std::ceil((float)strip / (nStrpZ)) - 1;

  return (moduleNumber - 1) * nCrystalsInModule + row * 2 + (col * nCellsXInStrip + (cell - 1) / 2) * nZ -
         (cell & 1 ? 1 : 0);
}

//____________________________________________________________________________
void Geometry::RelToAbsId(const int relid[2], Int_t& absId)
{
  // Converts the relative numbering into the absolute numbering
  // PHOS crystals:
  //  absId = from 1 to fNModules * fNPhi * fNZ
  absId = (relid[0] - 1) * fNZ // the offset along phi
          + relid[1];          // the offset along z
}

// //____________________________________________________________________________
// void Geometry::GetGlobalPHOS(const AliPHOSRecPoint* recPoint, TVector3 & gpos) const
// {
//   // Calculates the coordinates of a RecPoint and the error matrix in the ALICE global coordinate system

//   const AliPHOSRecPoint * tmpPHOS = recPoint ;
//   TVector3 localposition ;

//   tmpPHOS->GetLocalPosition(gpos) ;

//   if (!gGeoManager){
//     AliFatal("Geo manager not initialized\n");
//   }
//   //construct module name
//   TGeoHMatrix *m = 0x0;
//   char path[100] ;
//   Double_t dy ;
//   if(tmpPHOS->IsEmc()){
//     snprintf(path,100,"/ALIC_1/PHOS_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1",tmpPHOS->GetPHOSMod()) ;
//     if (!gGeoManager->CheckPath(path)){
//       snprintf(path,100,"/ALIC_1/PHOC_%d/PEMC_1/PCOL_1/PTIO_1/PCOR_1/PAGA_1/PTII_1",tmpPHOS->GetPHOSMod()) ;
//       if (!gGeoManager->CheckPath(path)){
//         snprintf(path,100,"/ALIC_1/PHOH_%d/PEMH_1/PCLH_1/PIOH_1/PCOH_1/PAGH_1/PTIH_1",tmpPHOS->GetPHOSMod()) ;
//         if(!gGeoManager->CheckPath(path)){
//           AliFatal("Geo manager can not find path \n");
// 	}
//       }
//     }
//     gGeoManager->cd(path) ;
//     m = gGeoManager->GetCurrentMatrix();
//     dy=fCrystalShift ;
//   }
//   else{
//     snprintf(path,100,"/ALIC_1/PHOC_%d/PCPV_1",tmpPHOS->GetPHOSMod());
//     if (!gGeoManager->CheckPath(path)){
//       snprintf(path,100,"/ALIC_1/PHOH_%d/PCPV_1",tmpPHOS->GetPHOSMod());
//       if (!gGeoManager->CheckPath(path))
//         AliFatal(Form("Geo manager can not find path /ALIC_1/PHOC(H)_%d/PCPV_1 \n",tmpPHOS->GetPHOSMod()));
//     }
//     gGeoManager->cd(path) ;
//     m = gGeoManager->GetCurrentMatrix();
//     dy= GetCPVBoxSize(1)/2. ; //center of CPV module
//   }
//   Double_t pos[3]={gpos.X(),gpos.Y()-dy,gpos.Z()} ;
//   if(tmpPHOS->IsEmc())
//     pos[2]=-pos[2] ; //Opposite z directions in EMC matrix and local frame!!!
//   Double_t posC[3] = {};
//   //now apply possible shifts and rotations
//   if (m){
//      m->LocalToMaster(pos,posC);
//   }
//   else{
//     AliFatal("Geo matrixes are not loaded \n") ;
//   }
//   gpos.SetXYZ(posC[0],posC[1],posC[2]) ;

// }
// //____________________________________________________________________________

// void Geometry::GetModuleCenter(TVector3& center,
// 				      const char *det,
// 				      Int_t module) const
// {
//   // Returns a position of the center of the CPV or EMC module
//   // in ideal (not misaligned) geometry
//   float rDet = 0.;
//   if      (strcmp(det,"CPV") == 0) rDet  = GetIPtoCPVDistance   ();
//   else if (strcmp(det,"EMC") == 0) rDet  = GetIPtoCrystalSurface();
//   else
//     AliFatal(Form("Wrong detector name %s",det));

//   float angle = GetPHOSAngle(module); // (40,20,0,-20,-40) degrees
//   angle *= TMath::Pi()/180;
//   angle += 3*TMath::Pi()/2.;
//   center.SetXYZ(rDet*TMath::Cos(angle), rDet*TMath::Sin(angle), 0.);
// }

void Geometry::Init()
{
  // Initializes the EMC parameters
  // Coordinate system chosen: x across beam, z along beam, y out of beam.
  // Reference point for all volumes incide module is
  // center of module in x,z on the upper surface of support beam

  // Distance from IP to surface of the crystals
  fIPtoCrystalSurface = 460.0;

  // CRYSTAL

  fCrystalHalfSize[0] = 2.2 / 2; // Half-Sizes of crystall
  fCrystalHalfSize[1] = 18.0 / 2;
  fCrystalHalfSize[2] = 2.2 / 2;

  // APD + preamplifier

  // fPinDiodeSize[0] = 1.71 ;   //Values of ame PIN diode
  // fPinDiodeSize[1] = 0.0280 ; // OHO 0.0280 is the depth of active layer
  // fPinDiodeSize[2] = 1.61 ;

  fPinDiodeHalfSize[0] = 0.5000 / 2; // APD 5 mm side
  fPinDiodeHalfSize[1] = 0.0100 / 2; // APD bulk thickness
  fPinDiodeHalfSize[2] = 0.5000 / 2; // APD 5 mm side

  fPreampHalfSize[0] = 1.5 / 2; // Preamplifier
  fPreampHalfSize[1] = 0.5 / 2;
  fPreampHalfSize[2] = 1.5 / 2;

  // Strip unit (8x2 crystals)

  fNCellsXInStrip = 8; // Number of crystals in strip unit along x-axis
  fNCellsZInStrip = 2; // Number of crystals in strip unit along z-axis
  fNStripX = 2;        // Number of strip units across along x-axis
  fNStripZ = 8;        // Number of strips along z-axis

  fStripWallWidthOut = 0.01; // Side to another strip
  fStripWallWidthIn = 0.02;  // Side betveen crystals in one strip

  fTyvecThickness = 0.0175; // Thickness of the tyvec

  fAirGapLed = 1.5 - 2 * fPreampHalfSize[1] - 2 * fPinDiodeHalfSize[1]; // Air gap before crystalls for LED system
                                                                        // Note, that Cell in Strip 1.5 longer then crystall

  //---Now calculate thechnical sizes for GEANT implementation

  fWrappedHalfSize[0] = (2 * fTyvecThickness + 2 * fCrystalHalfSize[0]) / 2; // This will be size of crystall
  fWrappedHalfSize[1] = fCrystalHalfSize[1];                                 // wrapped into tyvec
  fWrappedHalfSize[2] = (2 * fTyvecThickness + 2 * fCrystalHalfSize[2]) / 2; //

  fAirCellHalfSize[0] = fWrappedHalfSize[0] + 0.01;
  fAirCellHalfSize[1] = (fAirGapLed + 2 * fPreampHalfSize[1] +
                         2 * fPinDiodeHalfSize[1] + 2 * fWrappedHalfSize[1]) /
                        2; // in strip
  fAirCellHalfSize[2] = fWrappedHalfSize[2] + 0.01;

  //  fSupportPlateHalfSize[0] = ( (fNCellsXInStrip-1)*fStripWallWidthIn + 2*fStripWallWidthOut +
  //             fNCellsXInStrip * (2*fTyvecThickness + 2*fCrystalHalfSize[0]) )/2 ;
  fSupportPlateHalfSize[0] = 18.04 / 2;
  fSupportPlateHalfSize[1] = 6.0 / 2;
  //  fSupportPlateHalfSize[2] = ( (fNCellsZInStrip-1)*fStripWallWidthIn + 2*fStripWallWidthOut +
  //             fNCellsZInStrip * (2*fTyvecThickness + 2*fCrystalHalfSize[2]) )/2;
  fSupportPlateHalfSize[2] = 4.51 / 2;
  fSupportPlateThickness = 0.3;
  fSupportPlateInHalfSize[0] = fSupportPlateHalfSize[0];                          // Half-sizes of the air
  fSupportPlateInHalfSize[1] = fSupportPlateHalfSize[1] - fSupportPlateThickness; // box in the support plate
  fSupportPlateInHalfSize[2] = fSupportPlateHalfSize[2] - fSupportPlateThickness / 2;

  fStripHalfSize[0] = fSupportPlateHalfSize[0];
  fStripHalfSize[1] = (2 * fSupportPlateHalfSize[1] + 2 * fAirCellHalfSize[1]) / 2;
  fStripHalfSize[2] = fSupportPlateHalfSize[2];

  // ------- Inner hermoinsulation ---------------
  fInnerThermoWidthX = 2.0; // Width of the innerthermoinsulation across the beam
  fInnerThermoWidthY = 2.0; // Width of the upper cover of innerthermoinsulation
  fInnerThermoWidthZ = 2.0; // Width of the innerthermoinsulation along the beam

  fInnerThermoHalfSize[0] = (2 * fStripHalfSize[0] * fNStripX + 2 * fInnerThermoWidthX) / 2;
  fInnerThermoHalfSize[1] = (2 * fStripHalfSize[1] + fInnerThermoWidthY) / 2;
  fInnerThermoHalfSize[2] = (2 * fStripHalfSize[2] * fNStripZ + 2 * fInnerThermoWidthZ) / 2;

  // ------- Air gap between inner thermoinsulation and passive coller ---------

  fAirGapWidthX = 0.2; // Width of the air gap across the beam
  fAirGapWidthY = 0.2; // Width of the upper air gap
  fAirGapWidthZ = 0.2; // Width of the air gap along the beam

  fAirGapHalfSize[0] = (2 * fInnerThermoHalfSize[0] + 2 * fAirGapWidthX) / 2;
  fAirGapHalfSize[1] = (2 * fInnerThermoHalfSize[1] + fAirGapWidthY) / 2;
  fAirGapHalfSize[2] = (2 * fInnerThermoHalfSize[2] + 2 * fAirGapWidthZ) / 2;

  // ------- Passive Cooler ------------------------

  fCoolerWidthX = 2.0; // Width of the passive coller across the beam
  fCoolerWidthY = 0.3; // Width of the upper cover of cooler
  fCoolerWidthZ = 2.0; // Width of the passive cooler along the beam

  fCoolerHalfSize[0] = (2 * fAirGapHalfSize[0] + 2 * fCoolerWidthX) / 2;
  fCoolerHalfSize[1] = (2 * fAirGapHalfSize[1] + fCoolerWidthY) / 2;
  fCoolerHalfSize[2] = (2 * fAirGapHalfSize[2] + 2 * fCoolerWidthZ) / 2;

  // ------- Outer thermoinsulation and Al cover -------------------------------

  fAlCoverThickness = 0.1; // Thickness of the Al cover of the module

  fOuterThermoWidthXUp = 156.0 - fAlCoverThickness;
  // width of the upper surface of the PHOS module accross the beam
  fOuterThermoWidthY = 6.0; // with of the upper cover of outer thermoinsulation
  fOuterThermoWidthZ = 6.0; // width of the thermoinsulation along the beam

  fAlFrontCoverX = 6.0; // Width of Al strip around fiberglass window: across
  fAlFrontCoverZ = 6.0; // and along the beam

  // Calculate distance from IP to upper cover
  fIPtoOuterCoverDistance = fIPtoCrystalSurface - fAirGapLed - fInnerThermoWidthY - fAirGapWidthY -
                            fCoolerWidthY - fOuterThermoWidthY - fAlCoverThickness;

  float tanA = fOuterThermoWidthXUp / (2. * fIPtoOuterCoverDistance);
  // tan(a) where A = angle between IP to center and IP to side across beam

  fOuterThermoWidthXLow = fOuterThermoWidthXUp +
                          2 * (2 * fCoolerHalfSize[1] + fOuterThermoWidthY) * tanA - fAlCoverThickness;
  // width of the lower surface of the COOL section accross the beam

  fOuterThermoParams[0] = fOuterThermoWidthXUp / 2;  // half-length along x at the z surface positioned at -DZ;
  fOuterThermoParams[1] = fOuterThermoWidthXLow / 2; // half-length along x at the z surface positioned at +DZ;
  fOuterThermoParams[2] = (2 * fCoolerHalfSize[2] + 2 * fOuterThermoWidthZ) / 2;
  // `half-length along the y-axis' in out case this is z axis
  fOuterThermoParams[3] = (2 * fCoolerHalfSize[1] + fOuterThermoWidthY) / 2;
  // `half-length along the z-axis' in our case this is y axis

  fAlCoverParams[0] = fOuterThermoParams[0] + fAlCoverThickness;
  fAlCoverParams[1] = fOuterThermoParams[1] + fAlCoverThickness;
  fAlCoverParams[2] = fOuterThermoParams[2] + fAlCoverThickness;
  fAlCoverParams[3] = fOuterThermoParams[3] + fAlCoverThickness / 2;

  fFiberGlassHalfSize[0] = fAlCoverParams[0] - fAlFrontCoverX;
  fFiberGlassHalfSize[1] = fAlCoverParams[2] - fAlFrontCoverZ; // Note, here other ref. system
  fFiberGlassHalfSize[2] = fAlCoverThickness / 2;

  //============Now warm section======================
  // Al Cover
  fWarmAlCoverWidthX = 2 * fAlCoverParams[1]; // Across beam
  fWarmAlCoverWidthY = 159.0;                 // along beam

  // T-support
  fTSupport1Thickness = 3.5;
  fTSupport2Thickness = 5.0;
  fTSupport1Width = 10.6;
  fTSupport2Width = 3.1;
  fNTSupports = fNStripX + 1;
  fTSupportDist = 7.48;

  // Air space for FEE
  fAirSpaceFeeX = 148.6; // Across beam
  fAirSpaceFeeY = 135.0; // along beam
  fAirSpaceFeeZ = 19.0;  // out of beam

  // thermoinsulation
  fWarmBottomThickness = 4.0;
  fWarmUpperThickness = 4.0;

  // Frame
  fFrameThickness = 5.0;
  fFrameHeight = 15.0;

  // Fiberglass support
  fFiberGlassSup1X = 6.0;
  fFiberGlassSup1Y = 3.9 + fWarmUpperThickness;

  fFiberGlassSup2X = 3.0;
  fFiberGlassSup2Y = fFrameHeight;

  // Now calculate Half-sizes

  fWarmAlCoverWidthZ = fAirSpaceFeeZ + fWarmBottomThickness + fWarmUpperThickness +
                       fTSupport1Thickness + fTSupport2Thickness;

  fWarmAlCoverHalfSize[0] = fWarmAlCoverWidthX / 2;
  fWarmAlCoverHalfSize[1] = fWarmAlCoverWidthY / 2;
  fWarmAlCoverHalfSize[2] = fWarmAlCoverWidthZ / 2;

  fWarmThermoHalfSize[0] = fWarmAlCoverHalfSize[0] - fAlCoverThickness;
  fWarmThermoHalfSize[1] = fWarmAlCoverHalfSize[1] - fAlCoverThickness;
  fWarmThermoHalfSize[2] = fWarmAlCoverHalfSize[2] - fAlCoverThickness / 2;

  // T-support
  fTSupport1HalfSize[0] = fTSupport1Width / 2;                        // Across beam
  fTSupport1HalfSize[1] = (fAirSpaceFeeY + 2 * fFiberGlassSup1X) / 2; // along beam
  fTSupport1HalfSize[2] = fTSupport1Thickness / 2;                    // out of beam

  fTSupport2HalfSize[0] = fTSupport2Width / 2;     // Across beam
  fTSupport2HalfSize[1] = fTSupport1HalfSize[1];   // along beam
  fTSupport2HalfSize[2] = fTSupport2Thickness / 2; // out of beam

  // cables
  fTCables1HalfSize[0] = (2 * fTSupport1HalfSize[0] * fNTSupports + (fNTSupports - 1) * fTSupportDist) / 2; // Across beam
  fTCables1HalfSize[1] = fTSupport1HalfSize[1];                                                             // along beam
  fTCables1HalfSize[2] = fTSupport1HalfSize[2];                                                             // out of beam

  fTCables2HalfSize[0] = fTCables1HalfSize[0];  // Across beam
  fTCables2HalfSize[1] = fTSupport2HalfSize[1]; // along beam
  fTCables2HalfSize[2] = fTSupport2HalfSize[2]; // out of beam

  // frame: we define two frames along beam ...Z and across beam ...X
  fFrameXHalfSize[0] = (fAirSpaceFeeX + 2 * fFiberGlassSup2X + 2 * fFrameThickness) / 2;
  fFrameXHalfSize[1] = fFrameThickness / 2;
  fFrameXHalfSize[2] = fFrameHeight / 2;

  fFrameXPosition[0] = 0;
  fFrameXPosition[1] = fAirSpaceFeeY / 2 + fFiberGlassSup2X + fFrameXHalfSize[1];
  fFrameXPosition[2] = fWarmThermoHalfSize[2] - fFrameHeight / 2 - fWarmBottomThickness;

  fFrameZHalfSize[0] = fFrameThickness / 2;
  fFrameZHalfSize[1] = (fAirSpaceFeeY + 2 * fFiberGlassSup2X) / 2;
  fFrameZHalfSize[2] = fFrameHeight / 2;

  fFrameZPosition[0] = fAirSpaceFeeX / 2 + fFiberGlassSup2X + fFrameZHalfSize[0];
  fFrameZPosition[1] = 0;
  fFrameZPosition[2] = fWarmThermoHalfSize[2] - fFrameHeight / 2 - fWarmBottomThickness;

  // Fiberglass support define 4 fiber glass supports 2 along Z  and 2 along X

  fFGupXHalfSize[0] = fFrameXHalfSize[0];
  fFGupXHalfSize[1] = fFiberGlassSup1X / 2;
  fFGupXHalfSize[2] = fFiberGlassSup1Y / 2;

  fFGupXPosition[0] = 0;
  fFGupXPosition[1] = fAirSpaceFeeY / 2 + fFGupXHalfSize[1];
  fFGupXPosition[2] = fWarmThermoHalfSize[2] - fFrameHeight - fWarmBottomThickness - fFGupXHalfSize[2];

  fFGupZHalfSize[0] = fFiberGlassSup1X / 2;
  fFGupZHalfSize[1] = fAirSpaceFeeY / 2;
  fFGupZHalfSize[2] = fFiberGlassSup1Y / 2;

  fFGupZPosition[0] = fAirSpaceFeeX / 2 + fFGupZHalfSize[0];
  fFGupZPosition[1] = 0;
  fFGupZPosition[2] = fWarmThermoHalfSize[2] - fFrameHeight - fWarmBottomThickness - fFGupXHalfSize[2];

  fFGlowXHalfSize[0] = fFrameXHalfSize[0] - 2 * fFrameZHalfSize[0];
  fFGlowXHalfSize[1] = fFiberGlassSup2X / 2;
  fFGlowXHalfSize[2] = fFrameXHalfSize[2];

  fFGlowXPosition[0] = 0;
  fFGlowXPosition[1] = fAirSpaceFeeY / 2 + fFGlowXHalfSize[1];
  fFGlowXPosition[2] = fWarmThermoHalfSize[2] - fWarmBottomThickness - fFGlowXHalfSize[2];

  fFGlowZHalfSize[0] = fFiberGlassSup2X / 2;
  fFGlowZHalfSize[1] = fAirSpaceFeeY / 2;
  fFGlowZHalfSize[2] = fFrameZHalfSize[2];

  fFGlowZPosition[0] = fAirSpaceFeeX / 2 + fFGlowZHalfSize[0];
  fFGlowZPosition[1] = 0;
  fFGlowZPosition[2] = fWarmThermoHalfSize[2] - fWarmBottomThickness - fFGlowXHalfSize[2];

  // --- Air Gap for FEE ----

  fFEEAirHalfSize[0] = fAirSpaceFeeX / 2;
  fFEEAirHalfSize[1] = fAirSpaceFeeY / 2;
  fFEEAirHalfSize[2] = fAirSpaceFeeZ / 2;

  fFEEAirPosition[0] = 0;
  fFEEAirPosition[1] = 0;
  fFEEAirPosition[2] = fWarmThermoHalfSize[2] - fWarmBottomThickness - fFEEAirHalfSize[2];

  // --- Calculate the oveol dimentions of the EMC module

  fEMCParams[3] = fAlCoverParams[3] + fWarmAlCoverHalfSize[2];                                                     // Size out of beam
  fEMCParams[0] = fAlCoverParams[0];                                                                               // Upper size across the beam
  fEMCParams[1] = (fAlCoverParams[1] - fAlCoverParams[0]) * fEMCParams[3] / fAlCoverParams[3] + fAlCoverParams[0]; // Lower size across the beam
  fEMCParams[2] = fWarmAlCoverHalfSize[1];                                                                         // Size along the beam

  // fNPhi = fNStripX * fNCellsXInStrip; // number of crystals across the beam
  // fNZ = fNStripZ * fNCellsZInStrip;   // number of crystals along the beam

  // calculate offset to crystal surface
  fCrystalShift = -fInnerThermoHalfSize[1] + fStripHalfSize[1] + fSupportPlateHalfSize[1] +
                  fCrystalHalfSize[1] - fAirGapLed / 2. + fPinDiodeHalfSize[1] + fPreampHalfSize[1];

  // Geometry parameters are calculated

  const double kDEGRad = M_PI / 180.0;
  // angle of the module positioned so than normal vector is in horizontal plane and polar angle theta to z direction
  // A rotation matrix is described to GEANT by giving the polar and azimuthal angles of the axes of the DRS
  //  (x′, y′, z′) in the MRS
  fModuleAngle[0][0] = 90.;
  fModuleAngle[0][1] = 0.;
  fModuleAngle[1][0] = 90. - fModTheta;
  fModuleAngle[1][1] = 90.;
  fModuleAngle[2][0] = 180. - fModTheta;
  fModuleAngle[2][1] = -90.;

  fModuleCenter[0] = 0;
  fModuleCenter[1] = (GetIPtoOuterCoverDistance() + fEMCParams[3]) * std::sin(fModTheta * kDEGRad);
  fModuleCenter[2] = (GetIPtoOuterCoverDistance() + fEMCParams[3]) * std::cos(fModTheta * kDEGRad);
}