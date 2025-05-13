//_________________________________________________________________________
// Implementation of PHOS calorimeter

// --- Standard library ---

// --- ROOT system ---
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TParticle.h"
#include <TVirtualMC.h>
#include <TGeoElement.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>

// --- PHOS256 header files ---
#include "Phos.h"
#include "Geometry.h"
#include "SimParams.h"
#include "Stack.h"

//____________________________________________________________________________
Phos::Phos(float r, float theta) : fCurrentTrackID(-1),
                                   fCurrentCellID(-1),
                                   fCurentSuperParent(-1),
                                   fCurrentHit(nullptr)
{
  // Will create calorimeter at distance r to crystal surface and at polar ange theta
  // with normal pointing to interaction point
  Geometry::Instance(r, theta);
}

//____________________________________________________________________________
// void Phos::AddHit(int shunt, int primary, int Id, float * hits)
// {
// Add a hit to the hit list.
// A PHOS hit is the sum of all hits in a single crystal from one primary and within some time gate

// int hitCounter ;
// AliPHOSHit *newHit ;
// AliPHOSHit *curHit ;
// bool deja = kFALSE ;
// const Geometry * geom = Geometry::Instance() ;

// newHit = new AliPHOSHit(shunt, primary, Id, hits) ;

// for ( hitCounter = fNhits-1 ; hitCounter >= 0 && !deja ; hitCounter-- ) {
//   curHit = static_cast<AliPHOSHit*>((*fHits)[hitCounter]) ;
//   if(curHit->GetPrimary() != primary) break ;
//          // We add hits with the same primary, while GEANT treats primaries succesively
//   if( *curHit == *newHit ) {
//     *curHit + *newHit ;
//     deja = kTRUE ;
//   }
// }

// if ( !deja ) {
//   new((*fHits)[fNhits]) AliPHOSHit(*newHit) ;
//   // get the block Id number
//   int relid[4] ;
//   geom->AbsToRelNumbering(Id, relid) ;

//   fNhits++ ;
// }

// delete newHit;
// }

//____________________________________________________________________________
void Phos::FinishPrimary()
{
  // called at the end of each track (primary) by AliRun
  // hits are reset for each new track
  // accumulate the total hit-multiplicity
}

//____________________________________________________________________________
void Phos::FinishEvent()
{
  // Sort Hits
  // Add duplicates if any and remove them
  if (!fHits || fHits->size() == 0) {
    return;
  }

  auto first = fHits->begin();
  auto last = fHits->end();

  std::sort(first, last);

  first = fHits->begin();
  last = fHits->end();

  // this is copy of std::unique() method with addition: adding identical Hits
  auto itr = first;
  while (++first != last) {
    if (*itr == *first) {
      *itr += *first;
    } else {
      *(++itr) = *first;
    }
  }
  ++itr;

  fHits->erase(itr, fHits->end());

  // Apply Poisson smearing of light production
  first = fHits->begin();
  last = fHits->end();
  while (first != last) {
    float light = gRandom->Poisson(first->GetE() * SimParams::Instance()->fLightYieldPerGeV);
    first->SetEnergy(light / SimParams::Instance()->fLightYieldPerGeV);
    first++;
  }
}
void Phos::Reset()
{
  fSuperParents.clear();
  if (fHits) {
    fHits->clear();
  }
  fCurrentTrackID = -1;
  fCurrentCellID = -1;
  fCurentSuperParent = -1;
  fCurrentHit = nullptr;
}
bool Phos::ProcessHits()
{

  // 1. Remember all particles first entered PHOS (active medium)
  // 2. Collect all energy depositions in Cell by all secondaries from particle first entered PHOS

  // Check if this is first entered PHOS particle ("SuperParent")
  Stack* stack = static_cast<Stack*>(TVirtualMC::GetMC()->GetStack());
  const int partID = stack->GetCurrentTrackNumber();
  int superParent = -1;
  bool isNewPartile = false;       // Create Hit even if zero energy deposition
  if (partID != fCurrentTrackID) { // not same track as before, check: same SuperParent or new one?
    auto itTr = fSuperParents.find(partID);
    if (itTr == fSuperParents.end()) {
      // Search parent
      int parentID = stack->GetCurrentTrack()->GetMother(0);
      itTr = fSuperParents.find(parentID);
      if (itTr == fSuperParents.end()) { // Neither track or its parent found: new SuperParent
        fSuperParents[partID] = partID;
        superParent = partID;
        isNewPartile = true;
      } else { // parent found, this track - not
        superParent = itTr->second;
        fSuperParents[partID] = superParent;
        fCurrentTrackID = partID;
      }
    } else {
      superParent = itTr->second;
      fCurrentTrackID = partID;
    }
  } else {
    superParent = fCurentSuperParent;
  }

  if (isNewPartile) { // mark track to be kept by stack
    // stack->addHit(1); //TODO!!!! Mark track to be stored
  }

  double lostenergy = TVirtualMC::GetMC()->Edep();
  if (lostenergy < DBL_EPSILON && !isNewPartile) {
    return false; // do not create hits with zero energy deposition
  }

  int moduleNumber;
  TVirtualMC::GetMC()->CurrentVolOffID(
    11, moduleNumber); // 11: number of geom. levels between PXTL and PHOS module: get the PHOS module number ;
  int strip;
  TVirtualMC::GetMC()->CurrentVolOffID(3, strip); // 3: number of geom levels between PXTL and strip: get strip number in PHOS module
  int cell;
  TVirtualMC::GetMC()->CurrentVolOffID(2, cell); // 2: number of geom levels between PXTL and cell: get sell in strip number.
  int detID = Geometry::RelToAbsId(moduleNumber, strip, cell);
  if (superParent == fCurentSuperParent && detID == fCurrentCellID && fCurrentHit) {
    // continue with current hit
    fCurrentHit->AddEnergy(lostenergy);
    return true;
  }

  // try to find existing Hit
  if (!isNewPartile) {
    for (int itr = fHits->size() - 1; itr >= 0; itr--) {
      Hit* h = &(fHits->at(itr));
      if (h->GetLabel() != superParent) { // switched to another SuperParent, do not search further
        break;
      }
      if (h->GetCellID() == detID) { // found correct hit
        h->AddEnergy(lostenergy);
        fCurentSuperParent = superParent;
        fCurrentTrackID = partID;
        fCurrentCellID = detID;
        fCurrentHit = h;
        return true;
      }
    }
  }
  // Create new Hit
  float posX = 0., posY = 0., posZ = 0., momX = 0, momY = 0., momZ = 0., energy = 0.;
  TVirtualMC::GetMC()->TrackPosition(posX, posY, posZ);
  TVirtualMC::GetMC()->TrackMomentum(momX, momY, momZ, energy);
  double estart = TVirtualMC::GetMC()->Etot();
  double time = TVirtualMC::GetMC()->TrackTime(); // time in s

  fCurrentHit = &(fHits->emplace_back(detID, lostenergy, time, superParent));

  // fCurrentHit = fHits->back();
  fCurentSuperParent = superParent;
  fCurrentTrackID = partID;
  fCurrentCellID = detID;

  return true;
}

//____________________________________________________________________________
void Phos::CreateGeometryforEMC()
{
  // Create the PHOS-EMC geometry for GEANT
  // Author: Dmitri Peressounko August 2001
  // The used coordinate system:
  //   1. in Module: X along longer side, Y out of beam, Z along shorter side (along beam)
  //   2. In Strip the same: X along longer side, Y out of beam, Z along shorter side (along beam)

  // Get pointer to the array containing media indexes
  const Geometry* geom = Geometry::Instance();
  float par[4];
  int ipar;
  double* ubuf = 0;

  // ======= Define the strip ===============
  // still solid box
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetStripHalfSize() + ipar);
  int imed = gGeoManager->GetMedium("STEEL")->GetId();
  gGeoManager->Volume("PSTR", "BOX ", imed, par, 3); // Made of steel

  // --- define air hole in steel volume (cell of the strip unit)
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetAirCellHalfSize() + ipar);
  imed = gGeoManager->GetMedium("AIR")->GetId();
  gGeoManager->Volume("PCEL", "BOX ", imed, par, 3);

  // --- define wrapped crystal and put it into steel cell

  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetWrappedHalfSize() + ipar);
  imed = gGeoManager->GetMedium("TYVEK")->GetId();
  gGeoManager->Volume("PWRA", "BOX ", imed, par, 3);
  const float* pin = geom->GetAPDHalfSize();
  const float* preamp = geom->GetPreampHalfSize();
  float y = (geom->GetAirGapLed() - 2 * pin[1] - 2 * preamp[1]) / 2;
  gGeoManager->Node("PWRA", 1, "PCEL", 0.0, y, 0.0, 0, true, ubuf);

  // --- Define crystal and put it into wrapped crystall ---
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetCrystalHalfSize() + ipar);
  imed = gGeoManager->GetMedium("PWO")->GetId();
  gGeoManager->Volume("PXTL", "BOX ", imed, par, 3);
  gGeoManager->Node("PXTL", 1, "PWRA", 0.0, 0.0, 0.0, 0, true, ubuf);
  // void  Node (const char *name, Int_t nr, const char *mother, Double_t x, Double_t y, Double_t z, Int_t irot, Bool_t isOnly, Double_t *upar, Int_t npar=0)

  // --- define APD/PIN preamp and put it into AirCell

  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetAPDHalfSize() + ipar);
  imed = gGeoManager->GetMedium("SI")->GetId();
  gGeoManager->Volume("PPIN", "BOX ", imed, par, 3);
  const float* crystal = geom->GetCrystalHalfSize();
  y = crystal[1] + geom->GetAirGapLed() / 2 - preamp[1];
  gGeoManager->Node("PPIN", 1, "PCEL", 0.0, y, 0.0, 0, true, ubuf);
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetPreampHalfSize() + ipar);
  imed = gGeoManager->GetMedium("PCB")->GetId();
  gGeoManager->Volume("PREA", "BOX ", imed, par, 3); // Here I assumed preamp as a printed Circuit
  y = crystal[1] + geom->GetAirGapLed() / 2 + pin[1];
  gGeoManager->Node("PREA", 1, "PCEL", 0.0, y, 0.0, 0, true, ubuf);

  // --- Fill strip with wrapped cristals in steel cells
  const float* splate = geom->GetSupportPlateHalfSize();
  y = -splate[1];
  const float* acel = geom->GetAirCellHalfSize();

  for (int lev = 2, icel = 1;
       icel <= geom->GetNCellsXInStrip() * geom->GetNCellsZInStrip();
       icel += 2, lev += 2) {
    float x = (2 * (lev / 2) - 1 - geom->GetNCellsXInStrip()) * acel[0];
    float z = acel[2];
    gGeoManager->Node("PCEL", icel, "PSTR", x, y, +z, 0, true, ubuf);
    gGeoManager->Node("PCEL", icel + 1, "PSTR", x, y, -z, 0, true, ubuf);
  }

  // --- define the support plate, hole in it and position it in strip ----
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetSupportPlateHalfSize() + ipar);
  imed = gGeoManager->GetMedium("AL")->GetId();
  gGeoManager->Volume("PSUP", "BOX ", imed, par, 3);

  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetSupportPlateInHalfSize() + ipar);
  imed = gGeoManager->GetMedium("AIR")->GetId();
  gGeoManager->Volume("PSHO", "BOX ", imed, par, 3);
  float z = geom->GetSupportPlateThickness() / 2;
  gGeoManager->Node("PSHO", 1, "PSUP", 0.0, 0.0, z, 0, true, ubuf);

  y = acel[1];
  gGeoManager->Node("PSUP", 1, "PSTR", 0.0, y, 0.0, 0, true, ubuf);

  // ========== Fill module with strips and put them into inner thermoinsulation=============
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetInnerThermoHalfSize() + ipar);
  imed = gGeoManager->GetMedium("FOAM")->GetId();
  gGeoManager->Volume("PTII", "BOX ", imed, par, 3);

  const float* inthermo = geom->GetInnerThermoHalfSize();
  const float* strip = geom->GetStripHalfSize();
  y = inthermo[1] - strip[1];
  int irow;
  int nr = 1;
  int icol;

  for (irow = 0; irow < geom->GetNStripX(); irow++) {
    float x = (2 * irow + 1 - geom->GetNStripX()) * strip[0];
    for (icol = 0; icol < geom->GetNStripZ(); icol++) {
      z = (2 * icol + 1 - geom->GetNStripZ()) * strip[2];
      gGeoManager->Node("PSTR", nr, "PTII", x, y, z, 0, true, ubuf);
      nr++;
    }
  }

  // ------- define the air gap between thermoinsulation and cooler
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetAirGapHalfSize() + ipar);
  imed = gGeoManager->GetMedium("AIR")->GetId();
  gGeoManager->Volume("PAGA", "BOX ", imed, par, 3);
  const float* agap = geom->GetAirGapHalfSize();
  y = agap[1] - inthermo[1];

  gGeoManager->Node("PTII", 1, "PAGA", 0.0, y, 0.0, 0, true, ubuf);

  // ------- define the Al passive cooler
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetCoolerHalfSize() + ipar);
  imed = gGeoManager->GetMedium("AL")->GetId();
  gGeoManager->Volume("PCOR", "BOX ", imed, par, 3);

  const float* cooler = geom->GetCoolerHalfSize();
  y = cooler[1] - agap[1];

  gGeoManager->Node("PAGA", 1, "PCOR", 0.0, y, 0.0, 0, true, ubuf);

  // ------- define the outer thermoinsulating cover
  for (ipar = 0; ipar < 4; ipar++)
    par[ipar] = *(geom->GetOuterThermoParams() + ipar);
  imed = gGeoManager->GetMedium("FOAM")->GetId();
  gGeoManager->Volume("PTIO", "TRD1", imed, par, 4);
  const float* outparams = geom->GetOuterThermoParams();

  int idrotm = 0;
  TVirtualMC::GetMC()->Matrix(idrotm, 90.0, 0.0, 0.0, 0.0, 90.0, 270.0);
  // Frame in outer thermoinsulation and so on: z out of beam, y along beam, x across beam

  z = outparams[3] - cooler[1];
  gGeoManager->Node("PCOR", 1, "PTIO", 0., 0.0, z, idrotm, true, ubuf);

  // -------- Define the outer Aluminium cover -----
  for (ipar = 0; ipar < 4; ipar++)
    par[ipar] = *(geom->GetAlCoverParams() + ipar);
  imed = gGeoManager->GetMedium("AL")->GetId();
  gGeoManager->Volume("PCOL", "TRD1", imed, par, 4);

  const float* covparams = geom->GetAlCoverParams();
  z = covparams[3] - outparams[3];
  gGeoManager->Node("PTIO", 1, "PCOL", 0., 0.0, z, 0, true, ubuf);

  // --------- Define front fiberglass cover -----------
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetFiberGlassHalfSize() + ipar);
  imed = gGeoManager->GetMedium("FIBERGLASS")->GetId();
  gGeoManager->Volume("PFGC", "BOX ", imed, par, 3);
  z = -outparams[3];
  gGeoManager->Node("PFGC", 1, "PCOL", 0., 0.0, z, 0, true, ubuf);

  //=============This is all with cold section==============

  //------ Warm Section --------------
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetWarmAlCoverHalfSize() + ipar);
  imed = gGeoManager->GetMedium("AL")->GetId();
  gGeoManager->Volume("PWAR", "BOX ", imed, par, 3);
  const float* warmcov = geom->GetWarmAlCoverHalfSize();

  // --- Define the outer thermoinsulation ---
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetWarmThermoHalfSize() + ipar);
  imed = gGeoManager->GetMedium("FOAM")->GetId();
  gGeoManager->Volume("PWTI", "BOX ", imed, par, 3);
  const float* warmthermo = geom->GetWarmThermoHalfSize();
  z = -warmcov[2] + warmthermo[2];

  gGeoManager->Node("PWTI", 1, "PWAR", 0., 0.0, z, 0, true, ubuf);

  // --- Define cables area and put in it T-supports ----
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetTCables1HalfSize() + ipar);
  imed = gGeoManager->GetMedium("CABLES")->GetId();
  gGeoManager->Volume("PCA1", "BOX ", imed, par, 3);
  const float* cbox = geom->GetTCables1HalfSize();

  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetTSupport1HalfSize() + ipar);
  imed = gGeoManager->GetMedium("AL")->GetId();
  gGeoManager->Volume("PBE1", "BOX ", imed, par, 3);
  const float* beams = geom->GetTSupport1HalfSize();
  int isup;
  for (isup = 0; isup < geom->GetNTSuppots(); isup++) {
    float x = -cbox[0] + beams[0] + (2 * beams[0] + geom->GetTSupportDist()) * isup;
    gGeoManager->Node("PBE1", isup, "PCA1", x, 0.0, 0.0, 0, true, ubuf);
  }

  z = -warmthermo[2] + cbox[2];
  gGeoManager->Node("PCA1", 1, "PWTI", 0.0, 0.0, z, 0, true, ubuf);

  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetTCables2HalfSize() + ipar);
  imed = gGeoManager->GetMedium("CABLES")->GetId();
  gGeoManager->Volume("PCA2", "BOX ", imed, par, 3);
  const float* cbox2 = geom->GetTCables2HalfSize();

  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetTSupport2HalfSize() + ipar);
  imed = gGeoManager->GetMedium("AL")->GetId();
  gGeoManager->Volume("PBE2", "BOX ", imed, par, 3);
  for (isup = 0; isup < geom->GetNTSuppots(); isup++) {
    float x = -cbox[0] + beams[0] + (2 * beams[0] + geom->GetTSupportDist()) * isup;
    gGeoManager->Node("PBE2", isup, "PCA2", x, 0.0, 0.0, 0, true, ubuf);
  }

  z = -warmthermo[2] + 2 * cbox[2] + cbox2[2];
  gGeoManager->Node("PCA2", 1, "PWTI", 0.0, 0.0, z, 0, true, ubuf);

  // --- Define frame ---
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetFrameXHalfSize() + ipar);
  imed = gGeoManager->GetMedium("STEEL")->GetId();
  gGeoManager->Volume("PFRX", "BOX ", imed, par, 3);
  const float* posit1 = geom->GetFrameXPosition();
  gGeoManager->Node("PFRX", 1, "PWTI", posit1[0], posit1[1], posit1[2], 0, true, ubuf);
  gGeoManager->Node("PFRX", 2, "PWTI", posit1[0], -posit1[1], posit1[2], 0, true, ubuf);

  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetFrameZHalfSize() + ipar);
  imed = gGeoManager->GetMedium("STEEL")->GetId();
  gGeoManager->Volume("PFRZ", "BOX ", imed, par, 3);
  const float* posit2 = geom->GetFrameZPosition();
  gGeoManager->Node("PFRZ", 1, "PWTI", posit2[0], posit2[1], posit2[2], 0, true, ubuf);
  gGeoManager->Node("PFRZ", 2, "PWTI", -posit2[0], posit2[1], posit2[2], 0, true, ubuf);

  // --- Define Fiber Glass support ---
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetFGupXHalfSize() + ipar);
  imed = gGeoManager->GetMedium("FIBERGLASS")->GetId();
  gGeoManager->Volume("PFG1", "BOX ", imed, par, 3);
  const float* posit3 = geom->GetFGupXPosition();
  gGeoManager->Node("PFG1", 1, "PWTI", posit3[0], posit3[1], posit3[2], 0, true, ubuf);
  gGeoManager->Node("PFG1", 2, "PWTI", posit3[0], -posit3[1], posit3[2], 0, true, ubuf);

  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetFGupZHalfSize() + ipar);
  imed = gGeoManager->GetMedium("FIBERGLASS")->GetId();
  gGeoManager->Volume("PFG2", "BOX ", imed, par, 3);
  const float* posit4 = geom->GetFGupZPosition();
  gGeoManager->Node("PFG2", 1, "PWTI", posit4[0], posit4[1], posit4[2], 0, true, ubuf);
  gGeoManager->Node("PFG2", 2, "PWTI", -posit4[0], posit4[1], posit4[2], 0, true, ubuf);
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetFGlowXHalfSize() + ipar);
  imed = gGeoManager->GetMedium("FIBERGLASS")->GetId();
  gGeoManager->Volume("PFG3", "BOX ", imed, par, 3);
  const float* posit5 = geom->GetFGlowXPosition();
  gGeoManager->Node("PFG3", 1, "PWTI", posit5[0], posit5[1], posit5[2], 0, true, ubuf);
  gGeoManager->Node("PFG3", 2, "PWTI", posit5[0], -posit5[1], posit5[2], 0, true, ubuf);

  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetFGlowZHalfSize() + ipar);
  imed = gGeoManager->GetMedium("FIBERGLASS")->GetId();
  gGeoManager->Volume("PFG4", "BOX ", imed, par, 3);
  const float* posit6 = geom->GetFGlowZPosition();
  gGeoManager->Node("PFG4", 1, "PWTI", posit6[0], posit6[1], posit6[2], 0, true, ubuf);
  gGeoManager->Node("PFG4", 2, "PWTI", -posit6[0], posit6[1], posit6[2], 0, true, ubuf);

  // --- Define Air Gap for FEE electronics -----
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetFEEAirHalfSize() + ipar);
  imed = gGeoManager->GetMedium("AIR")->GetId();
  gGeoManager->Volume("PAFE", "BOX ", imed, par, 3);
  const float* posit7 = geom->GetFEEAirPosition();
  gGeoManager->Node("PAFE", 1, "PWTI", posit7[0], posit7[1], posit7[2], 0, true, ubuf);

  // Define the EMC module volume and combine Cool and Warm sections
  for (ipar = 0; ipar < 4; ipar++)
    par[ipar] = *(geom->GetPHOSParams() + ipar);
  imed = gGeoManager->GetMedium("AIR")->GetId();
  gGeoManager->Volume("PEMC", "TRD1", imed, par, 4);
  z = -warmcov[2];
  gGeoManager->Node("PCOL", 1, "PEMC", 0., 0., z, 0, true, ubuf);
  z = covparams[3];
  gGeoManager->Node("PWAR", 1, "PEMC", 0., 0., z, 0, true, ubuf);

  // Put created EMC geometry into PHOS volume
  gGeoManager->Node("PEMC", 1, "PHOS", 0., 0., 0., 0, true, ubuf);
}

//____________________________________________________________________________
void Phos::CreateMaterials()
{
  // Definitions of materials to build PHOS and associated tracking media.
  // media number in idtmed are 0 to 20.

  // Paremeter for tracking media
  double param[20] = {0};
  param[0] = 0;     // isvol  - Not used
  param[1] = 2;     // ifield - User defined magnetic field
  param[2] = 10.;   // fieldm - Maximum field value (in kiloGauss)
  param[3] = -20.;  // tmaxfd - Maximum angle due to field deflection
  param[4] = -0.01; // stemax - Maximum displacement for multiple scat
  param[5] = -.3;   // deemax - Maximum fractional energy loss, DLS
  param[6] = .001;  // epsil - Tracking precision
  param[7] = -.8;   // stmin

  double a, z, density;
  // --- The PbWO4 crystals ---
  TGeoElement* elPb = new TGeoElement("Lead", "Pb", z = 82, a = 207.19);
  TGeoElement* elW = new TGeoElement("Tangstate", "W", z = 74., a = 183.85);
  TGeoElement* elO = new TGeoElement("Oxygen", "O", z = 8., a = 16.00);
  TGeoElement* elC = new TGeoElement("Carbon", "C", z = 6., a = 12.011);
  TGeoElement* elN = new TGeoElement("Nitrogen", "N", z = 7., a = 14.01);
  TGeoElement* elH = new TGeoElement("Hydrogen", "H", z = 1., a = 1.00794);
  TGeoElement* elSi = new TGeoElement("Silicon", "Si", z = 14., a = 28.09);
  TGeoElement* elAl = new TGeoElement("Aluminium", "Al", z = 13., a = 26.98);
  TGeoElement* elFe = new TGeoElement("Ferrum", "Fe", z = 26., a = 55.845);
  TGeoElement* elCu = new TGeoElement("Cuprum", "Cu", z = 29., a = 63.546);

  // --- Air ---
  TGeoMixture* mixAir = new TGeoMixture("Air", 2, density = 1.29e-03);
  mixAir->AddElement(elN, 0.7);
  mixAir->AddElement(elO, 0.3);
  new TGeoMedium("AIR", kAir, mixAir, param);

  TGeoMixture* mixPWO = new TGeoMixture("PbWO4", 3, density = 8.28);
  mixPWO->AddElement(elPb, 1);
  mixPWO->AddElement(elW, 1);
  mixPWO->AddElement(elO, 4);
  new TGeoMedium("PWO", kPWO, mixPWO, param);

  // // --- Aluminium ---
  TGeoMixture* mixAl = new TGeoMixture("Al", 1, density = 2.7);
  new TGeoMedium("AL", kAl, mixAl, param);

  TGeoMixture* mixTyvek = new TGeoMixture("Tyvek", kTyvek, density = 0.331);
  mixTyvek->AddElement(elH, 1);
  mixTyvek->AddElement(elC, 2);
  new TGeoMedium("TYVEK", kTyvek, mixTyvek, param);

  TGeoMixture* mixFoam = new TGeoMixture("Foam", 2, density = 0.04);
  mixFoam->AddElement(elH, 1);
  mixFoam->AddElement(elC, 1);
  new TGeoMedium("FOAM", kFoam, mixFoam, param);

  // --- Silicon ---
  TGeoMixture* mixSi = new TGeoMixture("Si", 1, density = 2.33);
  new TGeoMedium("SI", kSi, mixSi, param);

  // --- PCB : Printed Circuit material ---
  TGeoMixture* mixPCB = new TGeoMixture("PCB", 4, density = 1.7);
  mixPCB->AddElement(elH, 1);
  mixPCB->AddElement(elC, 1);
  mixPCB->AddElement(elO, 1);
  mixPCB->AddElement(elN, 1);
  new TGeoMedium("PCB", kPCB, mixPCB, param);

  // --- Stainless steel (let it be pure iron) ---
  TGeoMixture* mixSteel = new TGeoMixture("STEEL", 1, density = 7.87);
  mixPCB->AddElement(elFe, 1);
  new TGeoMedium("STEEL", kSteel, mixSteel, param);

  // --- Fiberglass ---
  TGeoMixture* mixFiberglas = new TGeoMixture("Fiberglas", 4, density = 1.9);
  mixFiberglas->AddElement(elO, 4);
  mixFiberglas->AddElement(elN, 1);
  mixFiberglas->AddElement(elC, 7);
  mixFiberglas->AddElement(elH, 10);
  new TGeoMedium("FIBERGLASS", kFiberglas, mixFiberglas, param);

  // --- Cables in Air box  ---
  TGeoMixture* mixCables = new TGeoMixture("Cables", 4, density = 0.8);
  mixFiberglas->AddElement(elH, 1);
  mixFiberglas->AddElement(elC, 7);
  mixFiberglas->AddElement(elFe, 30);
  mixFiberglas->AddElement(elCu, 30);
  new TGeoMedium("CABLES", kCables, mixCables, param);
}

//____________________________________________________________________________
void Phos::CreateGeometry()
{
  // Create the PHOS geometry for Geant

  const Geometry* geom = Geometry::Instance();
  std::cout << "Phos::CreateGeometry, geom=" << geom << std::endl;

  // Create a PHOS module.
  // Gsvolu accepts non-const params. Avoid modification geometry
  float params[4] = {geom->GetPHOSParams()[0], geom->GetPHOSParams()[1], geom->GetPHOSParams()[2], geom->GetPHOSParams()[3]};
  int imed = gGeoManager->GetMedium("AIR")->GetId();
  gGeoManager->Volume("PHOS", "TRD1", imed, params, 4); // 20-> idtmed[20]

  CreateGeometryforEMC();

  //  CreateGeometryforSupport();

  // --- Position  PHOS mdules in Hall ---
  int iXYZ, iAngle;
  float angle[3][2];
  geom->GetModuleAngles(angle);
  int idrotm = 0;
  TVirtualMC::GetMC()->Matrix(idrotm, angle[0][0], angle[0][1],
                              angle[1][0], angle[1][1],
                              angle[2][0], angle[2][1]);

  float pos[3] = {0};
  // geom->GetModuleCenter(pos);
  gGeoManager->Node("PHOS", 0, "World", pos[0], pos[1], pos[2], idrotm, true, static_cast<double*>(nullptr));
}
