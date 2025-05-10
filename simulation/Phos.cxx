//_________________________________________________________________________
// Implementation of PHOS calorimeter

// --- Standard library ---

// --- ROOT system ---
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TParticle.h"
#include <TVirtualMC.h>

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

  // ======= Define the strip ===============
  // still solid box
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetStripHalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PSTR", "BOX ", fIdtmed[17], par, 3); // Made of steel

  // --- define air hole in steel volume (cell of the strip unit)
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetAirCellHalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PCEL", "BOX ", fIdtmed[20], par, 3);

  // --- define wrapped crystal and put it into steel cell

  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetWrappedHalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PWRA", "BOX ", fIdtmed[3], par, 3);
  const float* pin = geom->GetAPDHalfSize();
  const float* preamp = geom->GetPreampHalfSize();
  float y = (geom->GetAirGapLed() - 2 * pin[1] - 2 * preamp[1]) / 2;
  TVirtualMC::GetMC()->Gspos("PWRA", 1, "PCEL", 0.0, y, 0.0, 0, "ONLY");

  // --- Define crystal and put it into wrapped crystall ---
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetCrystalHalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PXTL", "BOX ", fIdtmed[0], par, 3);
  TVirtualMC::GetMC()->Gspos("PXTL", 1, "PWRA", 0.0, 0.0, 0.0, 0, "ONLY");

  // --- define APD/PIN preamp and put it into AirCell

  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetAPDHalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PPIN", "BOX ", fIdtmed[6], par, 3);
  const float* crystal = geom->GetCrystalHalfSize();
  y = crystal[1] + geom->GetAirGapLed() / 2 - preamp[1];
  TVirtualMC::GetMC()->Gspos("PPIN", 1, "PCEL", 0.0, y, 0.0, 0, "ONLY");
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetPreampHalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PREA", "BOX ", fIdtmed[12], par, 3);      // Here I assumed preamp as a printed Circuit
  y = crystal[1] + geom->GetAirGapLed() / 2 + pin[1];                    // May it should be changed
  TVirtualMC::GetMC()->Gspos("PREA", 1, "PCEL", 0.0, y, 0.0, 0, "ONLY"); // to ceramics?

  // --- Fill strip with wrapped cristals in steel cells
  const float* splate = geom->GetSupportPlateHalfSize();
  y = -splate[1];
  const float* acel = geom->GetAirCellHalfSize();

  for (int lev = 2, icel = 1;
       icel <= geom->GetNCellsXInStrip() * geom->GetNCellsZInStrip();
       icel += 2, lev += 2) {
    float x = (2 * (lev / 2) - 1 - geom->GetNCellsXInStrip()) * acel[0];
    float z = acel[2];
    TVirtualMC::GetMC()->Gspos("PCEL", icel, "PSTR", x, y, +z, 0, "ONLY");
    TVirtualMC::GetMC()->Gspos("PCEL", icel + 1, "PSTR", x, y, -z, 0, "ONLY");
  }

  // --- define the support plate, hole in it and position it in strip ----
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetSupportPlateHalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PSUP", "BOX ", fIdtmed[2], par, 3);

  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetSupportPlateInHalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PSHO", "BOX ", fIdtmed[20], par, 3);
  float z = geom->GetSupportPlateThickness() / 2;
  TVirtualMC::GetMC()->Gspos("PSHO", 1, "PSUP", 0.0, 0.0, z, 0, "ONLY");

  y = acel[1];
  TVirtualMC::GetMC()->Gspos("PSUP", 1, "PSTR", 0.0, y, 0.0, 0, "ONLY");

  // ========== Fill module with strips and put them into inner thermoinsulation=============
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetInnerThermoHalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PTII", "BOX ", fIdtmed[7], par, 3);

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
      TVirtualMC::GetMC()->Gspos("PSTR", nr, "PTII", x, y, z, 0, "ONLY");
      nr++;
    }
  }

  // ------- define the air gap between thermoinsulation and cooler
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetAirGapHalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PAGA", "BOX ", fIdtmed[20], par, 3);
  const float* agap = geom->GetAirGapHalfSize();
  y = agap[1] - inthermo[1];

  TVirtualMC::GetMC()->Gspos("PTII", 1, "PAGA", 0.0, y, 0.0, 0, "ONLY");

  // ------- define the Al passive cooler
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetCoolerHalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PCOR", "BOX ", fIdtmed[2], par, 3);

  const float* cooler = geom->GetCoolerHalfSize();
  y = cooler[1] - agap[1];

  TVirtualMC::GetMC()->Gspos("PAGA", 1, "PCOR", 0.0, y, 0.0, 0, "ONLY");

  // ------- define the outer thermoinsulating cover
  for (ipar = 0; ipar < 4; ipar++)
    par[ipar] = *(geom->GetOuterThermoParams() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PTIO", "TRD1", fIdtmed[7], par, 4);
  const float* outparams = geom->GetOuterThermoParams();

  int idrotm = 0;
  TVirtualMC::GetMC()->Matrix(idrotm, 90.0, 0.0, 0.0, 0.0, 90.0, 270.0);
  // Frame in outer thermoinsulation and so on: z out of beam, y along beam, x across beam

  z = outparams[3] - cooler[1];
  TVirtualMC::GetMC()->Gspos("PCOR", 1, "PTIO", 0., 0.0, z, idrotm, "ONLY");

  // -------- Define the outer Aluminium cover -----
  for (ipar = 0; ipar < 4; ipar++)
    par[ipar] = *(geom->GetAlCoverParams() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PCOL", "TRD1", fIdtmed[2], par, 4);

  const float* covparams = geom->GetAlCoverParams();
  z = covparams[3] - outparams[3];
  TVirtualMC::GetMC()->Gspos("PTIO", 1, "PCOL", 0., 0.0, z, 0, "ONLY");

  // --------- Define front fiberglass cover -----------
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetFiberGlassHalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PFGC", "BOX ", fIdtmed[18], par, 3);
  z = -outparams[3];
  TVirtualMC::GetMC()->Gspos("PFGC", 1, "PCOL", 0., 0.0, z, 0, "ONLY");

  //=============This is all with cold section==============

  //------ Warm Section --------------
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetWarmAlCoverHalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PWAR", "BOX ", fIdtmed[2], par, 3);
  const float* warmcov = geom->GetWarmAlCoverHalfSize();

  // --- Define the outer thermoinsulation ---
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetWarmThermoHalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PWTI", "BOX ", fIdtmed[7], par, 3);
  const float* warmthermo = geom->GetWarmThermoHalfSize();
  z = -warmcov[2] + warmthermo[2];

  TVirtualMC::GetMC()->Gspos("PWTI", 1, "PWAR", 0., 0.0, z, 0, "ONLY");

  // --- Define cables area and put in it T-supports ----
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetTCables1HalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PCA1", "BOX ", fIdtmed[718], par, 3);
  const float* cbox = geom->GetTCables1HalfSize();

  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetTSupport1HalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PBE1", "BOX ", fIdtmed[2], par, 3);
  const float* beams = geom->GetTSupport1HalfSize();
  int isup;
  for (isup = 0; isup < geom->GetNTSuppots(); isup++) {
    float x = -cbox[0] + beams[0] + (2 * beams[0] + geom->GetTSupportDist()) * isup;
    TVirtualMC::GetMC()->Gspos("PBE1", isup, "PCA1", x, 0.0, 0.0, 0, "ONLY");
  }

  z = -warmthermo[2] + cbox[2];
  TVirtualMC::GetMC()->Gspos("PCA1", 1, "PWTI", 0.0, 0.0, z, 0, "ONLY");

  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetTCables2HalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PCA2", "BOX ", fIdtmed[718], par, 3);
  const float* cbox2 = geom->GetTCables2HalfSize();

  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetTSupport2HalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PBE2", "BOX ", fIdtmed[2], par, 3);
  for (isup = 0; isup < geom->GetNTSuppots(); isup++) {
    float x = -cbox[0] + beams[0] + (2 * beams[0] + geom->GetTSupportDist()) * isup;
    TVirtualMC::GetMC()->Gspos("PBE2", isup, "PCA2", x, 0.0, 0.0, 0, "ONLY");
  }

  z = -warmthermo[2] + 2 * cbox[2] + cbox2[2];
  TVirtualMC::GetMC()->Gspos("PCA2", 1, "PWTI", 0.0, 0.0, z, 0, "ONLY");

  // --- Define frame ---
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetFrameXHalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PFRX", "BOX ", fIdtmed[17], par, 3);
  const float* posit1 = geom->GetFrameXPosition();
  TVirtualMC::GetMC()->Gspos("PFRX", 1, "PWTI", posit1[0], posit1[1], posit1[2], 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("PFRX", 2, "PWTI", posit1[0], -posit1[1], posit1[2], 0, "ONLY");

  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetFrameZHalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PFRZ", "BOX ", fIdtmed[17], par, 3);
  const float* posit2 = geom->GetFrameZPosition();
  TVirtualMC::GetMC()->Gspos("PFRZ", 1, "PWTI", posit2[0], posit2[1], posit2[2], 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("PFRZ", 2, "PWTI", -posit2[0], posit2[1], posit2[2], 0, "ONLY");

  // --- Define Fiber Glass support ---
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetFGupXHalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PFG1", "BOX ", fIdtmed[18], par, 3);
  const float* posit3 = geom->GetFGupXPosition();
  TVirtualMC::GetMC()->Gspos("PFG1", 1, "PWTI", posit3[0], posit3[1], posit3[2], 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("PFG1", 2, "PWTI", posit3[0], -posit3[1], posit3[2], 0, "ONLY");

  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetFGupZHalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PFG2", "BOX ", fIdtmed[18], par, 3);
  const float* posit4 = geom->GetFGupZPosition();
  TVirtualMC::GetMC()->Gspos("PFG2", 1, "PWTI", posit4[0], posit4[1], posit4[2], 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("PFG2", 2, "PWTI", -posit4[0], posit4[1], posit4[2], 0, "ONLY");
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetFGlowXHalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PFG3", "BOX ", fIdtmed[18], par, 3);
  const float* posit5 = geom->GetFGlowXPosition();
  TVirtualMC::GetMC()->Gspos("PFG3", 1, "PWTI", posit5[0], posit5[1], posit5[2], 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("PFG3", 2, "PWTI", posit5[0], -posit5[1], posit5[2], 0, "ONLY");

  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetFGlowZHalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PFG4", "BOX ", fIdtmed[18], par, 3);
  const float* posit6 = geom->GetFGlowZPosition();
  TVirtualMC::GetMC()->Gspos("PFG4", 1, "PWTI", posit6[0], posit6[1], posit6[2], 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("PFG4", 2, "PWTI", -posit6[0], posit6[1], posit6[2], 0, "ONLY");

  // --- Define Air Gap for FEE electronics -----
  for (ipar = 0; ipar < 3; ipar++)
    par[ipar] = *(geom->GetFEEAirHalfSize() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PAFE", "BOX ", fIdtmed[20], par, 3);
  const float* posit7 = geom->GetFEEAirPosition();
  TVirtualMC::GetMC()->Gspos("PAFE", 1, "PWTI", posit7[0], posit7[1], posit7[2], 0, "ONLY");

  // Define the EMC module volume and combine Cool and Warm sections
  for (ipar = 0; ipar < 4; ipar++)
    par[ipar] = *(geom->GetPHOSParams() + ipar);
  TVirtualMC::GetMC()->Gsvolu("PEMC", "TRD1", fIdtmed[20], par, 4);
  z = -warmcov[2];
  TVirtualMC::GetMC()->Gspos("PCOL", 1, "PEMC", 0., 0., z, 0, "ONLY");
  z = covparams[3];
  TVirtualMC::GetMC()->Gspos("PWAR", 1, "PEMC", 0., 0., z, 0, "ONLY");

  // Put created EMC geometry into PHOS volume
  TVirtualMC::GetMC()->Gspos("PEMC", 1, "PHOS", 0., 0., 0., 0, "ONLY");
}

//____________________________________________________________________________
void Phos::CreateMaterials()
{
  // Definitions of materials to build PHOS and associated tracking media.
  // media number in idtmed are 0 to 20.

  int kmat;
  // --- The PbWO4 crystals ---
  float aX[3] = {207.19, 183.85, 16.0};
  float zX[3] = {82.0, 74.0, 8.0};
  float wX[3] = {1.0, 1.0, 4.0};
  float dX = 8.28;

  TVirtualMC::GetMC()->Mixture(kmat, "PbWO4$", aX, zX, dX, -3, wX);
  fIdmate[0] = kmat;

  // --- The polysterene scintillator (CH) ---
  float aP[2] = {12.011, 1.00794};
  float zP[2] = {6.0, 1.0};
  float wP[2] = {1.0, 1.0};
  float dP = 1.032;

  TVirtualMC::GetMC()->Mixture(kmat, "Polystyrene$", aP, zP, dP, -2, wP);
  fIdmate[1] = kmat;

  // --- Aluminium ---
  TVirtualMC::GetMC()->Material(kmat, "Al$", 26.98, 13., 2.7, 8.9, 999., static_cast<float*>(nullptr), 0);
  // ---         Absorption length is ignored ^
  fIdmate[2] = kmat;

  // --- Tyvek (CnH2n) ---
  float aT[2] = {12.011, 1.00794};
  float zT[2] = {6.0, 1.0};
  float wT[2] = {1.0, 2.0};
  float dT = 0.331;

  TVirtualMC::GetMC()->Mixture(kmat, "Tyvek$", aT, zT, dT, -2, wT);
  fIdmate[3] = kmat;

  // --- Polystyrene foam ---
  float aF[2] = {12.011, 1.00794};
  float zF[2] = {6.0, 1.0};
  float wF[2] = {1.0, 1.0};
  float dF = 0.12;

  TVirtualMC::GetMC()->Mixture(kmat, "Foam$", aF, zF, dF, -2, wF);
  fIdmate[4] = kmat;

  // --- Titanium ---
  float aTIT[3] = {47.88, 26.98, 54.94};
  float zTIT[3] = {22.0, 13.0, 25.0};
  float wTIT[3] = {69.0, 6.0, 1.0};
  float dTIT = 4.5;

  TVirtualMC::GetMC()->Mixture(kmat, "Titanium$", aTIT, zTIT, dTIT, -3, wTIT);
  fIdmate[5] = kmat;

  // --- Silicon ---
  TVirtualMC::GetMC()->Material(kmat, "Si$", 28.0855, 14., 2.33, 9.36, 42.3, static_cast<float*>(nullptr), 0);
  fIdmate[6] = kmat;

  // --- Foam thermo insulation ---
  float aTI[2] = {12.011, 1.00794};
  float zTI[2] = {6.0, 1.0};
  float wTI[2] = {1.0, 1.0};
  float dTI = 0.04;

  TVirtualMC::GetMC()->Mixture(kmat, "Thermo Insul.$", aTI, zTI, dTI, -2, wTI);
  fIdmate[7] = kmat;

  // --- Textolith ---
  float aTX[4] = {16.0, 28.09, 12.011, 1.00794};
  float zTX[4] = {8.0, 14.0, 6.0, 1.0};
  float wTX[4] = {292.0, 68.0, 462.0, 736.0};
  float dTX = 1.75;

  TVirtualMC::GetMC()->Mixture(kmat, "Textolit$", aTX, zTX, dTX, -4, wTX);
  fIdmate[8] = kmat;

  //--- FR4  ---
  float aFR[4] = {16.0, 28.09, 12.011, 1.00794};
  float zFR[4] = {8.0, 14.0, 6.0, 1.0};
  float wFR[4] = {292.0, 68.0, 462.0, 736.0};
  float dFR = 1.8;

  TVirtualMC::GetMC()->Mixture(kmat, "FR4$", aFR, zFR, dFR, -4, wFR);
  fIdmate[9] = kmat;

  // --- The Composite Material for  micromegas (so far polyetylene) ---
  float aCM[2] = {12.01, 1.};
  float zCM[2] = {6., 1.};
  float wCM[2] = {1., 2.};
  float dCM = 0.935;

  TVirtualMC::GetMC()->Mixture(kmat, "Compo Mat$", aCM, zCM, dCM, -2, wCM);
  fIdmate[10] = kmat;

  // --- Copper ---
  TVirtualMC::GetMC()->Material(kmat, "Cu$", 63.546, 29, 8.96, 1.43, 14.8, static_cast<float*>(nullptr), 0);
  fIdmate[11] = kmat;

  // --- G10 : Printed Circuit material ---
  float aG10[4] = {12., 1., 16., 28.};
  float zG10[4] = {6., 1., 8., 14.};
  float wG10[4] = {.259, .288, .248, .205};
  float dG10 = 1.7;

  TVirtualMC::GetMC()->Mixture(kmat, "G10$", aG10, zG10, dG10, -4, wG10);
  fIdmate[12] = kmat;

  // --- Lead ---
  TVirtualMC::GetMC()->Material(kmat, "Pb$", 207.2, 82, 11.35, 0.56, 0., static_cast<float*>(nullptr), 0);
  fIdmate[13] = kmat;

  // --- The gas mixture ---
  // Co2
  float aCO[2] = {12.0, 16.0};
  float zCO[2] = {6.0, 8.0};
  float wCO[2] = {1.0, 2.0};
  float dCO = 0.001977;

  TVirtualMC::GetMC()->Mixture(kmat, "CO2$", aCO, zCO, dCO, -2, wCO);
  fIdmate[14] = kmat;

  // Ar
  float dAr = 0.001782;
  TVirtualMC::GetMC()->Material(kmat, "Ar$", 39.948, 18.0, dAr, 14.0, 0., static_cast<float*>(nullptr), 0);
  fIdmate[15] = kmat;

  // Ar+CO2 Mixture (80% / 20%)
  float arContent = 0.80; // Ar-content of the ArCO2-mixture
  float aArCO[3] = {39.948, 12.0, 16.0};
  float zArCO[3] = {18.0, 6.0, 8.0};
  float wArCO[3];
  wArCO[0] = arContent;
  wArCO[1] = (1 - arContent) * 1;
  wArCO[2] = (1 - arContent) * 2;
  float dArCO = arContent * dAr + (1 - arContent) * dCO;
  TVirtualMC::GetMC()->Mixture(kmat, "ArCO2$", aArCO, zArCO, dArCO, -3, wArCO);
  fIdmate[16] = kmat;

  // --- Stainless steel (let it be pure iron) ---
  TVirtualMC::GetMC()->Material(kmat, "Steel$", 55.845, 26, 7.87, 1.76, 0., static_cast<float*>(nullptr), 0);
  fIdmate[17] = kmat;

  // --- Fiberglass ---
  float aFG[4] = {16.0, 28.09, 12.011, 1.00794};
  float zFG[4] = {8.0, 14.0, 6.0, 1.0};
  float wFG[4] = {292.0, 68.0, 462.0, 736.0};
  float dFG = 1.9;

  TVirtualMC::GetMC()->Mixture(kmat, "Fibergla$", aFG, zFG, dFG, -4, wFG);
  fIdmate[18] = kmat;

  // --- Cables in Air box  ---
  // SERVICES

  float aCA[4] = {1., 12., 55.8, 63.5};
  float zCA[4] = {1., 6., 26., 29.};
  float wCA[4] = {.014, .086, .42, .48};
  float dCA = 0.8; // this density is raw estimation, if you know better - correct

  TVirtualMC::GetMC()->Mixture(kmat, "Cables  $", aCA, zCA, dCA, -4, wCA);
  fIdmate[19] = kmat;

  // --- Air ---
  float aAir[4] = {12.0107, 14.0067, 15.9994, 39.948};
  float zAir[4] = {6., 7., 8., 18.};
  float wAir[4] = {0.000124, 0.755267, 0.231781, 0.012827};
  float dAir = 1.20479E-3;

  TVirtualMC::GetMC()->Mixture(kmat, "Air$", aAir, zAir, dAir, 4, wAir);
  fIdmate[20] = kmat;

  // DEFINITION OF THE TRACKING MEDIA

  // for PHOS: idtmed[0->20] equivalent to fIdtmed[0->100]
  int isxfld = 0.;   //(TGeoGlobalMagField::Instance()->GetField())->Integ() ;
  float sxmgmx = 0.; //(TGeoGlobalMagField::Instance()->GetField())->Max() ;

  int kmed = 0;
  // The scintillator of the calorimeter made of PBW04                              -> idtmed[0]
  TVirtualMC::GetMC()->Medium(kmed, "PHOS Xtal    $", 0, 1,
                              isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, static_cast<float*>(nullptr), 0);
  fIdtmed[0] = kmed;

  // The scintillator of the CPV made of Polystyrene scintillator                   -> idtmed[700]
  TVirtualMC::GetMC()->Medium(kmed, "CPV scint.   $", 1, 1,
                              isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, static_cast<float*>(nullptr), 0);
  fIdtmed[1] = kmed;

  // Various Aluminium parts made of Al                                             -> idtmed[2]
  TVirtualMC::GetMC()->Medium(kmed, "Al parts     $", 2, 0,
                              isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.001, 0.001, static_cast<float*>(nullptr), 0);
  fIdtmed[2] = kmed;

  // The Tywek which wraps the calorimeter crystals                                 -> idtmed[3]
  TVirtualMC::GetMC()->Medium(kmed, "Tyvek wrapper$", 3, 0,
                              isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.001, 0.001, static_cast<float*>(nullptr), 0);
  fIdtmed[3] = kmed;

  // The Polystyrene foam around the calorimeter module                             -> idtmed[703]
  TVirtualMC::GetMC()->Medium(kmed, "Polyst. foam $", 4, 0,
                              isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, static_cast<float*>(nullptr), 0);
  fIdtmed[4] = kmed;

  // The Titanium around the calorimeter crystal                                    -> idtmed[704]
  TVirtualMC::GetMC()->Medium(kmed, "Titan. cover $", 5, 0,
                              isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.0001, 0.0001, static_cast<float*>(nullptr), 0);
  fIdtmed[5] = kmed;

  // The Silicon of the pin diode to read out the calorimeter crystal               -> idtmed[6]
  TVirtualMC::GetMC()->Medium(kmed, "Si PIN       $", 6, 0,
                              isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.01, 0.01, static_cast<float*>(nullptr), 0);
  fIdtmed[6] = kmed;

  // The thermo insulating material of the box which contains the calorimeter module -> idtmed[7]
  TVirtualMC::GetMC()->Medium(kmed, "Thermo Insul.$", 7, 0,
                              isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, static_cast<float*>(nullptr), 0);
  fIdtmed[7] = kmed;

  // The Textolit which makes up the box which contains the calorimeter module      -> idtmed[707]
  TVirtualMC::GetMC()->Medium(kmed, "Textolit     $", 8, 0,
                              isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, static_cast<float*>(nullptr), 0);
  fIdtmed[8] = kmed;

  // FR4: The Plastic which makes up the frame of micromegas                        -> idtmed[708]
  TVirtualMC::GetMC()->Medium(kmed, "FR4 $", 9, 0,
                              isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.0001, static_cast<float*>(nullptr), 0);
  fIdtmed[9] = kmed;

  // The Composite Material for  micromegas                                         -> idtmed[709]
  TVirtualMC::GetMC()->Medium(kmed, "CompoMat   $", 10, 0,
                              isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, static_cast<float*>(nullptr), 0);
  fIdtmed[10] = kmed;

  // Copper                                                                         -> idtmed[710]
  TVirtualMC::GetMC()->Medium(kmed, "Copper     $", 11, 0,
                              isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.0001, static_cast<float*>(nullptr), 0);
  fIdtmed[11] = kmed;

  // G10: Printed Circuit material                                                  -> idtmed[12]
  TVirtualMC::GetMC()->Medium(kmed, "G10        $", 12, 0,
                              isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.01, static_cast<float*>(nullptr), 0);
  fIdtmed[12] = kmed;

  // The Lead                                                                       -> idtmed[712]

  TVirtualMC::GetMC()->Medium(kmed, "Lead      $", 13, 0,
                              isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, static_cast<float*>(nullptr), 0);
  fIdtmed[13] = kmed;

  // The gas mixture: ArCo2                                                         -> idtmed[715]

  TVirtualMC::GetMC()->Medium(kmed, "ArCo2      $", 16, 1,
                              isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.01, static_cast<float*>(nullptr), 0);
  fIdtmed[16] = kmed;

  // Stainless steel                                                                -> idtmed[17]
  TVirtualMC::GetMC()->Medium(kmed, "Steel     $", 17, 0,
                              isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.0001, static_cast<float*>(nullptr), 0);
  fIdtmed[17] = kmed;

  // Fibergalss                                                                     -> idtmed[18]
  TVirtualMC::GetMC()->Medium(kmed, "Fiberglass$", 18, 0,
                              isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, static_cast<float*>(nullptr), 0);
  fIdtmed[18] = kmed;

  // Cables in air                                                                  -> idtmed[718]
  TVirtualMC::GetMC()->Medium(kmed, "Cables    $", 19, 0,
                              isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, static_cast<float*>(nullptr), 0);
  fIdtmed[19] = kmed;

  // Air                                                                            -> idtmed[20]
  TVirtualMC::GetMC()->Medium(kmed, "Air          $", 99, 0,
                              isxfld, sxmgmx, 10.0, 1.0, 0.1, 0.1, 10.0, static_cast<float*>(nullptr), 0);
  fIdtmed[20] = kmed;
}

//____________________________________________________________________________
void Phos::CreateGeometry()
{
  // Create the PHOS geometry for Geant

  const Geometry* geom = Geometry::Instance();

  // Create a PHOS module.
  // Gsvolu accepts non-const params. Avoid modification geometry
  float params[4] = {geom->GetPHOSParams()[0], geom->GetPHOSParams()[1], geom->GetPHOSParams()[2], geom->GetPHOSParams()[3]};
  TVirtualMC::GetMC()->Gsvolu("PHOS", "TRD1", 20, params, 4); // 20-> idtmed[20]

  CreateGeometryforEMC();

  CreateGeometryforSupport();

  // --- Position  PHOS mdules in Hall ---
  int iXYZ, iAngle;
  float angle[3][2];
  geom->GetModuleAngles(angle);
  int idrotm = 0;
  TVirtualMC::GetMC()->Matrix(idrotm, angle[0][0], angle[0][1],
                              angle[1][0], angle[1][1],
                              angle[2][0], angle[2][1]);

  float pos[3] = {0};
  geom->GetModuleCenter(pos);
  TVirtualMC::GetMC()->Gspos("PHOS", 0, "Hall", pos[0], pos[1], pos[2], idrotm, "ONLY");
}
