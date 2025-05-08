#ifndef HIT_H
#define HIT_H 1

#include "TObject.h"

class Hit : public TObject {
public:
   /** Default constructor **/
   MpdEmcHit() = default ;

   // /** Constructor with hit parameters (1)**/
   // MpdEmcHit(int detectorID, TVector3 pos, TVector3 dpos, int refIndex, int flag);

   // /** Constructor with hit parameters (2) [not the flag]**/
   // MpdEmcHit(int detectorID, TVector3 pos, TVector3 dpos, int refIndex);

   // MpdEmcHit(Uint sec, Uint row, Uint supMod, Uint mod, float e) {
   // } // TEMPORARY FIX: was not implemented -> undefined reference
   // MpdEmcHit(Uint sec, Uint row, Uint supMod, Uint mod, float e, float time);
   // MpdEmcHit(Uint detID, float e, float time);

   ~MpdEmcHit() = default;

   // void Print(const Option_t *opt = 0) const;

   // int GetFlag() const { return fFlag; };

   // int GetSec() const { return fDetectorID / (1000 * 12); };

   // int GetMod() const { return fDetectorID - (fDetectorID / 1000) * 1000; };

   // int GetSupMod() const { return -1; };

   // int GetRow() const
   // {
   //    return fDetectorID / 1000;
   //    ;
   // };

   // float GetE() const { return fE; };

   // float GetTime() const { return fTime; };

   // float GetRhoCenter() const { return sqrt(fX * fX + fY * fY); }

   // float GetZCenter() const { return fZ; };

   // float GetPhiCenter() const { return fPhiCenter; };

   // float GetThetaCenter() const { return fThetaCenter; };

   // int GetTrackId() const { return fTrackID; };

   // int GetPdg() const { return fPDG; };

   // int GetNumTracks() const { return fNumTracks; };

   // int GetDetectorID() const { return fDetectorID; };

   // float GetX() const { return fX; };

   // float GetY() const { return fY; };

   // float GetZ() const { return fZ; };

   // void SetFlag(int flag) { fFlag = flag; };

   // void SetEnergy(float e) { fE = e; };

   // void SetTime(float time) { fTime = time; };

   // void IncreaseEnergy(float e) { fE += e; };

   // void IncreaseEnergyTime(float timeEnergy) { fTime += timeEnergy; };

   // void SetTrackId(int id) { fTrackID = id; };

   // void SetPdg(int pdg) { fPDG = pdg; };

   // void SetNumTracks(int n) { fNumTracks = n; };

   // void SetPhiCenter(float phi) { fPhiCenter = phi; };

   // void SetZCenter(float z) { fZ = z; };

   // void SetThetaCenter(float theta) { fThetaCenter = theta; };

   // void SetX(float x) { fX = x; };

   // void SetY(float y) { fY = y; };

   // void SetZ(float z) { fZ = z; };

protected:
   float fE;       // energy
   float fTime;    // hit mean time
   int   fCellID;  // detector id of each hit
   int   fLabel;   // number of tracks, which have contribution in module
   ClassDef(Hit, 1);
};

#endif
