#ifndef GEOMETRY_H
#define GEOMETRY_H

// Geometry class  for PHOS : singleton

// --- ROOT system ---
#include "TObject.h"

// --- PHOS256 header files ---

class Geometry : public TObject
{

 public:
  // private constructor for singleton
  ~Geometry() = default;
  static const Geometry* Instance() { return fgGeom; }
  static const Geometry* Instance(float r, float theta);

  // local/master conversions
  static int RelToAbsId(int moduleNumber, int strip, int cell);
  static void RelToAbsId(const int relid[2], int& absId);

  // Return general PHOS parameters
  void GetModuleAngles(float angle[3][2]) const;
  void GetModuleCenter(float pos[3]) const;

  const float* GetPHOSParams() const { return fEMCParams; } // Half-sizes of PHOS trapecoid
  float GetOuterBoxSize(int index) const { return 2. * fEMCParams[index]; }
  float GetCrystalSize(int index) const { return 2. * fCrystalHalfSize[index]; }
  float GetCellStep() const { return 2. * fAirCellHalfSize[0]; }

  // Return ideal EMCA geometry parameters
  float GetIPtoCrystalSurface() const { return fIPtoCrystalSurface; }
  float GetIPtoOuterCoverDistance() const { return fIPtoOuterCoverDistance; }
  float GetAirGapLed() const { return fAirGapLed; }
  const Float_t* GetAirGapHalfSize() const { return fAirGapHalfSize; }
  const Float_t* GetAlCoverParams() const { return fAlCoverParams; }
  const float* GetCrystalHalfSize() const { return fCrystalHalfSize; }
  const Float_t* GetCoolerHalfSize() const { return fCoolerHalfSize; }
  const Float_t* GetFiberGlassHalfSize() const { return fFiberGlassHalfSize; }
  const float* GetStripHalfSize() const { return fStripHalfSize; }
  const float* GetAirCellHalfSize() const { return fAirCellHalfSize; }
  const float* GetWrappedHalfSize() const { return fWrappedHalfSize; }
  const float* GetAPDHalfSize() const { return fPinDiodeHalfSize; }
  const Float_t* GetOuterThermoParams() const { return fOuterThermoParams; }
  const float* GetPreampHalfSize() const { return fPreampHalfSize; }
  const Float_t* GetTCables1HalfSize() const { return fTCables1HalfSize; }
  const Float_t* GetTCables2HalfSize() const { return fTCables2HalfSize; }
  const Float_t* GetTSupport1HalfSize() const { return fTSupport1HalfSize; }
  const Float_t* GetTSupport2HalfSize() const { return fTSupport2HalfSize; }
  Float_t GetTSupportDist() const { return fTSupportDist; }
  const float* GetSupportPlateHalfSize() const { return fSupportPlateHalfSize; }
  const Float_t* GetSupportPlateInHalfSize() const { return fSupportPlateInHalfSize; }
  Float_t GetSupportPlateThickness() const { return fSupportPlateThickness; }
  const Float_t* GetFrameXHalfSize() const { return fFrameXHalfSize; }
  const Float_t* GetFrameZHalfSize() const { return fFrameZHalfSize; }
  const Float_t* GetFrameXPosition() const { return fFrameXPosition; }
  const Float_t* GetFrameZPosition() const { return fFrameZPosition; }
  const Float_t* GetFGupXHalfSize() const { return fFGupXHalfSize; }
  const Float_t* GetFGupXPosition() const { return fFGupXPosition; }
  const Float_t* GetFGupZHalfSize() const { return fFGupZHalfSize; }
  const Float_t* GetFGupZPosition() const { return fFGupZPosition; }
  const Float_t* GetFGlowXHalfSize() const { return fFGlowXHalfSize; }
  const Float_t* GetFGlowXPosition() const { return fFGlowXPosition; }
  const Float_t* GetFGlowZHalfSize() const { return fFGlowZHalfSize; }
  const Float_t* GetFGlowZPosition() const { return fFGlowZPosition; }
  const Float_t* GetFEEAirHalfSize() const { return fFEEAirHalfSize; }
  const Float_t* GetFEEAirPosition() const { return fFEEAirPosition; }
  const Float_t* GetInnerThermoHalfSize() const { return fInnerThermoHalfSize; }
  const Float_t* GetWarmAlCoverHalfSize() const { return fWarmAlCoverHalfSize; }
  const Float_t* GetWarmThermoHalfSize() const { return fWarmThermoHalfSize; }
  int GetNCellsXInStrip() const { return fNCellsXInStrip; }
  int GetNCellsZInStrip() const { return fNCellsZInStrip; }
  int GetNStripX() const { return fNStripX; }
  int GetNStripZ() const { return fNStripZ; }
  int GetNTSuppots() const { return fNTSupports; }
  int GetNPhi() const { return fNPhi; }
  int GetNZ() const { return fNZ; }
  int GetNCristalsInModule() const { return fNPhi * fNZ; }

  // void GetModuleCenter(TVector3& center, const char *det, int module) const;

 protected:
  Geometry(float r, float theta);
  void Init();

  //  void                     SetPHOSAngles();  // calculates the PHOS modules PHI angle

 private:
  static Geometry* fgGeom;

  static constexpr int fNPhi = 16; // Number of crystal units in X (phi) direction
  static constexpr int fNZ = 16;   // Number of crystal units in Z direction

  float fModR = 0.;         // Distance from IP to front surface of crystals
  float fModTheta = 0.;     // polar angle to module center
  float fCrystalShift;      // Distance from crystal center to front surface
  float fModuleCenter[3];   // xyz-position of the module center
  float fModuleAngle[3][2]; // polar and azymuth angles for 3 axes of modules

  // sizes
  float fStripHalfSize[3];          // Strip unit size/2
  float fAirCellHalfSize[3];        // geometry parameter
  float fWrappedHalfSize[3];        // geometry parameter
  float fSupportPlateHalfSize[3];   // geometry parameter
  float fSupportPlateInHalfSize[3]; // geometry parameter
  float fCrystalHalfSize[3];        // crystal size/2
  float fAirGapLed;                 // geometry parameter
  float fStripWallWidthOut;         // Side to another strip
  float fStripWallWidthIn;          // geometry parameter
  float fTyvecThickness;            // geometry parameter
  float fTSupport1HalfSize[3];      // geometry parameter
  float fTSupport2HalfSize[3];      // geometry parameter
  float fPreampHalfSize[3];         // geometry parameter
  float fPinDiodeHalfSize[3];       // Size of the PIN Diode

  float fOuterThermoParams[4];   // geometry parameter
  float fCoolerHalfSize[3];      // geometry parameter
  float fAirGapHalfSize[3];      // geometry parameter
  float fInnerThermoHalfSize[3]; // geometry parameter
  float fAlCoverParams[4];       // geometry parameter
  float fFiberGlassHalfSize[3];  // geometry parameter

  float fInnerThermoWidthX;      // geometry parameter
  float fInnerThermoWidthY;      // geometry parameter
  float fInnerThermoWidthZ;      // geometry parameter
  float fAirGapWidthX;           // geometry parameter
  float fAirGapWidthY;           // geometry parameter
  float fAirGapWidthZ;           // geometry parameter
  float fCoolerWidthX;           // geometry parameter
  float fCoolerWidthY;           // geometry parameter
  float fCoolerWidthZ;           // geometry parameter
  float fAlCoverThickness;       // geometry parameter
  float fOuterThermoWidthXUp;    // geometry parameter
  float fOuterThermoWidthXLow;   // geometry parameter
  float fOuterThermoWidthY;      // geometry parameter
  float fOuterThermoWidthZ;      // geometry parameter
  float fAlFrontCoverX;          // geometry parameter
  float fAlFrontCoverZ;          // geometry parameter
  float fFiberGlassSup2X;        // geometry parameter
  float fFiberGlassSup1X;        // geometry parameter
  float fFrameHeight;            // geometry parameter
  float fFrameThickness;         // geometry parameter
  float fAirSpaceFeeX;           // geometry parameter
  float fAirSpaceFeeZ;           // geometry parameter
  float fAirSpaceFeeY;           // geometry parameter
  float fTCables2HalfSize[3];    // geometry parameter
  float fTCables1HalfSize[3];    // geometry parameter
  float fWarmUpperThickness;     // geometry parameter
  float fWarmBottomThickness;    // geometry parameter
  float fWarmAlCoverWidthX;      // geometry parameter
  float fWarmAlCoverWidthY;      // geometry parameter
  float fWarmAlCoverWidthZ;      // geometry parameter
  float fWarmAlCoverHalfSize[3]; // geometry parameter
  float fWarmThermoHalfSize[3];  // geometry parameter
  float fFiberGlassSup1Y;        // geometry parameter
  float fFiberGlassSup2Y;        // geometry parameter
  float fTSupportDist;           // geometry parameter
  float fTSupport1Thickness;     // geometry parameter
  float fTSupport2Thickness;     // geometry parameter
  float fTSupport1Width;         // geometry parameter
  float fTSupport2Width;         // geometry parameter
  float fFrameXHalfSize[3];      // geometry parameter
  float fFrameZHalfSize[3];      // geometry parameter
  float fFrameXPosition[3];      // geometry parameter
  float fFrameZPosition[3];      // geometry parameter
  float fFGupXHalfSize[3];       // geometry parameter
  float fFGupXPosition[3];       // geometry parameter
  float fFGupZHalfSize[3];       // geometry parameter
  float fFGupZPosition[3];       // geometry parameter
  float fFGlowXHalfSize[3];      // geometry parameter
  float fFGlowXPosition[3];      // geometry parameter
  float fFGlowZHalfSize[3];      // geometry parameter
  float fFGlowZPosition[3];      // geometry parameter
  float fFEEAirHalfSize[3];      // geometry parameter
  float fFEEAirPosition[3];      // geometry parameter
  float fEMCParams[4];           // Half-sizes of PHOS trapecoid
  float fIPtoOuterCoverDistance; // Distances from interaction point to outer cover
  float fIPtoCrystalSurface;     // Distances from interaction point to Xtal surface

  float fSupportPlateThickness; // Thickness of the Aluminium support plate for Strip

  int fNCellsXInStrip = 8; // Number of cells in a strip unit in X
  int fNCellsZInStrip = 2; // Number of cells in a strip unit in Z
  int fNStripX = 2;        // Number of strip units in X
  int fNStripZ = 8;        // Number of strip units in Z
  int fNTSupports = 0;     // geometry parameter

  ClassDefNV(Geometry, 1) // PHOS geometry class
};

#endif // Geometry_H
