#ifndef SIMPARAMS_H
#define SIMPARAMS_H

#include "TObject.h"

class SimParams : public TObject
{

 public:
  ~SimParams() = default;
  static SimParams* Instance()
  {
    if (!fgSimParams)
      fgSimParams = new SimParams();
    return fgSimParams;
  }

 public:
  // Parameters used in conversion of deposited energy to APD response
  float fLightYieldPerGeV = 526.; ///< Average number of photoelectrons per GeV: 1.983 gamma/MeV * 0.2655 PDE eff of APD

  // Parameters used in electronic noise calculation and thresholds (Digitizer)
  float fReadoutTime = 5.;           ///< Read-out time in ns for default simulaionts
  float fDeadTime = 20.;             ///< PHOS dead time (includes Read-out time i.e. mDeadTime>=mReadoutTime)
  float fReadoutTimePU = 2000.;      ///< Read-out time in ns if pileup simulation on in DigitizerSpec
  float fDeadTimePU = 30000.;        ///< PHOS dead time if pileup simulation on in DigitizerSpec
  bool fSmearLightCollection = true; ///< Mimic light collection
  bool fApplyTimeResolution = false; ///< Apply time resolution in digitization
  bool fApplyNonLinearity = false;   ///< Apply energy non-linearity in digitization
  bool fSimulateNoise = true;        ///< Simulate noise in digitization
  bool fApplyDigitization = false;   ///< Apply energy digitization in digitization
  float fAPDNoise = 0.001;           ///< RMS of APD noise
  float fDigitThreshold = 2.;        ///< minimal energy to keep digit in ADC counts
  float fADCwidth = 0.005;           ///< width of ADC channel in GeV
  float fTOFa = 0.5e-9;              ///< constant term of TOF resolution
  float fTOFb = 1.e-9;               ///< stohastic term of TOF resolution
  float fCellNonLineaityA = 1.091;      ///< Amp of cel non-linearity
  float fCellNonLineaityB = 0.;   ///< Energy scale of cel non-linearity
  float fCellNonLineaityC = 1.;      ///< Overall calibration

  short fZSthreshold = 1;         ///< Zero Suppression threshold
  float fTimeResolutionA = 2.e-9; ///< Time resolution parameter A (in sec)
  float fTimeResolutionB = 2.e-9; ///< Time resolution parameter B (in sec/GeV)
  float fTimeResThreshold = 0.5;  ///< threshold for time resolution calculation (in GeV)
  float fMinNoiseTime = -200.e-9; ///< minimum time in noise channels (in sec)
  float fMaxNoiseTime = 2000.e-9; ///< minimum time in noise channels (in sec)

  // Parameters used in Raw simulation
  float fSampleDecayTime = 0.091; ///< Time parameter in Gamma2 function (1/tau, 100.e-9/2.1e-6)

  // //Parameters used in raw data reconstruction
  short fSpikeThreshold = 100;          ///< Single spike >100 ADC channels
  short fBaseLine = 0;                  ///<
  short fPreSamples = 2;                ///< number of pre-samples readout before sample (if no pedestal subtrauction)
  short fMCOverflow = 970;              ///< Overflow level for MC simulations: 1023-(pedestal~50)
  float fTimeTick = 100.;               ///< ns to PHOS digitization step conversion
  float fTRUTimeTick = 25.;             ///< ns to PHOS TRU digitization step
  float fSampleTimeFitAccuracy = 1.e-3; // Abs accuracy of time fit of saturated samples (in 100ns tick units)
  float fSampleAmpFitAccuracy = 1.e-2;  // Relative accuracy of amp. fit
  short fNIterations = 5;               ///< maximal number of iterations in oveflow sample fit

  // Parameters used in clusterization
  float fLogWeight = 4.5;              ///< Cutoff used in log. weight calculation
  float fDigitMinEnergy = 0.010;       ///< Minimal energy of digits to be used in cluster (GeV)
  float fClusteringThreshold = 0.050;  ///< Minimal energy of digit to start clustering (GeV)
  float fLocalMaximumCut = 0.05;      ///< Minimal height of local maximum over neighbours
  int fUnfoldMaxSize = 100;            ///< maximal number of cells in cluster to be unfolded
  bool fUnfoldClusters = true;         ///< To perform cluster unfolding
  float fUnfogingEAccuracy = 1.e-2;    ///< Accuracy of energy calculation in unfoding prosedure (GeV)
  float fUnfogingXZAccuracy = 1.e-1;   ///< Accuracy of position calculation in unfolding procedure (cm)
  float fUnfogingChi2Accuracy = 1.e-2; ///< critical chi2/NDF
  int nNMaxIterations = 10;            ///< Maximal number of iterations in unfolding procedure
  float fCoreR = 3.5;                  ///< Radius to caluclate core energy
  float fChi2radiusCut = 4.;           ///< Radius to calculate chi2
  float fSortingDelta = 1.;            ///< used in sorting clusters
  float fCluNonLineaityA = 1.3;         ///< Amp of cluster non-linearity
  float fCluNonLineaityB = 0.;         ///< Energy scale of cluster non-linearity
  float fCluNonLineaityC = 1.;         ///< Overall calibration
  int fNPrimMax = 5;                   ///< maximal number of primary particles per cluster
  int fNLMMax = 30;                    ///< maximal numner of local maxima
  int fNMaxIterations = 100;           ///< max number of iterations in unfolding

 protected:
  static SimParams* fgSimParams;
  SimParams() = default;

  ClassDefNV(SimParams, 1);
};
#endif /* SIMPARAMS_H */
