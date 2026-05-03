/// Class to extract Amlitude and time from one PHOS256 sample
/// Simplest: amplitude is the maximal value, time - linear fit to flont at K-level
///

#ifndef FITTER_H
#define FITTER_H

class Fitter
{

 public:
  enum FitStatus { kOK,
                   kNotEvaluated,
                   kEmptyBunch,
                   kOverflow,
                   kSpike,
                   kNoTime,
                   kFitFailed,
                   kBadPedestal,
                   kManyBunches };

 public:
  /// \brief Constructor
  Fitter() = default;

  /// \brief Destructor
  virtual ~Fitter() = default;

  /// \brief Evaluation Amplitude and TOF
  /// return status -1: not evaluated/empty bunch;
  ///                0: OK;
  ///                1: overflow;
  ///                4: single spikes
  ///                3: too large RMS;
  FitStatus Evaluate(const unsigned short* signal, int length);
  void Init();

  /// \brief estimate and subtract pedestals from pre-samples
  void SetPedSubtract(bool toSubtruct = false) { mPedSubtract = toSubtruct; }

  void SetNPresamples(int n) { mPreSamples = n; }

  /// \brief amplitude in last fitted sample
  float GetAmp() const { return mAmp; }

  /// \brief amplitude in last fitted sample
  float GetPedestal() const { return mPedestal; }

  /// \brief Chi2/NDF of last performed fit
  float GetChi2() const { return mChi2; }

  /// \brief time in last fitted sample
  float GetTime() const { return mTime; }

  /// \brief is last fitted sample has overflow
  bool IsOverflow() const { return mOverflow; }

 public:
  static constexpr int NMAXSAMPLES = 40; ///< maximal expected number of samples per bunch

 protected:
  bool mPedSubtract = false;         ///< should one evaluate and subtract pedestals
  bool mOverflow;                    ///< is last sample saturated
  FitStatus mStatus = kNotEvaluated; ///< status of last evaluated sample: -1: not yet evaluated; 0: OK; 1: overflow; 2: too large RMS; 3: single spikes
  short mPreSamples = 0;             ///< number of presamples
  short mMaxSample = 0;              ///< maximal sample
  float mPedestal = 0.;              ///< Pedestal
  float mAmp;                        ///< amplitude of last processed sample
  float mTime = 0;                   ///< time of last processed sample
  float mChi2 = 1e10;                ///< chi2 calculated in last fit
  short mSpikeThreshold = 100;
  short mBaseLine = 0;
  float mSampleInt[NMAXSAMPLES]; ///< Integral of first i samples of a signal
  // ClassDef(Fitter, 1)
}; // End of Fitter

#endif