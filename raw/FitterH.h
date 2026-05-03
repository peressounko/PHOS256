/// Class to extract Amlitude and time from one PHOS256 sample
/// Simplest: amplitude is the maximal value, time - linear fit to flont at K-level
///

#ifndef FITTERH_H
#define FITTERH_H
#include "TH1F.h"
#include "TF1.h"

#include "Fitter.h"

class FitterH final : public Fitter
{

 public:
  static constexpr int NMAXSAMPLES = 40; ///< maximal expected number of samples per bunch
  /// \brief Constructor
  FitterH() = default;

  /// \brief Destructor
  ~FitterH() = default;

  /// \brief Evaluation Amplitude and TOF
  Fitter::FitStatus Evaluate(const unsigned short* signal, int length);
  void Init();

 protected:
  void init();

 private:
  TH1F* mh = nullptr;
  TF1* mfff = nullptr;

  ClassDef(FitterH, 1);
}; // End of FitterH

#endif