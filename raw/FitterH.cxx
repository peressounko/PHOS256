/// Class to extract Amlitude and time from one PHOS256 sample
/// Simplest: amplitude is the maximal value, time - linear fit to flont at K-level
///
#include "FitterH.h"

//==================Implementation
double SampleShape(double* x, double* p)
{
  // parameters
  const double dectime = 0.091; ///< Time parameter in Gamma2 function (1/tau, 100.e-9/2.1e-6)
  const double shapenorm = 0.25 * dectime * dectime * exp(2.);
  double ped = p[0];
  double amp = p[1];
  double tim = p[2];
  if (x[0] < tim) {
    return ped;
  }
  double z = x[0] - tim;
  return ped + shapenorm * amp * z * z * exp(-z * dectime);
}

void FitterH::Init()
{
  mh = new TH1F("hSample", "", NMAXSAMPLES, 0., NMAXSAMPLES);
  mh->Sumw2();
  for (int i = 1; i <= NMAXSAMPLES; i++) {
    mh->SetBinError(i, 1.);
  }
  mfff = new TF1("SampleShape", SampleShape, 0., float(NMAXSAMPLES), 3);
}

Fitter::FitStatus FitterH::Evaluate(const unsigned short* signal, int sigLength)
{
  // Calculate signal parameters (energy, time, quality) from array of samples
  // Energy is a maximum sample minus pedestal 9
  // Time is the first time bin
  // Signal overflows is there are at least 3 samples of the same amplitude above 900

  const double minAmpToFit = 3.; // do not fit samples with max amplitude smaller or eq. than this

  if (sigLength == 0) {
    return kEmptyBunch;
  }
  if (sigLength < 2) {
    if (sigLength == 1) {
      mAmp = signal[0];
      mTime = 0;
    } else {
      mAmp = std::max(signal[0], signal[1]);
      mTime = 0;
    }
    return kOK;
  }

  float maxAmp = 0, minAmp = 1024;
  ushort previous = 0; // to exclude overflow bins
  int maxAt = 0, minAt = 0.;
  for (int i = 0; i < sigLength; i++) {
    mh->SetBinContent(sigLength - i, signal[i]);
    mh->SetBinError(sigLength - i, 1.);
    if (signal[i] > 900 && abs(signal[i] - previous) < 3) {
      mh->SetBinError(sigLength - i, 0.);     // To remove point from the fit.
      mh->SetBinError(sigLength - i + 1, 0.); // To remove point from the fit.
    }
    previous = signal[i];
    if (maxAmp < static_cast<float>(signal[i])) {
      maxAmp = signal[i];
      maxAt = i;
    }
    if (minAmp > static_cast<float>(signal[i])) {
      minAmp = signal[i];
      minAt = i;
    }
  }

  mPedestal = 0.;
  if (mPedSubtract) {
    for (int i = 1; i <= mPreSamples; i++) {
      mPedestal += mh->GetBinContent(i);
    }
    if (mPreSamples > 0)
      mPedestal /= mPreSamples;
  }

  // If too small signal, do not fit.
  if (maxAmp - mPedestal < minAmpToFit &&
      mPedestal - minAmp < minAmpToFit) { // Too small signal, no sense to fit
    mAmp = double(maxAmp) - mPedestal;
    mTime = sigLength - maxAt - 3; // 3 - approximate distance between start and maximum.
    return kOK;
  }

  mfff->SetParLimits(0, mPedestal - 2, mPedestal + 2.);
  mfff->SetParLimits(1, -20., 5000.);
  mfff->SetParLimits(2, -3., float(sigLength));
  if (maxAmp - mPedestal > mPedestal - minAmp) { // Positive signal
    mfff->SetParameters(mPedestal, maxAmp - mPedestal, maxAt - 2.);
  } else {
    mfff->SetParameters(mPedestal, minAmp - mPedestal, minAt - 2.);
  }
  // mfff->SetParameters(mPedestal,0.,2.);
  int status = mh->Fit(mfff, "QRL", "", 0., float(sigLength));
  if (status == 3) {
    status = mh->Fit(mfff, "QRL", "", 0., float(sigLength));
  }
  mAmp = mfff->GetParameter(1);
  mPedestal = mfff->GetParameter(0);
  mTime = mfff->GetParameter(2);
  mChi2 = mfff->GetChisquare();

  // TCanvas * c = (TCanvas*)gROOT->FindObject("FitCanvas");
  // c->cd();
  // mh->SetStats(0);
  // mh->Draw();
  // c->Update();
  // printf("Amp = %f, sigLength=%d, status=%d, maxAmp=%d\n",mAmp,sigLength, status,maxAmp);

  // getchar();

  if (mAmp > 2. * (maxAmp - mfff->GetParameter(0))) { // Fit failed
    mAmp = maxAmp - mfff->GetParameter(0);
    // printf("mAmp=%f, maxAmp=%d, ped=%f\n",mAmp,maxAmp,mfff->GetParameter(0));
    return kOK;
  }

  // Minuit2Minimizer::Minimize: migrad return status
  // 0 Valid minimum: The minimization converged successfully, and the covariance matrix is valid and positive-definite.
  // 1 Valid minimum, but covariance not accurate: The minimization converged, but there were issues with the covariance matrix calculation (e.g., forced positive-definite).
  // 2 Hesse is invalid: The Hesse matrix calculation failed or is invalid.
  // 3 Edm is above max: The Estimated Distance to Minimum (Edm) is larger than the specified tolerance, meaning convergence criteria were not met.
  // 4 Reached call limit: The maximum number of function calls was reached before convergence.
  // 5 Covariance is not positive defined: The calculated covariance matrix is not positive-definite, which implies issues with the minimum found (e.g., it might be a saddle point or a very flat region).

  // if(status>1){
  //  // if(mAmp>5){
  //    TCanvas * c = (TCanvas*)gROOT->FindObject("FitCanvas");
  //    c->cd();
  //    mh->SetStats(0);
  //    mh->Draw();
  //    c->Update();
  //    printf("Amp = %f, sigLength=%d, status=%d, maxAmp=%d\n",mAmp,sigLength, status,maxAmp);

  //    getchar();
  //  }

  if (status == 0 || status == 1) { // 0: OK, 1: covariance inaccurate
    return kOK;
  } else {
    mAmp = double(maxAmp) - mPedestal;
    mTime = sigLength - maxAt - 3; // 3 - approximate distance between start and maximum.
                                   // printf(" mAmp = %f, MaxAmp=%d, ped=%f \n",mAmp,maxAmp,ped);
    return kFitFailed;
  }

  // if (mOverflow) {
  //   return kOverflow;
  // } else {
  //   return kOK;
  // }
}
