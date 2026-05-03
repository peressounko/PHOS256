/// Class to extract Amlitude and time from one PHOS256 sample
/// Simplest: amplitude is the maximal value, time - linear fit to flont at K-level
///
#include <cmath>
#include "Fitter.h"

//==================Implementation
void Fitter::Init()
{
  // Calculate integral of a sample with unit height
  // parameters
  const double dectime = 0.091; ///< Time parameter in Gamma2 function (1/tau, 100.e-9/2.1e-6)
  const double shapenorm = 0.25 * dectime * dectime * std::exp(2.);
  double z = 0.9;
  mSampleInt[0] = shapenorm * z * z * std::exp(-z * dectime);
  for (int i = 1; i < NMAXSAMPLES; i++) {
    z = i + 0.9;
    mSampleInt[i] = mSampleInt[i - 1] + shapenorm * z * z * std::exp(-z * dectime);
  }
}

Fitter::FitStatus Fitter::Evaluate(const unsigned short* signal, int sigLength)
{
  // Calculate signal parameters (energy, time, quality) from array of samples
  // Energy is a maximum sample minus pedestal 9
  // Time is the first time bin
  // Signal overflows is there are at least 3 samples of the same amplitude above 900

  // TCanvas * c = (TCanvas*)gROOT->FindObject("FtiCanvas");
  // c->cd();
  // mmh->Reset();
  // for(int i=0; i<sigLength; i++){
  //   mmh->SetBinContent(sigLength-i,signal[i]);
  //   mmh->SetBinError(sigLength-i,0);
  // }

  if (sigLength == 0) {
    return kEmptyBunch;
  }

  float pedMean = 0;
  int nPed = 0;
  mMaxSample = 0;
  int nMax = 0; // number of consequitive maximal samples
  bool spike = false;
  mOverflow = false;
  int integral = 0;

  int ap = -1, app = -1; // remember previous values to evaluate spikes
  for (int it = 0; it < sigLength; ++it) {
    unsigned short a = signal[sigLength - it - 1]; // inverse signal time order
    if (mPedSubtract) {
      if (nPed < mPreSamples) {
        nPed++;
        pedMean += a;
      } else {
        mAmp += a;
      }
    } else {
      mAmp += a;
    }
    if (a > mMaxSample) {
      mMaxSample = a;
      nMax = 0;
    }
    if (a == mMaxSample) {
      nMax++;
    }
    // check if there is a spike
    if (app >= 0 && ap >= 0) {
      spike |= (2 * ap - (a + app) > 2 * mSpikeThreshold);
    }
    app = ap;
    ap = a;
  }
  mAmp = (float)mMaxSample; // maximal amplitude

  float pedestal = 0;
  if (mPedSubtract) {
    if (nPed > 0) {
      pedMean /= nPed;
    } else {
      mAmp = 0.;
      mTime = -2;
      mOverflow = false;
      return kBadPedestal;
    }
  } else {
    pedMean = 0.;
  }

  //   if(mPedSubtract){
  //     // int l = sigLength - mPreSamples;
  //     // if(l>0){
  //     //   mAmp-=l*pedMean;
  //     //   mAmp/= mSampleInt[l-1];
  //     // }
  //     // else{
  //     //   mAmp = 0.;
  //     // }
  //   }
  //   else{
  //     mAmp /= mSampleInt[sigLength-1];
  //     // mAmp /= mSampleInt[39];
  //   }

  // printf("ped = %f, Amp=%f, sum=%f, sigLength=%d, mPreSamples=%d, mPedSubtract=%d\n",pedMean,mAmp,pedMean+mAmp,sigLength,mPreSamples,mPedSubtract);

  //   mmh->Draw();
  //   c->Update();
  //   getchar();

  // if (spike) {
  //   mTime = -2;
  //   mOverflow = false;
  //   return kSpike;
  // }

  if (mMaxSample > 900 && nMax > 2) {
    mOverflow = true;
  }

  if (mAmp < mBaseLine) {
    mAmp = 0;
  }

  // Evaluate time
  mTime = -2;
  const int nLine = 6;       // Parameters of fitting
  const float eMinTOF = 10.; // Choosed from beam-test and cosmic analyis
  const float kAmp = 0.35;   // Result slightly depends on them, so no getters
  // Avoid too low peak:
  if (mAmp < eMinTOF) {
    return kOK; // use estimated time
  }

  // Find index posK (kLevel is a level of "timestamp" point Tk):
  int posK = sigLength - 1; // last point before crossing k-level
  float levelK = pedestal + kAmp * mAmp;
  while (posK >= 0 && signal[posK] <= levelK) {
    posK--;
  }
  posK++;

  if (posK == 0 || posK == sigLength - 1) {
    return kNoTime; //
  }

  // Find crossing point by solving linear equation (least squares method)
  int np = 0;
  int iup = posK - 1;
  int idn = posK;
  double sx = 0., sy = 0., sxx = 0., sxy = 0.;
  double x, y;

  while (np < nLine) {
    // point above crossing point
    if (iup >= 0) {
      x = sigLength - iup - 1;
      y = signal[iup];
      sx += x;
      sy += y;
      sxx += (x * x);
      sxy += (x * y);
      np++;
      iup--;
    }
    // Point below crossing point
    if (idn < sigLength) {
      if (signal[idn] < pedestal) {
        idn = sigLength - 1; // do not scan further
        idn++;
        continue;
      }
      x = sigLength - idn - 1;
      y = signal[idn];
      sx += x;
      sy += y;
      sxx += (x * x);
      sxy += (x * y);
      np++;
      idn++;
    }
    if (idn >= sigLength && iup < 0) {
      break; // can not fit futher
    }
  }

  double det = np * sxx - sx * sx;
  if (det == 0) {
    return kNoTime;
  }
  if (np == 0) {
    return kEmptyBunch;
  }
  double c1 = (np * sxy - sx * sy) / det; // slope
  double c0 = (sy - c1 * sx) / np;        // offset
  if (c1 == 0) {
    return kNoTime;
  }

  // Find where the line cross kLevel:
  mTime += (levelK - c0) / c1 - 5.; // 5: mean offset between k-Level and start times

  if (mOverflow) {
    return kOverflow;
  } else {
    return kOK;
  }
}
