#ifndef CALIBPARAMS_H
#define CALIBPARAMS_H
#include <array>
#include "TObject.h"

class TH2;

class CalibParams
{
 public:
  /// \brief Constructor
  CalibParams() = default;

  /// \brief Constructor for tests
  CalibParams(int test);

  /// \brief Constructor for tests
  CalibParams(CalibParams& a) = default;

  CalibParams& operator=(const CalibParams& other) = default;

  /// \brief Destructor
  ~CalibParams() = default;

  /// \brief Get High Gain energy calibration coefficients
  /// \param cellID Absolute ID of cell
  /// \return high gain energy calibration coefficient of the cell
  float GetGain(int cellID) const { return mGainCalib.at(cellID - OFFSET); }

  /// \brief Set High Gain energy calibration coefficient
  /// \param cellID Absolute ID of cell
  /// \param c is the calibration coefficient
  void SetGain(int cellID, float c) { mGainCalib.at(cellID - OFFSET) = c; }

  /// \brief Set High Gain energy calibration coefficients for one module in the form of 2D histogram
  /// \param 2D(64,56) histogram with calibration coefficients
  /// \param module number
  /// \return Is successful
  bool SetGain(TH2* h);

  /// \brief Get High Gain to Low Gain ratio calibration coefficients
  /// \param cellID Absolute ID of cell
  /// \return High Gain to Low Gain ratio of the cell
  float GetHGLGRatio(int cellID) const { return mHGLGRatio.at(cellID - OFFSET); }

  /// \brief Set High Gain to Low Gain ratio
  /// \param cellID Absolute ID of cell
  /// \param r is the calibration coefficient
  void SetHGLGRatio(int cellID, float r) { mHGLGRatio.at(cellID - OFFSET) = r; }

  /// \brief Set High Gain to Low Gain ratio for one module in the form of 2D histogram
  /// \param 2D(64,56) histogram with High Gain to Low Gain ratio
  /// \param module number
  /// \return Is successful
  bool SetHGLGRatio(TH2* h);

  /// \brief Get High Gain time calibration coefficients
  /// \param cellID Absolute ID of cell
  /// \return high gain time calibration coefficient of the cell
  float GetHGTimeCalib(int cellID) const { return mHGTimeCalib.at(cellID - OFFSET); }

  /// \brief Set High Gain time calibration coefficient
  /// \param cellID Absolute ID of cell
  /// \param t is the calibration coefficient
  void SetHGTimeCalib(int cellID, float t) { mHGTimeCalib.at(cellID - OFFSET) = t; }

  /// \brief Set High Gain time calibration coefficients for one module in the form of 2D histogram
  /// \param 2D(64,56) histogram with calibration coefficients
  /// \param module number
  /// \return Is successful
  bool SetHGTimeCalib(TH2* h);

  /// \brief Get Low Gain time calibration coefficient
  /// \param cellID Absolute ID of cell
  /// \return low gain time calibration coefficient of the cell
  float GetLGTimeCalib(int cellID) const { return mLGTimeCalib.at(cellID - OFFSET); }

  /// \brief Set time calibration coefficient
  /// \param cellID Absolute ID of cell
  /// \param t is the calibration coefficient
  void SetLGTimeCalib(int cellID, float t) { mLGTimeCalib.at(cellID - OFFSET) = t; }

  /// \brief Set Low Gain time calibration coefficients for one module in the form of 2D histogram
  /// \param 2D(64,56) histogram with calibration coefficients
  /// \param module number
  /// \return Is successful
  bool SetLGTimeCalib(TH2* h);

 private:
  static constexpr int NCHANNELS = 256;      ///< Number of channels = 14336-1792
  static constexpr int OFFSET = 1;           ///< Non-existing channels 56*64*0.5+1
  std::array<float, NCHANNELS> mGainCalib;   ///< Container for the gain calibration coefficients
  std::array<float, NCHANNELS> mHGLGRatio;   ///< Container for the High Gain to Low Gain ratios
  std::array<float, NCHANNELS> mHGTimeCalib; ///< Container for the High Gain time calibration coefficients
  std::array<float, NCHANNELS> mLGTimeCalib; ///< Container for the Low Gain time calibration coefficients

  ClassDefNV(CalibParams, 1);
};

#endif
