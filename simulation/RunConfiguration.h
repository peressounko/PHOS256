#ifndef RUN_CONFIGURATION_H
#define RUN_CONFIGURATION_H

#include "TG4RunConfiguration.h"

class RunConfiguration : public TG4RunConfiguration
{
 public:
  RunConfiguration(const TString& physicsList = "emStandard",
                   const TString& specialProcess = "stepLimiter");
  virtual ~RunConfiguration() = default;

  // methods
  virtual G4VUserDetectorConstruction* CreateDetectorConstruction() { return static_cast<G4VUserDetectorConstruction*>(nullptr); }

  // set methods
  void SetUseLocalMagField(Bool_t localMagField)
  {
    fUseLocalMagField = localMagField;
  }

 private:
  /// Option to use local magnetic field
  Bool_t fUseLocalMagField;
};

#endif // RUN_CONFIGURATION1_H
