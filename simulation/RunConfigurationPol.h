#ifndef RUN_CONFIGURATIONPOL_H
#define RUN_CONFIGURATIONPOL_H

#include "TG4RunConfiguration.h"

class RunConfigurationPol : public TG4RunConfiguration
{
 public:
  RunConfigurationPol(const TString& physicsList = "emStandard",
                      const TString& specialProcess = "stepLimiter");
  virtual ~RunConfigurationPol() = default;

  // methods
  G4VUserPhysicsList* CreatePhysicsList();

 private:
  /// Option to use local magnetic field
};

#endif // RUN_CONFIGURATION_POL_H
