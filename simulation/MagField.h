#ifndef MAGFIELD_H
#define MAGFIELD_H
#include <TVirtualMagField.h>

/// \brief Definition of a uniform magnetic field within a given region
///
class MagField : public TVirtualMagField
{
 public:
  MagField(Double_t Bx, Double_t By, Double_t Bz);
  MagField() = default;
  virtual ~MagField() = default;

  virtual void Field(const Double_t* /*x*/, Double_t* B);

 private:
  MagField(const MagField&);
  MagField& operator=(const MagField&);

  Double_t fB[3] = {0.2, 0., 0.}; ///< Magnetic field vector

  ClassDef(MagField, 1) // Uniform magnetic field
};

#endif // MAGFIELD_H
