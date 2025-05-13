#ifndef HALL_H
#define HALL_H

class Hall
{

 public:
  Hall() = default;
  virtual ~Hall() = default;
  void CreateGeometry();
  void Init() {}

 private:
  int fIdmix[2] = {0};  // material/mixtures
  int fIdtmed[2] = {0}; // media

  ClassDef(Hall, 1) // Class for ALICE experimental hall
};

#endif
