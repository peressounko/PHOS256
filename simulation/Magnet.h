#ifndef MAGNET_H
#define MAGNET_H

class Magnet
{

 public:
  Magnet() = default;
  virtual ~Magnet() = default;
  void CreateGeometry();
  void Init() {}

 private:
 
  ClassDef(Magnet, 1) // Class for ALICE experimental hall
};

#endif
