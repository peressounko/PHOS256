void RunRaw(){
  TGeoManager::Import("geometry_50_90.root");
  Geometry::Instance(120,30);
  SimParams * sp = SimParams::Instance();
  sp->fClusteringThreshold = 0.05;
  sp->fDigitMinEnergy = 0.03;
  sp->fADCwidth = 0.012;
  Reconstruction r;
//  r.SetRawFileName("../Muon/PHOS_run2998.raw");
  r.SetRawFileName("../Muon/PHOS_run2993.raw");
  r.SetFitter(1);
  r.SetNPreSamples(10);
  r.Run();

}
