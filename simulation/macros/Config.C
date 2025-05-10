//#include <TPGCode.h>

void Config()
{

   // libraries required by geant321
    // gSystem->Load("libgeant321");

    new     TGeant4("C++ Interface to Geant4");

    TGeant4 *geant = dynamic_cast<TGeant4*>gMC;

    // Set External decayer
//    TVirtualMCDecayer *decayer = new AliDecayerPythia();
    // decayer->SetForceDecay(kAll);
    // decayer->Init();
    // gMC->SetExternalDecayer(decayer);

    //
    //
    //=======================================================================
    // ******* GEANT STEERING parameters FOR ALICE SIMULATION *******
    // geant3->SetTRIG(1);         //Number of events to be processed 
    // geant3->SetSWIT(4, 10);
    // geant3->SetDEBU(0, 0, 1);
    // //geant3->SetSWIT(2,2);
    // geant3->SetDCAY(1);
    // geant3->SetPAIR(1);
    // geant3->SetCOMP(1);
    // geant3->SetPHOT(1);
    // geant3->SetPFIS(0);
    // geant3->SetDRAY(0);
    // geant3->SetANNI(1);
    // geant3->SetBREM(1);
    // geant3->SetMUNU(1);
    // geant3->SetCKOV(1);
    // geant3->SetHADR(1);         //Select pure GEANH (HADR 1) or GEANH/NUCRIN (HADR 3)
    // geant3->SetLOSS(2);
    // geant3->SetMULS(1);
    // geant3->SetRAYL(1);
    // geant3->SetAUTO(1);         //Select automatic STMIN etc... calc. (AUTO 1) or manual (AUTO 0)
    // geant3->SetABAN(0);         //Restore 3.16 behaviour for abandoned tracks
    // geant3->SetOPTI(2);         //Select optimisation level for GEANT geometry searches (0,1,2)
    // geant3->SetERAN(5.e-7);

    // Float_t cut = 1.e-3;        // 1MeV cut by default
    // Float_t tofmax = 1.e10;

    //             GAM ELEC NHAD CHAD MUON EBREM MUHAB EDEL MUDEL MUPA TOFMAX
    // geant3->SetCUTS(cut, cut, cut, cut, cut, cut, cut, cut, cut, cut,
                    // tofmax);
    
    //Configure generator
    GenBox *gener = new GenBox(1);
    gener->SetPart(kGamma);
    gener->SetMomentumRange(10,11.);
    gener->SetPhiRange(270.5,270.7);
    gener->SetThetaRange(90.5,90.7);
    gener->Init();
 
    // Magnetic Field 
    //TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 1., 1., AliMagF::k5kG));

    //=================== Create detector(s) =============================
    BODY *Hall = new Hall("Experimental hall");

    Phos *ph = new Phos(0.75,20.); //Position: radius in m, polar angle, degrees

}
