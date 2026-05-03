void WriteMapping()
{
  // Create the RCU mapping files for PHOS256
  // The geometrical mapping within one FEE card is read
  // from the file CSP2ALTRO_new.dat prepared manually beforehand.
  //
  // The hardware address of the FEE channels is a 12-bit word formed 
  // from the branch number, FEE number, ALTRO chip number and ALTRO channel
  // as follows:
  // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  // |0|0|0|0| |       |     |       |
  // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  //          ^   FEE    chip  channel
  //       branch
  //
  // Max address: 1 1110 100 1111 = 3919
  // 
  // Author: Yuri Kharlov + Dmitri Peresunko
  // Date  : 24 April 2026
  // $Id$


  char string[128];
  uint xcell,zcell,csp,altro,chanHG,chanLG;
  uint hwAddress;
  const char *FEEmapping[2] = {"CSP2ALTRO_new.dat", "CSP2ALTRO_new.dat"};
  int map=0;

  FILE *fd = fopen("CSP2ALTRO_new.dat","r");
    
  FILE *fRCU = fopen("Mapping.data.unsorted","w");
  uint maxHWAddress=0;
  uint nHWaddress=0;
    
  while (fgets(string,128,fd)) {
    if (string[0]=='*') {
	    continue;
    }
    sscanf(string,"%d %d %d %d %d %d", &xcell,&zcell,&csp,&altro,&chanHG,&chanLG);
	  int iBranch=0; // iBranch==1; 
	  for (int iFEE=1; iFEE<=8; iFEE++) {
	    // High gain
	    hwAddress = chanHG | (altro<<4) | (iFEE<<7) | (iBranch<<11);
	    if (hwAddress > maxHWAddress){ 
        maxHWAddress=hwAddress;
      }
	    nHWaddress++;
	    fprintf(fRCU,"%4d %4d %4d %4d\n",
		    hwAddress,
		    xcell,
		    zcell + (iFEE-1)*2,1);
	    // Low gain
	    hwAddress = chanLG | (altro<<4) | (iFEE<<7) | (iBranch<<11);
	    if (hwAddress > maxHWAddress) maxHWAddress=hwAddress;
	    nHWaddress++;
	    fprintf(fRCU,"%4d %4d %4d %4d\n",
		    hwAddress,
		    xcell,
		    zcell + (iFEE-1)*2,0);
	  }
	}
  fclose(fd);

  fclose(fRCU);
    
  // Post-process the RCU mapping files
          
  // Add the number of channels and maximum HW address
      
  fRCU = fopen("Mapping.data","w");
  fprintf(fRCU,"%d\n%d\n",nHWaddress,maxHWAddress);
  fclose(fRCU);
      
  // Sort HW addresses      
  gSystem->Exec("sort -n Mapping.data.unsorted >> Mapping.data");
  gSystem->Exec("rm -f Mapping.data.unsorted");
}
