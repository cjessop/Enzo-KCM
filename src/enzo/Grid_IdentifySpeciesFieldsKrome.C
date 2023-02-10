/***********************************************************************
/
/  GRID CLASS (IDENTIFY THE MULTI-SPECIES FIELDS)
/
/  written by: KROME developers, based on the original ENZO routine
/  date:       2013
/
************************************************************************/
 
#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
 
 
int grid::IdentifySpeciesFieldsKrome(
 int &DeNum, int &HMNum, int &HINum, int &HeINum,
 int &H2INum, int &OINum, int &OHINum, int &H2OINum,
 int &O2INum, int &CHINum, int &COINum, int &CINum,
 int &CO2INum, int &CH2INum, int &CH3INum,
 int &CH4INum, int &HIINum, int &HeIINum,
 int &H2IINum, int &OIINum, int &OHIINum,
 int &H2OIINum, int &H3OIINum, int &O2IINum,
 int &CIINum, int &COIINum, int &HeIIINum
){

 DeNum = HMNum = HINum = HeINum = H2INum =
 OINum = OHINum = H2OINum = O2INum = CHINum =
 COINum = CINum = CO2INum = CH2INum = CH3INum =
 CH4INum = HIINum = HeIINum = H2IINum = OIINum =
 OHIINum = H2OIINum = H3OIINum = O2IINum =
 CIINum = COIINum = HeIIINum = 0;
 
  /* Find Fields for the KROME species */
 DeNum = FindField(ElectronDensity, FieldType, NumberOfBaryonFields);
 HMNum = FindField(HMDensity, FieldType, NumberOfBaryonFields);
 HINum = FindField(HIDensity, FieldType, NumberOfBaryonFields);
 HeINum = FindField(HeIDensity, FieldType, NumberOfBaryonFields);
 H2INum = FindField(H2IDensity, FieldType, NumberOfBaryonFields);
 OINum = FindField(OIDensity, FieldType, NumberOfBaryonFields);
 OHINum = FindField(OHIDensity, FieldType, NumberOfBaryonFields);
 H2OINum = FindField(H2OIDensity, FieldType, NumberOfBaryonFields);
 O2INum = FindField(O2IDensity, FieldType, NumberOfBaryonFields);
 CHINum = FindField(CHIDensity, FieldType, NumberOfBaryonFields);
 COINum = FindField(COIDensity, FieldType, NumberOfBaryonFields);
 CINum = FindField(CIDensity, FieldType, NumberOfBaryonFields);
 CO2INum = FindField(CO2IDensity, FieldType, NumberOfBaryonFields);
 CH2INum = FindField(CH2IDensity, FieldType, NumberOfBaryonFields);
 CH3INum = FindField(CH3IDensity, FieldType, NumberOfBaryonFields);
 CH4INum = FindField(CH4IDensity, FieldType, NumberOfBaryonFields);
 HIINum = FindField(HIIDensity, FieldType, NumberOfBaryonFields);
 HeIINum = FindField(HeIIDensity, FieldType, NumberOfBaryonFields);
 H2IINum = FindField(H2IIDensity, FieldType, NumberOfBaryonFields);
 OIINum = FindField(OIIDensity, FieldType, NumberOfBaryonFields);
 OHIINum = FindField(OHIIDensity, FieldType, NumberOfBaryonFields);
 H2OIINum = FindField(H2OIIDensity, FieldType, NumberOfBaryonFields);
 H3OIINum = FindField(H3OIIDensity, FieldType, NumberOfBaryonFields);
 O2IINum = FindField(O2IIDensity, FieldType, NumberOfBaryonFields);
 CIINum = FindField(CIIDensity, FieldType, NumberOfBaryonFields);
 COIINum = FindField(COIIDensity, FieldType, NumberOfBaryonFields);
 HeIIINum = FindField(HeIIIDensity, FieldType, NumberOfBaryonFields);


  /* Error if any not found. */
  if (DeNum<0 || HMNum<0 || HINum<0 || HeINum<0 ||
 H2INum<0 || OINum<0 || OHINum<0 || H2OINum<0 ||
 O2INum<0 || CHINum<0 || COINum<0 || CINum<0 ||
 CO2INum<0 || CH2INum<0 || CH3INum<0 || CH4INum<0 ||
 HIINum<0 || HeIINum<0 || H2IINum<0 || OIINum<0 ||
 OHIINum<0 || H2OIINum<0 || H3OIINum<0 ||
 O2IINum<0 || CIINum<0 || COIINum<0 || HeIIINum<0) {
    ENZO_VFAIL("De=%"ISYM", HM=%"ISYM", HI=%"ISYM", HeI=%"ISYM", H2I=%"ISYM", OI=%"ISYM", OHI=%"ISYM", H2OI=%"ISYM", O2I=%"ISYM", CHI=%"ISYM", COI=%"ISYM", CI=%"ISYM", CO2I=%"ISYM", CH2I=%"ISYM", CH3I=%"ISYM", CH4I=%"ISYM", HII=%"ISYM", HeII=%"ISYM", H2II=%"ISYM", OII=%"ISYM", OHII=%"ISYM", H2OII=%"ISYM", H3OII=%"ISYM", O2II=%"ISYM", CII=%"ISYM", COII=%"ISYM", HeIII=%"ISYM"\n",
	    DeNum, HMNum, HINum, HeINum, H2INum, OINum,
 OHINum, H2OINum, O2INum, CHINum, COINum,
 CINum, CO2INum, CH2INum, CH3INum, CH4INum,
 HIINum, HeIINum, H2IINum, OIINum, OHIINum,
 H2OIINum, H3OIINum, O2IINum, CIINum, COIINum,
 HeIIINum)
  }
 
  return SUCCESS;
}
