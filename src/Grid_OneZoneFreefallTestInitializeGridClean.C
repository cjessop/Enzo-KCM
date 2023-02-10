/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID TO A UNIFORM POOL OF GAS)
/
/  written by: Britton Smith
/  date:       February, 2008
/  modified1:  Chris Jessop
/  date:       August 2020
/
/  PURPOSE: Only different from InitializeUniformGrid in how 
/           species fractions are initialized.
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::OneZoneFreefallTestInitializeGrid(float InitialDensity,
					    float MinimumEnergy,
					    float MaximumEnergy,
					    float MinimumMetallicity,
					    float MaximumMetallicity)
{
  /* declarations */

  int dim, i, j, k, size, GCM, index;

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum, CINum, CIINum, OINum, OIINum, CHINum, CH2INum, COINum,
      OHINum, H2OINum, O2INum, COIINum, O2IINum, OHIINum, H2OIINum, H3OIINum,
      CH3INum, CH4INum, CO2INum ,MetalNum;

  int ExtraField[2];

  /* initialize density and pressure history storage */
  freefall_density = new float*[3];
  freefall_pressure = new float*[3];

  for (i = 0;i < 3;i++) {
    freefall_density[i] = NULL;
    freefall_pressure[i] = NULL;
  }

  /* create fields */
 
  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  int vel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1)
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;

  int colorfields = NumberOfBaryonFields;

  // Enzo's standard multispecies (primordial chemistry - H, D, He)
  if (MultiSpecies) {
    FieldType[DeNum     = NumberOfBaryonFields++] = ElectronDensity;
    FieldType[HINum     = NumberOfBaryonFields++] = HIDensity;
    FieldType[HIINum    = NumberOfBaryonFields++] = HIIDensity;
    FieldType[HeINum    = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum   = NumberOfBaryonFields++] = HeIIDensity;
    FieldType[HeIIINum  = NumberOfBaryonFields++] = HeIIIDensity;
    if (MultiSpecies > 1) {
      FieldType[HMNum   = NumberOfBaryonFields++] = HMDensity;
      FieldType[H2INum  = NumberOfBaryonFields++] = H2IDensity;
      FieldType[H2IINum = NumberOfBaryonFields++] = H2IIDensity;
    }
    if (MultiSpecies > 2) {
      FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
      FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
      FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
    }
  }

  if (KCM) {
    FieldType[DeNum   = NumberOfBaryonFields++] = ElectronDensity;
    FieldType[HINum   = NumberOfBaryonFields++] = HIDensity;
    FieldType[HIINum  = NumberOfBaryonFields++] = HIIDensity;
    FieldType[HeINum  = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum = NumberOfBaryonFields++] = HeIIDensity;
    FieldType[HeIIINum= NumberOfBaryonFields++] = HeIIIDensity;
    FieldType[HMNum   = NumberOfBaryonFields++] = HMDensity;
    FieldType[H2INum  = NumberOfBaryonFields++] = H2IDensity;
    FieldType[H2IINum = NumberOfBaryonFields++] = H2IIDensity;
    FieldType[CINum   = NumberOfBaryonFields++] = CIDensity;
    FieldType[CIINum  = NumberOfBaryonFields++] = CIIDensity;
    FieldType[OINum   = NumberOfBaryonFields++] = OIDensity;
    FieldType[OIINum  = NumberOfBaryonFields++] = OIIDensity;
    FieldType[CHINum  = NumberOfBaryonFields++] = CHIDensity;
    FieldType[CH2INum = NumberOfBaryonFields++] = CH2IDensity;
    FieldType[COINum  = NumberOfBaryonFields++] = COIDensity;
    FieldType[OHINum  = NumberOfBaryonFields++] = OHIDensity;
    FieldType[H2OINum = NumberOfBaryonFields++] = H2OIDensity;
    FieldType[O2INum  = NumberOfBaryonFields++] = O2IDensity;
    FieldType[COIINum = NumberOfBaryonFields++] = COIIDensity;
    FieldType[O2IINum = NumberOfBaryonFields++] = O2IIDensity;
    FieldType[OHIINum = NumberOfBaryonFields++] = OHIIDensity;
    FieldType[H2OIINum= NumberOfBaryonFields++] = H2OIIDensity;
    FieldType[H3OIINum= NumberOfBaryonFields++] = H3OIIDensity;
    FieldType[CH3INum = NumberOfBaryonFields++] = CH3IDensity;
    FieldType[CH4INum = NumberOfBaryonFields++] = CH4IDensity;
    FieldType[CO2INum = NumberOfBaryonFields++] = CO2IDensity;
  }
  

  //  Metal fields, including the standard 'metallicity' as well 
  // as two extra fields
  if (TestProblemData.UseMetallicityField) {
    FieldType[MetalNum = NumberOfBaryonFields++] = Metallicity;

    if (MultiMetals){
      FieldType[ExtraField[0] = NumberOfBaryonFields++] = ExtraType0;
      FieldType[ExtraField[1] = NumberOfBaryonFields++] = ExtraType1;
    }
  }

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (GridRank > 2) {
    ENZO_FAIL("GridRank must be 1 or 2 for OneZoneFreefallTest.\n");
  }

  // Get the units so we can convert temperature later
 
  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1;
 
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits,
	       Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  /* compute size of fields */
  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* allocate fields */

  
  this->AllocateGrids();

  /* Set up a 1d grid that varies in energy or a 2d grid that also varies in 
     metallicity. */

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) { // nothing
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) { // metallicity
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) { // energy

	index = i + j*GridDimension[0] + k*GridDimension[0]*GridDimension[1];

	BaryonField[0][index] = InitialDensity;

	/* set energy */
	BaryonField[1][index] = MinimumEnergy * 
	  POW((MaximumEnergy / MinimumEnergy), (float(i - GridStartIndex[0]) / 
						(GridEndIndex[0] - GridStartIndex[0])));
	if (DualEnergyFormalism) 
	  BaryonField[2][index] = BaryonField[1][index];

	/* set metallicity */
	if (TestProblemData.UseMetallicityField && (GridRank > 1)) {
	  BaryonField[MetalNum][index] = BaryonField[0][index] * 
	    CoolData.SolarMetalFractionByMass *
	    MinimumMetallicity * POW((MaximumMetallicity / MinimumMetallicity), 
				     (float(j - GridStartIndex[1]) / 
				      (GridEndIndex[1] - GridStartIndex[1])));
	  if (MultiMetals) {
	    BaryonField[MetalNum+1][index] = BaryonField[MetalNum][index];
	    BaryonField[MetalNum+2][index] = BaryonField[MetalNum][index];
	  }
	}

      }
    }
  }
 
  /* set velocities */
 
  for (dim = 0; dim < GridRank; dim++)
    for (i = 0; i < size; i++)
      BaryonField[vel+dim][i] = 0.0;
 
  /* set density of color fields to user-specified values (if user doesn't specify, 
     the defaults are set in SetDefaultGlobalValues.  Do some minimal amount of error
     checking to try to ensure charge conservation when appropriate */
  for (i = 0; i < size; i++){

    // Set multispecies fields!
    // this attempts to set them such that species conservation is maintained,
    // using the method in CosmologySimulationInitializeGrid.C
    if(MultiSpecies){

      BaryonField[HINum][i] = TestProblemData.HI_Fraction * 
	TestProblemData.HydrogenFractionByMass * BaryonField[0][i];
 
      BaryonField[HeINum][i] =  TestProblemData.HeII_Fraction *
	BaryonField[0][i] * (1.0-TestProblemData.HydrogenFractionByMass);

      BaryonField[HeIINum][i] = TestProblemData.HeIII_Fraction *
	BaryonField[0][i] * (1.0-TestProblemData.HydrogenFractionByMass);

      BaryonField[HeIIINum][i] =
	(1.0 - TestProblemData.HydrogenFractionByMass) * BaryonField[0][i] -
	BaryonField[HeINum][i] - BaryonField[HeIINum][i];

      if(MultiSpecies > 1){
	BaryonField[HMNum][i] = TestProblemData.HM_Fraction *
	  TestProblemData.HydrogenFractionByMass * BaryonField[0][i];

	BaryonField[H2INum][i] = 2 * TestProblemData.H2I_Fraction *
	  TestProblemData.HydrogenFractionByMass * BaryonField[0][i];

	BaryonField[H2IINum][i] = 2 * TestProblemData.H2II_Fraction * 
	  TestProblemData.HydrogenFractionByMass * BaryonField[0][i];
      }

      // HII density is calculated by subtracting off the various ionized fractions
      // from the total
      BaryonField[HIINum][i] = TestProblemData.HydrogenFractionByMass * BaryonField[0][i]
	- BaryonField[HINum][i];
      if (MultiSpecies > 1)
	BaryonField[HIINum][i] -= (BaryonField[HMNum][i] + BaryonField[H2IINum][i]
				  + BaryonField[H2INum][i]);

      // Electron "density" (remember, this is a factor of m_p/m_e scaled from the 'normal'
      // density for convenience) is calculated by summing up all of the ionized species.
      // The factors of 0.25 and 0.5 in front of HeII and HeIII are to fix the fact that we're
      // calculating mass density, not number density (because the BaryonField values are 4x as
      // heavy for helium for a single electron)
      BaryonField[DeNum][i] = BaryonField[HIINum][i] +
	0.25*BaryonField[HeIINum][i] + 0.5*BaryonField[HeIIINum][i];
      if (MultiSpecies > 1)
	BaryonField[DeNum][i] += 0.5*BaryonField[H2IINum][i] -
	  BaryonField[HMNum][i];

      // Set deuterium species (assumed to be a negligible fraction of the total, so not
      // counted in the conservation)
      if(MultiSpecies > 2){
	BaryonField[DINum ][i]  = CoolData.DeuteriumToHydrogenRatio * BaryonField[HINum][i];
	BaryonField[DIINum][i] = CoolData.DeuteriumToHydrogenRatio * BaryonField[HIINum][i];
	BaryonField[HDINum][i] = 0.75 * CoolData.DeuteriumToHydrogenRatio * BaryonField[H2INum][i];
      }

    } // if(MultiSpecies)

    if (KCM) {

      BaryonField[HINum][i] = TestProblemData.HI_Fraction *
        TestProblemData.HydrogenFractionByMass * BaryonField[0][i];
 
      BaryonField[HeINum][i] =  TestProblemData.HeII_Fraction *
        BaryonField[0][i] * (1.0-TestProblemData.HydrogenFractionByMass);

      BaryonField[HeIINum][i] = TestProblemData.HeIII_Fraction *
        BaryonField[0][i] * (1.0-TestProblemData.HydrogenFractionByMass);

      BaryonField[HeIIINum][i] =
        (1.0 - TestProblemData.HydrogenFractionByMass) * BaryonField[0][i] -
        BaryonField[HeINum][i] - BaryonField[HeIINum][i];

      BaryonField[HMNum][i] = TestProblemData.HM_Fraction *
        TestProblemData.HydrogenFractionByMass * BaryonField[0][i];

      BaryonField[H2INum][i] = 2 * TestProblemData.H2I_Fraction *
        TestProblemData.HydrogenFractionByMass * BaryonField[0][i];

      BaryonField[H2IINum][i] = 2 * TestProblemData.H2II_Fraction *
        TestProblemData.HydrogenFractionByMass * BaryonField[0][i];

      BaryonField[HIINum][i] = TestProblemData.HydrogenFractionByMass * BaryonField[0][i]
        - BaryonField[HINum][i];

      BaryonField[DeNum][i] = BaryonField[HIINum][i] +
        0.25*BaryonField[HeIINum][i] + 0.5*BaryonField[HeIIINum][i];
      if (MultiSpecies > 1)
        BaryonField[DeNum][i] += 0.5*BaryonField[H2IINum][i] -
          BaryonField[HMNum][i];

      BaryonField[CINum][i] = TestProblemData.CI_Fraction *
        TestProblemData.CarbonFractionByMass * BaryonField[0][i];

      BaryonField[CIINum][i] = TestProblemData.CII_Fraction * 
        TestProblemData.CarbonFractionByMass * BaryonField[0][i];

      BaryonField[OINum][i] = TestProblemData.OI_Fraction * 
        TestProblemData.OxygenFractionByMass * BaryonField[0][i];

      BaryonField[OIINum][i] = TestProblemData.OII_Fraction * 
        TestProblemData.OxygenFractionByMass * BaryonField[0][i];

      BaryonField[CHINum][i] = TestProblemData.CHI_Fraction * 
        TestProblemData.MethineFractionByMass * BaryonField[0][i];

      BaryonField[CH2INum][i] = TestProblemData.CH2I_Fraction * 
        TestProblemData.Methine2FractionByMass * BaryonField[0][i];

      BaryonField[COINum][i] = TestProblemData.COI_Fraction * 
        TestProblemData.CarbonMonoxideFractionByMass * BaryonField[0][i];

      BaryonField[OHINum][i] = TestProblemData.OHI_Fraction * 
        TestProblemData.HydroxyFractionByMass * BaryonField[0][i];

      BaryonField[H2OINum][i] = TestProblemData.H2OI_Fraction * 
        TestProblemData.WaterFractionByMass * BaryonField[0][i];

      BaryonField[O2INum][i] = TestProblemData.O2I_Fraction * 
        TestProblemData.Oxygen2FractionByMass * BaryonField[0][i];

      BaryonField[COIINum][i] = TestProblemData.COII_Fraction * 
        TestProblemData.CarbonMonoxideFractionByMass * BaryonField[0][i];

      BaryonField[O2IINum][i] = TestProblemData.O2II_Fraction * 
        TestProblemData.Oxygen2FractionByMass * BaryonField[0][i];

      BaryonField[OHIINum][i] = TestProblemData.OHII_Fraction * 
        TestProblemData.HydroxyFractionByMass * BaryonField[0][i];

      BaryonField[H2OIINum][i] = TestProblemData.H2OII_Fraction * 
        TestProblemData.WaterFractionByMass * BaryonField[0][i];

      BaryonField[H3OIINum][i] = TestProblemData.H3OII_Fraction * 
        TestProblemData.HydroniumFractionByMass * BaryonField[0][i];

      BaryonField[CH3INum][i] = TestProblemData.CH3I_Fraction * 
        TestProblemData.Methine3FractionByMass * BaryonField[0][i];

      BaryonField[CH4INum][i] = TestProblemData.CH4I_Fraction * 
        TestProblemData.Methine4FractionByMass * BaryonField[0][i];

      BaryonField[CO2INum][i] = TestProblemData.CO2I_Fraction * 
        TestProblemData.CarbonDioxideFractionByMass * BaryonField[0][i];

    }

  } // for (i = 0; i < size; i++)

  return SUCCESS;
}
