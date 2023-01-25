//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

// Application includes
#include "custom_constitutive/joint_bilinear_cohesive_2D_law.hpp"

namespace Kratos
{

void JointBilinearCohesive2DLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
	rFeatures.mOptions.Set( PLANE_STRAIN_LAW );
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
	//rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the spacedimension
	rFeatures.mSpaceDimension = 2;

	//Set the strain size
	rFeatures.mStrainSize = 2;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void JointBilinearCohesive2DLaw::ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables,
                                                         Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
    {
        rVariables.EquivalentStrain = std::sqrt(StrainVector[0]*StrainVector[0]+
                                                StrainVector[1]*StrainVector[1])/rVariables.CriticalDisplacement;
	}
	else // Contact between interfaces
	{
        rVariables.EquivalentStrain = std::abs(StrainVector[0])/rVariables.CriticalDisplacement;
	}
}

//----------------------------------------------------------------------------------------

void JointBilinearCohesive2DLaw::ComputeConstitutiveMatrix(Matrix& rConstitutiveMatrix,
                                                           ConstitutiveLawVariables& rVariables,
                                                           Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
    {
        if(rVariables.LoadingFlag) // Loading
        {
            rConstitutiveMatrix(0,0) = rVariables.YieldStress/((1.0-rVariables.DamageThreshold)*rVariables.CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                        StrainVector[0]*StrainVector[0]/(rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
            rConstitutiveMatrix(1,1) = rVariables.YieldStress/((1.0-rVariables.DamageThreshold)*rVariables.CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                        StrainVector[1]*StrainVector[1]/(rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );

            rConstitutiveMatrix(0,1) = -rVariables.YieldStress*StrainVector[0]*StrainVector[1]/( (1.0-rVariables.DamageThreshold)*
                                        rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable );
            rConstitutiveMatrix(1,0) = rConstitutiveMatrix(0,1);
        }
        else  // Unloading
        {
            rConstitutiveMatrix(0,0) = rVariables.YieldStress/(rVariables.CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-rVariables.DamageThreshold);
            rConstitutiveMatrix(1,1) = rConstitutiveMatrix(0,0);

            rConstitutiveMatrix(0,1) = 0.0;
            rConstitutiveMatrix(1,0) = 0.0;
        }
    }
    else // Contact between interfaces
    {
        if(rVariables.LoadingFlag) // Loading
        {
            rConstitutiveMatrix(0,0) = rVariables.YieldStress/((1.0-rVariables.DamageThreshold)*rVariables.CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                        StrainVector[0]*StrainVector[0]/(rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
            rConstitutiveMatrix(1,1) = rVariables.YoungModulus/(rVariables.DamageThreshold*rVariables.CriticalDisplacement);

            if(StrainVector[0] > 1.0e-20)
            {
                rConstitutiveMatrix(0,1) = -rVariables.YieldStress*StrainVector[0]*StrainVector[1]/( (1.0-rVariables.DamageThreshold)*
                                            rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable );
            }
            else if(StrainVector[0] < -1.0e-20)
            {
                rConstitutiveMatrix(0,1) = -rVariables.YieldStress*StrainVector[0]*StrainVector[1]/( (1.0-rVariables.DamageThreshold)*
                                            rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*rVariables.CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable );
            }
            else
            {
                rConstitutiveMatrix(0,1) = 0.0;
            }

            rConstitutiveMatrix(1,0) = 0.0;
        }
        else  // Unloading
        {
            rConstitutiveMatrix(0,0) = rVariables.YieldStress/(rVariables.CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-rVariables.DamageThreshold);
            rConstitutiveMatrix(1,1) = rVariables.YoungModulus/(rVariables.DamageThreshold*rVariables.CriticalDisplacement);
            rConstitutiveMatrix(0,1) = 0.0;
            rConstitutiveMatrix(1,0) = 0.0;
        }
    }
}

//----------------------------------------------------------------------------------------

void JointBilinearCohesive2DLaw::ComputeStressVector(Vector& rStressVector,
                                                     ConstitutiveLawVariables& rVariables,
                                                     Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
    {
        rStressVector[0] = rVariables.YieldStress/(rVariables.CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-rVariables.DamageThreshold) * StrainVector[0];
        rStressVector[1] = rVariables.YieldStress/(rVariables.CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-rVariables.DamageThreshold) * StrainVector[1];
    }
    else // Contact between interfaces
    {
        // Note: StrainVector[1] < 0.0
        rStressVector[1] = rVariables.YoungModulus/(rVariables.DamageThreshold*rVariables.CriticalDisplacement)*StrainVector[1];

        if(StrainVector[0] > 0.0)
        {
            rStressVector[0] = rVariables.YieldStress/(rVariables.CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-rVariables.DamageThreshold)*StrainVector[0];
        }
        else if(StrainVector[0] < 0.0)
        {
            rStressVector[0] = rVariables.YieldStress/(rVariables.CriticalDisplacement*mStateVariable)*(1.0-mStateVariable)/(1.0-rVariables.DamageThreshold)*StrainVector[0];
        }
        else
        {
            rStressVector[0] = 0.0;
        }
    }
    // Add Uplift Pressure
    rStressVector[1] += mUpliftPressure;
}

} // Namespace Kratos
