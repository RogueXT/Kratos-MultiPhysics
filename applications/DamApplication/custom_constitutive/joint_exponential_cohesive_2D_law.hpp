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

#if !defined (KRATOS_JOINT_EXPONENTIAL_COHESIVE_2D_LAW_H_INCLUDED)
#define  KRATOS_JOINT_EXPONENTIAL_COHESIVE_2D_LAW_H_INCLUDED

// System includes

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/joint_exponential_cohesive_3D_law.hpp"
#include "dam_application_variables.h"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(DAM_APPLICATION) JointExponentialCohesive2DLaw : public JointExponentialCohesive3DLaw
{

public:

    /// Definition of the base class
    typedef JointExponentialCohesive3DLaw BaseType;

    KRATOS_CLASS_POINTER_DEFINITION(JointExponentialCohesive2DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    JointExponentialCohesive2DLaw()
    {
    }

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<JointExponentialCohesive2DLaw>(JointExponentialCohesive2DLaw(*this));
    }

    // Copy Constructor
    JointExponentialCohesive2DLaw (const JointExponentialCohesive2DLaw& rOther) : JointExponentialCohesive3DLaw(rOther)
    {
    }

    // Destructor
    ~JointExponentialCohesive2DLaw() override
    {
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeConstitutiveLawVariables(ConstitutiveLawVariables& rVariables, Parameters& rValues) override;

    void ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables, Parameters& rValues) override;

    void ComputeCriticalDisplacement(ConstitutiveLawVariables& rVariables, Parameters& rValues) override;

    void ComputeConstitutiveMatrix(Matrix& rConstitutiveMatrix,
                                    ConstitutiveLawVariables& rVariables,
                                    Parameters& rValues) override;

    void ComputeStressVector(Vector& rStressVector,
                                ConstitutiveLawVariables& rVariables,
                                Parameters& rValues) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

}; // Class JointExponentialCohesive2DLaw
}  // namespace Kratos.
#endif // KRATOS_JOINT_EXPONENTIAL_COHESIVE_2D_LAW_H_INCLUDED  defined
