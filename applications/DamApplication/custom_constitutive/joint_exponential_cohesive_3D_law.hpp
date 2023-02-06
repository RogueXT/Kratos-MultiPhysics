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

#if !defined (KRATOS_JOINT_EXPONENTIAL_COHESIVE_3D_LAW_H_INCLUDED)
#define  KRATOS_JOINT_EXPONENTIAL_COHESIVE_3D_LAW_H_INCLUDED

// System includes

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/joint_bilinear_cohesive_3D_law.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(DAM_APPLICATION) JointExponentialCohesive3DLaw : public JointBilinearCohesive3DLaw
{

public:

    /// Definition of the base class
    typedef JointBilinearCohesive3DLaw BaseType;

    KRATOS_CLASS_POINTER_DEFINITION(JointExponentialCohesive3DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    JointExponentialCohesive3DLaw()
    {
    }

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<JointExponentialCohesive3DLaw>(JointExponentialCohesive3DLaw(*this));
    }

    // Copy Constructor
    JointExponentialCohesive3DLaw (const JointExponentialCohesive3DLaw& rOther) : JointBilinearCohesive3DLaw(rOther)
    {
    }

    // Destructor
    ~JointExponentialCohesive3DLaw() override
    {
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const override;

    void InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues ) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void FinalizeMaterialResponseCauchy (Parameters & rValues) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double& GetValue( const Variable<double>& rThisVariable, double& rValue ) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables
    double mDamageVariable;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeConstitutiveLawVariables(ConstitutiveLawVariables& rVariables, Parameters& rValues) override;

    void ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables, Parameters& rValues) override;

    void ComputeDamageVariable(ConstitutiveLawVariables& rVariables, Parameters& rValues);

    virtual void ComputeCriticalDisplacement(ConstitutiveLawVariables& rVariables, Parameters& rValues);

    void ComputeConstitutiveMatrix(Matrix& rConstitutiveMatrix,
                                    ConstitutiveLawVariables& rVariables,
                                    Parameters& rValues) override;

    void ComputeStressVector(Vector& rStressVector,
                                ConstitutiveLawVariables& rVariables,
                                Parameters& rValues) override;

    double MacaulayBrackets(const double& Value);
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

}; // Class JointExponentialCohesive3DLaw
}  // namespace Kratos.
#endif // KRATOS_JOINT_EXPONENTIAL_COHESIVE_3D_LAW_H_INCLUDED  defined
