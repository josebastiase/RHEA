#pragma once

#include "ComputeElasticityTensorBase.h"

/**
 * Defines an elastic tensor for the model.
 * The elastic tensor is assumed isotropic, so it is defined by its bulk modulus and shear modulus.
 * However, these moduli may vary throughout the model (they are MOOSE AuxVariables) so the resulting
 * elastic tensor is heterogeneous - different at each point in the model.
 * The bulk modulus and shear modulus may be defined using RHEA's python workflow
 */

class IsotropicElasticModulusFromVar : public ComputeElasticityTensorBase
{
public:
  IsotropicElasticModulusFromVar(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void computeQpElasticityTensor() override;

  /// Elastic constants
  const VariableValue & _bulk_modulus;
  const VariableValue & _shear_modulus;

  /// Individual elasticity tensor
  RankFourTensor _Cijkl;
};
