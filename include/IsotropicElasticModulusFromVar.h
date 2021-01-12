#pragma once

#include "ComputeElasticityTensorBase.h"

class IsotropicElasticModulusFromVar;

template <>
InputParameters validParams<IsotropicElasticModulusFromVar>();

/**
 * Compute isotropic elastic tensor defines an elasticity tensor material for
 * isotropic materials from data file. This may be coupled with the python workflow to
 * allocate data at qp of the mesh. From bulk modulus and shear modulus.
 */

class IsotropicElasticModulusFromVar : public ComputeElasticityTensorBase
{
public:
  IsotropicElasticModulusFromVar(const InputParameters & parameters);

protected:
  virtual void computeQpElasticityTensor() override;

  /// Elastic constants
  const VariableValue & _bulk_modulus;
  const VariableValue & _shear_modulus;

  /// Individual elasticity tensor
  RankFourTensor _Cijkl;
};
