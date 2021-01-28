#include "IsotropicElasticModulusFromVar.h"

registerMooseObject("TensorMechanicsApp", IsotropicElasticModulusFromVar);

template <>
InputParameters
validParams<IsotropicElasticModulusFromVar>()
{
  InputParameters params = validParams<ComputeElasticityTensorBase>();
  params.addClassDescription("Compute a isotropic, heterogeneous elasticity tensor using the "
			     "spatially-varying bulk_modulus and shear_modulus values.");
  params.addRequiredCoupledVar("bulk_modulus", "The bulk modulus for the material, which can "
			       "vary throughout the mesh using RHEA's python workflow");
  params.addRequiredCoupledVar("shear_modulus", "The shear modulus of the material, which can "
			       "vary throughout the mesh using RHEA's python workflow");
  return params;
}

IsotropicElasticModulusFromVar::IsotropicElasticModulusFromVar(const InputParameters & parameters)
  : ComputeElasticityTensorBase(parameters),
    _bulk_modulus(coupledValue("bulk_modulus")),
    _shear_modulus(coupledValue("shear_modulus"))
{
   issueGuarantee(_elasticity_tensor_name, Guarantee::ISOTROPIC);
}

void
IsotropicElasticModulusFromVar::computeQpElasticityTensor()
{
 // Fill elasticity tensor

   const Real E = (9 * _bulk_modulus[_qp] * _shear_modulus[_qp]) / (3 * _bulk_modulus[_qp] + _shear_modulus[_qp]);
   const Real nu = (3 * _bulk_modulus[_qp] - 2 * _shear_modulus[_qp]) / (2 * (3 * _bulk_modulus[_qp] + _shear_modulus[_qp]));

  _Cijkl.fillSymmetricIsotropicEandNu(E, nu);
  _elasticity_tensor[_qp] = _Cijkl;
}
