#include "BiotModulusFromVar.h"

registerMooseObject("PorousFlowApp", BiotModulusFromVar);

template <>
InputParameters
validParams<BiotModulusFromVar>()
{
  InputParameters params = validParams<PorousFlowMaterialVectorBase>();
  params.addRangeCheckedParam<Real>(
      "biot_coefficient", 1.0, "biot_coefficient>=0 & biot_coefficient<=1", "Biot coefficient");
  params.addRangeCheckedParam<Real>(
      "fluid_bulk_modulus", 2.0E9, "fluid_bulk_modulus>0", "Fluid bulk modulus");
  params.addRequiredCoupledVar("bulk_modulus", "The bulk modulus for the material.");
  params.addPrivateParam<std::string>("pf_material_type", "biot_modulus");
  params.addClassDescription("Computes the Biot Modulus, which is assumed to be constant for all "
                             "time.  Sometimes 1 / BiotModulus is called storativity");
  return params;
}

BiotModulusFromVar::BiotModulusFromVar(const InputParameters & parameters)
  : PorousFlowMaterialVectorBase(parameters),
    _biot_coefficient(getParam<Real>("biot_coefficient")),
    _fluid_bulk_modulus(getParam<Real>("fluid_bulk_modulus")),
    _bulk_modulus(coupledValue("bulk_modulus")),
    _porosity(_nodal_material ? getMaterialProperty<Real>("PorousFlow_porosity_nodal")
                              : getMaterialProperty<Real>("PorousFlow_porosity_qp")),
    _biot_modulus(_nodal_material ? declareProperty<Real>("PorousFlow_constant_biot_modulus_nodal")
                                  : declareProperty<Real>("PorousFlow_constant_biot_modulus_qp")),
    _biot_modulus_old(_nodal_material
                          ? getMaterialPropertyOld<Real>("PorousFlow_constant_biot_modulus_nodal")
                          : getMaterialPropertyOld<Real>("PorousFlow_constant_biot_modulus_qp"))
{
}

void
BiotModulusFromVar::initQpStatefulProperties()
{
  _biot_modulus[_qp] = 1.0 / ((1.0 - _biot_coefficient) * (_biot_coefficient - _porosity[_qp]) /
                                  _bulk_modulus[_qp] +
                              _porosity[_qp] / _fluid_bulk_modulus);
}

void
BiotModulusFromVar::computeQpProperties()
{
  _biot_modulus[_qp] = _biot_modulus_old[_qp];
}
