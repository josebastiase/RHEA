#pragma once

#include "PorousFlowMaterialVectorBase.h"

/**
 * Material designed to provide a time-invariant
 * Biot Modulus, M, where
 * 1 / M = (1 - alpha) * (alpha - phi) * C + phi / Kf .
 * Here
 * alpha = Biot coefficient (assumed constant)
 * phi = initial value of porosity
 * K = drained bulk modulus of the porous material from var
 * Kf = fluid bulk modulus (assumed constant)
 */
class BiotModulusFromVar : public PorousFlowMaterialVectorBase
{
public:
  static InputParameters validParams();

  BiotModulusFromVar(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// Biot coefficient
  const Real _biot_coefficient;

  /// Fluid bulk modulus
  const Real _fluid_bulk_modulus;

  /// Solid bulk modulus
  const VariableValue & _bulk_modulus;

  /// porosity at the nodes or quadpoints.  Only the initial value is ever used
  const MaterialProperty<Real> & _porosity;

  /// Computed Biot modulus
  MaterialProperty<Real> & _biot_modulus;

  /// Old value of Biot modulus.  This variable is necessary in order to keep Biot modulus constant even if porosity is changing.
  const MaterialProperty<Real> & _biot_modulus_old;
};
