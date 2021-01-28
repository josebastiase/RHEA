#pragma once

#include "PorousFlowMaterialVectorBase.h"

/**
 * Material designed to provide a time-invariant
 * Biot Modulus, M, where
 * 1 / M = (1 - alpha) * (alpha - phi) / K + phi / Kf .
 * Here
 * alpha = Biot coefficient (assumed time and space independent)
 * phi = initial value of porosity (can vary throughout the mesh)
 * K = drained bulk modulus of the porous skeleton (can vary throughout the mesh)
 * Kf = fluid bulk modulus (assumed time and space independent)
 *
 * Note that RHEA allows the result to be heterogeneous, ie, different at each point in the mesh.
 * This is because the drained bulk modulus and the porosity may be defined using RHEA's python
 * workflow.
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

  /// Drained bulk modulus of the porous skeleton
  const VariableValue & _bulk_modulus;

  /// porosity at the nodes or quadpoints.  Only the initial value is ever used
  const MaterialProperty<Real> & _porosity;

  /// Computed Biot modulus
  MaterialProperty<Real> & _biot_modulus;

  /// Old value of Biot modulus.  This variable is necessary in order to keep Biot modulus constant even if porosity is changing.
  const MaterialProperty<Real> & _biot_modulus_old;
};
