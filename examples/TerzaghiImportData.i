#Model to try import spatially distribuited data to MOOSE. Imports .data files for hydraulic conductivity, porosity, bulk modulus and shear modulus.

# K: Hydraulic conductivity
# P: Porosity
# L: Bulk modulus
# G: Shear modulus

# Generate quadratic mesh
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 99 # Discretization -1
  xmin = -1
  xmax = 1
  ymin = 0
  ymax = 100
[]

# The easiest way to use PorousFlow, is using actions/dictator
[GlobalParams]
 displacements = 'disp_x disp_y'
 PorousFlowDictator = dictator
 multiply_by_density = false
 biot_coefficient = 0.9
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [porepressure]
    scaling = 1E11
  []
[]

# Import data
[Functions]
  # Hydraulic conductivity
  [K_fcn]
    type = PiecewiseMulticonstant
    direction = 'left right'
    data_file = Data/K.data
  []
  # Porosity
  [P_fcn]
    type = PiecewiseMulticonstant
    direction = 'left right'
    data_file = Data/p.data
  []
  # Bulk modulus
  [L_fcn]
    type = PiecewiseMulticonstant
    direction = 'left right'
    data_file = Data/L.data
  []
  # Shear modulus
  [G_fcn]
    type = PiecewiseMulticonstant
    direction = 'left right'
    data_file = Data/G.data
  []
[]

# Storage the imported data in a variable
[AuxVariables]
  [K]
    family = MONOMIAL
    order = CONSTANT
  []
  [p]
    family = MONOMIAL
    order = CONSTANT
  []
  [L]
    order = CONSTANT
    family = MONOMIAL
  []
  [G]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [K]
    type = FunctionAux
    function = K_fcn
    variable = K
    execute_on = initial
  []
  [p]
    type = FunctionAux
    function = P_fcn
    variable = p
    execute_on = initial
  [../]
  [L]
    type = FunctionAux
    function = L_fcn
    variable = L
    execute_on = initial
  []
  [G]
    type = FunctionAux
    function = G_fcn
    variable = G
    execute_on = initial
  []
[]

# Boundary conditions. System open only at the top
[BCs]
  [no_x_disp]
    type = DirichletBC
    preset = true
    variable = disp_x
    value = 0
    boundary = 'left right'
  []
  [no_y_disp]
    type = DirichletBC
    preset = true
    variable = disp_y
    value = 0
    boundary = 'bottom'
  []
  [top_drained]
    type = DirichletBC
    variable = porepressure
    value = 0
    boundary = 'top'
  []
  [./top_load]
    type = NeumannBC
    variable = disp_y
    value = -1000
    boundary = 'top'
  []
[]

[Modules]
  [FluidProperties]
    [the_simple_fluid]
      type = SimpleFluidProperties
      thermal_expansion = 0.0
      bulk_modulus = 2.2E9
      viscosity = 1E-3
      density0 = 1000.0
    []
  []
[]

[PorousFlowBasicTHM]
  coupling_type = HydroMechanical
  displacements = 'disp_x disp_y'
  porepressure = porepressure
  gravity = '0 0 0'
  fp = the_simple_fluid
[]

[Materials]
  [elasticity_tesnors]
    type = IsotropicElasticModulusFromVar
    bulk_modulus = L
    shear_modulus = G
  []
  [strain]
    type = ComputeSmallStrain
  []
  [stress]
    type = ComputeLinearElasticStress
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = p
  []
  [biot_modulus]
    type = BiotModulusFromVar
    bulk_modulus = L
  []
  [permeability]
    type = PorousFlowPermeabilityConstFromVar
    perm_xx = K
    perm_yy = K
    perm_zz = K
  []
[]


[Postprocessors]
  [p0]
    type = PointValue
    outputs = csv
    point = '0 0 0'
    variable = porepressure
    use_displaced_mesh = false
  []
  [p1]
    type = PointValue
    outputs = csv
    point = '0 10 0'
    variable = porepressure
    use_displaced_mesh = false
  []
  [./p2]
    type = PointValue
    outputs = csv
    point = '0 20 0'
    variable = porepressure
    use_displaced_mesh = false
  []
  [p3]
    type = PointValue
    outputs = csv
    point = '0 30 0'
    variable = porepressure
    use_displaced_mesh = false
  []
  [p4]
    type = PointValue
    outputs = csv
    point = '0 40 0'
    variable = porepressure
    use_displaced_mesh = false
  []
  [p5]
    type = PointValue
    outputs = csv
    point = '0 50 0'
    variable = porepressure
    use_displaced_mesh = false
  []
  [p6]
    type = PointValue
    outputs = csv
    point = '0 60 0'
    variable = porepressure
    use_displaced_mesh = false
  []
  [p7]
    type = PointValue
    outputs = csv
    point = '0 70 0'
    variable = porepressure
    use_displaced_mesh = false
  []
  [p8]
    type = PointValue
    outputs = csv
    point = '0 80 0'
    variable = porepressure
    use_displaced_mesh = false
  []
  [p9]
    type = PointValue
    outputs = csv
    point = '0 90 0'
    variable = porepressure
    use_displaced_mesh = false
  []
  [p99]
    type = PointValue
    outputs = csv
    point = '0 100 0'
    variable = porepressure
    use_displaced_mesh = false
  []
  [dt]
    type = FunctionValuePostprocessor
    outputs = console
    function = if(t<0,1,10)
  []
[]

[Preconditioning]
  [mumps_is_best_for_parallel_jobs]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  []
[]

[Executioner]
  type = Transient
  solve_type = Newton
  start_time = -1
  end_time = 1000
  [./TimeStepper]
    type = PostprocessorDT
    postprocessor = dt
    dt = 1
  [../]
[]

[Outputs]
  execute_on = 'timestep_end'
  file_base = gold/TerzaghisImportData
  [csv]
    type = CSV
  []
[]
