# Terzaghi model using spatially-distribuited data for hydraulic conductivity, porosity, bulk modulus and shear modulus.
# It is assumed that the notebook Workflow_TerzaghiImportDataLayers.ipynb has created files containing the spatial distribution of the hydraulic conductivity, porosity, bulk modulus and shear modulus.
# In this case the hydraulic conductivity, porosity, bulk modulus and shear modulus correspond to a layered medium with spatially-varying hydraulic conductivity

# K: Hydraulic conductivity
# P: Porosity
# L: Bulk modulus
# G: Shear modulus

# Smaller time-step sizes yield better agreement with the analytical solution

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
    data_file = K.data
  []
  # Porosity
  [P_fcn]
    type = PiecewiseMulticonstant
    direction = 'left right'
    data_file = p.data
  []
  # Bulk modulus
  [L_fcn]
    type = PiecewiseMulticonstant
    direction = 'left right'
    data_file = L.data
  []
  # Shear modulus
  [G_fcn]
    type = PiecewiseMulticonstant
    direction = 'left right'
    data_file = G.data
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

[VectorPostprocessors]
  [pp]
    type = LineValueSampler
    sort_by = y
    variable = porepressure
    start_point = '0 100 0'
    end_point = '0 0 0'
    num_points = 21
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Postprocessors]
  [dt]
    type = FunctionValuePostprocessor
    outputs = console
    function = 'min(1E6, max(1E4, 0.1 * t))' # note: smaller values yields better agreement with the analytical solution
  []
[]

[Executioner]
  type = Transient
  solve_type = Newton
  end_time = 7E6
  [./TimeStepper]
    type = PostprocessorDT
    postprocessor = dt
    dt = 1E4
  [../]
[]

[Outputs]
  [csv]
    type = CSV
    sync_times = '0 1E4 1E5 3E5 1E6 2E6 3E6 5E6 7E6'
    sync_only = true
  []
[]
