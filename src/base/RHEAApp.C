#include "RHEAApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
RHEAApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy DirichletBC, that is, set DirichletBC default for preset = true
  params.set<bool>("use_legacy_dirichlet_bc") = false;

  // Can output sensible material properties on INITIAL
  params.set<bool>("use_legacy_material_output") = false;
 
  return params;
}

RHEAApp::RHEAApp(InputParameters parameters) : MooseApp(parameters)
{
  RHEAApp::registerAll(_factory, _action_factory, _syntax);
}

RHEAApp::~RHEAApp() {}

void
RHEAApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"RHEAApp"});
  Registry::registerActionsTo(af, {"RHEAApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
RHEAApp::registerApps()
{
  registerApp(RHEAApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
RHEAApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  RHEAApp::registerAll(f, af, s);
}
extern "C" void
RHEAApp__registerApps()
{
  RHEAApp::registerApps();
}
