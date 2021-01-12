//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "RHEATestApp.h"
#include "RHEAApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
RHEATestApp::validParams()
{
  InputParameters params = RHEAApp::validParams();
  return params;
}

RHEATestApp::RHEATestApp(InputParameters parameters) : MooseApp(parameters)
{
  RHEATestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

RHEATestApp::~RHEATestApp() {}

void
RHEATestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  RHEAApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"RHEATestApp"});
    Registry::registerActionsTo(af, {"RHEATestApp"});
  }
}

void
RHEATestApp::registerApps()
{
  registerApp(RHEAApp);
  registerApp(RHEATestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
RHEATestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  RHEATestApp::registerAll(f, af, s);
}
extern "C" void
RHEATestApp__registerApps()
{
  RHEATestApp::registerApps();
}
