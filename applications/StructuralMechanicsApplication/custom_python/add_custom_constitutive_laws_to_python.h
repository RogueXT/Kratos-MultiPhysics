// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Riccardo Rossi
//

#pragma once

// System includes
#include <pybind11/pybind11.h>
// External includes

// Project includes
#include "includes/define_python.h"

namespace Kratos
{

  namespace Python
  {

    void  AddCustomConstitutiveLawsToPython(pybind11::module& m);

  }  // namespace Python.

}  // namespace Kratos.
