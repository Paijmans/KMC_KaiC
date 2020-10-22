/* 
  Propagate function propagates the system one reaction.
*/

#ifndef _PROPAGATE_H_
#define _PROPAGATE_H_

#include "main.hpp"
#include "Monomer.hpp"
#include "PropensityContainer.hpp"

bool propagate(SystemVariables *sys, Hexamer *hexamers, PropensityContainer *prop_cont);

#endif
