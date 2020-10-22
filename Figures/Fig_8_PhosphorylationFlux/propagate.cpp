/* propagates the simulator one event,
   Returns TRUE if simulation has not reached required end time yet.
*/

#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>
#include <math.h>

#include "main.hpp"
#include "Monomer.hpp"
#include "Hexamer.hpp"
#include "PropensityContainer.hpp"
#include "propagate.hpp"
#include "random.hpp"

using namespace std;

bool propagate(SystemVariables *sys, Hexamer *hexamers, PropensityContainer *prop_cont)
{
  //Draw two random numbers.
  double rnd1( ran_get() ), rnd2( ran_get() );

  ///update tsim, choose monomer and propagate it.
  sys->tsim -= log(rnd1) / prop_cont->get_qtot();

  /* For time dependent flux registration.   */
  sys->register_flux=1;
  
  int fired_hex_idx( prop_cont->choose_index(rnd2) ); 

  /* When KaiA binds to hexamer CI or CII domain,
     propensities of all hexamers change.
  */
  if( hexamers[fired_hex_idx].propagate() )
  {   
    //Change Afree dependent reactions in all hexamers.
    prop_cont->set_qext_all(hexamers);
  }
 
  sys->step_cntr++;
  return sys->tsim < sys->tend;
}
