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

  sys->register_flux=1;
  /* For time dependent flux registration.   
  if(!sys->register_flux) 
  { 
    if(sys->tsim >= sys->treg[sys->index_reg])
    {
      sys->register_flux = 1;
      for(int i(0); i<sys->N_hexamers; i++)
      {
        hexamers[i].tlastswitch = sys->treg[sys->index_reg];
      }     
      cerr << "\n## Temporal flux integration #" << sys->index_reg 
           << ", ENGAGED@t=" << sys->tsim << endl;
    }
  } 
  else
  {
    if(sys->tsim > sys->treg[sys->index_reg] + sys->interval_reg)
    {
      for(int i(0); i<sys->N_hexamers; i++)
      {
        hexamers[i].get_Monomer(0)->fire_reaction(12);
      }
      cerr << "\n## Temporal flux integration DISENGAGED@t=" << sys->tsim << endl;      
      sys->register_flux = 0;
      sys->index_reg++;      
    }  
  }
  */

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
