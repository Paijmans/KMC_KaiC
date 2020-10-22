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
#include <iomanip>

using namespace std;

bool propagate(SystemVariables *sys, Hexamer *hexamers, PropensityContainer *prop_cont)
{
  //Draw two random numbers.
  double rnd1( ran_get() ), rnd2( ran_get() );

  ///update tsim, choose monomer and propagate it.
  sys->tsim -= log(rnd1) / prop_cont->get_qtot();

  /* For time dependent flux registration. */
  if(sys->tsim > sys->treg[sys->index_reg+1])
  {
    //Reset CI/CII*KaiA bound times; print data for all hexamers.
    for(int i(0); i<sys->N_hexamers; i++)
    {
      if(hexamers[i].get_CIKaiA_bound()>0)
      {
        hexamers[i].tboundCI += sys->treg[sys->index_reg+1] - hexamers[i].tlastbindCI;
        hexamers[i].tlastbindCI =sys->treg[sys->index_reg+1];
      }
      if(hexamers[i].get_CIIKaiA_bound())
      {
        hexamers[i].tboundCII += sys->treg[sys->index_reg+1] - hexamers[i].tlastbindCII;
        hexamers[i].tlastbindCII = sys->treg[sys->index_reg+1];
      }
  
      cout << setprecision(10) << hexamers[i].tboundCI << "\t" << hexamers[i].tboundCII << endl;
    
      hexamers[i].tboundCI=0;
      hexamers[i].tboundCII=0;
    }
    cerr << "\n## KaiA bound times reset @" << sys->tsim << endl;      
    sys->index_reg++;      
  }  
    
  int fired_hex_idx( prop_cont->choose_index(rnd2) ); 

  /* When KaiA binds to hexamer CI or CII domain,
     propensities of all hexamers change. */
  if( hexamers[fired_hex_idx].propagate() )
  {   
    //Change Afree dependent reactions in all hexamers.
    prop_cont->set_qext_all(hexamers);
  }
 
  sys->step_cntr++;
  return sys->tsim < sys->tend;
}
