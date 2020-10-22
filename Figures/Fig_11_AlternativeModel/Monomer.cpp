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
#include "random.hpp"

using namespace std; 

/* Set state variables of monomer. */
void Monomer::initialize_state( int pindex, bool pT, bool pS, bool pnI, bool pnII )
{
  this->index = pindex;

  this->T = pT;
  this->S = pS;
  this->nI = pnI;
  this->nII = pnII;
}

/* Propagate monomer */
bool Monomer::propagate()
{
  double rnd( ran_get() );
  int reaction_idx( this->choose_reaction(rnd) );
    
  this->fire_reaction( reaction_idx );
  
  //update_prop( reaction_idx );
  set_propensities();
  
  return 1;
}

/* Give random number rnd [0,1), function returns index of fired reaction.
   Index to reaction mapping is defined in calc_propensities. */
int Monomer::choose_reaction(double rnd)
{
  double qacc(this->prop_list[0]);
  int i(0);
  
  rnd *= SCALE;
  rnd *= mon_prop_cont->get_qtot_el(this->index);
  while( rnd > qacc )
  {
    qacc += this->prop_list[++i];
  }
  
  return i;
}

/* Calculate int and ext (that depend on hexamer state) propensities, 
   and update prop container. */
void Monomer::update_mon_prop_cont()
{
  double sum(0.0);
  int i(4);
    
  while(i < MONOMER_N_REACTS)
  {
    sum += this->prop_list[i++];
  }

  mon_prop_cont->update_qint_qext_el(sum, 
    prop_list[0] + prop_list[1] + prop_list[2] + prop_list[3], index);
}


/* Given fired reaction channel index, execute state change in monomer. */
void Monomer::fire_reaction(int reaction_channel)
{
  if( sys->register_flux && reaction_channel > 3 )
  {
    sys->hist_pTVSpS[hexamer->get_hex_T()][hexamer->get_hex_S()] += 
      sys->tsim - hexamer->tlastswitch;
    hexamer->tlastswitch = sys->tsim;
    
    switch (reaction_channel)
    {
      case 4:  sys->rate_pTp1[hexamer->get_hex_T()][hexamer->get_hex_S()]++;
      break;
      case 5:  sys->rate_pTm1[hexamer->get_hex_T()][hexamer->get_hex_S()]++;
      break;
      case 6:  sys->rate_pSp1[hexamer->get_hex_T()][hexamer->get_hex_S()]++;
      break;
      case 7:  sys->rate_pSm1[hexamer->get_hex_T()][hexamer->get_hex_S()]++;
      break;
      case 8:  sys->rate_pTp1[hexamer->get_hex_T()][hexamer->get_hex_S()]++;
      break;
      case 9:  sys->rate_pTm1[hexamer->get_hex_T()][hexamer->get_hex_S()]++;
      break;
      case 10: sys->rate_pSp1[hexamer->get_hex_T()][hexamer->get_hex_S()]++;
      break;
      case 11: sys->rate_pSm1[hexamer->get_hex_T()][hexamer->get_hex_S()]++;
      break;     
    }
  }

  switch (reaction_channel)
  {
    case 0:  nI = 0; sys->CIATPcons++; //CI-ATP -> CI-ADP + Pi
    break;
    case 1:  nI = 1; //CI-ADP + ATP -> CI-ATP + ADP
    break;
    case 2:  nII = 0; //sys->CIIATPcons++; //CII-ATP -> CII-ADP + Pi
    break;
    case 3:  nII = 1; //CII-ADP + ATP -> CII-ATP + ADP
    break;
    case 4:  nII = 0; T = 1; S = 0; sys->CIIATPcons++; // U <-> T
    break;
    case 5:  nII = 1; T = 0; S = 0; sys->CIIATPcons--; 
    break;
    case 6:  nII = 0; T = 1; S = 1; sys->CIIATPcons++; // T <-> D
    break;
    case 7:  nII = 1; T = 1; S = 0; sys->CIIATPcons--;    
    break;
    case 8:  nII = 0; T = 1; S = 1; sys->CIIATPcons++; // S <-> D
    break;
    case 9:  nII = 1; T = 0; S = 1; sys->CIIATPcons--;
    break;
    case 10: nII = 0; T = 0; S = 1; sys->CIIATPcons++; // U <-> S
    break;
    case 11: nII = 1; T = 0; S = 0; sys->CIIATPcons--;             
    break;  
  }
}

/* CII nucleotide exchange rate, depends on KaiA bound state */
double Monomer::kCIInucloff()
{
  double kCIInucloff0( reaction_consts->kCIInucloff0 );
  double kCIInucloffA( reaction_consts->kCIInucloffA );
  double CIIAbound( hexamer->get_CIIKaiA_bound() );
  
  double retval(0.);
  
  if( hexamer->get_active() )
    retval = (double) kCIInucloffA * CIIAbound;
  else
    retval = kCIInucloff0;
  
  return retval;
}




/* Calculate all reaction propensities possible in this monomer.
   Monomer has NREACT reaction channels, which propensities are
   stored in prop_list.
*/
void Monomer::set_propensities()
{ 
  //System variables
  double ATPfrac = sys->ATPfrac;
  
  //Nucleotide related
  double kCIhyd    ( reaction_consts->kCIhyd );
  double kCIoffATP ( reaction_consts->kCIATPoff );
  double kCIoffADP ( hexamer->kCIADPoff() );
  
  double kCIIhyd ( reaction_consts->kCIIhyd0 );
  double kCIIoffnucl ( this->kCIInucloff() ); 
  double KATPKADP ( reaction_consts->KATPoKADP );
  
  /* Relative Nucleotide Affinity depends of conformation */
  if( 1 - hexamer->get_active() )
    KATPKADP = 1.0/KATPKADP;
          
  //Bare phosphotransfer rates 
  double kUT( reaction_consts->kUT );
  double kTU( reaction_consts->kTU );
  double kTD( reaction_consts->kTD );
  double kDT( reaction_consts->kDT );
  double kSD( reaction_consts->kSD );
  double kDS( reaction_consts->kDS );
  double kUS( reaction_consts->kUS );
  double kSU( reaction_consts->kSU );
  
  /* Phosphotransfer is modified by CII.KaiA differential affinity. */
  if(hexamer->get_CIIKaiA_bound())
  {  
    double dgUTbind( reaction_consts->dgACIIT - reaction_consts->dgACIIU );
    double dgTDbind( reaction_consts->dgACIID - reaction_consts->dgACIIT );
    double dgUSbind( reaction_consts->dgACIIS - reaction_consts->dgACIIU );
    double dgSDbind( reaction_consts->dgACIID - reaction_consts->dgACIIS );  
    
    kUT *= exp(-.5 * dgUTbind);
    kTU *= exp( .5 * dgUTbind);
    kTD *= exp(-.5 * dgTDbind);
    kDT *= exp( .5 * dgTDbind);
    kSD *= exp(-.5 * dgSDbind);
    kDS *= exp( .5 * dgSDbind);
    kUS *= exp(-.5 * dgUSbind);
    kSU *= exp( .5 * dgUSbind);
    
    kCIIhyd = reaction_consts->kCIIhydA;
  } 
    
  /* Calculate reactionchannel propensities */  
  // CI:  T <-> D (Reactions depend on hexamer state)
  prop_list[0]  = (double) nI * kCIhyd;
  prop_list[1]  = (double) (1-nI) * kCIoffADP;

  // CII: T <-> D (Reactions depend on hexamer state, affinity is set by OFF-rates)
  prop_list[2]  = (double) nII * (kCIIhyd + (1. - ATPfrac) * kCIIoffnucl);
  prop_list[3]  = (double) (1-nII) * ATPfrac * kCIIoffnucl / KATPKADP;
  
  // CII: T <-> D (Reactions depend on hexamer state, affinity is set by ON-rates)
  // KATPKADP = kDon/kTon (knucloff = 6.0, KATPKADP=0.25
  //double temp( (1.0 - ATPfrac)/ATPfrac * KATPKADP );
  //prop_list[2]  = (double) nII * (kCIIhyd + kCIIoffnucl * temp / ( 1.0 + temp ) );
  //prop_list[3]  = (double) (1-nII) * kCIIoffnucl / ( 1.0 + temp );
  
  // CII: U <-> T
  prop_list[4]  = (double) (1-T) * (1-S) * nII * kUT;
  prop_list[5]  = (double) T * (1-S) * (1-nII) * kTU;
  // CII: T <-> D
  prop_list[6]  = (double) T * (1-S) * nII * kTD;
  prop_list[7]  = (double) T * S * (1-nII) * kDT;
  // CII: S <-> D
  prop_list[8]  = (double) (1-T) * S * nII * kSD;
  prop_list[9]  = (double) T * S * (1-nII) * kDS;
  // CII: U <-> S    
  prop_list[10] = (double) (1-T) * (1-S) * nII * kUS;
  prop_list[11] = (double) (1-T) * S * (1-nII) * kSU;
     
  /* Calculate total, external and internal reaction propensities */
  update_mon_prop_cont();
}


/* Recalculate reactions that change due to state changes in the hexamer. */
double Monomer::update_ext_prop()
{
  set_propensities();

  return 0;
}

//Print reaction propensities of this monomer.
void Monomer::print_mon_props()
{
  cout << endl;
  cout << "### Propensities of monomer idx " << hexamer->get_index() * 6 + index << endl;
  cout << "State - nI: " << nI << ", nII: " << nII << ", T: " << T 
       << ", S: " << S << endl;
  cout << "CI  : T->D    (0)  " << prop_list[0] << endl;
  cout << "CI  : D->T    (1)  " << prop_list[1] << endl;
  cout << "CII : T->D    (2)  " << prop_list[2] << endl;
  cout << "CII : D->T    (3)  " << prop_list[3] << endl;
  cout << "CII : pU->pT  (4)  " << prop_list[4] << endl; 
  cout << "CII : pT->pU  (5)  " << prop_list[5] << endl;  
  cout << "CII : pT->pD  (6)  " << prop_list[6] << endl;
  cout << "CII : pD->pT  (7)  " << prop_list[7] << endl;  
  cout << "CII : pS->pD  (8)  " << prop_list[8] << endl;
  cout << "CII : pD->pS  (9)  " << prop_list[9] << endl;
  cout << "CII : pU->pS  (10) " << prop_list[10] << endl;
  cout << "CII : pS->pU  (11) " << prop_list[11] << endl;
  cout << "Props - tot: " << mon_prop_cont->get_qtot_el(index) 
       << ", int: "       << mon_prop_cont->get_qint_el(index) 
       << ", ext: "       << mon_prop_cont->get_qext_el(index) << endl;
}

