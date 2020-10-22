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

/* Set state variables of hexamer */
void Hexamer::initialize_state( int pindex, bool pactive, 
  bool pCIIKaiA_bound, int pCIKaiA_bound, int pCIKaiB_bound )
{
  this->index = pindex;
  this->active = pactive;
  this->CIIKaiA_bound = pCIIKaiA_bound;
  this->CIKaiA_bound = pCIKaiA_bound;
  this->CIKaiB_bound = pCIKaiB_bound;  
  this->tlastbindCI = 0;
  this->tlastbindCII = 0;  
  this->tboundCI = 0;
  this->tboundCII = 0; 
}

/* Set which monomers are in the hexamer, and vice versa. */
void Hexamer::set_sextet(Monomer *monomers)
{ 
  for(int i(6*index); i < (index + 1) * 6; i++)
  {
    sextet[i-6*index] = &monomers[i];
    monomers[i].set_hexamer( this );
    monomers[i].set_prop_container( &(this->mon_prop_cont) );
    monomers[i].set_reaction_consts(this->reaction_consts);
    monomers[i].set_sysvars(this->sys);
    
    if(this->sys->start_phosphorylated)
    {   
      double rnd( ran_get() );
      bool S(0),T(0);      
      
      //double Tfrac(0.12), Sfrac(0.13), Dfrac(0.08);
      //double Tfrac(0.06), Sfrac(0.25), Dfrac(0.10);      
      //double Tfrac(0.07), Sfrac(0.33), Dfrac(0.18);
      double Tfrac(0.10), Sfrac(0.15), Dfrac(0.40);
      //double Tfrac(0.00), Sfrac(0.0), Dfrac(1.0);
      
      
      if(rnd < Dfrac) 
      {
        T=1; S=1; //D-state
      }    
      if(rnd > Dfrac && rnd < Dfrac + Tfrac) 
      {
        T=1; S=0; //T-state
      }
      if(rnd > Dfrac + Tfrac && rnd < Dfrac + Tfrac + Sfrac) 
      {
        T=0; S=1; //S-state
      }
           
      monomers[i].initialize_state(i - 6*index,T,S,1,ran_bool(sys->ATPfrac));
    }
    else
    {
      monomers[i].initialize_state(i - 6*index,0,0,1,ran_bool(sys->ATPfrac));
    }
  }
  
  /* Monomer propensities depend on hexamer state,
     so they can only be calculated after all monomer states are set.
  */
  for(int i(0); i < 6; i++)
  {
    sextet[i]->set_propensities();
    //sextet[i]->print_mon_props();
  }
}

/* Propagate hexamer, returns TRUE if reaction is fired
   that induces a change in external variables */
bool Hexamer::propagate()
{
  double rnd( ran_get() );

  int reaction_idx( this->choose_reaction(rnd) );   
  bool ext_change( this->fire_reaction( reaction_idx ) );
   
  /* Update propensities due to state change */
  set_propensities();
   
  //Reactions 0,1,3,4,5 change Afree.   
  return ext_change;
}

/* Give random number rnd [0,1), function returns index of fired reaction.
   Index to reaction mapping is defined in calc_propensities.
*/
int Hexamer::choose_reaction(double rnd)
{
  double qacc(this->prop_list[0]);
  int i(0);
  
  rnd*=SCALE;
  rnd*=hex_prop_cont->get_qtot_el(this->index);
  while( rnd > qacc )
  {
    qacc += this->prop_list[++i];
  }
    
  return i;
}

/* Given fired reaction channel index, execute state change in hexamer or monomer.
   Returns TRUE if a change in a monomers alters the hexamer propensities. */
/* Given fired reaction channel index, execute state change in hexamer or monomer.
   Returns TRUE if a change in a monomers alters the hexamer propensities. */
bool Hexamer::fire_reaction(int reaction_channel)
{
  bool ext_change(0), mon_ext_change(0);

  switch (reaction_channel)
  { 
    case 0: CIIKaiA_bound = 1; 
    sys->Afree--; 
    sys->cAfree = (double) sys->Afree / sys->volume;
    ext_change = 1;
    
    tlastbindCII = sys->tsim;
    break;
    case 1: CIIKaiA_bound = 0; 
    sys->Afree++;
    sys->cAfree = (double) sys->Afree / sys->volume;   
    ext_change = 1; 
       
    tboundCII += sys->tsim - this->tlastbindCII;
    break;
    case 2: CIKaiB_bound++;
    break;
    case 3: CIKaiB_bound--;
    break;
    case 4: CIKaiA_bound++; 
    sys->Afree--;
    sys->cAfree = (double) sys->Afree / sys->volume;
    ext_change = 1;
      
    if( CIKaiA_bound == 1 )
      tlastbindCI = sys->tsim;  
    break;
    case 5: CIKaiA_bound--; 
    sys->Afree++;
    sys->cAfree = (double) sys->Afree / sys->volume; 
    ext_change = 1;  
    
    if( CIKaiA_bound == 0 )
      tboundCI += sys->tsim - this->tlastbindCI;     
    break;
    case 6: active = 0; 
    break;
    case 7: active = 1;
    break;
    case 8: mon_ext_change = sextet[0]->propagate();
    break;
    case 9: mon_ext_change = sextet[1]->propagate();
    break;
    case 10: mon_ext_change = sextet[2]->propagate();
    break;
    case 11: mon_ext_change = sextet[3]->propagate();
    break;
    case 12: mon_ext_change = sextet[4]->propagate();
    break;
    case 13: mon_ext_change = sextet[5]->propagate();
    break;
  }
  
  return ext_change;
}

/* Add propensities of all reaction channels, and return total. */
void Hexamer::update_prop_cont()
{
  double hex_qint( prop_list[1] + prop_list[2] + prop_list[3]);
  double hex_qext( prop_list[0] + prop_list[4]);
  
  for(int i(5); i<14; i++)
  {
    hex_qint += prop_list[i];
  }

  hex_prop_cont->update_qint_qext_el( hex_qint, hex_qext, this->index);
}

/* Calculate all reaction propensities in this hexamer
   that depend on the external state variable Afree. */
double Hexamer::update_ext_prop()
{
  /* Sequestration by CI and CII domain depend on Afree */
  // CII: CII + A -> ACII
  prop_list[0] = prefactor_r0 * sys->cAfree;
   
  // CI: CIina + A -> ACIina
  prop_list[4] = prefactor_r4 * sys->cAfree;  

  return prop_list[0] + prop_list[4];
}

/* Calculate all reaction propensities possible in this Hexamer.
   Hexamer has HEXAMER_N_REACT + 6 reaction channels, 
   which propensities are stored in prop_list. */
void Hexamer::set_propensities()
{  
  /* Calculate reactionchannel propensities */   
  // CII: CII + A <-> ACII
  prefactor_r0 = (1-CIIKaiA_bound) * kCIIAon();
  prop_list[0] = prefactor_r0 * sys->cAfree;
  prop_list[1] = CIIKaiA_bound * kCIIAoff();
  
  // CI: CI + B <-> CI*B (NOT [KaiB] dependent)
  prop_list[2] = ( CIKaiB_bound < reaction_consts->nBseq ) ? kCIBon() : 0.0;
  //prop_list[2] = ( reaction_consts->nBseq - CIKaiB_bound ) * kCIBon();
  prop_list[3] = CIKaiB_bound * kCIBoff();
  //prop_list[3] = ( CIKaiA_bound == 0 ) ? CIKaiB_bound * kCIBoff() : 0.0;  
  //prop_list[3] = ( CIKaiB_bound > 0 ) ? kCIBoff() : 0.0;
  
  // CI: CI + A <-> CI*A
  prefactor_r4 = ( CIKaiB_bound == reaction_consts->nBseq ) ? 
    (reaction_consts->nAseq - CIKaiA_bound) * kCIAon() : 0.0;
  //prefactor_r4 = ( CIKaiB_bound > CIKaiA_bound ) ? 
  //  (reaction_consts->nAseq - CIKaiA_bound) * kCIAon() : 0.0;  
  prop_list[4] = prefactor_r4 * sys->cAfree;
  prop_list[5] = CIKaiA_bound * kCIAoff();
  
  // CI: active <--> inactive
  prop_list[6] = (double) active * kconf_f();
  prop_list[7] = (double) (1-active) * kconf_b();
  
  // Recalculate monomer propensities due to changes in hexamer state variables.
  for( int i(0); i < 6; i++ )
  {
    sextet[i]->set_propensities();  
  }
  
  // Monomer propensities
  prop_list[8]  = mon_prop_cont.get_qtot_el(0);
  prop_list[9]  = mon_prop_cont.get_qtot_el(1);
  prop_list[10] = mon_prop_cont.get_qtot_el(2);
  prop_list[11] = mon_prop_cont.get_qtot_el(3);
  prop_list[12] = mon_prop_cont.get_qtot_el(4);
  prop_list[13] = mon_prop_cont.get_qtot_el(5);  
      
  /* Calculate total, external and internal reaction propensities */
  update_prop_cont();
}

/* Conformation switch rate from Active to Inactive state. */
double Hexamer::kconf_f()
{
  double beta_energy( dGconfADP() + dGconfKaiA() + dGconfKaiB() ); 

  return reaction_consts->kconf0 * exp(-0.5 * beta_energy);
}

double Hexamer::kconf_b()
{
  double beta_energy( dGconfADP() + dGconfKaiA() + dGconfKaiB() ); 

  return reaction_consts->kconf0 * exp(0.5 * beta_energy);
}


/* Free energy difference between Active and Inactive state,
   dependent on ADP in CI domain. */
double Hexamer::dGconfADP()
{
  int    nIADPtot( 6 - this->get_hex_CIATP() );
   
  return reaction_consts->ddGTDconf * (reaction_consts->nIADPconfref - nIADPtot);
}


/* Free energy difference between Active and Inactive state,
   dependent on KaiB bound to the CI domain. */
double Hexamer::dGconfKaiB()
{
  const double KdAct( reaction_consts->kACIBoff/reaction_consts->kACIBon );
  const double KdIna( reaction_consts->kICIBoff/reaction_consts->kICIBon );  
    
  double beta_energy( (double)CIKaiB_bound * log( KdIna/KdAct ) );
  
  return beta_energy;
}


/* Free energy difference between Active and Inactive state,
   dependent on KaiA bound to the CI and/or CI domain. */
double Hexamer::dGconfKaiA()
{
  const double KdAct( reaction_consts->kACIAoff/reaction_consts->kACIAon );
  const double KdIna( reaction_consts->kICIAoff/reaction_consts->kICIAon );  
    
  double beta_energy( (double)CIKaiA_bound * log( KdIna/KdAct )  );

  beta_energy += CIIKaiA_bound * reaction_consts->dGICII;
  
  return beta_energy;
}


/* Binding free energy difference between CII and CII.A, in Active and Inactive state */
double Hexamer::dGACIIbind()
{
  double dgnU((double) get_hex_pU() * reaction_consts->dgACIIU);
  double dgnT((double) get_hex_pT() * reaction_consts->dgACIIT);
  double dgnD((double) get_hex_pD() * reaction_consts->dgACIID);
  double dgnS((double) get_hex_pS() * reaction_consts->dgACIIS);
  
  double beta_energy( dgnU + dgnT + dgnD + dgnS + (1-active) * reaction_consts->dGICII );

  return beta_energy;
}


/* Activation energy between ADP bound state ADP unbound state,
   in Active state dependent, on CII phos state of whole hexamer. */
double Hexamer::dGACIActADPoff()
{
  double dgnU((double) get_hex_pU() * reaction_consts->dgACIActU);
  double dgnT((double) get_hex_pT() * reaction_consts->dgACIActT);
  double dgnD((double) get_hex_pD() * reaction_consts->dgACIActD);
  double dgnS((double) get_hex_pS() * reaction_consts->dgACIActS);
  
  return dgnU + dgnT + dgnD + dgnS;
}


/* Activation energy between ADP bound state ADP unbound state,
   in Inactive state, dependent on CII phos state of whole hexamer. */
double Hexamer::dGICIActADPoff()
{
  double dgnU((double) get_hex_pU() * reaction_consts->dgICIActU);
  double dgnT((double) get_hex_pT() * reaction_consts->dgICIActT);
  double dgnD((double) get_hex_pD() * reaction_consts->dgICIActD);
  double dgnS((double) get_hex_pS() * reaction_consts->dgICIActS);
  
  return dgnU + dgnT + dgnD + dgnS;
}


double Hexamer::kCIBon()
{  
  return (double) active * reaction_consts->kACIBon + (1-active) * reaction_consts->kICIBon;
}


double Hexamer::kCIBoff()
{
  return (double) active * reaction_consts->kACIBoff + (1-active) * reaction_consts->kICIBoff;
}


/* KaiA on-rate to CI domain, Active/Inactive -conformation dependent. */
double Hexamer::kCIAon()
{
  return (double) active * reaction_consts->kACIAon + (1-active) * reaction_consts->kICIAon;
}


/* KaiA off-rate to CI domain, Active/Inactive -conformation dependent. */
double Hexamer::kCIAoff()
{
  return (get_CIKaiB_bound() == reaction_consts->nBseq) ? reaction_consts->kICIAoff : reaction_consts->kACIAoff;
  //return active * reaction_consts->kACIAoff + (1-active) * reaction_consts->kICIAoff;
  //return (get_CIKaiA_bound() > get_CIKaiB_bound()) ? reaction_consts->kACIAoff : reaction_consts->kICIAoff;
}


/*  KaiA on-rate to CII domain */
double Hexamer::kCIIAon()
{
  return reaction_consts->kCIIAon * exp( -0.5*dGACIIbind() );
}


/*  KaiA off-rate to CII domain */
double Hexamer::kCIIAoff()
{
  return reaction_consts->kCIIAoff * exp( 0.5*dGACIIbind() );
}


/* ADP off rate in CI domain */
double Hexamer::kCIADPoff()
{
  double kCIADPoff(reaction_consts->kCIADPoff0), 
         beta_actenergy(0);

  //Energy of transition state.
  if(get_CIKaiB_bound() == reaction_consts->nBseq && reaction_consts->nBseq > 0)
    beta_actenergy += this->dGICIActADPoff();
  else
    beta_actenergy += this->dGACIActADPoff();    

  // KaiB lowers ADP off rate by stabilization of CI-ADP state.
  beta_actenergy += get_CIKaiB_bound() * reaction_consts->kAkIDoff;
  
  return kCIADPoff * exp( -beta_actenergy );
}


//Print reaction propensities of this Hexamer.
void Hexamer::print_hex_props()
{
  cout << endl;
  cout << "### Propensities of Hexamer - idx " << index << endl;
  cout << "Hex state - active: " << active << ", CIIA: " << CIIKaiA_bound 
       << ", CIB " << CIKaiB_bound << ", CIA " << CIKaiA_bound << endl;
  cout << "Mon state - S/T tot: " << get_hex_S() << "/" << get_hex_T() 
       << ", CI/CII-ADP tot: " << 6-get_hex_CIATP() << "/" << 6-get_hex_CIIATP() << endl;
  cout << "CII + A -> CIIA    (0) " << prop_list[0] << endl;
  cout << "CIIA -> A + CII    (1) " << prop_list[1] << endl;
  cout << "CIB(m) + B -> CIB(m+1) (2) " << prop_list[2] << endl;
  cout << "CIB(m) -> B + CIB(m-1) (3) " << prop_list[3] << endl;  
  cout << "CI + A -> CIA      (4) " << prop_list[4] << endl; 
  cout << "CIA -> A + CI      (5) " << prop_list[5] << endl;  
  cout << "Active -> Inactive (6) " << prop_list[6] << endl;
  cout << "Inactive -> Active (7) " << prop_list[7] << endl;  
  cout << "Monomer qtot 1     (8) " << prop_list[8] << endl;
  cout << "Monomer qtot 2     (9) " << prop_list[9] << endl;
  cout << "Monomer qtot 3     (10) " << prop_list[10] << endl;
  cout << "Monomer qtot 4     (11) " << prop_list[11] << endl;
  cout << "Monomer qtot 5     (12) " << prop_list[12] << endl;      
  cout << "Monomer qtot 6     (13) " << prop_list[13] << endl;
  cout << "Props - tot: " << hex_prop_cont->get_qtot_el(index) 
       << ", int: " << hex_prop_cont->get_qint_el(index) 
       << ", ext: " << hex_prop_cont->get_qext_el(index) << endl;
}

//Total #nucleotide binding pockets with ATP in CI domain.
int Hexamer::get_hex_CIATP()
{
  int nItot(0);
  
  for(int i(0); i<6; i++)
  {
    nItot += sextet[i]->get_nI();
  }

  return nItot;
}

int Hexamer::get_hex_CIIATP()
{
  int nIItot(0);
  
  for(int i(0); i<6; i++)
  {
    nIItot += sextet[i]->get_nII();
  }

  return nIItot;
}

// Return #monomers in U-state.
int Hexamer::get_hex_pU()
{
  int Utot(0);
  
  for(int i(0); i<6; i++)
  {
    Utot += (1-sextet[i]->get_T()) * (1-sextet[i]->get_S());
  }

  return Utot;
}

// Return #monomers in T-state.
int Hexamer::get_hex_pT()
{
  int Ttot(0);
  
  for(int i(0); i<6; i++)
  {
    Ttot += (1-sextet[i]->get_S()) * sextet[i]->get_T();
  }

  return Ttot;
}

// Return #monomers in S-state.
int Hexamer::get_hex_pS()
{
  int Stot(0);
  
  for(int i(0); i<6; i++)
  {
    Stot += (1-sextet[i]->get_T()) * sextet[i]->get_S();
  }

  return Stot;
}

// Return #monomers in D-state.
int Hexamer::get_hex_pD()
{
  int Dtot(0);
  
  for(int i(0); i<6; i++)
  {
    Dtot += sextet[i]->get_T() * sextet[i]->get_S();
  }

  return Dtot;
}

// Return #monomers with T site phosphorylated.
int Hexamer::get_hex_T()
{
  int Ttot(0);
  
  for(int i(0); i<6; i++)
  {
    Ttot += sextet[i]->get_T();
  }

  return Ttot;
}

// Return #monomers with S site phosphorylated.
int Hexamer::get_hex_S()
{
  int Stot(0);
  
  for(int i(0); i<6; i++)
  {
    Stot += sextet[i]->get_S();
  }

  return Stot;
}


