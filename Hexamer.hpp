/*
  Definition of the Hexamer class, which includes all state variables 
  describing a hexamer, and all functions acting on these variables.
*/

#ifndef _HEXAMER_H_
#define _HEXAMER_H_

#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>
#include <math.h>

#include "main.hpp"
#include "PropensityContainer.hpp"
#include "Monomer.hpp"
#include "Hexamer.hpp"

class Hexamer
{
private: // Constant definitions
    
public: //Constructors
    Hexamer()
    : index(0), active(1), CIIKaiA_bound(0), CIKaiA_bound(0), 
      CIKaiB_bound(0), mon_prop_cont(6) {}
    
public: // Public function definitions.
  
    /* Overhead functions */
    //Initialize hexamer state.
    void initialize_state( int pindex, bool pactive, 
                           bool pCIIKaiA_bound, int pCIKaiA_bound, int pCIKaiB_bound );
    //Propagate hexamer one reaction.
    bool propagate();
    //Set reaction firing propensities given hexamer state.
    void set_propensities();
    //Update reaction that change propensity when KaiA concentration changes.
    double update_ext_prop();
    //Choose reaction to fire given random number.
    int choose_reaction(double rnd);
    //Fire reaction with index reaction_channel.
    bool fire_reaction(int reaction_channel);
    //Update reaction propensities in reaction container.
    void update_prop_cont();
    //Print all hexamer propensities.
    void print_hex_props();

    
    /* Reaction rate functions */
    //Conformational forward and backward flipping rates.
    double kconf_f();
    double kconf_b();
    
    //KaiB on and off rates for the CI domain.
    double kCIBon();
    double kCIBoff();
    
    //KaiA on and off rates for the CI domain.
    double kCIAon(); 
    double kCIAoff();
    
    //ADP dissociation rate from CI domain.
    double kCIADPoff();
    
    //KaiA on and off rates for the CII domain.
    double kCIIAon();
    double kCIIAoff();

   
    /* State energies */
    //KaiA CII binding energy.
    double dGACIIbind();
    //Conformation state energy difference due to nucleotides in CI.
    double dGconfADP();
    //Conformational state energy difffernce due to bound KaiA.
    double dGconfKaiA();
    //Conformational state energy difffernce due to bound KaiB.
    double dGconfKaiB();
    //ADP binding on CI activation energy in Active state.
    double dGACIActADPoff();
    //ADP binding on CI activation energy in Inactive state.          
    double dGICIActADPoff();    
      
    
    /* Getters */
    int get_hex_CIATP();
    int get_hex_CIIATP();
    
    int get_hex_T();
    int get_hex_S();
    int get_hex_pU();
    int get_hex_pT();
    int get_hex_pS();
    int get_hex_pD();
    
    bool get_CIIKaiA_bound()
    {
      return this->CIIKaiA_bound;
    }
    
    int get_CIKaiA_bound()
    {
      return this->CIKaiA_bound;
    }
    
    int get_CIKaiB_bound()
    {
      return this->CIKaiB_bound;
    }    
       
    bool get_active()
    {
      return this->active;
    }
    
    int get_index()
    {
      return this->index;
    }
    
    Monomer* get_Monomer(int i)
    {
      if( i < 0 || i > 6)
      {
        std::cerr << "Invalid Monomer index: " << i << std::endl;
        exit(1); 
      } 
       
      return sextet[i];
    }    
    
    
    /* Setters */
    void set_sextet(Monomer *monomers);
    
    void set_reaction_consts(ReactionConstants *preaction_consts)
    {
      this->reaction_consts = preaction_consts;
    }
    
    void set_prop_container(PropensityContainer *pprop_cont)
    {
      this->hex_prop_cont = pprop_cont;
    }
    
    void set_sysvars(SystemVariables *psys)
    {
      this->sys = psys;
    }

public: // Public member variables.


private: // Private member variables.

    /* Monomer state variables. */
    // Hexamer unique index.
    int index;
    
    //Conformational state of hexamer: 1=active, 0=inactive.
    bool active;
    
    // Is KaiA bound to CII domain.
    bool CIIKaiA_bound;
    
    // Number of KaiA dimers sequestered by CI domain.
    int CIKaiA_bound;

    // Number of KaiB monomers sequestered by CI domain.
    int CIKaiB_bound;    
       
    // Pointers to the six monomers in this hexamer.
    Monomer *sextet[6];

    // Reaction propensities of this hexamer + monomers.
    double prop_list[HEXAMER_N_REACTS + 6];  

    // Pointers to lists
    ReactionConstants   *reaction_consts;
    PropensityContainer *hex_prop_cont;
    PropensityContainer  mon_prop_cont;
    SystemVariables *sys;
    
    //Prefactors for external propensity calculations.
    double prefactor_r0, prefactor_r4;   
};

#endif

