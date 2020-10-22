/* 
  Definition of the Monomer class. In here, al state variables that make a monomer
  and the fuctions acting on them are defined.
*/
#ifndef _MONOMER_H_
#define _MONOMER_H_

#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>
#include <math.h>

#include "main.hpp"
#include "PropensityContainer.hpp"


class Monomer
{
private: // Constant definitions
    // Number of reactions in monomer;
    //static const int NREACTS = ; 
    
public: //Constructors
    Monomer()
    : index(0), T(0), S(0), nI(1), nII(1) {}
    
    //Monomer(ReactionConstants *rc, PropensityContainer *pc)
    //: index(0), T(0), S(0), nI(1), nII(1), active(1), reaction_consts(rc), prop_cont(pc) {}

    Monomer(int pindex, bool pT, bool pS, bool pnI, bool pnII)
    : index(pindex), T(pT), S(pS), nI(pnI), nII(pnII) {}    

public: // Public function definitions.

    /* Overhead functions */
    // Initialize monomer state
    void initialize_state( int pindex, bool pT, bool pS, bool pnI, bool pnII );
    //Propagate monomer by one reaction.
    bool propagate();
    //Set reaction firing propensities given monomer state.
    void set_propensities();
    //Update reaction propensities that depend on hexamer state.
    double update_ext_prop();
    //Update external and internal reactions.
    void update_mon_prop_cont();   
    //Choose monomer reacion to fire given random number.
    int choose_reaction(double rnd);
    //Fire reaction; change mononer state accordingly.
    void fire_reaction(int reaction_channel);       
    //Print all reaction firing propensities.
    void print_mon_props();

    
    /* Reaction rates */
    double kCIInucloff();
    
    
    /* Getters */
    bool get_T()
    {
      return this->T;
    }

    bool get_S()
    {
      return this->S;
    }
    
    bool get_nI()
    {
      return this->nI;
    }

    bool get_nII()
    {
      return this->nII;
    }

    double get_propList(int i)
    {
      if( i < 0 || i >= MONOMER_N_REACTS)
        return 0;
      
      return this->prop_list[i];
    }
    
    
    /* Setters */    
    void set_hexamer(Hexamer *phexamer)
    {
      this->hexamer = phexamer;
    }
    
    void set_reaction_consts(ReactionConstants *preaction_consts)
    {
      this->reaction_consts = preaction_consts;
    }
    
    void set_prop_container(PropensityContainer *pprop_cont)
    {
      this->mon_prop_cont = pprop_cont;
    }
    
    void set_sysvars(SystemVariables *psys)
    {
      this->sys = psys;
    }

private: // Private member variables.

    /* Monomer state variables. */
    //Monomer unique index, and a ptr to the hexamer object it belongs to.
    int index;
    Hexamer *hexamer;
    // Phosphorylation state T-site (1: phos, 0: unphos)
    bool T;
    // Phosphorylation state S-site 
    bool S;
    // Nucleotide binding pocket of CI domain (1: ATP, 0: ADP)
    bool nI;
    // Nucleotide binding pocket of CII domain.
    bool nII;  
        
    //Reaction propensities of this monomer.
    double prop_list[MONOMER_N_REACTS];

    // Pointer to lists
    ReactionConstants *reaction_consts;
    PropensityContainer *mon_prop_cont;
    SystemVariables *sys;
};

#endif

