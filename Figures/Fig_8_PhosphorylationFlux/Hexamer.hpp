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
  
    void initialize_state( int pindex, bool pactive, 
                           bool pCIIKaiA_bound, int pCIKaiA_bound, int pCIKaiB_bound );
    bool propagate();

    void set_propensities();
    double update_ext_prop();
    int choose_reaction(double rnd);
    bool fire_reaction(int reaction_channel);
    //void update_propensities();  
    void update_prop_cont();
    void register_state(int idx);        
    
    
    //Hexamer-state dependent reactions
    double kCIBon();
    double kCIBoff();    
    double kCIAon(); 
    double kCIAoff();
    double kCIADPoff();
    
    //double kCIInucloff();
    double kCIIAon();
    double kCIIAoff();

    //Conformational state flipping rates.
    double kconf_f();
    double kconf_b();
    
    //KaiA CII binding energy.
    double dGACIIbind();
    //Conformation state energy difference due to nucleotides in CI.
    double dGconfADP();
    //Conformational state energy difffernce due to bound KaiA.
    double dGconfKaiA();
    //Conformational state energy difffernce due to bound KaiA.
    double dGconfKaiB();
    //ADP binding on CI activation energy in Active state.
    double dGACIActADPoff();
    //ADP binding on CI activation energy in Inactive state.          
    double dGICIActADPoff();    
    
    void print_hex_props();
    
    
    
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

public:
    //Some data gathering.
    bool Acheck, Icheck, Bcheck, Scheck;
    double state[7][4];    
    
    //Flux plot for ind. hexamers.
    double tlastswitch;

private: // Private member variables.

    /* Monomer state variables. */
    // Hexamer unique index.
    int index;
    
    //Conformational state of hexamer: 1=active, 0=inactive.
    bool active;
    
    // Is KaiA bound to CII domain.
    bool CIIKaiA_bound;
    
    // Number of KaiA sequestered by CI domain.
    int CIKaiA_bound;

    // Number of KaiB sequestered by CI domain.
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

