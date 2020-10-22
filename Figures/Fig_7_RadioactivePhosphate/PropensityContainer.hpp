/*
  Data container for internal, external and total reaction propensities.
  Total is always the sum of int and ext propensities.
  
  Int: Propensity that only depends on internal state.
  Ext: Propensity that depends on both internal and external states.
       If external state changes, all propensities change.
*/

#ifndef _PROPENSITYCONTAINER_H_
#define _PROPENSITYCONTAINER_H_

#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>
#include <math.h>

#include "main.hpp"
#include "Monomer.hpp"

class PropensityContainer
{
private:
    
public: //Constructors
    //Create new propensity container with length plen.
    PropensityContainer(int plen)
    {
      len = plen;
      qint = 0.;
      qext = 0.;
      qtot = 0.;
      
      qint_list = new double[len];
      qext_list = new double[len];
      qtot_list = new double[len];
      
      for(int i(0); i<plen; i++)
      {
        qint_list[i] = 0;
        qext_list[i] = 0;
        qtot_list[i] = 0;
      }
    }
   

public: // Public function definitions.   

    int choose_index(double rnd);
   
    double sum_qint();
    double sum_qext();
    double sum_qtot();

    
    /* Setters */
    void update_qint_qext_el(double, double, int);
    void set_qint_qext_el(double, double, int);
      
    void update_qint_el(double, int);
    void update_qext_el(double, int);

    void set_qint_el(double, int);
    void set_qext_el(double, int);
    void set_qext_all(Hexamer *hexamers);
    
          
    /* Getters */
    double get_qint()
    {
      return this->qint;
    }

    double get_qext()
    {
      return this->qext;
    }
    
    double get_qtot()
    {
      return this->qtot;
    }
    
    double get_qint_el(int index)
    {
      return this->qint_list[index];
    }

    double get_qext_el(int index)
    {
      return this->qext_list[index];
    }
    
    double get_qtot_el(int index)
    {
      return this->qtot_list[index];
    }
    
    int get_len()
    {
      return this->len;
    }

private: // Private member variables.
    double *qint_list;
    double *qext_list;
    double *qtot_list;
    
    double qint;
    double qext;
    double qtot;
    
    int len;
};

#endif
