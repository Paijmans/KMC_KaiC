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


void PropensityContainer::update_qint_qext_el(double qint_el, double qext_el, int index)
{
  qint -= qint_list[index];
  qint += qint_el;
  
  qext -= qext_list[index];
  qext += qext_el;
  
  qtot = qint + qext;
 
  qint_list[index] = qint_el;
  qext_list[index] = qext_el;
  qtot_list[index] = qint_list[index] + qext_list[index];
}


void PropensityContainer::set_qint_qext_el(double qint_el, double qext_el, int index)
{ 
  qint_list[index] = qint_el;
  sum_qint();
  
  qext_list[index] = qext_el;
  sum_qext();
  
  qtot_list[index] = qint_list[index] + qext_list[index];  
  sum_qtot();
}


void PropensityContainer::update_qint_el(double qint_el, int index)
{
  qint -= qint_list[index];
  qint += qint_el;
  qtot = qint + qext;
    
  qint_list[index] = qint_el;
  qtot_list[index] = qint_list[index] + qext_list[index];
}


void PropensityContainer::set_qint_el(double qint_el, int index)
{      
  qint_list[index] = qint_el;
  sum_qint();
 
  qtot_list[index] = qint_list[index] + qext_list[index];
  sum_qtot();
}
  
    
void PropensityContainer::update_qext_el(double qext_el, int index)
{     
  qext -= qext_list[index];
  qext += qext_el;
  qtot = qint + qext;
    
  qext_list[index] = qext_el;
  qtot_list[index] = qint_list[index] + qext_list[index];    
}


void PropensityContainer::set_qext_el(double qext_el, int index)
{     
  qext_list[index] = qext_el;
  sum_qext();
 
  qtot_list[index] = qint_list[index] + qext_list[index];
  sum_qtot();     
}


void PropensityContainer::set_qext_all(Hexamer *hexamers)
{     
  for(int i(0); i<len; i++)
  {
    qext_list[i] = hexamers[i].update_ext_prop();
    qtot_list[i] = qint_list[i] + qext_list[i];
  }
      
  sum_qext();
  sum_qtot();       
}


double PropensityContainer::sum_qint()
{
  double sum(0.0);
    
  for(int i(0); i<len; i++)
  {
    sum += qint_list[i];
  }
      
  return this->qint = sum;
}
    
    
double PropensityContainer::sum_qext()
{
  double sum(0.0);
    
  for(int i(0); i<len; i++)
  {
    sum += qext_list[i];
  }
      
  return this->qext = sum;
}
    
    
double PropensityContainer::sum_qtot()
{
  double sum(0.0);
    
  for(int i(0); i<len; i++)
  {
    sum += qtot_list[i];
  }
      
  return this->qtot = sum;
}


int PropensityContainer::choose_index(double rnd)
{
  double qacc(qtot_list[0]);
  int i(0);
  
  rnd*=SCALE;  
  rnd*=qtot;
  while( rnd > qacc )
  {
    qacc += qtot_list[++i];
  }
  
  return i;
}
