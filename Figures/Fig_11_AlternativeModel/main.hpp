#ifndef _GKaiC_H_
#define _GKaiC_H_

/*---------------------------INCLUDE'S---------------------------------------*/
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>
#include <math.h>

/* SIMULATION WIDE FUNCTIONS */
#define F_COMP(a,b) (fabs ( a-b ) < 1e-10)
#define MIN(a,b)    ( a < b ? a : b )
#define MAX(a,b)    ( a > b ? a : b )

/* SIMULATON WIDE CONSTANTS */
#define EPSILON     1e-9
#define SCALE       1 - EPSILON
#define PI          3.141592653589793
#define INF_POS     1.0/0.0

#define CANCELLED   -1
#define FALSE 0
#define TRUE  1

#define MONOMER_N_REACTS 12
#define HEXAMER_N_REACTS 8
#define DATA_NBR 3

#define HIST_RES 1

/* Structure definitions */
class Monomer;
class Hexamer;
class PropensityContainer;

/* Array of integers or doubles */
typedef struct int_array 
{
  int len,  ///< The length of the array
      *val; ///< A pointer to the first element
} IntArray;

typedef struct double_array 
{
  int len;  ///< The length of the array
  double *val; ///< A pointer to the first element
} DoubleArray;


/* Structure for all reaction constants. */
typedef struct reaction_constants
{ 
  //CI-nucleotide exchange
  double kCIhyd;
  double kCIADPoff0;
  double kCIADPoffA;
  double kCIATPoff; //0

  double kAkIDoff;
  double kconf0;
  double kBswitch_f;
  double ddGTDconf;
  double ddGTDAbind;  
  int    nIADPAref;
  int    nIADPconfref;  
  double dGIbindB;
  
  
  //CI-KaiA sequestration.
  double kACIAon;
  double kACIAoff;
  double kICIAon;
  double kICIAoff;
  double kICIBon;
  double kICIBoff;  
  double kACIBon;
  double kACIBoff;  
  double nAseq;
  double nBseq;  
  
  //Act-Ina dG due to phosph.
  double dgICIU;
  double dgICIT;
  double dgICID;
  double dgICIS;      

  double dgACIIU;
  double dgACIIT;
  double dgACIID;
  double dgACIIS;      
  double dGICII;
  
  //Activation energy of ADP off rate in CI domain, Active state
  double dgACIActU;
  double dgACIActT;
  double dgACIActD;
  double dgACIActS;

  //Activation energy of ADP off rate in CI domain, Inactive state
  double dgICIActU;
  double dgICIActT;
  double dgICIActD;
  double dgICIActS;
  
  //CII-nucleotide exchange
  double kCIIAon;
  double kCIIAoff;
  double kCIIhyd0;
  double kCIIhydA;  
  double kCIInucloff0;
  double kCIInucloffA;
  double KATPoKADP;

  //CII Phosphotransfer rates.
  double kUT;
  double kTU;
  double kTD;
  double kDT;
  double kSD;
  double kDS;
  double kUS;
  double kSU;  
  //Spontanious (Cp -> C + Pi) decay rates
  double kTUs;
  double kDTs;
  double kDSs;
  double kSUs;    
  
} ReactionConstants;

/* Data container for output data */
typedef struct hexamer_avr
{
  double t;
  double p;
  double Afree;
  double ACI;
  double BCI;  
  double ACII;
  double Atot;
  double CIATP;
  double CIIATP;
  double Stot;
  double Ttot;  
  double pU;
  double pT;
  double pS;
  double pD;
  double CIATPcons;
  double CIIATPcons;  
  double dCIATPcons;
  double dCIIATPcons;  
  double active;

  double dGACIIbind;
  double dGICIbind;
  double dGconfADP;
  double kADPoff;
  
} HexamerAvr;

/* Global state variables */
typedef struct system_variables
{
  //Static vars.
  double ATPfrac;
  double KaiA0; //Initial conditions in uM.
  double KaiB0;
  double KaiC0;
  double volume; //Volume is 1 cubic micron -> 600.
  
  double tend;
  double tsim;
  double tequ;
  
  bool start_phosphorylated;
  
  //Number of hexamers in system.
  int N_hexamers;
  
  //File names
  char output_filename[128];
  char input_filename[128];
  char param_filename[128];
  
  //Dynamic vars
  int Afree; //Free KaiA dimers; 
  double cAfree; //Free KaiA concentration 
  int CIATPcons; //Number of consumed ATP molecules in CI domain.
  int CIIATPcons; //Number of consumed ATP molecules in CII domain.
  
  //Radioactive phosphate counters;
  int radPi;
  int radST;
  int radATP;
  
  //Vars related to ADP shift experiment.
  bool nowisADPphase;
  double tstartADP;
  double tendADP;
  double ATPfrac_ADPphase;
  
  //Counters
  int step_cntr;
  double t_sample_incr;
  int sample_cnt; //Max sample number.
  int sample_cntr; //Current sample number.
  
  int rnd_seed;
  
  //For flux diagram; track all changes in S,T space.
  double hist_pTVSpS[6*HIST_RES+1][6*HIST_RES+1];
  uint rate_pTp1[6*HIST_RES+1][6*HIST_RES+1];
  uint rate_pTm1[6*HIST_RES+1][6*HIST_RES+1];
  uint rate_pSp1[6*HIST_RES+1][6*HIST_RES+1];
  uint rate_pSm1[6*HIST_RES+1][6*HIST_RES+1]; 
  
  //Flux measured at certain intervals only.
  bool register_flux;
  uint index_reg;
  double interval_reg;
  double treg[100000];
  double start_frac_reg;
  
  //State flips
  uint FlipCntrFw;
  uint FlipCntrBw;
  
  /* Pointer to data struct for output data, 
     split to active and inactive hexamers*/
  HexamerAvr *Aoutput_data;
  HexamerAvr *Ioutput_data;  
} SystemVariables;
#endif
