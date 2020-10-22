/*
  Program to simulate the in vitro Kai system using Gillespie.
  We will use the model by van Zon et al, combined with the 
  new model which includes nucleotide binding pockets and KaiA
  acting as an nucleotide exchange factor.

  Hexamer consist of 6 monomers. 
  Each monomer calculates its propensities separately.
  
  USAGE:
  
  > ./GKaiC <parameter_filename> 
*/

#include <sstream>
#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>
#include <math.h>
#include <stdio.h>

#include "main.hpp"
#include "Monomer.hpp"
#include "Hexamer.hpp"
#include "PropensityContainer.hpp"
#include "propagate.hpp"
#include "random.hpp"

using namespace std;

/* local function definitions */
ReactionConstants initialize_reaction_consts(SystemVariables *sys);
void initialize_system_vars(SystemVariables *sys);
void write_outputfile(SystemVariables *sys);
void allocate_data_array(SystemVariables *sys);
HexamerAvr calc_hex_averages( Hexamer *hexamers, SystemVariables *sys, ReactionConstants *react_const );
string delUnnecessary(string &str);

/* Let's roll */
int main(int argc, char *argv[])
{ 
  /* Set-up propensities, monomers, system variables 
     and reaction constants - data structures */
  
  // First read in parameter filename;
  double temp;
  SystemVariables sys;
  if(argc > 1)
  {
    strncpy(sys.param_filename, argv[1], 128);
  }
  else
  {
    strcpy(sys.param_filename,"default.par");
  }
  cerr << "Parameter file set to: " << sys.param_filename << endl; 
  initialize_system_vars(&sys);
  ReactionConstants reaction_consts( initialize_reaction_consts(&sys) );

  //Initialize random number generator with (time-dependent, 0) seed.
  ran_init(sys.rnd_seed);

  PropensityContainer prop_cont(sys.N_hexamers); 
  Hexamer *hexamers = new Hexamer[sys.N_hexamers];
  Monomer *monomers = new Monomer[6*sys.N_hexamers];
  
  allocate_data_array(&sys);

  //Initialize hexamer states and calculate monomer and hexamer propensities.
  for(int i(0); i<sys.N_hexamers; i++)
  {
    hexamers[i].set_prop_container(&prop_cont);
    hexamers[i].set_reaction_consts(&reaction_consts); 
    hexamers[i].set_sysvars(&sys);
    hexamers[i].initialize_state(i,1,0,0,0); 
    hexamers[i].set_sextet(monomers); //sets monomers state and propensities.
    hexamers[i].set_propensities();
    //hexamers[i].print_hex_props();
  }

  double t_sample(sys.tequ);
  cerr << "Afree: " << sys.Afree << ", tsim: " << sys.tsim << ", tend: " << sys.tend << endl;   

  //Propagate simulation until tend is reached.
  //sys.register_flux=1;
  while( propagate(&sys, hexamers, &prop_cont) ) 
  { 
    /* Record samples for output data */
    if(sys.tsim > t_sample)
    {   
      //Calculate quantities of interest.
      HexamerAvr hex_avr( calc_hex_averages( hexamers, &sys, &reaction_consts ) );
          
      //Print to screen.
      if(sys.sample_cntr%100 == 0)
      {     
        cerr << endl << "### SIMULATOR STATE ###" << endl;
        cerr << "t: " << sys.tsim << ", p(t): " << hex_avr.p 
             << ", Act(t): " << hex_avr.active << endl;
        cerr << "Afree: " << hex_avr.Afree << ", Atot: " << hex_avr.Atot 
             << ", ACII: " << hex_avr.ACII << ", ACI: " << hex_avr.ACI << endl;
        cerr << "CIATP: " << hex_avr.CIATP << ", CIIATP: " << hex_avr.CIIATP << endl;
        cerr << "pU: " << hex_avr.pU << ", pT: " << hex_avr.pT 
             << ", pD: " << hex_avr.pD << ", pS: " << hex_avr.pS << endl;
      }

      t_sample += sys.t_sample_incr;
      sys.sample_cntr++;
    }
    
    /* Change ATP fraction in bulk to simulate phase-shift experiment */
    if(sys.tsim >= sys.tstartADP && sys.tsim < sys.tendADP && !sys.nowisADPphase)
    {   
      sys.nowisADPphase=1;
      sys.ATPfrac = sys.ATPfrac_ADPphase;
      
      for(int i(0); i<sys.N_hexamers; i++)
      {
        hexamers[i].set_propensities();
      }
    
      cerr << "ADP phase: ENGAGED, tsim: " << sys.tsim << endl;      
    }
    else if(sys.tsim > sys.tendADP && sys.nowisADPphase)
    {
      sys.nowisADPphase = 0;
      sys.ATPfrac = 1.0;
      
      for(int i(0); i<sys.N_hexamers; i++)
      {
        hexamers[i].set_propensities();
      }  
      
      cerr << "ADP phase: DISENGAGED, tsim: " << sys.tsim << endl;      
    }

 
  }
  
  /* ### For normal output ### */
  //write_outputfile(&sys);  
  
  /* ### For flux histogram output ### */
  temp = 0.;
  for( int nT=0; nT<HIST_RES*6+1; nT++ )
  {
    for( int nS=0; nS<HIST_RES*6+1; nS++ )
    {
      temp += sys.hist_pTVSpS[nT][nS];
      cout << nT << "\t" << nS << "\t" 
           << (double)sys.hist_pTVSpS[nT][nS]/(sys.tsim * sys.N_hexamers) << "\t"
           << (double)sys.rate_pTp1[nT][nS]/(sys.tsim * sys.N_hexamers) << "\t"
           << (double)sys.rate_pTm1[nT][nS]/(sys.tsim * sys.N_hexamers) << "\t"
           << (double)sys.rate_pSp1[nT][nS]/(sys.tsim * sys.N_hexamers) << "\t"
           << (double)sys.rate_pSm1[nT][nS]/(sys.tsim * sys.N_hexamers) << endl;
    }
  }  
  cerr << "Total acc time: " << temp/sys.N_hexamers << endl;

  cerr << "System Ready" << endl;
}
  
  
//Loops over all hexamers and calculates avergaes which are store in hex_avr_data;
HexamerAvr calc_hex_averages( Hexamer *hexamers, SystemVariables *sys, ReactionConstants *reaction_consts )
{
  HexamerAvr Ahex_avr_data, Ihex_avr_data;
  
  //State dependent variables.
  int ICIBA(0), ACIBA(0);
  int ACIB(0), ICIB(0); 
  int IACII(0), AACII(0);
  int ICIATP(0), ACIATP(0), ICIATPstd(0), ACIATPstd(0);
  int ICIIATP(0), ACIIATP(0);
  int IStot(0), AStot(0);
  int ITtot(0), ATtot(0);
  int IpU(0), ApU(0);
  int IpT(0), ApT(0);
  int IpS(0), ApS(0);
  int IpD(0), ApD(0);
  int Ip(0), Ap(0);

  //State independent variables;
  int prev_CIATPcons(0), prev_CIIATPcons(0);
  int active(0), inactive(0);
  int Atot = 0;    
  
  //Energies
  double AdGconfADP(0.), IdGconfADP(0.), 
         IdGACIIbind(0.), AdGACIIbind(0.), 
         AADPoff(0.), IADPoff(0.);
  
  double N_hexamers(sys->N_hexamers);
  //double N_hexamers(1);  
  
  if(sys->sample_cntr > 0)
  {
    prev_CIATPcons = sys->Aoutput_data[sys->sample_cntr-1].CIATPcons * (6*N_hexamers);
    prev_CIIATPcons = sys->Aoutput_data[sys->sample_cntr-1].CIIATPcons * (6*N_hexamers);    
  }
      
  for(int i(0); i<sys->N_hexamers; i++)
  {  
    if( hexamers[i].get_active() )
    {
      Ap += hexamers[i].get_hex_pT() + hexamers[i].get_hex_pS() + hexamers[i].get_hex_pD();
      if( hexamers[i].get_CIKaiB_bound() == reaction_consts->nBseq && 
          hexamers[i].get_CIKaiA_bound() == 0 ) ACIB++;
      if( hexamers[i].get_CIKaiA_bound() > 0 ) ACIBA += 1;
      AACII += hexamers[i].get_CIIKaiA_bound();
      ACIATP += hexamers[i].get_hex_CIATP();
      ACIATPstd += hexamers[i].get_hex_CIATP()*hexamers[i].get_hex_CIATP();     
      ACIIATP += hexamers[i].get_hex_CIIATP();
      ATtot += hexamers[i].get_hex_T();    
      AStot += hexamers[i].get_hex_S();
      ApU += hexamers[i].get_hex_pU();
      ApT += hexamers[i].get_hex_pT();
      ApD += hexamers[i].get_hex_pD();
      ApS += hexamers[i].get_hex_pS();
      AdGconfADP += hexamers[i].kconf_f();
      AdGACIIbind += hexamers[i].kCIIAon();
      AADPoff += hexamers[i].kCIADPoff();      
            
      active++;
    }
    else
    {  
      Ip += hexamers[i].get_hex_pT() + hexamers[i].get_hex_pS() + hexamers[i].get_hex_pD();     
      if( hexamers[i].get_CIKaiB_bound() == reaction_consts->nBseq && 
          hexamers[i].get_CIKaiA_bound() == 0 ) ICIB++;
      if( hexamers[i].get_CIKaiA_bound() > 0 ) ICIBA += 1;
      IACII += hexamers[i].get_CIIKaiA_bound();
      ICIATP += hexamers[i].get_hex_CIATP();
      ICIATPstd += hexamers[i].get_hex_CIATP()*hexamers[i].get_hex_CIATP();
      ICIIATP += hexamers[i].get_hex_CIIATP();
      ITtot += hexamers[i].get_hex_T();    
      IStot += hexamers[i].get_hex_S();
      IpU += hexamers[i].get_hex_pU();
      IpT += hexamers[i].get_hex_pT();
      IpD += hexamers[i].get_hex_pD();
      IpS += hexamers[i].get_hex_pS();
      IdGconfADP += hexamers[i].kconf_b();            
      IADPoff += hexamers[i].kCIADPoff();
      IdGACIIbind += hexamers[i].kCIIAon();
    }
    Atot += hexamers[i].get_CIKaiA_bound() + hexamers[i].get_CIIKaiA_bound();      
  }
  Atot += sys->Afree;
  
  inactive = sys->N_hexamers - active;  
 
  double KaiA0(sys->KaiA0);
  double KaiC0(sys->KaiC0);
  if(KaiA0 == 0)
    KaiA0 = 1;
  if(KaiC0 == 0)
    KaiC0 = 1;

  Ahex_avr_data.CIATPcons  = (double) sys->CIATPcons/(6*N_hexamers);
  Ahex_avr_data.dCIATPcons = 
    (double) 24. * (sys->CIATPcons - prev_CIATPcons)/(6*N_hexamers*sys->t_sample_incr);
  Ahex_avr_data.CIIATPcons  = (double) sys->CIIATPcons/(6*N_hexamers);
  Ahex_avr_data.dCIIATPcons = 
    (double) 24. * (sys->CIIATPcons - prev_CIIATPcons)/(6*N_hexamers*sys->t_sample_incr);

  Ahex_avr_data.active   = (double) active;  
  Ihex_avr_data.active   = (double) inactive;

  Ahex_avr_data.Atot     = (double) Atot / (KaiA0 * sys->volume);    
  Ahex_avr_data.t        = (double) sys->tsim;
  Ahex_avr_data.Afree    = (double) sys->Afree / (KaiA0 * sys->volume);  
  Ihex_avr_data.Atot     = (double) Atot / (KaiA0 * sys->volume);    
  Ihex_avr_data.t        = (double) sys->tsim;
  Ihex_avr_data.Afree    = (double) sys->Afree / (KaiA0 * sys->volume);
    
  Ahex_avr_data.p        = (double) Ap/(6*N_hexamers);
  Ahex_avr_data.ACI      = (double) ACIBA/(KaiC0 * sys->volume);
  Ahex_avr_data.BCI      = (double) ACIB/(KaiC0 * sys->volume);    
  Ahex_avr_data.ACII     = (double) AACII/(KaiC0 * sys->volume);
  Ahex_avr_data.CIATP    = (double) ACIATP/N_hexamers;
  Ahex_avr_data.CIIATP   = (double) ACIIATP/N_hexamers;
  Ahex_avr_data.Ttot     = (double) ATtot/N_hexamers;
  Ahex_avr_data.Stot     = (double) AStot/N_hexamers;  
  Ahex_avr_data.pU       = (double) ApU/N_hexamers;
  Ahex_avr_data.pT       = (double) ApT/N_hexamers;
  Ahex_avr_data.pD       = (double) ApD/N_hexamers;
  Ahex_avr_data.pS       = (double) ApS/N_hexamers;

  Ahex_avr_data.dGconfADP = (double) AdGconfADP/N_hexamers;
  Ahex_avr_data.kADPoff = (double) AADPoff/N_hexamers;  
  Ahex_avr_data.dGACIIbind = (double) AdGACIIbind/N_hexamers;   

  Ihex_avr_data.p        = (double) Ip/(6*N_hexamers);
  Ihex_avr_data.ACI      = (double) ICIBA/(KaiC0 * sys->volume);
  Ihex_avr_data.BCI      = (double) ICIB/(KaiC0 * sys->volume);
  Ihex_avr_data.ACII     = (double) IACII/(KaiC0 * sys->volume);
  Ihex_avr_data.CIATP    = (double) ICIATP/N_hexamers;
  Ihex_avr_data.CIIATP   = (double) ICIIATP/N_hexamers;
  Ihex_avr_data.Ttot     = (double) ITtot/N_hexamers;
  Ihex_avr_data.Stot     = (double) IStot/N_hexamers;  
  Ihex_avr_data.pU       = (double) IpU/N_hexamers;
  Ihex_avr_data.pT       = (double) IpT/N_hexamers;
  Ihex_avr_data.pD       = (double) IpD/N_hexamers;
  Ihex_avr_data.pS       = (double) IpS/N_hexamers;
    
  Ihex_avr_data.dGconfADP = (double) IdGconfADP/N_hexamers;   
  Ihex_avr_data.kADPoff = (double) IADPoff/N_hexamers;
  Ihex_avr_data.dGACIIbind = (double) IdGACIIbind/N_hexamers;


   

  //Put data in array for writing to file.
  sys->Aoutput_data[sys->sample_cntr] = Ahex_avr_data;
  sys->Ioutput_data[sys->sample_cntr] = Ihex_avr_data;  

  //Calculate some cummulatives for printing to screen.
  HexamerAvr hex_avr;

  hex_avr.Afree = (double) sys->Afree / (KaiA0 * sys->volume);
  hex_avr.Atot = Atot / (KaiA0 * sys->volume);

  hex_avr.active = (double)Ahex_avr_data.active / N_hexamers;
  hex_avr.p = (Ahex_avr_data.p + Ihex_avr_data.p);
  hex_avr.ACI = Ahex_avr_data.ACI + Ihex_avr_data.ACI;
  hex_avr.ACII = Ahex_avr_data.ACII + Ihex_avr_data.ACII;
  hex_avr.CIATP = Ahex_avr_data.CIATP + Ihex_avr_data.CIATP;
  hex_avr.CIIATP = Ahex_avr_data.CIIATP + Ihex_avr_data.CIIATP;
  hex_avr.pU = Ahex_avr_data.pU + Ihex_avr_data.pU;
  hex_avr.pT = Ahex_avr_data.pT + Ihex_avr_data.pT;
  hex_avr.pD = Ahex_avr_data.pD + Ihex_avr_data.pD;
  hex_avr.pS = Ahex_avr_data.pS + Ihex_avr_data.pS;
  
  return hex_avr;  
}


//Allocate data for writing datparam_name.
void allocate_data_array(SystemVariables *sys)
{
  /* set-up data arrays for timetraces */
  sys->sample_cnt = (int) ((sys->tend - sys->tequ) / sys->t_sample_incr + 0.5);
  //Number of timetraces;
  sys->Aoutput_data = (HexamerAvr*) calloc( sys->sample_cnt ,sizeof(HexamerAvr) );
  sys->Ioutput_data = (HexamerAvr*) calloc( sys->sample_cnt ,sizeof(HexamerAvr) );    
}


/* Write timetraces to file
   3 files are produced: Averages for Active state hexamers
                         Averages for Inactive state hexamers
                         Averages for All hexamers 
*/
void write_outputfile(SystemVariables *sys)
{

  
  char Tfilename[128]; memset(&Tfilename[0], 0, sizeof(Tfilename));
  char Afilename[128]; memset(&Afilename[0], 0, sizeof(Afilename));
  char Ifilename[128]; memset(&Ifilename[0], 0, sizeof(Ifilename));
 
  strcat(Tfilename, sys->output_filename);  
  strcat(Tfilename, ".dat");  

  strcpy(Afilename,"A");  
  strcat(Afilename, Tfilename);
  strcpy(Ifilename,"I");
  strcat(Ifilename, Tfilename);
  
  FILE *Afp = fopen( Afilename, "w" );
  FILE *Ifp = fopen( Ifilename, "w" );
  FILE *Tfp = fopen( Tfilename, "w" );      
  
  cerr << "Writing output of Active hexamers to file: " << Afilename << endl;

  // Write file for Active hexamers.
  fprintf(Afp,"#time\tphos_frac\tAfree\tACI\tACII\tAtot\tCIATP\tCIIATP\tTtot\tStot\tpU\tpT\tpD\tpS\tCIATPcons\tdCIATPcons\tCIIATPcons\tdCIIATPcons\tdGconfADP\tdGACIIbind\tkADPoff\tactive\tCIBB\n");  

  // iterate through all samples
  for(int j(0); j < sys->sample_cnt; j++ ) 
  {
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].t );       //1
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].p );       //2
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].Afree );   //3
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].ACI );     //4
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].ACII );    //5
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].Atot );    //6
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].CIATP );   //7                      
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].CIIATP );  //8
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].Ttot );    //9
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].Stot );    //10
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].pU );      //11
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].pT );      //12
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].pD );      //13
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].pS );      //14
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].CIATPcons ); //15
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].dCIATPcons );//16                      
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].CIIATPcons ); //17
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].dCIIATPcons );//18    
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].dGconfADP ); //19     
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].dGACIIbind ); //20
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].kADPoff ); //21   
    fprintf( Afp, "%e\t", sys->Aoutput_data[j].active/sys->N_hexamers );  //22
    fprintf( Afp, "%e\n", sys->Aoutput_data[j].BCI );  //23
  }
    
  // close file
  fclose(Afp);

  cerr << "Writing output of Inactive hexamers to file: " << Ifilename << endl;
 
  // Write file for Inactive hexamers.
  fprintf(Ifp,"#time\tphos_frac\tAfree\tACI\tACII\tAtot\tCIATP\tCIIATP\tTtot\tStot\tpU\tpT\tpD\tpS\tCIATPcons\tdCIATPcons\tCIIATPcons\tdCIIATPcons\tdGconfADP\tdGACIIbind\tkADPoff\tinactive\tCIBB\n");

  // iterate through all samples
  for(int j(0); j < sys->sample_cnt; j++ ) 
  {
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].t );       //1
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].p );       //2
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].Afree );   //3
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].ACI );     //4
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].ACII );    //5
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].Atot );    //6
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].CIATP );   //7                      
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].CIIATP );  //8
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].Ttot );    //9
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].Stot );    //10
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].pU );      //11
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].pT );      //12
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].pD );      //13
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].pS );      //14
    fprintf( Ifp, "%e\t", sys->Aoutput_data[j].CIATPcons ); //15
    fprintf( Ifp, "%e\t", sys->Aoutput_data[j].dCIATPcons );//16                      
    fprintf( Ifp, "%e\t", sys->Aoutput_data[j].CIIATPcons ); //17
    fprintf( Ifp, "%e\t", sys->Aoutput_data[j].dCIIATPcons );//18                     
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].dGconfADP ); //19     
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].dGACIIbind ); //20
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].kADPoff ); //21
    fprintf( Ifp, "%e\t", sys->Ioutput_data[j].active/sys->N_hexamers ); //22            
    fprintf( Ifp, "%e\n", sys->Ioutput_data[j].BCI );  //23    
  }
    
  // close file
  fclose(Ifp);
  
  cerr << "Writing output of all hexamers to file: " << Tfilename << endl;
  
  // Write file for All hexamers.
    fprintf(Tfp,"#time\tphos_frac\tAfree\tACI\tACII\tAtot\tCIATP\tCIIATP\tTtot\tStot\tpU\tpT\tpD\tpS\tCIATPcons\tdCIATPcons\tCIIATPcons\tdCIIATPcons\tdGconfADP\tdGACIIbind\tkADPoff\tactive\tCIBB\n");  

  // iterate through all samples
  for(int j(0); j < sys->sample_cnt; j++ ) 
  { 
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].t );       //1
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].p + sys->Ioutput_data[j].p );       //2
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].Afree );   //3
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].ACI + sys->Ioutput_data[j].ACI );     //4
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].ACII + sys->Ioutput_data[j].ACII );    //5
    fprintf( Tfp, "%e\t", sys->KaiA0 );    //6
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].CIATP + sys->Ioutput_data[j].CIATP ) ; //7                      
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].CIIATP + sys->Ioutput_data[j].CIIATP ); //8
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].Ttot + sys->Ioutput_data[j].Ttot );    //9
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].Stot + sys->Ioutput_data[j].Stot );    //10
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].pU + sys->Ioutput_data[j].pU );      //11
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].pT + sys->Ioutput_data[j].pT );      //12
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].pD + sys->Ioutput_data[j].pD );      //13
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].pS + sys->Ioutput_data[j].pS );      //14
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].CIATPcons ); //15
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].dCIATPcons );//16                      
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].CIIATPcons ); //17
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].dCIIATPcons );//18                     
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].dGconfADP + sys->Ioutput_data[j].dGconfADP ); //19     
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].dGACIIbind + sys->Ioutput_data[j].dGACIIbind ); //20
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].kADPoff + sys->Ioutput_data[j].kADPoff ); //21
    fprintf( Tfp, "%e\t", sys->Aoutput_data[j].active/sys->N_hexamers );  //22   
    fprintf( Tfp, "%e\n", sys->Aoutput_data[j].BCI + sys->Ioutput_data[j].BCI ); //23
  }  
  
  // close file
  fclose(Tfp);  
}


/* Load system variables from file defined in command line argument. */
void initialize_system_vars(SystemVariables *sys)
{ 
  ifstream infile(sys->param_filename);
  string line, param_name, param_string, temp;
  double param_value;
  bool paramset(0);
  
  //Set parameter default values (can be overwritten by .par file)
  sys->volume                    = 1200;
  sys->KaiA0                     = 0.6;
  sys->KaiB0                     = 0.6;
  sys->KaiC0                     = 0.6;
  sys->tend                      = 120;
  sys->tequ                      = 0;
  sys->ATPfrac                   = 1.0;
  sys->t_sample_incr             = 0.1;
  sys->start_phosphorylated      = 0;
  sys->rnd_seed                  = 42;
  
  //cerr << endl << "SYSTEM VARIABLES USED IN THIS SIMULATION: " << endl;
  while (getline(infile, line))
  {
    //remove unnecessary spaces in line.
    line = delUnnecessary( line ); 
    istringstream iss(line);
    paramset = 0;
    
    // Check if line has formatting: '<string> <number>'
    if (iss >> param_name >> param_value) 
    {        
      /* Loading system variables */   
      if( param_name.compare("volume") == 0 )
      {
        sys->volume = param_value;
        paramset = 1;
      }  

      if( param_name.compare("KaiA0") == 0 )
      {
        sys->KaiA0 = param_value;
        paramset = 1;
      }  

      if( param_name.compare("KaiB0") == 0 )
      {
        sys->KaiB0 = param_value;
        paramset = 1;
      }  
     
      if( param_name.compare("KaiC0") == 0 )
      {
        sys->KaiC0 = param_value;
        paramset = 1;
      }  
      
      if( param_name.compare("tend") == 0 )
      {
        sys->tend = param_value;
        paramset = 1;
      }  
      
      if( param_name.compare("tequ") == 0 )
      {
        sys->tequ = param_value;
        paramset = 1;
      }  
      
      if( param_name.compare("ATPfrac") == 0 )
      {
        sys->ATPfrac = param_value;
        paramset = 1;
      }        
      
      if( param_name.compare("t_sample_incr") == 0 )
      {
        sys->t_sample_incr = param_value;
        paramset = 1;
      }       
      
      if( param_name.compare("start_phosphorylated") == 0 )
      {
        sys->start_phosphorylated = param_value;
        paramset = 1;
      }
                    
      if( param_name.compare("tstartADP") == 0 )
      {
        sys->tstartADP = param_value;
        paramset = 1;
      }
                
      if( param_name.compare("tendADP") == 0 )
      {
        sys->tendADP = param_value;
        paramset = 1;
      } 
                     
      if( param_name.compare("ATPfrac_ADPphase") == 0 )
      {
        sys->ATPfrac_ADPphase = param_value;
        paramset = 1;
      }   
      
      if( param_name.compare("rnd_seed") == 0 )
      {
        sys->rnd_seed = param_value;
        paramset = 1;
      }               

      if( param_name.compare("start_frac_reg") == 0 )
      {
        sys->start_frac_reg = param_value;
        paramset = 1;
      }   
      
      if( param_name.compare("interval_reg") == 0 )
      {
        sys->interval_reg = param_value;
        paramset = 1;
      }   
             
      if(!paramset)
      {
        //cerr << param_name << " - is not a valid parameter name." << endl;
      }
      else
      {
        //cerr << param_name << " = " << param_value << endl;
        //fprintf( Pfp, "%s\t%e\n", param_name.c_str(), param_value );     
      }
          
    }
    else
    {
      // Check if line has formatting: '<string> <string>'
      iss.clear();
      iss.str( line );
      if (iss >> param_name >> param_string)
      { 
         if( param_name.compare("input_filename") == 0 )
         {
           strcpy(sys->input_filename, param_string.c_str());
           paramset = 1;
         }
         
         if( param_name.compare("output_filename") == 0 )
         {
           strcpy(sys->output_filename, param_string.c_str());
           paramset = 1;
         }
                     
         if(!paramset)
         {
           //cerr << param_name << " - is not a valid parameter name." << endl;
         }
         else
         {
           //cerr << param_name << " = " << param_string << endl;
           //fprintf( Pfp, "%s\t%s\n", param_name.c_str(), param_string.c_str() );     
         }  
      }
      else
      {       
        //Check for comment and empty lines
        if(!(param_name.compare(0,1,"#") == 0) && !line.empty())
        {
           cerr << "Error - Line with incorrect formatting: -" << param_name << 
           "-" << param_string << "." << endl; 
           exit(1); 
        }   
        continue;
      }
    } //if not variable_name variable_value
  }//while 

  //Set system variables
  sys->tsim = 0.0;
  sys->N_hexamers = sys->KaiC0 * sys->volume;
  sys->Afree = sys->KaiA0 * sys->volume;
  sys->cAfree = sys->Afree / sys->volume;
  sys->CIATPcons = 0; sys->CIIATPcons = 0;
  
  sys->step_cntr = 0;
  sys->sample_cntr = 0;
  sys->radPi = 0;
  sys->radST = 0;
  sys->radATP = 0;  
  
  sys->nowisADPphase = 0;
 

  //RELATED TO fLUX HISTOGRAM. 
  sys->FlipCntrFw = 0;
  sys->FlipCntrBw = 0;
  sys->index_reg=0;
  sys->register_flux=0;
  for(int i=0;i<6*HIST_RES+1;i++)
  { for(int j=0; j<6*HIST_RES+1; j++)
    {
      sys->rate_pTp1[i][j]=0; 
      sys->rate_pTm1[i][j]=0;
      sys->rate_pSp1[i][j]=0;
      sys->rate_pSm1[i][j]=0; 
      sys->hist_pTVSpS[i][j]=0.;                 
    }
  }  
}


/* Function reads in reaction constants from file. */
ReactionConstants initialize_reaction_consts(SystemVariables *sys)
{
  ReactionConstants reaction_consts;

  ifstream infile(sys->param_filename);
  string line, param_name;
  double param_value;
  string param_string;
  bool paramset(0);

  //Write read-in parameters to file output_filename.parameters.
  //char Pfilename[128];
  //strcat(Pfilename, "NuEvenNiet");  
  //strcat(Pfilename, ".parameters"); 
  //FILE *Pfp = fopen( Pfilename, "a" );

  //cerr << endl << "PARAMETER VALUES USED IN THIS SIMULATION: " << endl;
  while (getline(infile, line))
  {
    //remove unnecessary spaces in line.
    line = delUnnecessary( line );
  
    istringstream iss(line);
    paramset = 0;
    
    // Check if line has formatting: '<string> <number>'
    if (iss >> param_name >> param_value) 
    {
      /* Loading model parameters */              
      if( param_name.compare("kCIhyd") == 0 )
      {
        reaction_consts.kCIhyd = param_value; 
        paramset = 1;
      }            

      if( param_name.compare("kCIADPoff0") == 0 )
      {
        reaction_consts.kCIADPoff0 = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("kCIADPoffA") == 0 )
      {
        reaction_consts.kCIADPoffA = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("kCIATPoff") == 0 )
      {
        reaction_consts.kCIATPoff = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("kAkIDoff") == 0 )
      {
        reaction_consts.kAkIDoff = param_value; 
        paramset = 1;
      }

      if( param_name.compare("ddGTDconf") == 0 )
      {
        reaction_consts.ddGTDconf = param_value; 
        paramset = 1;
      }
            
      if( param_name.compare("ddGTDAbind") == 0 )
      {
        reaction_consts.ddGTDAbind = param_value; 
        paramset = 1;
      }
            
            
      if( param_name.compare("kconf0") == 0 )
      {
        reaction_consts.kconf0 = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("kBswitch_f") == 0 )
      {
        reaction_consts.kBswitch_f = param_value; 
        paramset = 1;
      }

      if( param_name.compare("nIADPconfref") == 0 )
      {
        reaction_consts.nIADPconfref = param_value; 
        paramset = 1;
      }

      if( param_name.compare("nIADPAref") == 0 )
      {
        reaction_consts.nIADPAref = param_value; 
        paramset = 1;
      }

      if( param_name.compare("dGIbindB") == 0 )
      {
        reaction_consts.dGIbindB = param_value; 
        paramset = 1;
      }

      if( param_name.compare("dgICIU") == 0 )
      {
        reaction_consts.dgICIU = param_value; 
        paramset = 1;
      }

      if( param_name.compare("dgICIT") == 0 )
      {
        reaction_consts.dgICIT = param_value; 
        paramset = 1;
      }

      if( param_name.compare("dgICID") == 0 )
      {
        reaction_consts.dgICID = param_value; 
        paramset = 1;
      }

      if( param_name.compare("dgICIS") == 0 )
      {
        reaction_consts.dgICIS = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("dgACIIU") == 0 )
      {
        reaction_consts.dgACIIU = param_value; 
        paramset = 1;
      }

      if( param_name.compare("dgACIIT") == 0 )
      {
        reaction_consts.dgACIIT = param_value; 
        paramset = 1;
      }

      if( param_name.compare("dgACIID") == 0 )
      {
        reaction_consts.dgACIID = param_value; 
        paramset = 1;
      }

      if( param_name.compare("dgACIIS") == 0 )
      {
        reaction_consts.dgACIIS = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("dGICII") == 0 )
      {
        reaction_consts.dGICII = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("dgACIActU") == 0 )
      {
        reaction_consts.dgACIActU = param_value; 
        paramset = 1;
      }   
      
      if( param_name.compare("dgACIActT") == 0 )
      {
        reaction_consts.dgACIActT = param_value; 
        paramset = 1;
      }  
      
      if( param_name.compare("dgACIActD") == 0 )
      {
        reaction_consts.dgACIActD = param_value; 
        paramset = 1;
      }  
      
      if( param_name.compare("dgACIActS") == 0 )
      {
        reaction_consts.dgACIActS = param_value; 
        paramset = 1;
      }  
      
      if( param_name.compare("dgICIActU") == 0 )
      {
        reaction_consts.dgICIActU = param_value; 
        paramset = 1;
      }   
      
      if( param_name.compare("dgICIActT") == 0 )
      {
        reaction_consts.dgICIActT = param_value; 
        paramset = 1;
      }  
      
      if( param_name.compare("dgICIActD") == 0 )
      {
        reaction_consts.dgICIActD = param_value; 
        paramset = 1;
      }  
      
      if( param_name.compare("dgICIActS") == 0 )
      {
        reaction_consts.dgICIActS = param_value; 
        paramset = 1;
      }                                    

      if( param_name.compare("kICIBon") == 0 )
      {
        reaction_consts.kICIBon = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("kICIBoff") == 0 )
      {
        reaction_consts.kICIBoff = param_value; 
        paramset = 1;
      }         

      if( param_name.compare("kACIBon") == 0 )
      {
        reaction_consts.kACIBon = param_value; 
        paramset = 1;
      }

      if( param_name.compare("kACIBoff") == 0 )
      {
        reaction_consts.kACIBoff = param_value; 
        paramset = 1;
      }
                
      if( param_name.compare("kACIAon") == 0 )
      {
        reaction_consts.kACIAon = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("kACIAoff") == 0 )
      {
        reaction_consts.kACIAoff = param_value; 
        paramset = 1;
      }

      if( param_name.compare("kICIAon") == 0 )
      {
        reaction_consts.kICIAon = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("kICIAoff") == 0 )
      {
        reaction_consts.kICIAoff = param_value; 
        paramset = 1;
      }

      if( param_name.compare("nAseq") == 0 )
      {
        reaction_consts.nAseq = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("nBseq") == 0 )
      {
        reaction_consts.nBseq = param_value; 
        paramset = 1;
      }      
      
      if( param_name.compare("kCIIAon") == 0 )
      {
        reaction_consts.kCIIAon = param_value; 
        paramset = 1;
      }                                                      

      if( param_name.compare("kCIIAoff") == 0 )
      {
        reaction_consts.kCIIAoff = param_value; 
        paramset = 1;
      }                                                      

      if( param_name.compare("kCIIhyd0") == 0 )
      {
        reaction_consts.kCIIhyd0 = param_value; 
        paramset = 1;
      }
      
      if( param_name.compare("kCIIhydA") == 0 )
      {
        reaction_consts.kCIIhydA = param_value; 
        paramset = 1;
      }                                                             

      if( param_name.compare("kCIInucloff0") == 0 )
      {
        reaction_consts.kCIInucloff0 = param_value; 
        paramset = 1;
      }                                                      

      if( param_name.compare("kCIInucloffA") == 0 )
      {
        reaction_consts.kCIInucloffA = param_value; 
        paramset = 1;
      }  
      
      if( param_name.compare("KATPoKADP") == 0 )
      {
        reaction_consts.KATPoKADP = param_value; 
        paramset = 1;
      }      

      if( param_name.compare("kUT") == 0 )
      {
        reaction_consts.kUT = param_value; 
        paramset = 1;
      }  

      if( param_name.compare("kTU") == 0 )
      {
        reaction_consts.kTU = param_value;
        paramset = 1;
      }  

      if( param_name.compare("kTD") == 0 )
      {
        reaction_consts.kTD = param_value;
        paramset = 1;
      }  
      
      if( param_name.compare("kDT") == 0 )
      {
        reaction_consts.kDT = param_value;
        paramset = 1;
      }  
      
      if( param_name.compare("kSD") == 0 )
      {
        reaction_consts.kSD = param_value;
        paramset = 1;
      }  
      
      if( param_name.compare("kDS") == 0 )
      {
        reaction_consts.kDS = param_value;
        paramset = 1;
      }  
      
      if( param_name.compare("kUS") == 0 )
      {
        reaction_consts.kUS = param_value;
        paramset = 1;
      }  
      
      if( param_name.compare("kSU") == 0 )
      {
        reaction_consts.kSU = param_value;
        paramset = 1;
      }  
      
      if( param_name.compare("kTUs") == 0 )
      {
        reaction_consts.kTUs = param_value;
        paramset = 1;
      }  

      if( param_name.compare("kDTs") == 0 )
      {
        reaction_consts.kDTs = param_value;
        paramset = 1;
      }  
           
      if( param_name.compare("kDSs") == 0 )
      {
        reaction_consts.kDSs = param_value;
        paramset = 1;
      }  
      
      if( param_name.compare("kSUs") == 0 )
      {
        reaction_consts.kSUs = param_value;
        paramset = 1;
      }        
      
      if(!paramset)
      {
        //cerr << param_name << " - is not a valid parameter name." << endl;
      }
      else
      {
        //cerr << param_name << " = " << param_value << endl;
        //fprintf( Pfp, "%s\t%e\n", param_name.c_str(), param_value );             
      }
          
    }
    else
    {
      // Check if line has formatting: '<string> <string>'
      iss.clear();
      iss.str( line );
      if (iss >> param_name >> param_string)
      {                
         if(!paramset)
         {
           //cerr << param_name << " - is not a valid parameter name." << endl;
         }
         else
         {
           //cerr << param_name << " = " << param_string << endl;
           //fprintf( Pfp, "%s\t%s\n", param_name.c_str(), param_string.c_str() );                
         }              
      }
      else
      { 
        //Check for comment and empty lines
        if(!(param_name.compare(0,1,"#") == 0) && !line.empty())
        {
           cerr << "Error - Line with incorrect formatting: " << param_name << endl; 
           exit(1); 
        }   
        continue;
      }
      

    }
  }//while 

  //fclose(Pfp);    
  return reaction_consts;  
}

string delUnnecessary(string &str)
{
    int size = str.length();
    for(int j = 0; j<=size; j++)
    {
        for(int i = 0; i <=j; i++)
        {
            if(str[i] == ' ' && str[i+1] == ' ')
            {
                str.erase(str.begin() + i);
            }
            else if(str[0]== ' ')
            {
                str.erase(str.begin());
            }
            else if(str[i] == '\0' && str[i-1]== ' ')
            {
                str.erase(str.end() - 1);
            }
        }
    }
    return str;
}
