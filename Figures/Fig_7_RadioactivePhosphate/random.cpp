/* Random Number Generator as per Numerical Recipes */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "random.hpp"
#include "main.hpp"

int ran_iff = FALSE;  ///< global variable stating if the random generator is initialized
int ran_seed = MSEED; ///< global variable holding the seed for the random generator


/**
     Returns a seed for any random number generator.
     @return Number to be used as a seed
     @author
     Robert G. Brown                        http://www.phy.duke.edu/~rgb/
     Duke University Dept. of Physics, Box 90305
     Durham, N.C. 27708-0305
     Phone: 1-919-660-2567  Fax: 919-660-2525     email:rgb@phy.duke.edu
 */
unsigned long int random_seed()
{
  unsigned int seed;
  struct timeval tv;
  FILE *devrandom;

  // can't use /dev/random all of the times since a server can't collect enought entropy and the program would halt in this case
  if( (devrandom = fopen( "/dev/urandom", "r" )) == NULL )
  {
    gettimeofday( &tv, 0 );
    seed = tv.tv_sec + tv.tv_usec;
//     fprintf( stderr, "Got seed %u from gettimeofday()\n", seed );
  } 
  else 
  {
    fread( &seed, sizeof(seed), 1, devrandom );
    fclose( devrandom );
//     fprintf( stderr, "Got seed %u from /dev/random\n", seed );
  }

  return seed;
}


/**
  Initialize the random number generator with a given seed.
  If the seed is set to 0, it is determined by using the current time.
  @param[in]  seed  Integer used as a seed for the generator
*/
void ran_init( long seed )
{
  
  if( 0 == seed )
  {
    ran_seed = random_seed();

#ifdef __alpha__    /* correct for 64bit time values */
    ran_seed %= MSEED;
#endif
  }
  else
  {
    ran_seed = seed;
  }
  
  fprintf( stderr, "The seed for the random number generator has been set to: %d\n", ran_seed );
  
  ran_iff = FALSE;
}


#define MBIG 1000000000
#define FAC  (1.0/MBIG)
#define MZ   0
/**
  Function returns a random number of a uniform distribution of the range (0,1].
  The generator needs no initialization, but will always return the same sequence, if
  it is not initialized by ran_init().
  The algorithm is the ran3 algorithm from numerical recipies.
  @return The random number out of (0,1]
*/
double ran_get()
{
  static int inext, inextp;
  static long ma[56];
  int i, ii, k;
  long mj, mk;
  
  if( ran_iff == FALSE )
  {
    // random number generator must be initialized
    mj = labs(ran_seed); // seed has to be positive to ensure ran_get() being in [0,1)
    mj = mj % MBIG;
    ma[55] = mj;
    mk = 1;
    for( i=1; i<55; i++ )
    {
      ii = (21*i) % 55;
      ma[ii] = mk;
      mk = mj - mk;
      if( mk < MZ ) 
        mk = mk + MBIG;
      mj = ma[ii];
    }
    for( k=1; k<=4; k++ )
    {
      for( i=1; i<=55; i++ )
      {
        ma[i] = ma[i] - ma[1 + ((i+30)%55)];
        if( ma[i] < MZ )
          ma[i] = ma[i] + MBIG;
      }
    }
    inext=0;
    inextp=31;
    ran_iff = TRUE; // random number generator is now initialized
  }
  
  // algorithm for creating random numbers
  if( ++inext == 56 )
    inext = 1;
  if( ++inextp == 56 )
    inextp = 1;
  mj = ma[inext] - ma[inextp];
  if( mj < MZ )
    mj = mj + MBIG;
  ma[inext] = mj;
  return 1. - mj*FAC;
}
#undef MBIG
#undef FAC
#undef MZ


/* Returns 1 with probability p, and 0 with prob 1-p. */
bool ran_bool(double p)
{
  if( ran_get() <= p )
    return 1;
  else
    return 0;
}


/**
  Generate gaussian distributed random numbers using the Box-Muller algorithm.
  @param[in]  sigma   The widht of the Gaussian
  @return the random number
 */
double ran_gaussian( const double sigma )
{
  // use static variables to speed up process, since always two random numb
  static int stored=FALSE;  // state variable keeping track of the number of stored random numbers in memory
  static double y;
  double x1, x2, w;
 
  // check, if there is a number in store
  if( ! stored )
  {
    // create two uniformly distributed random numbers in a circle of radius 1
    do 
    {
      x1 = 2.0 * ran_get() - 1.0;
      x2 = 2.0 * ran_get() - 1.0;
      w = x1*x1 + x2*x2;
    } while ( w >= 1.0 );
    
    w = sqrt( -2.0*log(w)/w );
    
    // store the first number
    y = x1 * w;
    stored = TRUE;
    
    // return the second number
    return x2 * w * sigma;
  }
  else
  {
    // return the number, which has been stored
    stored = FALSE;
    return y * sigma;
  }
}


/**
  Returns a random number of a binomial distribution
  @param[in]  pp  Success probability parameter 0 <= p <=1
  @param[in]  n   Number of trials parameter
  @return The random number
*/
int ran_binomial( double pp, int n )
{
  int j;
  static int nold=(-1);
  double am,em,g,angle,p,bnl,sq,t,y;
  static double pold=(-1.0),pc,plog,pclog,en,oldg;

  p=(pp <= 0.5 ? pp : 1.0-pp);
  am=n*p;
  if( n < 25 ) 
  {
    bnl=0.0;
    for( j=1; j<=n; j++ )
      if( ran_get() < p ) 
        ++bnl;
  } 
  else if( am < 1.0 ) 
  {
    g=exp(-am);
    t=1.0;
    for( j=0; j<=n; j++ ) 
    {
      t *= ran_get();
      if (t < g) 
        break;
    }
    bnl=(j <= n ? j : n);
  } 
  else 
  {
    if( n != nold ) 
    {
      en=n;
      oldg=gammln(en+1.0);
      nold=n;
    } 
    if( p != pold ) 
    {
      pc=1.0-p;
      plog=log(p);
      pclog=log(pc);
      pold=p;
    }
    sq=sqrt(2.0*am*pc);
    do 
    {
      do 
      {
        angle=M_PI*ran_get();
        y=tan(angle);
        em=sq*y+am;
      } while( em < 0.0 || em >= (en+1.0) );
      em=floor(em);
      t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
          -gammln(en-em+1.0)+em*plog+(en-em)*pclog);
    } while( ran_get() > t );
    bnl=em;
  }
  if (p != pp) 
    bnl=n-bnl;
  return (int) bnl;
}


/**
  Computes ln(gamma(xx)).
  @param[in]  xx   The argument to be evaluated
  @return ln(gamma(xx))
*/
double gammln( double xx )
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
    24.01409824083091,-1.231739572450155,
    0.1208650973866179e-2,-0.5395239384953e-5};
  int j;

  y = x = xx;
  tmp = x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser = 1.000000000190015;
  for( j=0; j<=5; j++ )
    ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}
