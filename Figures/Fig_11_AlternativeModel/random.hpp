#ifndef _RANDOM_H_
#define _RANDOM_H_
// Header file for the random number generator

// default seed, if random number generator is not initialized
#define MSEED 161803398

// functions
unsigned long int random_seed();
void ran_init( long seed );
double ran_get();
double ran_gaussian( const double sigma );
int ran_binomial( double pp, int n );
double gammln( double xx );
bool ran_bool(double p);

#endif
