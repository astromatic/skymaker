/*
                                  random.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*       Part of:        a program that uses randomly generated numbers
*
*       Author:         E.BERTIN (IAP)
*
*       Contents:       functions returning random numbers.
*
*       Last modify:    23/08/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "define.h"
#include "globals.h"
#include "random.h"

#ifndef	THREADS_NMAX
#define	THREAD_NMAX 16
#endif

#ifdef USE_THREADS
static unsigned int	glob_seed[THREADS_NMAX];
#endif

/****** random_gauss ********************************************************
PROTO   double random_gauss(double sigma)
PURPOSE Generate a random number with a (centered) Gaussian pdf.
INPUT   Standard deviation.
OUTPUT  Gauss-distributed random number.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 17/08/2006
*/
double	random_gauss(double sigma, int p)
  {
   double	x,y,z, r;

  while((z=pow(x=random_double(p)-0.5,2.0) + pow(y=random_double(p)-0.5,2.0))
	> 0.25);
  while ((r=random_double(p)) <= 0.0);

  return sigma*sqrt(-2.0*log(r)/z)*x;
  }


/****** random_int **********************************************************
PROTO   int random_int(void)
PURPOSE Generate a random integer over the range [0,RAND_MAX].
INPUT   Process index (ignored for single threads).
OUTPUT  Random integer number with uniform distribution.
NOTES   The actual upper bound of the range is implementation-dependent.
AUTHOR  E. Bertin (IAP)
VERSION 17/08/2006
*/
int	random_int(int p)
  {
#ifdef USE_THREADS
  return (int)rand_r(&glob_seed[p]);
#else
  return (int)rand();
#endif
  }


/****** random_double ********************************************************
PROTO   double random_double(void)
PURPOSE Generate a random number with uniform distribution over [0.0,1.0[
INPUT   Process index (ignored for single threads).
OUTPUT  Random double with uniform distribution.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 17/08/2006
*/
double	random_double(int p)
  {
#ifdef USE_THREADS
  return (double)rand_r(&glob_seed[p]) / RAND_MAX;
#else
  return (double)rand() / RAND_MAX;
#endif
  }


/****** init_random **********************************************************
PROTO	void init_random(int seed)
PURPOSE	Initialize the random number generator.
INPUT	Seed.
OUTPUT	-.
NOTES	The seed is used to initialize the random sequence at a particular
	position, which is implementation-dependent. If seed = 0, then the
	actual seed is taken from the time() function (which varies each
	second). For multithreaded code, a series of secondary seeds is
	generated.
AUTHOR	E. Bertin (IAP)
VERSION	17/08/2006
*/
void	init_random(int seed)
  {
#ifdef USE_THREADS
   int	p;
#endif

  if (seed)
    srand((unsigned int)seed);
  else
    srand((unsigned int)time(NULL));
#ifdef USE_THREADS
  for (p=0; p<THREADS_NMAX; p++)
    glob_seed[p] = (unsigned int)rand();
#endif

  return;
  }


/****i* gammln ***************************************************************
PROTO   double gammln(double xx)
PURPOSE Returns the log of the Gamma function (from Num. Recipes in C, p.168).
INPUT   A double.
OUTPUT  Log of the Gamma function.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 29/10/97
*/
static double	gammln(double xx)

  {
   double		x,tmp,ser;
   static double	cof[6]={76.18009173,-86.50532033,24.01409822,
			-1.231739516,0.120858003e-2,-0.536382e-5};
   int			j;

  tmp=(x=xx-1.0)+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.0;
  for (j=0;j<6;j++)
    ser += cof[j]/(x+=1.0);

  return log(2.50662827465*ser)-tmp;
  }


/****** random_poisson *******************************************************
PROTO   double random_poisson(double xm, int p)
PURPOSE Returns a random number with Poisson deviate (from Num. Recipes in C.,
        p.222) centered on xm.
INPUT   Mean of the Poisson distribution,
	process index (ignored for single threads).
OUTPUT  A double containing the integer (!) variable with Poisson deviate.
NOTES   I am still searching for a faster algorithm!!
AUTHOR  E. Bertin (IAP)
VERSION 23/08/2006
*/
double	random_poisson(double xm, int p)
  {
   double		gammln();
   double		sq,alxm,g,oldm,em,t,y;

  sq = alxm = g = 0.0;
  oldm = -1.0;
  if (xm < 12.0)
    {
    if (xm != oldm)
      {
      oldm=xm;
      g=exp(-xm);
      }
    em = -1.0;
    t=1.0;
    do
      {
      em += 1.0;
      t *= random_double(p);
      } while (t > g);
    }
  else
    {
    if (xm != oldm)
      {
      oldm=xm;
      sq=sqrt(2.0*xm);
      alxm=log(xm);
      g=xm*alxm-gammln(xm+1.0);
      }
    do
      {
      do
        {
        y=tan(PI*random_double(p));
        em=sq*y+xm;
        } while (em < 0.0);
      em=floor(em);
      t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
      } while (random_double(p) > t);
    }

  return em;
  }

