/*
 				random.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	a program that uses randomly generated numbers
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Definitions related to the generation of random numbers
*
*	Last modify:	17/08/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*--------------------------------- constants -------------------------------*/

#ifndef	RAND_MAX
#define	RAND_MAX	0x7fffffffL	/* default dynamic range of rand() */
#endif

/*-------------------------------- protos -----------------------------------*/

double	random_double(int p),
	random_gauss(double sigma, int p),
	random_poisson(double xm, int p);

int	random_int(int p);

void	init_random(int seed);


