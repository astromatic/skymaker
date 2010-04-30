/*
				stars.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SkyMaker
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Routines for creating star fields.
*
*	Last modify:	30/04/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "globals.h"
#include "fft.h"
#include "image.h"
#include "list.h"
#include "prefs.h"
#include "psf.h"
#include "random.h"
#include "simul.h"
#include "stars.h"

#ifdef USE_THREADS
#include "threads.h"

static void	pthread_makestarfield(simstruct *sim, int nstar),
		*pthread_makeobj(void *arg);
static int	pthread_nextobj(int objindex, int proc);

static pthread_t		*thread;
static pthread_mutex_t		objmutex, imagemutex;
static pthread_cond_t		*objcond;

static simstruct	*pthread_sim;
static objstruct	*pthread_obj;
static char		pthread_str[MAXCHAR];
static double		pthread_xrange,pthread_yrange;
static int		*pthread_addobjflag, *pthread_objqueue,
			pthread_nstar,pthread_nstarmax,
			pthread_addobji,pthread_procobji, pthread_nobj,
			pthread_step, pthread_gridoffsetx,pthread_gridoffsety,
			pthread_ngridx, pthread_gridstep, pthread_gridindex,
			pthread_ngrid, pthread_gridflag;
#endif

#define	EPS_RANDFLUX	1e-8	/* = dynamic range of 20 mag in stellar flux */

/******************************** makestarfield ******************************/
/*
Render a star field.
*/
void	makestarfield(simstruct *sim)

  {
   char			str[MAXCHAR];
   double		xrange,yrange, dnstars;
   int			nstars, step, gridflag;
#ifndef USE_THREADS
   static objstruct	obj={0};
   int			i, gridindex, ngrid;
#endif

  init_random(sim->starposseed);
  gridflag = (sim->imatype==GRID || sim->imatype==GRID_NONOISE);
  xrange = sim->fimasize[0]-1.0;
  yrange = sim->fimasize[1]-1.0;
  dnstars = sim->pixscale[0]*sim->pixscale[1] / (sim->mscan[0]*sim->mscan[1])
		*xrange*yrange*sim->scountdens*(ARCSEC*ARCSEC)/(DEG*DEG);
  if (dnstars<2e12)
    nstars = (int)random_poisson(dnstars, 0);
  else
    {
    nstars = 0;			/* To avoid gcc -Wall warnings */
    error(EXIT_FAILURE, "Sorry, too many stars: it would take DAYS to ",
	"render, even on a super-fast computer!");
    }
  sprintf(str, "Adding %d stars...", nstars);
  NFPRINTF(OUTPUT, str);
/*- rough estimate of computing time necessary between each announcement */
  step = (int)(5000000.0*sim->psfoversamp*sim->psfoversamp)
	/(sim->psfsize[0]*sim->psfsize[1]*sim->mscan[0]*sim->mscan[1]);

#ifdef USE_THREADS
  pthread_makestarfield(sim, nstars);
#else
/* Let's initialize some object parameters */
  obj.flux = 0.0;
  obj.type = 100;

  if (gridflag)
    {
    gridindex = sim->gridindex;
    ngrid = sim->ngrid[0]*sim->ngrid[1];
    }
  else
    gridindex = ngrid = 0;		/* To avoid gcc -Wall warnings */

  for (i=0; i<nstars; i++)
    {
/*-- Stop (in grid mode) once objects fill the grid */
    if (gridflag && gridindex>=ngrid)
      break;

    if (!(i%step))
      {
      sprintf(str, "Adding %d stars... (%2.0f%% done)",nstars, 100.0*i/nstars);
      NFPRINTF(OUTPUT, str);
      }
    if (!gridflag)
      {
      obj.x = xrange*random_double(0)-sim->margin[0];
      obj.y = yrange*random_double(0)-sim->margin[1];
      }
    obj.mag = sim->maglim[1] + log10(random_double(0)+EPS_RANDFLUX)
		/ sim->scountslope;
    if (obj.mag > sim->maglim[0])
      {
      if (gridflag)
        {
        obj.x = gridindex%sim->ngrid[0] + sim->gridoffset[0]+random_double(0)-0.5;
        obj.y = gridindex/sim->ngrid[0] + sim->gridoffset[1]+random_double(0)-0.5;
        gridindex++;
        }
      make_star(sim, &obj);
      add_image(obj.subimage, obj.subsize[0], obj.subsize[1],
	sim->image, sim->fimasize[0], sim->fimasize[1],
	obj.subx, obj.suby, (float)obj.subfactor);
       obj.flux = 0.0;
      writeobj(sim, &obj);
      }
    }

  endobj(&obj);
  sim->gridindex = gridindex;

#endif

  return;
  }


/****** make_star *******************************************************
PROTO	void make_star(simstruct *sim, objstruct *obj)
PURPOSE	Generate a stellar (integrated) profile.
INPUT	Pointer to the sim structure,
	pointer to the object structure.
OUTPUT	-.
NOTES	Writes to an allocated image buffer, not directly to the image to
	allow multithreading.
AUTHOR	E. Bertin (IAP)
VERSION	30/04/2010
 ***/
void	make_star(simstruct *sim, objstruct *obj)

  {
   PIXTYPE	*subt, *psf;
   double	dx,dy,flux, osamp;
   int		i, nsub2, nsubo;

/* Convert magnitude to linear units */ 
  if (obj->flux)
    flux = (float)obj->flux;
  else
    obj->flux = (float)(flux = DEXP(0.4*(sim->magzero2-obj->mag)));
  osamp = sim->psfoversamp;
  dx = (obj->x + sim->margin[0]-sim->dpsfc[0])/sim->mscan[0];
  dy = (obj->y + sim->margin[1]-sim->dpsfc[1])/sim->mscan[1];
  dx -= (double)(obj->subx = (int)(dx+0.49999));
  dy -= (double)(obj->suby = (int)(dy+0.49999));
/* Resample to lower resolution */
  nsubo = obj->subsize[0]*obj->subsize[1];
  obj->subsize[0] = obj->subsize[1] = (int)(sim->psfsize[0]/osamp+1.0);
  nsub2 = obj->subsize[0]*obj->subsize[1];
  if (!nsubo)
    {
    QMALLOC(obj->subimage, PIXTYPE, nsub2);
    }
  else if (nsub2>nsubo)
    {
    QREALLOC(obj->subimage, PIXTYPE, nsub2);
    }
  psf = interp_psf(sim, obj->x, obj->y);
  resample_image(psf, sim->psfsize[0], sim->psfsize[1], obj,
	-dx*osamp, -dy*osamp, osamp);
  if (sim->npsf>1)
    free(psf);
  flux = 0.0;
  for (i=nsub2,subt=obj->subimage; i--;)
    flux += (double)*(subt++);

  obj->subfactor = obj->flux/flux;

  return;
  }

#ifdef USE_THREADS

/****** pthread_makeobj ******************************************************
PROTO	void *pthread_makeobj(void *arg)
PURPOSE	Object generation thread.
INPUT	Pointer to the thread number.
OUTPUT	NULL void pointer.
NOTES	Relies on global variables.
AUTHOR	E. Bertin (IAP)
VERSION	17/08/2006
 ***/
static void	*pthread_makeobj(void *arg)
  {
   objstruct	*obj;
   int		obji , proc;

  proc = *((int *)arg);
  obji = -1;
/* Exit if the end of file is reached by any thread, including this one */
  while ((obji=pthread_nextobj(obji, proc))!=-1)
    {
    obj = &pthread_obj[obji];
    if (obj->mag > pthread_sim->maglim[0])
      make_star(pthread_sim, obj);
    }

  pthread_exit(NULL);
  return (void *)NULL;
  }


/****** pthread_nextobj ******************************************************
PROTO	int pthread_nextobj(int objindex, int proc)
PURPOSE	Manage object list generation and thread synchronisation.
INPUT	Previous object index in list,
	Process index.
OUTPUT	Next object index in list.
NOTES	Relies on global variables.
AUTHOR	E. Bertin (IAP)
VERSION	28/09/2007
 ***/
static int	pthread_nextobj(int obji, int proc)
  {
   objstruct	*obj;
   int		q;

  QPTHREAD_MUTEX_LOCK(&objmutex);
/* The newly processed object is ready to be added to the image */
  if (obji>=0)
    pthread_addobjflag[obji] = 2;
/* If we just finished the "right" object, add it to image! */
  if (obji == pthread_addobji)
    {
    while (pthread_addobjflag[pthread_addobji]==2)
      {
      obj = &pthread_obj[pthread_addobji];
/*---- Add the object to the image */
      if (obj->subimage && obj->mag > pthread_sim->maglim[0])
        {
        QPTHREAD_MUTEX_LOCK(&imagemutex);
        add_image(obj->subimage, obj->subsize[0], obj->subsize[1],
		pthread_sim->image, pthread_sim->fimasize[0],
		pthread_sim->fimasize[1],
		obj->subx, obj->suby, (float)obj->subfactor);
        QPTHREAD_MUTEX_UNLOCK(&imagemutex);
/*------ Add the object to the output list */
        obj->flux = 0.0;
        writeobj(pthread_sim, obj);
        }
      pthread_addobjflag[pthread_addobji] = 0;
      QPTHREAD_COND_BROADCAST(&objcond[pthread_addobji]);
      pthread_addobji = (pthread_addobji+1)%pthread_nobj;
      }
    }
/* If no more object to process, return a "-1" (meaning exit thread) */
  if (pthread_nstar++ < pthread_nstarmax
	&& (!pthread_gridflag || pthread_gridindex<pthread_ngrid))
    {
    if (!(pthread_nstar%pthread_step))
      {
      sprintf(pthread_str, "Adding %d stars... (%2.0f%% done)",
		pthread_nstarmax, 100.0*pthread_nstar/pthread_nstarmax);
      NFPRINTF(OUTPUT, pthread_str);
      }
    obji = pthread_procobji;
    pthread_procobji = (pthread_procobji+1)%pthread_nobj;
/*-- If the next available buffer has not been flushed yet, wait */
    q=++pthread_objqueue[obji];
    while (pthread_addobjflag[obji] || --q)
      QPTHREAD_COND_WAIT(&objcond[obji], &objmutex);
/*-- Set content */
    obj = &pthread_obj[obji];
    obj->type = 100;
    obj->mag = pthread_sim->maglim[1] + log10(random_double(proc)+EPS_RANDFLUX)
		/ pthread_sim->scountslope;
    if (obj->mag > pthread_sim->maglim[0])
      {
      if (pthread_gridflag)
        {
        obj->x = (pthread_gridindex%pthread_ngridx)*pthread_gridstep
		+ pthread_gridoffsetx + random_double(proc)-0.5;
        obj->y = (pthread_gridindex/pthread_ngridx)*pthread_gridstep
		+ pthread_gridoffsety + random_double(proc)-0.5;
        pthread_gridindex++;
        }
      else
        {
        obj->x = pthread_xrange*random_double(proc)-pthread_sim->margin[0];
        obj->y = pthread_yrange*random_double(proc)-pthread_sim->margin[1];
        }
      }
    if (pthread_objqueue[obji])
      pthread_objqueue[obji]--;
    pthread_addobjflag[obji] = 1;
    }
  else
    obji=-1;

  QPTHREAD_MUTEX_UNLOCK(&objmutex);

  return obji;
  }


/****** pthread_makestarfield ***************************************************
PROTO	void pthread_makestarfield(simstruct *sim, int nstars)
PURPOSE	Generate a star field using multithreads.
INPUT	Pointer to the sim structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	28/09/2007
 ***/
static void	pthread_makestarfield(simstruct *sim, int nstars)
  {
   static pthread_attr_t        pthread_attr;
   int                          *proc,
                                o, p;

  pthread_step = (int)(5000000.0*sim->psfoversamp*sim->psfoversamp)
	/(sim->psfsize[0]*sim->psfsize[1]*sim->mscan[0]*sim->mscan[1]);
/* Number of active threads */
  nproc = prefs.nthreads;
  pthread_nobj = 2*nproc;	/* A margin of 2X for better efficiency */
  QCALLOC(pthread_obj, objstruct, pthread_nobj);
  QCALLOC(pthread_addobjflag, int, pthread_nobj);
  QCALLOC(pthread_objqueue, int, pthread_nobj);
/* Set up multi-threading stuff */
  QMALLOC(objcond, pthread_cond_t, pthread_nobj);
  for (o=0; o<pthread_nobj; o++)
    {
    QPTHREAD_COND_INIT(&objcond[o], NULL);
    }
  QPTHREAD_MUTEX_INIT(&objmutex, NULL);
  QPTHREAD_MUTEX_INIT(&imagemutex, NULL);
  QMALLOC(proc, int, nproc);
  QMALLOC(thread, pthread_t, nproc);
  QPTHREAD_ATTR_INIT(&pthread_attr);
  QPTHREAD_ATTR_SETDETACHSTATE(&pthread_attr, PTHREAD_CREATE_JOINABLE);
  pthread_sim = sim;
  pthread_nstar = 0;
  pthread_nstarmax = nstars;
  pthread_gridflag = (sim->imatype==GRID || sim->imatype==GRID_NONOISE);
  if (pthread_gridflag)
    {
    pthread_ngridx = sim->ngrid[0];
    pthread_ngrid = sim->ngrid[0]*sim->ngrid[1];
    pthread_gridoffsetx = sim->gridoffset[0];
    pthread_gridoffsety = sim->gridoffset[1];
    pthread_gridstep = sim->gridstep;
    pthread_gridindex = sim->gridindex;
    }
  else
    {
    pthread_xrange = sim->fimasize[0]-1.0;
    pthread_yrange = sim->fimasize[1]-1.0;
    }
  pthread_addobji = pthread_procobji = 0;
/* Start the reading/generation threads */
  for (p=0; p<nproc; p++)
    {
    proc[p] = p;
    QPTHREAD_CREATE(&thread[p], &pthread_attr, &pthread_makeobj, &proc[p]);
    }
  for (p=0; p<nproc; p++)
    QPTHREAD_JOIN(thread[p], NULL);
/* Clean up multi-threading stuff */
  QPTHREAD_MUTEX_DESTROY(&objmutex);
  QPTHREAD_MUTEX_DESTROY(&imagemutex);
  QPTHREAD_ATTR_DESTROY(&pthread_attr);
  for (o=0; o<pthread_nobj; o++)
    {
    QPTHREAD_COND_DESTROY(&objcond[o]);
    endobj(&pthread_obj[o]);
    }

  if (pthread_gridflag)
    sim->gridindex = pthread_gridindex;

  free(objcond);
  free(pthread_obj);
  free(pthread_addobjflag);
  free(pthread_objqueue);
  free(proc);
  free(thread);

  return;
  }

#endif

