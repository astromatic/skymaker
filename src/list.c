/*
*				list.c
*
* Manage source lists.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 2003-2020 IAP/CNRS/SorbonneU
*
*	License:		GNU General Public License
*
*	SkyMaker is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*	SkyMaker is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with SkyMaker. If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		01/12/2020
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "globals.h"
#include "fft.h"
#include "fitswcs.h"
#include "galaxies.h"
#include "image.h"
#include "list.h"
#include "prefs.h"
#include "rasters.h"
#include "simul.h"
#include "stars.h"
#include "random.h"

#ifdef USE_THREADS
#include "threads.h"

static void		pthread_list_read(simstruct *sim),
			*pthread_list_readobj(void *arg);
static int		pthread_list_nextobj(int objindex, char *str, int nproc);

static pthread_t	*thread;
static pthread_mutex_t	objmutex, imagemutex;
static pthread_cond_t	*objcond;

static simstruct	*pthread_sim;
static objstruct	*pthread_obj;
static char		pthread_str[MAXCHAR];
static int		*pthread_addobjflag, *pthread_objqueue, pthread_endflag,
			pthread_addobji,pthread_procobji, pthread_nobj,
			pthread_objn, pthread_objnp,
			pthread_gridoffsetx,pthread_gridoffsety,
			pthread_ngridx, pthread_gridstep, pthread_gridindex,
			pthread_ngrid, pthread_gridflag;
#endif

static double		list_strtofloat(char *str, char **buf, double def);
static int		list_strtoint(char *str, char **buf, int def);

/****** list_read ************************************************************
PROTO	void readlist(simstruct *sim)
PURPOSE	Read the object list.
INPUT	Pointer to the sim structure.
OUTPUT	-.
NOTES	Global prefs variables are used.
AUTHOR	E. Bertin (IAP)
VERSION	09/03/2020
 ***/
void    list_read(simstruct *sim) {

#ifndef USE_THREADS
   objstruct		obj = {0};
   char			str[MAXCHAR], msg[MAXCHAR];
   int			i, gridindex, gridflag, ngrid;
#endif

/* First test if a list-filename has been provided in the command line */
  if (!(*(sim->inlistname)))
    return;

  NFPRINTF(OUTPUT, "Reading input list...");

  if ((sim->inlistfile = fopen(sim->inlistname,"r")) == NULL)
    error(EXIT_FAILURE,"*ERROR*: can't read ", sim->inlistname);

  fft_init(1);

#ifdef USE_THREADS
  pthread_list_read(sim);
#else
  obj.flux = 0.0;
  gridflag = (sim->imatype==GRID || sim->imatype==GRID_NONOISE);
  if (gridflag) {
    gridindex = sim->gridindex;
    ngrid = sim->ngrid[0]*sim->ngrid[1];
  } else
    ngrid = gridindex = 0;		/* To avoid gcc -Wall warnings */

  for (i=0; fgets(str, MAXCHAR, sim->inlistfile); i++) {
    if (!(i%READLIST_DISPSTEP)) {
      sprintf(msg, "Painting input list... (%d objects painted / %d read)",
		prefs.nobj ,i);
      NFPRINTF(OUTPUT, msg);
    }
    if (gridflag && gridindex>=ngrid)
      break;
    if (list_readobj(sim, &obj, str, 0) == RETURN_ERROR)
      continue;
    if (obj.mag < prefs.listmag_limits[0] || obj.mag > prefs.listmag_limits[1])
      continue;
    if (gridflag) {
      obj.pos[0] = gridindex%sim->ngrid[0] + sim->gridoffset[0]
			+ random_double(0) - 0.5;
      obj.pos[1] = gridindex/sim->ngrid[0] + sim->gridoffset[1]
			+ random_double(0) - 0.5;
      gridindex++;
    }
    list_make_obj(sim, &obj);
/*-- Try adding the object to the image (and continue if out of frame) */
    if ((add_image(obj.subimage, obj.subsize[0], obj.subsize[1],
	sim->image, sim->fimasize[0], sim->fimasize[1],
	obj.subpos[0], obj.subpos[1], (float)obj.subfactor))
      continue;
/*-- Add the object to the output list */
    list_writeobj(sim, &obj);
  }

  sim->gridindex = gridindex;
#endif

  fclose(sim->inlistfile);
  fft_end(1);

  return;
}


/****** list_readobj ********************************************************
PROTO	int list_readobj(simstruct *sim, objstruct *obj, char *str, int proc)
PURPOSE	Read the data an object.
INPUT	Pointer to the sim structure,
	pointer to the obj structure,
	character string,
	process index (ignored if not from a thread).
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	17/05/2020
 ***/
int	list_readobj(simstruct *sim, objstruct *obj, char *str, int proc) {
   char		*cptr, *cptr2, *strtokbuf;
   int		gridflag;

  gridflag = (sim->imatype==GRID || sim->imatype==GRID_NONOISE);

/*-- Examine current input line (discard empty and comment lines) */
  if (!*str || strchr("#\t\n",*str))
    return RETURN_ERROR;
  if (!(cptr=strtok_r(str, " \t", &strtokbuf)))
    return RETURN_ERROR;
  obj->type = atoi(cptr);
  obj->flux = 0.0;
/* Note: in FITS, the center of the first pixel is at position (1.0,1.0) */
  if (!(cptr=strtok_r(NULL, " \t", &strtokbuf)))
    return RETURN_ERROR;
  if (!(cptr2=strtok_r(NULL, " \t", &strtokbuf)))
    return RETURN_ERROR;
  if (!gridflag) {
    obj->wcspos[0] = obj->pos[0]  = atof(cptr);
    obj->wcspos[1] = obj->pos[1]  = atof(cptr2);
    if (sim->wcsflag && sim->wcs)
/*---- Convert World coordinates to pixel coordinates */
      wcs_to_raw(sim->wcs, obj->wcspos, obj->pos);
  }

  if (!(cptr=strtok_r(NULL, " \t", &strtokbuf)))
    return RETURN_ERROR;
  obj->mag = atof(cptr);
  if (obj->mag < prefs.listmag_limits[0] || obj->mag > prefs.listmag_limits[1])
    return RETURN_ERROR;
  switch(obj->type) {
    case 200:
    case 210:
      obj->bulge_ratio = list_strtofloat(NULL, &strtokbuf, 1.0);
      obj->bulge_sersicn = (obj->type == 210) ?
	list_strtofloat(NULL, &strtokbuf, 4.0) : 4.0;
      obj->bulge_req = list_strtofloat(NULL, &strtokbuf, 1.0);
      obj->bulge_aspect = list_strtofloat(NULL, &strtokbuf, 1.0);
      obj->bulge_posang = list_strtofloat(NULL, &strtokbuf,
	random_double(proc)*360.0 - 180.0);
      obj->disk_scale = list_strtofloat(NULL, &strtokbuf, 1.0);
      obj->disk_aspect = list_strtofloat(NULL, &strtokbuf, 1.0);
      obj->disk_posang = list_strtofloat(NULL, &strtokbuf,
	random_double(proc)*360.0 - 180.0);
      obj->z = list_strtofloat(NULL, &strtokbuf, 0.0);
      obj->hubble_type =  list_strtofloat(NULL, &strtokbuf, 0.0);
      break;
    case 300:
      obj->raster_index =  list_strtoint(NULL, &strtokbuf, 0);
      obj->raster_size = list_strtofloat(NULL, &strtokbuf, 1.0);
      obj->raster_aspect = list_strtofloat(NULL, &strtokbuf, 1.0);
      obj->raster_posang = list_strtofloat(NULL, &strtokbuf,
				random_double(proc)*360.0 - 180.0);
      obj->z = list_strtofloat(NULL, &strtokbuf, 0.0);
      break;
    default:
      break;
  }

  return RETURN_OK;
}


/****** list_strtofloat *******************************************************
PROTO	double list_strtofloat(char *str, char *buf, double default)
PURPOSE	Write the data an object.
INPUT	Pointer to the sim structure,
	pointer to the obj structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	09/03/2020
 ***/
static double    list_strtofloat(char *str, char **buf, double def) {

   char	*cptr;

  return (cptr=strtok_r(str, " \t", buf))? atof(cptr) : def;
}


/****** list_strtoint ********************************************************
PROTO	int list_strtoint(char *str, char *buf, int default)
PURPOSE	Write the data an object.
INPUT	Pointer to the sim structure,
	pointer to the obj structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	18/05/2020
 ***/
static int    list_strtoint(char *str, char **buf, int def) {

   char	*cptr;

  return (cptr=strtok_r(str, " \t", buf))? atoi(cptr) : def;
}


/****** list_makeobj *********************************************************
PROTO	int list_makeobj(simstruct *sim, objstruct *obj)
PURPOSE	Switch to the right renderer based on object type.
INPUT	Pointer to the sim structure,
	pointer to the obj structure.
OUTPUT	RETURN_OK if the object was succesfully rasterized or RETURN_ERROR
	otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	05/05/2020
 ***/
int	list_makeobj(simstruct *sim, objstruct *obj) {

  int ret;

  switch(obj->type) {
    case 100:
      ret = make_star(sim, obj);
      break;
    case 200:
    case 210:
      ret = make_galaxy(sim, obj);
      break;
    case 300:
    case 310:
      ret = make_raster(sim, obj);
      break;
    default:
      return RETURN_ERROR;
  }

  return ret;
}


/****** list_writeobj ********************************************************
PROTO	void list_writeobj(simstruct *sim, objstruct *obj)
PURPOSE	Write the data an object.
INPUT	Pointer to the sim structure,
	pointer to the obj structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	01/12/2020
 ***/
void    list_writeobj(simstruct *sim, objstruct *obj) {
   char			str[MAXCHAR];

 /*-- The format depends on object type */
  switch(obj->type) {
    case 100:
      fprintf(sim->outlistfile, "%3d %11.4f %11.4f %8.4f\n",
	obj->type, obj->pos[0], obj->pos[1], obj->mag);
      break;
    case 200:
      fprintf(sim->outlistfile, "%3d %11.4f %11.4f %8.4f %5.3f %9.3f %5.3f "
		"%+7.2f %9.3f %5.3f %+7.2f %8.5f %+4.1f %11.2f\n",
	obj->type, obj->pos[0], obj->pos[1], obj->mag,
	obj->bulge_ratio, obj->bulge_req, obj->bulge_aspect, obj->bulge_posang,
	obj->disk_scale, obj->disk_aspect, obj->disk_posang, obj->z,
	obj->hubble_type, obj->noiseqarea);
      break;
    case 210:
      fprintf(sim->outlistfile, "%3d %11.4f %11.4f %8.4f %5.3f %9.3f %5.3f "
		"%4.2f %+7.2f %9.3f %5.3f %+7.2f %8.5f %+4.1f %11.2f\n",
	obj->type, obj->pos[0], obj->pos[1], obj->mag,
	obj->bulge_sersicn, obj->bulge_ratio, obj->bulge_req, obj->bulge_aspect,
	obj->bulge_posang, obj->disk_scale, obj->disk_aspect, obj->disk_posang,
	obj->z, obj->hubble_type, obj->noiseqarea);
      break;
    case 300:
      fprintf(sim->outlistfile, "%3d %11.4f %11.4f %8.4f "
		"%010d %9.3f %5.3f %+7.2f %8.5f %11.2f\n",
	obj->type, obj->pos[0], obj->pos[1], obj->mag,
	obj->raster_index, obj->raster_size, obj->raster_aspect,
	obj->raster_posang, obj->z, obj->noiseqarea);
      break;
    default:
      sprintf(str, "%d", obj->type);
      error(EXIT_FAILURE, "*Error*: Unknown object type in input list: ", str);
  }
  prefs.nobj++;

  return;
}


/****** list_openout *********************************************************
PROTO	void list_openout(simstruct *sim)
PURPOSE	Open the output list for subsequent writings.
INPUT	Pointer to the sim structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	09/03/2020
 ***/
void    list_openout(simstruct *sim)

  {
  if ((sim->outlistfile = fopen(sim->outlistname, "w")) == NULL)
    error(-1,"*ERROR*: can't create ", sim->outlistname);

  return;
  }


/****** list_closeout ********************************************************
PROTO	void list_closeout(simstruct *sim)
PURPOSE	Close the output list.
INPUT	Pointer to the sim structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	09/03/2020
 ***/
void    list_closeout(simstruct *sim)

  {
  fclose(sim->outlistfile);

  return;
  }


/****** list_endobj **********************************************************
PROTO	void list_endobj(objstruct *obj)
PURPOSE	Free memory allocated for an object.
INPUT	Pointer to the obj structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	09/03/2020
 ***/
void	list_endobj(objstruct *obj)
  {
  if (obj->subsize[0] && obj->subsize[1])
    free(obj->subimage);

  if (obj->buf1size)
    {
    free(obj->maskbuf);
    free(obj->nmaskbuf);
    free(obj->startbuf);
    }
  if (obj->buf2size)
    {
    free(obj->buf2);
    }

  return;
  }

#ifdef USE_THREADS

/****** pthread_list_readobj *************************************************
PROTO	void *pthread_list_readobj(void *arg)
PURPOSE	Object generation thread.
INPUT	Pointer to the thread number.
OUTPUT	NULL void pointer.
NOTES	Relies on global variables.
AUTHOR	E. Bertin (IAP)
VERSION	10/03/2020
 ***/
static void	*pthread_list_readobj(void *arg)
  {
   objstruct	*obj;
   char		str[MAXCHAR];
   int		obji, proc;

  proc = *((int *)arg);
  obji = -1;
/* Exit if the end of file is reached by any thread, including this one */
  while ((obji=pthread_list_nextobj(obji, str, proc)) != -1)
    {
    obj = &pthread_obj[obji];
    obj->ok = 0;
    if (list_readobj(pthread_sim, obj, str, proc) != RETURN_OK)
      continue;
    list_makeobj(pthread_sim, obj);
    obj->ok = 1;
    }

  pthread_exit(NULL);
  return (void *)NULL;
  }


/****** pthread_list_nextobj *************************************************
PROTO	int pthread_list_nextobj(int objindex, char *str)
PURPOSE	Manage object list reading and thread synchronisation.
INPUT	Previous object index in list,
	string buffer.
OUTPUT	Next object index in list.
NOTES	Relies on global variables.
AUTHOR	E. Bertin (IAP)
VERSION	19/05/2020
 ***/
static int	pthread_list_nextobj(int obji, char *str, int proc) {

   objstruct	*obj;
   int		q, retcode;

  QPTHREAD_MUTEX_LOCK(&objmutex);
/* The newly processed object is ready to be added to the image */
  if (obji>=0)
    pthread_addobjflag[obji] = 2;
/* If we just finished the "right" object, add it to image! */
  if (obji == pthread_addobji) {
    while (pthread_addobjflag[pthread_addobji]==2) {
      obj = &pthread_obj[pthread_addobji];
/*---- Add the object to the image */
      if (obj->subimage && obj->ok) {
        QPTHREAD_MUTEX_LOCK(&imagemutex);
        retcode = image_add(obj->subimage, obj->subsize[0], obj->subsize[1],
		pthread_sim->image, pthread_sim->fimasize[0],
		pthread_sim->fimasize[1],
		obj->subpos[0], obj->subpos[1], (float)obj->subfactor);
        obj->ok = 0;
        QPTHREAD_MUTEX_UNLOCK(&imagemutex);
/*------ Add the object to the output list */
        obj->flux = 0.0;
        if (!retcode) {
          list_writeobj(pthread_sim, obj);
          pthread_objnp++;
        }
      }
      pthread_addobjflag[pthread_addobji] = 0;
      QPTHREAD_COND_BROADCAST(&objcond[pthread_addobji]);
      pthread_addobji = (pthread_addobji+1)%pthread_nobj;
    }
  }
/* If no more object to process, return a "-1" (meaning exit thread) */
  if (!pthread_endflag && fgets(str, MAXCHAR, pthread_sim->inlistfile)
	&& (!pthread_gridflag || pthread_gridindex<pthread_ngrid)) {
    if (!(pthread_objn%READLIST_DISPSTEP)) {
      sprintf(pthread_str,
		"Painting input list... (%d objects painted / %d read)",
		pthread_objnp, pthread_objn);
      NFPRINTF(OUTPUT, pthread_str);
    }
    pthread_objn++;
    obji = pthread_procobji;
    pthread_procobji = (pthread_procobji+1)%pthread_nobj;
/*-- If the next available buffer has not been flushed yet, wait */
    q=++pthread_objqueue[obji];
    while (pthread_addobjflag[obji] || --q)
      QPTHREAD_COND_WAIT(&objcond[obji], &objmutex);
/*-- Set content */
    if (pthread_objqueue[obji])
      pthread_objqueue[obji]--;
    pthread_addobjflag[obji] = 1;
    if (pthread_gridflag) {
      pthread_obj[obji].pos[0] = (pthread_gridindex%pthread_ngridx)
		* pthread_gridstep
		+ pthread_gridoffsetx + random_double(proc)-0.5;
      pthread_obj[obji].pos[1] = (pthread_gridindex/pthread_ngridx)
		* pthread_gridstep
		+ pthread_gridoffsety + random_double(proc)-0.5;
      pthread_gridindex++;
    }
  } else {
    pthread_endflag = 1;
    obji=-1;
  }

  QPTHREAD_MUTEX_UNLOCK(&objmutex);

  return obji;
}


/****** pthread_list_read ****************************************************
PROTO	void pthread_list_read(simstruct *sim)
PURPOSE	Read the object list using multithreads.
INPUT	Pointer to the sim structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	09/03/2020
 ***/
static void	pthread_list_read(simstruct *sim)
  {
   static pthread_attr_t	pthread_attr;
   int				*proc,
				o, p;

/* Number of active threads */
  nproc = prefs.nthreads;
  pthread_nobj = 2*nproc;	/* A margin of 2X for better efficiency */
  QCALLOC(pthread_obj, objstruct, pthread_nobj);
  QCALLOC(pthread_addobjflag, int, pthread_nobj);
  QCALLOC(pthread_objqueue, int, pthread_nobj);
/* Set up multi-threading stuff */
  QMALLOC(objcond, pthread_cond_t, pthread_nobj);
  for (o=0; o<pthread_nobj; o++) {
    QPTHREAD_COND_INIT(&objcond[o], NULL);
  }
  QPTHREAD_MUTEX_INIT(&objmutex, NULL);
  QPTHREAD_MUTEX_INIT(&imagemutex, NULL);
  QMALLOC(proc, int, nproc);
  QMALLOC(thread, pthread_t, nproc);
  QPTHREAD_ATTR_INIT(&pthread_attr);
  QPTHREAD_ATTR_SETDETACHSTATE(&pthread_attr, PTHREAD_CREATE_JOINABLE);
  pthread_sim = sim;
  pthread_endflag = 0;
  pthread_addobji = pthread_procobji = pthread_objn = pthread_objnp = 0;
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
/* Start the reading/generation threads */
  for (p=0; p<nproc; p++)
    {
    proc[p] = p;
    QPTHREAD_CREATE(&thread[p], &pthread_attr, &pthread_list_readobj, &proc[p]);
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
    list_endobj(&pthread_obj[o]);
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

