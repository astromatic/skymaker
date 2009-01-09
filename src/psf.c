/*
				psf.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SkyMaker
*
*	Author:		E.BERTIN, (IAP)
*
*	Contents:	Routines dealing with the PSF.
*
*	Last modify:	19/03/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fftw3.h>

#include "define.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "fft.h"
#include "image.h"
#include "imaout.h"
#include "prefs.h"
#include "psf.h"
#include "simul.h"

/*********************************** makepsf *********************************/
/*
Create the PSF pixmap of the system Atmosphere+Optics+Target.
*/
void	makepsf(simstruct *sim)
  {
   PIXTYPE	*bmp,*bmp2, *pix,*pix2;
   double	pixscale[2], fsize[2], fscale[2], redscale[2],step[2],pstep[2],
		*carms, *carmst,
		redscale2, redscalem, xcf,ycf, r, r1,r2,h2, theta,ct,st, dx,dy,
		sig2, motfact, ft, phi, 
                phidefoc, phispher, phicomax, phicomay, phiast00, phiast45,
                phitri00, phitri30, phiqua00, phiqua22,
		grey, smallest, step2, rmax, norm, dsum,
		a2,b2, cxx,cyy,cxy, bandpass, rmax2, rsig,invrsig2, rmin,rmin2;
   float	*fbmp, *fbmpt;
   int		psfsize[2], hpsfsize[2], nstep[2],
		i, psfwm1, psfw2,psfh2,
		nbpix, nbpix2, x,y, x2,y2, x22,y22, xc,yc, osamp,osampm1,
		bstep, narms,dx0,dy0,r0, mask, sx,sy, trackflag;

/*-- Prepare the main PSF parameters */
  NFPRINTF(OUTPUT, "Preparing the PSF...");
/*-- Subsampling of the image pixel */
  osamp = (int)(sim->psfoversamp+0.5);
  pixscale[0] = sim->pixscale[0]/osamp;
  pixscale[1] = sim->pixscale[1]/osamp;
  nbpix = (psfsize[0]=sim->psfsize[0])*(psfsize[1]=sim->psfsize[1]);
  fsize[0] = (sim->lambdaeq*MICRON)/(pixscale[0]*ARCSEC);
  fsize[1] = (sim->lambdaeq*MICRON)/(pixscale[1]*ARCSEC);
/* Subpixel size in the pupil plane (in m) */
  fscale[0] = fsize[0]/psfsize[0];
  fscale[1] = fsize[1]/psfsize[1];
  xc = (int)(xcf = psfsize[0]/2.0);
  yc = (int)(ycf = psfsize[1]/2.0);

/* Reduced scale for computing phase factor: it includes fwidth/D, which
   is 1/(fu), where f is the spatial frequency D/lambda and u the subpixel size
   in radian, and the factor 2/psfw to normalize dx0 and dy0 to 1 */
  redscale[0] = 2.*fscale[0]/sim->psfdm1;
  redscale[1] = 2.*fscale[1]/sim->psfdm1;
  redscale2 = redscale[0]*redscale[1];
  redscalem = sqrt(redscale2);
/* The smallest element scales define the subsampling */
  smallest = sim->psfdm1/1000.0;	/* default: 1/1000 of M1 diameter */
  if (smallest>=fscale[0])
    smallest = fscale[0];
  if (smallest>=fscale[1])
    smallest = fscale[1];
  if (sim->psfdm2>sim->lambdaeq*MICRON && sim->psfdm2/4.0<smallest)
    smallest = sim->psfdm2/4.0;
  if (sim->psfarmw>sim->lambdaeq*MM && sim->psfarmw*MM/4.0<smallest)
    smallest = sim->psfarmw*MM/4.0;

/* Convert to a subsampling unit of the pupil */
  nstep[0]=(int)(fscale[0]/smallest);
  if (nstep[0]>10)
    nstep[0] = 10;
  step[0] = 1.0/nstep[0];
  nstep[1]=(int)(fscale[1]/smallest);
  if (nstep[1]>10)
    nstep[1] = 10;
  step[1] = 1.0/nstep[1];
  step2 = step[0]*step[1];
  pstep[0] = 0.5*(step[0]-1);
  pstep[1] = 0.5*(step[1]-1);

/* Now compute the ``phase factor'' to create defocused or aberrated images.
   Set-up parameters in the configuration file corresponds to d80%, which is
   roughly FWHM for a gaussian. */
/* Original formula (EB)
   phidefoc   = PI*sim->psfdefocfwhm*ARCSEC*fwscale*fwscale
	       /(sim->psfdm1*sim->lambdaeq*MICRON);
*/
  phidefoc = PI*sim->psfdm1/(sim->lambdaeq*MICRON)*sim->psfd80defoc*ARCSEC
             *0.279*redscale2;
  phispher = PI*sim->psfdm1/(sim->lambdaeq*MICRON)*sim->psfd80spher*ARCSEC
             *0.179*redscale2*redscale2;
  phicomax = PI*sim->psfdm1/(sim->lambdaeq*MICRON)*sim->psfd80comax*ARCSEC
             *0.534*redscale2*redscalem;
  phicomay = PI*sim->psfdm1/(sim->lambdaeq*MICRON)*sim->psfd80comay*ARCSEC
             *0.534*redscale2*redscalem;
  phiast00 = PI*sim->psfdm1/(sim->lambdaeq*MICRON)*sim->psfd80ast00*ARCSEC
             *0.279*redscale2;
  phiast45 = PI*sim->psfdm1/(sim->lambdaeq*MICRON)*sim->psfd80ast45*ARCSEC
             *0.279*redscale2;
  phitri00 = PI*sim->psfdm1/(sim->lambdaeq*MICRON)*sim->psfd80tri00*ARCSEC
             *0.208*redscale2*redscalem;
  phitri30 = PI*sim->psfdm1/(sim->lambdaeq*MICRON)*sim->psfd80tri30*ARCSEC
             *0.208*redscale2*redscalem;
  phiqua00 = PI*sim->psfdm1/(sim->lambdaeq*MICRON)*sim->psfd80qua00*ARCSEC
             *0.172*redscale2*redscale2;
  phiqua22 = PI*sim->psfdm1/(sim->lambdaeq*MICRON)*sim->psfd80qua22*ARCSEC
             *0.172*redscale2*redscale2;
  fprintf(OUTPUT, "Pupil range: %.1f-%.1fmm / Subsampling: %dX / "
		"ro=%.1fcm / motion=%.2f'' rms\n",
	fscale[0]/MM, fsize[0]/MM, nstep[0], sim->ro/CM, sim->psfmotion);
  if ((sim->psfseeingtype!=LONG_EXPOSURE || fsize[0]<sim->psfdm1+sim->ro)
	&& fsize[0]<2*sim->psfdm1)
    warning("PSF mask is too coarse (important aliasing): ",
	"increase oversampling");
  narms = sim->psfnarms;
  if (narms && fscale[0]>sim->psfarmw)
    warning("PSF mask is too small to yield accurate diffraction spikes: ",
	"increase map width");
  if (narms && fscale[1]>sim->psfarmw)
    warning("PSF mask is too small to yield accurate diffraction spikes: ",
	"increase map height");

/* Generate the PSF in the Fourier plane */
  hpsfsize[0] = psfsize[0]/2;
  hpsfsize[1] = psfsize[1]/2;
  r1 = sim->psfdm1/fscale[0]/2;
  r1 *= r1;
  r2 = sim->psfdm2/fscale[0]/2;
  r2 *= r2;
  h2 = sim->psfarmw*MM/fscale[0]/2;

  NFPRINTF(OUTPUT, "Creating the pupil mask...");
/* Allocate memory */
  QFFTWMALLOC(fbmp, float, nbpix*2);
  if (narms)
    QCALLOC(carms, double, narms*2)
  else
    carms = NULL;
  carmst = carms;
  for (i=0; i<narms; i++)
    {
    theta = fmod((sim->psfarmang+i*360.0/narms), 360.0)*DEG;
    *(carmst++) = sin(theta);
    *(carmst++) = cos(theta);
    }
/* Create the aperture mask */
  fbmpt = fbmp;
  for (y=0; y<psfsize[1]; y++)
    {
    dy0 = y - hpsfsize[1] + pstep[0];
    for (x=0; x<psfsize[0]; x++)
      {
      dx0 = x - hpsfsize[0] + pstep[1];
      grey = 0.0;
      for (dy=dy0,sy=0; sy<nstep[1]; sy++,dy+=step[1])
        for (dx=dx0,sx=0; sx<nstep[0]; sx++,dx+=step[0])
          {
          r = dx*dx+dy*dy;
/* Put light between central hole and border, except on spider arms */
          if ((mask = r>=r2 && r<r1))
            {
            carmst = carms;
            for (i=0; i<narms; i++)
              {
              st = *(carmst++);
              ct = *(carmst++);
              r = st*dx-ct*dy;
              if (st<0.0)
                r = -r;
              if (r>-h2 && r<h2 && ct*dx+st*dy>0.0)
                mask = 0;
              }
            }
          if (mask)
            grey += step2;
          }
      r0 = dx0*dx0+dy0*dy0;
/* Original formula (EB): note that r0 is a square
      phi = phidefoc*r0;
*/
      phi = phidefoc*r0 + phispher*r0*r0 + (phicomax*dx0+phicomay*dy0)*r0
	+ phiast00*(dx0*dx0-dy0*dy0)+phiast45*2.*dx0*dy0
	+ phitri00*dx0*(dx0*dx0-3.*dy0*dy0)+phitri30*dy0*(3.*dx0*dx0-dy0*dy0)
	+ phiqua00*(dx0*dx0*dx0*dx0+dy0*dy0*dy0*dy0-6.*dx0*dx0*dy0*dy0)
	+ phiqua22*4.*dx0*dy0*(dx0*dx0-dy0*dy0);

      *(fbmpt++) = grey*cos(phi);	/* real part */
      *(fbmpt++) = grey*sin(phi);	/* imaginary part */
      }
    }

  free(carms);

/* Save the image of the pupil ("no-return" option!) */
  if  (sim->imatype == PUPIL_REAL
	|| sim->imatype == PUPIL_IMAGINARY
	|| sim->imatype == PUPIL_MODULUS
	|| sim->imatype == PUPIL_PHASE)
    {
    QMALLOC(bmp, float, nbpix);
    pix = bmp;
    fbmpt = fbmp;
    switch(sim->imatype)
      {
      case PUPIL_REAL:
        for (i=0; i<nbpix; i++,fbmpt+=2)
          *(pix++) = *fbmpt;
        break;
      case PUPIL_IMAGINARY:
        fbmpt++;
        for (i=0; i<nbpix; i++,fbmpt+=2)
          *(pix++) = *fbmpt;
        break;
      case PUPIL_MODULUS:
        for (i=0; i<nbpix; i++,fbmpt++)
          {
          *pix = *fbmpt**fbmpt;
          fbmpt++;
          *pix = sqrt(*pix+*fbmpt**fbmpt);
          pix++;
          }
        break;
      case PUPIL_PHASE:
        for (i=0; i<nbpix; i++,fbmpt+=2)
          {
          *(pix++) = atan2(*(fbmpt+1), *fbmpt?*fbmpt:1e-6)/DEG;
          }
        break;
      default:
        error(EXIT_FAILURE, "*Internal Error*: no such imatype in ",
		"makepsf()" );
      }
    sim->image = bmp;
    sim->imasize[0] = psfsize[0];
    sim->imasize[1] = psfsize[1];
    NFPRINTF(OUTPUT, "Writing image...");
    writeima(sim);
    NFPRINTF(OUTPUT, "Freeing memory...");
    QFFTWFREE(bmp);
    NFPRINTF(OUTPUT, "All done");
    QPRINTF(OUTPUT, "\n");
    exit(EXIT_SUCCESS);
    }

/* Let's Fourier transform to get the amplitude field */

  NFPRINTF(OUTPUT, "1st FFT...");
  fft_init(prefs.nthreads);
  fft_ctf(fbmp, psfsize[0], psfsize[1], FFTW_FORWARD);

/* Now convert to Intensity */

  NFPRINTF(OUTPUT, "Converting to intensity...");
  fbmpt = fbmp+1;
  for (i=0; i<nbpix; i++)
    {
    ft = *fbmpt;
    *(fbmpt--) = 0.0;
    *fbmpt = (*fbmpt * *fbmpt + ft*ft)/nbpix;
    fbmpt += 3;
    }

/* Go back to frequency domain */

  NFPRINTF(OUTPUT, "2nd FFT...");
  fft_ctf(fbmp, psfsize[0], psfsize[1], FFTW_BACKWARD);

/* Save the MTF of the pupil ("no-return" option!) */
  if  (sim->imatype == PUPIL_MTF)
    {
    NFPRINTF(OUTPUT, "Descrambling pixels...");
    QMALLOC(bmp, float, nbpix);
    pix = bmp;
    dsum = 0.0;
    for (y=0; y<psfsize[1]; y++)
      {
      y2 = y-yc;
      if (y2<0)
        y2 += psfsize[1];
      y2 *= psfsize[0];
      for (x=0; x<psfsize[0]; x++)
        {
        x2 = x-xc;
        if (x2<0)
          x2 += psfsize[0];
        dsum += *(pix++) = *(fbmp+2*(x2+y2));
        }
      }
    pix = bmp;
    for (i=0; i<nbpix; i++)
      *(pix++) /= (PIXTYPE)dsum;
    sim->image = bmp;
    sim->imasize[0] = psfsize[0];
    sim->imasize[1] = psfsize[1];
    NFPRINTF(OUTPUT, "Writing image...");
    writeima(sim);
    NFPRINTF(OUTPUT, "Freeing memory...");
    QFFTWFREE(bmp);
    fft_end(prefs.nthreads);
    NFPRINTF(OUTPUT, "All done");
    NPRINTF(OUTPUT, "\n");
    exit(EXIT_SUCCESS);
    }

/* Preparing tracking parameters */
  sig2 = sim->ro/fscale[0];
  sig2 *= sig2;
  motfact = fscale[0]*fscale[1]/(sim->psfdm1*sim->psfdm1);
  rmax = 1.0/motfact;
  ct=st=a2=b2=cxx=cyy=cxy = 0.0;	/* To avoid gcc -Wall warnings */
  if (sim->psftracktype!=NO_TRACKERR && (sim->psftrackmaj || sim->psftrackmin))
    {
    trackflag = 1;
    a2 = 2*PI*sim->psftrackmaj/(psfsize[0]*pixscale[0]);
    b2 = 2*PI*sim->psftrackmin/(psfsize[0]*pixscale[0]);
    ct = cos(sim->psftrackang*DEG);
    st = sin(sim->psftrackang*DEG);
    switch(sim->psftracktype)
      {
      case JITTER:
        a2 *= a2;
        b2 *= b2;
        cxx = ct*ct*a2+st*st*b2;
        cyy = st*st*a2+ct*ct*b2;
        cxy = 2*ct*st*(a2-b2);
        break;
      case DRIFT:
        a2 /= 2.0;
        b2 *= b2;
        break;
      default:
        error(EXIT_FAILURE, "*Internal ERROR*: Track-Error Type Unknown",
                        " in makepsf()");
        break;
      }
    }
  else
    trackflag = 0;

/* Convolve by seeing */
  if (sim->psfseeingtype != NO_SEEING)
    {
    NFPRINTF(OUTPUT, "Convolving by seeing...");
    fbmpt = fbmp;
    for (y=0; y<psfsize[1]; y++)
      {
      y2 = y<hpsfsize[1]?y:y-psfsize[1];
      y22 = y2*y2;
      for (x=0; x<psfsize[0]; x++)
        {
        x2 = x<hpsfsize[0]?x:x-psfsize[0];
        x22 = x2*x2;
        r = x22+y22;
        if (sim->psfseeingtype == SHORT_EXPOSURE)
          bandpass = r<rmax? exp(-3.4419*pow(r/sig2,5./6.)
		*(1-pow(motfact*r,1./6.))) : 0.0;
        else
          bandpass = exp(-3.4419*pow(r/sig2,5./6.));
        *(fbmpt++) *= bandpass;
        *(fbmpt++) *= bandpass;
        }
      }
    }

/* Convolve by tracking errors */
  if (trackflag)
    {
    NFPRINTF(OUTPUT, "Convolving by tracking errors...");
    fbmpt = fbmp;
    for (y=0; y<psfsize[1]; y++)
      {
      y2 = y<hpsfsize[1]?y:y-psfsize[1];
      y22 = y2*y2;
      for (x=0; x<psfsize[0]; x++)
        {
        x2 = x<hpsfsize[0]?x:x-psfsize[0];
        x22 = x2*x2;
        if (sim->psftracktype == JITTER)
          bandpass = exp(-(cxx*x22+cyy*y22+cxy*x2*y2)/2);
        else
          {
          r = a2*(ct*x2+st*y2);
          bandpass = (r?sin(r)/r:1.0);
          if (b2)
            {
            r = x22+y22;
            bandpass *= exp(-b2*r/2);
            }
          }
        *(fbmpt++) *= bandpass;
        *(fbmpt++) *= bandpass;
        }
      }
    }

/* Get the normalization from zero-frequency Fourier component */
  norm = *fbmp;

/* Save the MTF of the PSF ("no-return" option!) */
  if  (sim->imatype == PSF_MTF)
    {
    NFPRINTF(OUTPUT, "Descrambling pixels...");
    QMALLOC(bmp, float, nbpix);
    pix = bmp;
    for (y=0; y<psfsize[1]; y++)
      {
      y2 = y-yc;
      if (y2<0)
        y2 += psfsize[1];
      y2 *= psfsize[0];
      for (x=0; x<psfsize[0]; x++)
        {
        x2 = x-xc;
        if (x2<0)
          x2 += psfsize[0];
        *(pix++) = *(fbmp+2*(x2+y2));
        }
      }
    pix = bmp;
    for (i=0; i<nbpix; i++)
      *(pix++) /= (PIXTYPE)norm;
    sim->image = bmp;
    sim->imasize[0] = psfsize[0];
    sim->imasize[1] = psfsize[1];
    NFPRINTF(OUTPUT, "Writing image...");
    writeima(sim);
    NFPRINTF(OUTPUT, "Freeing memory...");
    QFFTWFREE(bmp);
    fft_end(prefs.nthreads);
    NFPRINTF(OUTPUT, "All done");
    QPRINTF(OUTPUT, "\n");
    exit(EXIT_SUCCESS);
    }

  NFPRINTF(OUTPUT, "3rd FFT...");
  fft_ctf(fbmp, psfsize[0], psfsize[1], FFTW_FORWARD);
  fft_end(prefs.nthreads);
  norm *= (PIXTYPE)nbpix;

/* Copy to the psf pixmap with the correct arrangement (and normalize) */
  NFPRINTF(OUTPUT, "Descrambling pixels...");
  QMALLOC(bmp, float, nbpix);
  pix = bmp;
  for (y=0; y<psfsize[1]; y++)
    {
    y2 = y-yc;
    y22 = y2*y2;
    if (y2<0)
      y2 += psfsize[1];
    y2 *= psfsize[0];
    for (x=0; x<psfsize[0]; x++, pix++)
      {
      x2 = x-xc;
      x22 = x2*x2;
      if (x2<0)
        x2 += psfsize[0];
      *pix = *(fbmp+2*(x2+y2))/norm;
      }
    }

  free(fbmp);

/* Truncate to a disk that has diameter = (box width - 2*PSF "radius")*/
  NFPRINTF(OUTPUT, "Truncating...");
  rmax = psfsize[0] - 1.0 - xc;
  if (rmax > (r=xc))
    rmax = r;
  if (rmax > (r=psfsize[1]-1.0-yc))
    rmax = r;
  if (rmax > (r=yc))
    rmax = r;
  if (rmax<1.0)
    rmax = 1.0;
  rmax2 = rmax*rmax;
  rsig = sqrt(sim->psfarea/PI)*osamp;
  invrsig2 = 1/(2*rsig*rsig);
  rmin = rmax - (3*rsig);	/* 3 sigma annulus (almost no aliasing) */
  rmin2 = rmin*rmin;

  pix = bmp;
  dy = -yc;
  for (y=psfsize[1]; y--; dy+=1.0)
    {
    dx = -xc;
    for (x=psfsize[0]; x--; dx+=1.0, pix++)
      if ((r2=dx*dx+dy*dy)>rmin2)
        *pix *= (r2>rmax2)?0.0:exp((2*rmin*sqrt(r2)-r2-rmin2)*invrsig2);
    }

/* Save the full resolution image of the PSF ("no-return" option!) */
  if  (sim->imatype == PSF_FULLRES)
    {
    sim->image = bmp;
    sim->imasize[0] = psfsize[0];
    sim->imasize[1] = psfsize[1];
    NFPRINTF(OUTPUT, "Writing image...");
    writeima(sim);
    NFPRINTF(OUTPUT, "Freeing memory...");
    QFFTWFREE(bmp);
    NFPRINTF(OUTPUT, "All done");
    QPRINTF(OUTPUT, "\n");
    exit(EXIT_SUCCESS);
    }

  osampm1 = osamp-1;
  nbpix2 = (psfw2=psfsize[0]-osampm1)*(psfh2=psfsize[1]-osampm1);
  QCALLOC(bmp2, PIXTYPE, nbpix2);

/* Convolve by the pixel */
  NFPRINTF(OUTPUT, "Convolving by detector pixel...");

  bstep = psfsize[0]*osamp - 1;
  psfwm1 = psfw2 - 1;
  pix = bmp;
  pix2 = bmp2;
  dsum = 0.0;
  for (y=0; y<psfh2; y++, pix+=osampm1)
    for (x=0; x<psfw2; x++)
      {
      for (y2=0; y2<osamp; y2++, pix+=psfwm1)
        for (x2=0; x2<osamp; x2++)
          *pix2 += (PIXTYPE)*(pix++);
      dsum += *(pix2++);
      pix -= bstep;
      }

  QFFTWFREE(bmp);

/* Update PSF parameters in the simulation structure */
  sim->psf=bmp2;
  sim->psfsize[0] = psfw2;
  sim->psfsize[1] = psfh2;
  sim->psfc[0] = xcf - osampm1/2.0;
  sim->psfc[1] = ycf - osampm1/2.0;

/* Normalize the PSF */
  NFPRINTF(OUTPUT, "Normalizing PSF...");
  dsum /= osamp*osamp;
  pix = bmp2;
  for (i=0; i<nbpix2; i++)
    *(pix++) /= (PIXTYPE)dsum;

/* Save the final resolution image of the PSF ("no-return" option!) */
  if  (sim->imatype == PSF_FINALRES)
    {
    sim->image = bmp2;
    sim->imasize[0] = psfw2;
    sim->imasize[1] = psfh2;
    NFPRINTF(OUTPUT, "Writing image...");
    writeima(sim);
    NFPRINTF(OUTPUT, "Freeing memory...");
    free(bmp2);
    NFPRINTF(OUTPUT, "All done");
    QPRINTF(OUTPUT, "\n");
    exit(EXIT_SUCCESS);
    }

  return;
  }


/********************************* makeaureole *******************************/
/*
Create the pixmap (in final resolution) of the aureole generated by diffusion.
*/
void	makeaureole(simstruct *sim)
  {
   PIXTYPE	*pix, amp;
   double	sum;
   int		i,j, x,y, y2,r,o, w,h,haw,hah, rmax, rmed, rmedax;

/* First compute the mask-size that should yield the max. processing speed */
  if ((i=4*sim->margin[0])>sim->fimasize[0])
    i = sim->fimasize[0];
  i /= sim->mscan[0];
  for (i--, j=1; i!=0; i/=2, j*=2);
  w = sim->aursize[0] = j;
  if ((i=4*sim->margin[1])>sim->fimasize[1])
    i = sim->fimasize[1];
  i /= sim->mscan[1];
  for (i--, j=1; i!=0; i/=2, j*=2);
  h = sim->aursize[1] = j;
  QFFTWCALLOC(sim->aureole, PIXTYPE, w*h);
  haw = sim->aurange/sim->mscan[0];
  hah = sim->aurange/sim->mscan[1];
  amp = 60.0*60.0*DEXP(-0.4*sim->psfhalosb);	/* 60 arcsec = 1' */
  rmax = haw*hah;
  rmed = rmax*81.0/100.0;
  rmedax = rmax-rmed;
/* Lay down the aureole and normalize it */
  o = 0;
  do
    {
    sum = 0.0;
    pix = sim->aureole + w*(h-hah+1) - haw;
    for (y=-hah; y<hah; y++, pix+=2*(w-haw))
      {
      if (!y)
        pix -= w*h;
      y2 = y*y;
      for (x=-haw; x<haw; x++)
        {
        if (!x)
          pix -= w;
        if ((r=x*x+y2) && r<rmax)
          sum += (double)(*(pix++) = (r<rmed?amp:amp*(rmax-r)/rmedax)/(r+o*o));
        else
          pix++;
        }
      }
    }  while ((*sim->aureole=1.0-(PIXTYPE)sum)<0.0 && o++<8);

/* In case even the Lorentz law did not succeed... */
  if (o>10)
    error(EXIT_FAILURE, "*Internal ERROR*: bad aureole normalization in ",
	"makeaureole()");

  return;
  }


/********************************** makedft **********************************/
/*
Create the Fourier transform of the PSF components. Sizes are assumed to be
even numbers.
*/
void	makedft(simstruct *sim, int order)
  {
   PIXTYPE	*mask,*maskt, *psf,*psft;
   double	dx,dy, r,r2,rmin,rmin2,rmax,rmax2,rsig,invrsig2;
   int		width,height,npix,offset, psfwidth,psfheight,psfnpix,
		cpwidth, cpheight,hcpwidth,hcpheight, i,j,p,x,y;

  psfwidth = sim->psfsize[0];
  psfheight = sim->psfsize[1];
  psfnpix = psfwidth*psfheight;
  width = height = 1<<order;
  npix = width*height;
  cpwidth = (width>psfwidth)?psfwidth:width;
  hcpwidth = cpwidth>>1;
  cpwidth = hcpwidth<<1;
  offset = width - cpwidth;
  cpheight = (height>psfheight)?psfheight:height;
  hcpheight = cpheight>>1;
  cpheight = hcpheight<<1;
  QCALLOC(mask, PIXTYPE, npix);

  for (p=0; p<sim->npsf; p++)
    {
    psf = sim->psf + p*psfnpix;
/*-- Frame and descramble the PSF data */
    psft = psf + (psfheight/2)*psfwidth + psfwidth/2;
    maskt = mask;
    for (j=hcpheight; j--; psft+=psfwidth)
      {
      for (i=hcpwidth; i--;)
        *(maskt++) = *(psft++);      
      psft -= cpwidth;
      maskt += offset;
      for (i=hcpwidth; i--;)
        *(maskt++) = *(psft++);      
      }

    psft = psf + ((psfheight/2)-hcpheight)*psfwidth + psfwidth/2;
    maskt += width*(height-cpheight);
    for (j=hcpheight; j--; psft+=psfwidth)
      {
      for (i=hcpwidth; i--;)
        *(maskt++) = *(psft++);      
      psft -= cpwidth;
      maskt += offset;
      for (i=hcpwidth; i--;)
        *(maskt++) = *(psft++);      
      }

/* Truncate to a disk that has diameter = (box width) */
    rmax = cpwidth - 1.0 - hcpwidth;
    if (rmax > (r=hcpwidth))
      rmax = r;
    if (rmax > (r=cpheight-1.0-hcpheight))
      rmax = r;
    if (rmax > (r=hcpheight))
      rmax = r;
    if (rmax<1.0)
      rmax = 1.0;
    rmax2 = rmax*rmax;
    rsig = sqrt(sim->psfarea/PI)*sim->psfoversamp;
    invrsig2 = 1/(2*rsig*rsig);
    rmin = rmax - (3*rsig);	/* 3 sigma annulus (almost no aliasing) */
    rmin2 = rmin*rmin;

    maskt = mask;
    dy = 0.0;
    for (y=hcpheight; y--; dy+=1.0)
      {
      dx = 0.0;
      for (x=hcpwidth; x--; dx+=1.0, maskt++)
        if ((r2=dx*dx+dy*dy)>rmin2)
          *maskt *= (r2>rmax2)?0.0:exp((2*rmin*sqrt(r2)-r2-rmin2)*invrsig2);
      dx = -hcpwidth;
      maskt += offset;
      for (x=hcpwidth; x--; dx+=1.0, maskt++)
        if ((r2=dx*dx+dy*dy)>rmin2)
          *maskt *= (r2>rmax2)?0.0:exp((2*rmin*sqrt(r2)-r2-rmin2)*invrsig2);
      }
    dy = -hcpheight;
    maskt += width*(height-cpheight);
    for (y=hcpheight; y--; dy+=1.0)
      {
      dx = 0.0;
      for (x=hcpwidth; x--; dx+=1.0, maskt++)
        if ((r2=dx*dx+dy*dy)>rmin2)
          *maskt *= (r2>rmax2)?0.0:exp((2*rmin*sqrt(r2)-r2-rmin2)*invrsig2);
      dx = -hcpwidth;
      maskt += offset;
      for (x=hcpwidth; x--; dx+=1.0, maskt++)
        if ((r2=dx*dx+dy*dy)>rmin2)
          *maskt *= (r2>rmax2)?0.0:exp((2*rmin*sqrt(r2)-r2-rmin2)*invrsig2);
      }

/* Finally move to Fourier space */
    sim->psfdft[order] = fft_rtf(mask, width,height);
    }

  free(mask);

  return;
  }


/******************************** freepsf ***********************************/
/*
Free allocated memory for the PSF.
*/
void    freepsf(simstruct *sim)

  {
   int	i;

  free(sim->psf);
  for (i=0; i<PSF_NORDER; i++)
    {
    QFFTWFREE(sim->psfdft[i]);
    }

  return;
  }


/******************************** readpsf ***********************************/
/*
Read an image containing the centered PSF.
*/
void	readpsf(simstruct *sim)
  {
   catstruct	*cat;
   tabstruct	*tab;
   double	dsum;
   PIXTYPE	*pix,
		sum;
   char		lstr[MAXCHAR],
		*filename,*rfilename, *str, *str2;
   int		i, ext, nbpix;

  filename = sim->psfname;
  if ((str = strrchr(filename, '[')))
    {
    *str = '\0';
    if ((str2 = strrchr(str, ']')))
      *str2 = '\0';
    }

/* A short, "relative" version of the filename */
  if (!(rfilename = strrchr(filename, '/')))
    rfilename = filename;
  else
    rfilename++;

  sprintf(lstr,"Examining %s", rfilename);
  NFPRINTF(OUTPUT, lstr);
  if (!(cat=read_cat(filename)))
    error(EXIT_FAILURE, "*Error*: no FITS data found in ", filename);

  if (str)
    {
    if (!(tab = name_to_tab(cat, str+1, 0)))
      {
      ext = atoi(str+1);
      tab = pos_to_tab(cat, ext, 0);
      }
    }
  else
    tab = cat->tab;

/*-- Force the data to be at least 2D */
  if (tab->naxis<2)
    error(EXIT_FAILURE, "*Error*: no 2D FITS data in ", filename);
  sim->psfsize[0] = tab->naxisn[0];
  sim->psfsize[1] = tab->naxisn[1];
  sim->psfc[0] = (sim->psfsize[0]-1.0)/2.0;
  sim->psfc[1] = (sim->psfsize[1]-1.0)/2.0;
  sprintf(lstr,"Loading %s", rfilename);
  NFPRINTF(OUTPUT, lstr);
  sim->psf = alloc_body(tab, NULL);
  tab->bodybuf = NULL;
  free_cat(&cat, 1);

/* Normalize the PSF */
  NFPRINTF(OUTPUT, "Normalizing PSF...");
  nbpix = sim->psfsize[0]*sim->psfsize[0];
  pix = sim->psf;
  dsum = 0.0;
  for (i=nbpix; i--;)
    dsum += *(pix++);
  sum = (PIXTYPE) (dsum / (sim->psfoversamp*sim->psfoversamp));
  pix = sim->psf;
  for (i=nbpix; i--;)
    *(pix++) /= sum;

  return;
  }
