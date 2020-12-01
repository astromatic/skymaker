/*
*				main.c
*
* Command-line parsing
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 1998-2020 IAP/CNRS/SorbonneU
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
*	Last modified:		13/11/2020
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "globals.h"
#include "prefs.h"

#define		SYNTAX \
EXECUTABLE " [<list_file>] [-c <configuration_file>] [-<keyword> <value>]\n"\
"> to dump a default configuration file: " EXECUTABLE " -d \n" \
"> to dump a default extended configuration file: " EXECUTABLE " -dd \n"

extern const char	notokstr[];
char			prefsname[MAXCHAR];

//------------------- global variables for multithreading -------------------

int			nproc;	// Number of child threads

/********************************** main ************************************/
int	main(int argc, char *argv[])

  {
   char		**argkey, **argval;
   double	tdiff, nobj;
   int		a, narg, opt,opt2;

  if (argc<1)
    {
    fprintf(OUTPUT, "\n         %s  version %s (%s)\n", BANNER,MYVERSION,DATE);
    fprintf(OUTPUT, "\nWritten by %s\n", AUTHORS);
    fprintf(OUTPUT, "Copyright %s\n", COPYRIGHT);
    fprintf(OUTPUT, "\nvisit %s\n", WEBSITE);
    fprintf(OUTPUT, "\n%s\n", DISCLAIMER);
    error(EXIT_SUCCESS, "SYNTAX: ", SYNTAX);
    }

  QMALLOC(argkey, char *, argc);
  QMALLOC(argval, char *, argc);

/*default parameters */
  strcpy(prefsname, "sky.conf");
  *(prefs.inlistname) = 0;
  narg = 0;

  for (a=1; a<argc; a++)
    {
    if (*(argv[a]) == '-')
      {
      opt = (int)argv[a][1];
      if (strlen(argv[a])<4 || opt == '-')
        {
        opt2 = (int)tolower((int)argv[a][2]);
        if (opt == '-')
          {
          opt = opt2;
          opt2 = (int)tolower((int)argv[a][3]);
          }
        switch(opt)
          {
          case 'c':
            if (a<(argc-1))
              strcpy(prefsname, argv[++a]);
            break;
          case 'd':
            dumpprefs(opt2=='d' ? 1 : 0);
            exit(EXIT_SUCCESS);
            break;
          case 'v':
            printf("%s version %s (%s)\n", BANNER,MYVERSION,DATE);
            exit(EXIT_SUCCESS);
            break;
          case 'h':
            fprintf(OUTPUT, "\nSYNTAX: %s", SYNTAX);
            exit(EXIT_SUCCESS);
            break;
          default:
            error(EXIT_SUCCESS,"SYNTAX: ", SYNTAX);
          }
        }
      else
        {
        argkey[narg] = &argv[a][1];
        argval[narg++] = argv[++a];
        }       
      }
    else
      strcpy(prefs.inlistname, argv[a]);
    }

  readprefs(prefsname, argkey, argval, narg);
  useprefs();

  free(argkey);
  free(argval);

  makeit();

  endprefs();

  NFPRINTF(OUTPUT, "");
  QPRINTF(OUTPUT, "%d sources added\n", prefs.nobj);
  tdiff = prefs.time_diff>0.0? prefs.time_diff : 0.001;
  nobj = (double)prefs.nobj/tdiff;
  NPRINTF(OUTPUT, "> All done (in %.1f s: %.1f source%s/s)\n",
        prefs.time_diff, nobj, nobj>1.0? "s":"");

  exit(EXIT_SUCCESS);
  }
