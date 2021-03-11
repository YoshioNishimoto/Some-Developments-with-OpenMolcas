/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Per-Olof Widmark                                       *
***********************************************************************/
/**************************************************************************/
/*                                                                        */
/* This routine does xml dump of strings.                                 */
/*                                                                        */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* Author:  Per-Olof Widmark                                              */
/*          Lund University, Sweden                                       */
/*                                                                        */
/**************************************************************************/
#include <stdio.h>
#include "molcastype.h"
#include "xmlapi.h"
#ifdef _CAPITALS_
#define xml_cdumps XML_CDUMPS
#else
#ifndef ADD_
#define xml_cdumps xml_cdumps_
#endif
#endif
void xml_cdumps(char *name, int *nx_name) {
   FILE *f;
   char  line[256];
   int   n_name;
   int   k;

   n_name=*nx_name;

   if((f=fopen(XMLDUMP,"a"))==NULL) return;
   for(k=0; k<n_name; k++) { line[k]=name[k]; if(line[k]==' ') line[k]=0; }; line[n_name]=0;
   fprintf(f,"<%s>\n",line);

   fclose(f);
}
