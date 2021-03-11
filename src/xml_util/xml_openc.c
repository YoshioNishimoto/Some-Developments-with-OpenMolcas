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
/* This routine opens an xml container.                                   */
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
#define xml_openc XML_OPENC
#else
#ifndef ADD_
#define xml_openc xml_openc_
#endif
#endif
void xml_openc(char *name,        INT *nx_name,
               char *appear,      INT *nx_appear,
               char *units, INT *nx_units, INT *Level,
               char *value,       INT *nx_value) {
   FILE *f;
   char  line[256];
   int   n_name;
   int   n_appear;
   int   n_units;
   int   n_value;
   int   k;
   int   level;

   n_name=*nx_name;
   n_appear=*nx_appear;
   n_units=*nx_units;
   n_value=*nx_value;
   level=*Level;

   if((f=fopen(XMLDUMP,"a"))==NULL) return;
   for(k=0; k<n_name; k++) { line[k]=name[k]; if(line[k]==' ') line[k]=0; }; line[n_name]=0;
   fprintf(f,"<%s",line);
   xml_prspec(f,"appear",appear,n_appear);
   xml_prspec(f,"units",units,n_units);
   if(level>0) fprintf(f," level=\"%i\"",level);
   xml_prspec(f,"value",value,n_value);
   fprintf(f,">\n");
   fclose(f);
}
