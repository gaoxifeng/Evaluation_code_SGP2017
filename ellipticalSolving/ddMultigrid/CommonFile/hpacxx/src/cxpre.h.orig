/*
   Copyright (C)  2000    Daniel A. Atkinson  <DanAtk@aol.com>
   Copyright (C)  2004    Ivano Primi  <ivprimi@libero.it>    

   This file is part of the HPA Library.

   The HPA Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The HPA Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the HPA Library; if not, write to the Free
   Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
   02110-1301 USA.
*/

#ifndef _CXPRE_H_
#define _CXPRE_H_

#include<stdio.h>
#include"xpre.h"

#ifdef __cplusplus
extern "C"
{
#endif

  struct cxpr
  {
    struct xpr re, im;
  };

  struct cxprcmp_res
  {
    int re, im;
  };

  extern const struct cxpr cxZero;
  extern const struct cxpr cxOne;
  extern const struct cxpr cxIU;

  struct cxpr cxreset (struct xpr re, struct xpr im);
  struct cxpr cxconv (struct xpr x);
  struct xpr cxre (struct cxpr z);
  struct xpr cxim (struct cxpr z);
  struct cxpr cxswap (struct cxpr z);

  struct xpr cxabs (struct cxpr z);
  struct xpr cxarg (struct cxpr z);
  int cxrec (struct cxpr z, struct cxpr *w);

  struct cxpr cxadd (struct cxpr z1, struct cxpr z2, int k);
  struct cxpr cxsum (struct cxpr z1, struct cxpr z2);
  struct cxpr cxsub (struct cxpr z1, struct cxpr z2);
  struct cxpr cxmul (struct cxpr z1, struct cxpr z2);
  /* Multiplication by a real number */
  struct cxpr cxrmul (struct xpr c, struct cxpr z);

  /* Multiplication by +i */
  struct cxpr cxdrot (struct cxpr z);

  /* Multiplication by -i */
  struct cxpr cxrrot (struct cxpr z);
  struct cxpr cxdiv (struct cxpr z1, struct cxpr z2);

  struct cxpr cxgdiv (struct cxpr z1, struct cxpr z2);
  struct cxpr cxidiv (struct cxpr z1, struct cxpr z2);
  struct cxpr cxgmod (struct cxpr z1, struct cxpr z2);
  struct cxpr cxmod (struct cxpr z1, struct cxpr z2);
  struct cxpr cxpwr (struct cxpr z, int n);
  struct cxpr cxsqr (struct cxpr z);
  struct cxpr cxpow (struct cxpr z1, struct cxpr z2);
  struct cxpr cxroot (struct cxpr z, int i, int n);
  struct cxpr cxsqrt (struct cxpr z);

  struct cxprcmp_res cxprcmp (const struct cxpr* z1, const struct cxpr* z2);
  int cxis0 (const struct cxpr* z);
  int cxnot0 (const struct cxpr* z);
  int cxeq (struct cxpr z1, struct cxpr z2);
  int cxneq (struct cxpr z1, struct cxpr z2);
  int cxgt (struct cxpr z1, struct cxpr z2);
  int cxge (struct cxpr z1, struct cxpr z2);
  int cxlt (struct cxpr z1, struct cxpr z2);
  int cxle (struct cxpr z1, struct cxpr z2);

  struct cxpr cxconj (struct cxpr z);
  struct cxpr cxneg (struct cxpr z);
  struct cxpr cxinv (struct cxpr z);

  struct cxpr cxexp (struct cxpr z);
  struct cxpr cxexp10 (struct cxpr z);
  struct cxpr cxexp2 (struct cxpr z);
  struct cxpr cxlog (struct cxpr z);
  struct cxpr cxlog10 (struct cxpr z);
  struct cxpr cxlog2 (struct cxpr z);
  struct cxpr cxlog_sqrt (struct cxpr z);
  struct cxpr cxsin (struct cxpr z);
  struct cxpr cxcos (struct cxpr z);
  struct cxpr cxtan (struct cxpr z);
  struct cxpr cxsinh (struct cxpr z);
  struct cxpr cxcosh (struct cxpr z);
  struct cxpr cxtanh (struct cxpr z);
  struct cxpr cxasin (struct cxpr z);
  struct cxpr cxacos (struct cxpr z);
  struct cxpr cxatan (struct cxpr z);
  struct cxpr cxasinh (struct cxpr z);
  struct cxpr cxacosh (struct cxpr z);
  struct cxpr cxatanh (struct cxpr z);

  struct cxpr cxfloor (struct cxpr z);
  struct cxpr cxceil (struct cxpr z);
  struct cxpr cxround (struct cxpr z);
  struct cxpr cxtrunc (struct cxpr z);
  struct cxpr cxfrac (struct cxpr z);
  struct cxpr cxfix (struct cxpr z);

/* Conversion's functions */
  struct cxpr strtocx (const char *q, char **endptr);
  struct cxpr atocx (const char *s);
  char *cxpr_asprint (struct cxpr z, int sc_not, int sign, int lim);
  char *cxtoa (struct cxpr z, int lim);
  struct cxpr dctocx (double re, double im);
  struct cxpr fctocx (float re, float im);
  struct cxpr ictocx (long re, long im);
  struct cxpr uctocx (unsigned long re, unsigned long im);
  void cxtodc (const struct cxpr *z, double *re, double *im);
  void cxtofc (const struct cxpr *z, float *re, float *im);

/* Output functions */

#define CX1I_CHAR 'i'
#define CX1I_STR  "i"

  void cxpr_print (FILE * stream, struct cxpr z, int sc_not, int sign,
		   int lim);
  void cxprcxpr (struct cxpr z, int m);
  void cxprint (FILE * stream, struct cxpr z);

/* Special output functions and related macros */

#define XFMT_STD       0
#define XFMT_RAW       1
#define XFMT_ALT       2

#define CXDEF_LDEL   '('
#define CXDEF_RDEL   ')'
#define CX_SEPARATOR ", "	/* TO BE USED WITH THE ALT FORMAT */
#define CX_EMPTY_SEP "  "	/* TO BE USED WITH THE RAW FORMAT */
#define CX_SEP_L     2		/* LENGTH OF THE SEPARATOR        */

  int cxfout (FILE * pf, struct xoutflags ofs, struct cxpr z);
  int cxout (struct xoutflags ofs, struct cxpr z);
  unsigned long
    cxsout (char *s, unsigned long n, struct xoutflags ofs, struct cxpr z);

#define CXRESET(re, im) (struct cxpr){re, im}
#define CXCONV(x) (struct cxpr){x, xZero}
#define CXRE(z) (z).re
#define CXIM(z) (z).im
#define CXSWAP(z) (struct cxpr){(z).im, (z).re}

#define cxconvert cxconv
#define cxdiff cxsub
#define cxprod cxmul
#define cxipow cxpwr

#ifdef __cplusplus
}
#endif
#endif				/* _CXPRE_H_ */
