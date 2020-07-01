// Header file for ab.c
#ifndef _AB_H_
#define _AB_H_

#include "nsode.h"
#include <string.h>

void ab1 (  real *                  t,
            real *                  x,
            real *                  xout,
            real *                  dxdt,
            const real              h,
            const uint              numvars,
            DERIVATIVE_FUNCTION     f);

void ab2 (  real *                  t,
            real *                  x,
            real *                  xout,
            real *                  dxdt,
            real *                  dxdt_1,
            const real              h,
            const uint              numvars,
            DERIVATIVE_FUNCTION     f);

void ab3 (  real *                  t,
            real *                  x,
            real *                  xout,
            real *                  dxdt,
            real *              	dxdt_1,
            real *                  dxdt_2,
            const real              h,
            const uint              numvars,
            DERIVATIVE_FUNCTION     f);

void ab4 (  real *                  t,
            real *                  x,
            real *                  xout,
            real *                  dxdt,
            real *              	dxdt_1,
            real *                  dxdt_2,
            real *              	dxdt_3,
            const real              h,
            const uint              numvars,
            DERIVATIVE_FUNCTION     f);

#endif
