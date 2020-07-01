// Header file for rk.c
#ifndef _RK_H_
#define _RK_H_

#include "nsode.h"

void rk1 (  real *                  t,
            real *                  x,
            real *                  xout,
            const real              h,
            const uint              numvars,
            DERIVATIVE_FUNCTION     f);

void rk2 (  real *                  t,
            real *                  x,
            real *                  xout,
            real *                  k,
            const real              h,
            const uint              numvars,
            DERIVATIVE_FUNCTION     f);

void rk4 (  real *                  t,
            real *                  x,
            real *                  xout,
            real *              	k13,
            real *                  k24,
            const real              h,
            const uint              numvars,
            DERIVATIVE_FUNCTION     f);

#endif
