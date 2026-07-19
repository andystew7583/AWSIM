/**
 * pressure.h
 *
 * Rigid-lid pressure correction interface for AWSIM.
 *
 */
#ifndef _AWSIM_PRESSURE_H_
#define _AWSIM_PRESSURE_H_

#include "defs.h"

bool initPressureSolvers ( uint    Nx,
                           uint    Ny,
                           real    Lx,
                           real    Ly,
                           real    dx,
                           real    dy,
                           real ** Hc,
                           real ** Hw,
                           real ** Hs,
                           bool    useMG,
                           bool    useFullMG,
                           bool    useWallEW,
                           bool    useWallNS,
                           real    tol,
                           uint    maxIters,
                           real    omegaMG,
                           bool    debugOutput );
uint surfPressure (real *** uu, real *** vv, real *** hh, real dt, real ** pi, real rp, bool update_diags);
real optimizeSOR (real *** uu, real *** vv, real *** hh, real ** pi, real *** uu_buf, real *** vv_buf, real ** pi_buf, real dt, real rp);

#endif
