/**
 * multigrid.h
 *
 * MultiGrid pressure solver interface for AWSIM.
 *
 */
#ifndef _AWSIM_MULTIGRID_H_
#define _AWSIM_MULTIGRID_H_

#include "defs.h"

bool init_MG ( uint    Nx,
               uint    Ny,
               real    Lx,
               real    Ly,
               real ** Hc,
               bool    useWallEW,
               bool    useWallNS,
               bool    useFullMG,
               real    tol,
               uint    maxIters,
               real    omega,
               bool    debugOutput );
uint solve_MG (real ** pi, real ** pi_rhs);

#endif
