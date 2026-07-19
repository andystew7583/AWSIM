/**
 * sor.h
 *
 * Successive over-relaxation pressure solver interface for AWSIM.
 *
 */
#ifndef _AWSIM_SOR_H_
#define _AWSIM_SOR_H_

#include "defs.h"

bool init_SOR ( uint    Nx,
                uint    Ny,
                real    dx,
                real    dy,
                real ** Hw,
                real ** Hs,
                bool    useWallEW,
                bool    useWallNS,
                real    tol,
                uint    maxIters,
                bool    debugOutput );
uint solve_SOR (real ** pi, real ** pi_rhs, real rp);

#endif
