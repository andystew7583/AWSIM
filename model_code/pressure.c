/**
 * pressure.c
 *
 * Rigid-lid pressure correction workflow for AWSIM.
 *
 */
#include <math.h>

#include "pressure.h"
#include "multigrid.h"
#include "sor.h"

// Avoids memory errors associated with usual definition of bool
typedef int mybool;

extern const mybool debug;
extern uint Nx;
extern uint Ny;
extern uint Nlay;
extern uint N;
extern real dx;
extern real dy;
extern mybool useWallEW;
extern mybool useWallNS;
extern mybool use_MG;
extern real pi_tol;
extern uint maxiters;
extern real rp_acc_max;
extern real rp_opt_max;
extern real rp_opt_min;
extern uint N_rp;

extern real ** pi_rhs;
extern real *** h_west;
extern real *** h_south;
extern real avg_fac_hu;
extern real avg_fac_hv;
extern real avg_fac_e;
extern real dt_avg_hu;
extern real dt_avg_hv;
extern real dt_avg_e;
extern real *** hu_tend_gradM;
extern real *** hv_tend_gradM;
extern real *** e_tend_gradM;
extern real *** e_flux_uP;
extern real *** e_flux_vP;

void calcFaceThickness (real ** hh, real ** h_west, real ** h_south, mybool use_ghost, uint Nx, uint Ny, real ** uu, real ** vv);

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
                           bool    debugOutput )
{
  if (useMG)
  {
    return init_MG(Nx,Ny,Lx,Ly,Hc,useWallEW,useWallNS,useFullMG,tol,maxIters,omegaMG,debugOutput);
  }

  return init_SOR(Nx,Ny,dx,dy,Hw,Hs,useWallEW,useWallNS,tol,maxIters,debugOutput);
}

/**
 * surfPressure
 *
 * Calculates the surface pressure from the intermediate velocity. Or, equivalently,
 * computes a barotropic correction to the velocity field that renders it non-divergent
 * in a depth-integral sense. This correction is then added to the 3D velocity field,
 * i.e. uu and vv are modified by this function.
 *
 * The derivatives are approximated using central differences, and MultiGrid or Successive
 * Over-Relaxation is used to solve the Poisson equation
 *   div (H grad pi) = div (Hu)
 * where H is the ocean depth and u is the depth-averaged velocity. The boundary 
 * condition on lateral boundaries is 
 *   grad pi.n = u.n = 0.
 * See Rempfer (2006) for further details on the boundary condition.
 *
 * The matrix pi is used as a prior for the iterative procedure, and is modified to
 * store the updated value of pi when this function returns.
 *
 * The flag update_diags tells this function whether to update momentum and energy
 * diagnostics once the surface pressure has been calculated.
 *
 * Returns the number of iterations required to achieve convergence,
 * or 0 if the method did not converge.
 *
 */
uint surfPressure (real *** uu, real *** vv, real *** hh, real dt, real ** pi, real rp, bool update_diags)
{
  // Volume fluxes at cell faces
  real hu = 0;
  real hv = 0;
  real rhs_u = 0;
  real rhs_v = 0;
  
  // Number of iterations required to converge
  uint iters = 0;
  
  // Looping variables
  int i,j,k,im1,imin,jm1,jmin;
  
  // Set right-hand side of Poisson equation
  memset(*pi_rhs,0,Nx*Ny*sizeof(real));
  
  for (k = 0; k < Nlay; k ++)
  {

    
    // Calculate layer thickness on cell faces
    // N.B. Here we use h_west and h_south, which are Nlay x Nx+2*Ng x Ny+2*Ng matrices,
    // as Nlay x Nx x Ny matrices.
    // NOTE: If using advection schemes in which the face thickensses depend on the velocities,
    // following the pressure solve the velocities will be adjusted, which may in turn modify
    // the layer thicknesses when they are next calculated. Thus, the resulting volume fluxes
    // may not be perfectly nondivergent. Deviations in the water column thickness will be
    // corrected using correctThickness after each time step.
    calcFaceThickness(hh[k],h_west[k],h_south[k],false,Nx,Ny,uu[k],vv[k]);
    
    // Add contribution due to x-volume fluxes
    imin = useWallEW ? 1 : 0;
    for (i = imin; i < Nx; i ++)
    {
      im1 = (i+Nx-1) % Nx;
      
      for (j = 0; j < Ny; j ++)
      {
        hu = uu[k][i][j]*h_west[k][i][j];
        pi_rhs[i][j] -= hu / (dx*dt);
        pi_rhs[im1][j] += hu / (dx*dt);
      }
    }

    
    // Add contribution due to y-volume fluxes
    jmin = useWallNS ? 1 : 0;
    for (j = jmin; j < Ny; j ++)
    {
      jm1 = (j+Ny-1) % Ny;
      for (i = 0; i < Nx; i ++)
      {
        hv = vv[k][i][j]*h_south[k][i][j];
        pi_rhs[i][j] -= hv / (dy*dt);
        pi_rhs[i][jm1] += hv / (dy*dt);
      }
    }
    
  }

  // Solve for surface pressure using selected scheme
  if (use_MG)
  {
    iters = solve_MG(pi,pi_rhs);
  }
  else
  {
    iters = solve_SOR(pi,pi_rhs,rp);
  }
  
  // Correct u-velocity
  for (k = 0; k < Nlay; k ++)
  {
    
    
    imin = useWallEW ? 1 : 0;
    for (i = imin; i < Nx; i ++)
    {
      im1 = (i+Nx-1) % Nx;
      for (j = 0; j < Ny; j ++)
      {
        rhs_u = - dt*(pi[i][j]-pi[im1][j])/dx;
        uu[k][i][j] += rhs_u;
        if (update_diags && (dt_avg_hu > 0))
        {
          hu_tend_gradM[k][i][j] += h_west[k][i][j] * avg_fac_hu * rhs_u;
        }
        if (update_diags && (dt_avg_e > 0))
        {
          e_tend_gradM[k][i][j] += rhs_u*h_west[k][i][j]*uu[k][i][j]*avg_fac_e;
        }
      }
    }
    
  }

  // Correct v-velocity
  for (k = 0; k < Nlay; k ++)
  {
    
    
    jmin = useWallNS ? 1 : 0;
    for (j = jmin; j < Ny; j ++)
    {
      jm1 = (j+Ny-1) % Ny;
      for (i = 0; i < Nx; i ++)
      {
        rhs_v = - dt*(pi[i][j]-pi[i][jm1])/dy;
        vv[k][i][j] += rhs_v;
        if (update_diags && (dt_avg_hv > 0))
        {
          hv_tend_gradM[k][i][j] += h_south[k][i][j] * avg_fac_hv * rhs_v;
        }
        if (update_diags && (dt_avg_e > 0))
        {
          e_tend_gradM[k][i][j] += rhs_v*h_south[k][i][j]*vv[k][i][j]*avg_fac_e;
        }
      }
    }
    
  }
  
  // Additional energy budget terms
  if (update_diags && (dt_avg_e > 0))
  {
    for (k = 0; k < Nlay; k ++)
    {
      for (j = 0; j < Ny; j ++)
      {
        jm1 = (j+Ny-1) % Ny;
        
        for (i = 0; i < Nx; i ++)
        {
          im1 = (i+Nx-1) % Nx;
          
          // Barotropic contribution to pressure flux
          e_flux_uP[k][i][j] += h_west[k][i][j] * uu[k][i][j] * 0.5*(pi[i][j]+pi[im1][j]) * avg_fac_e*dt;
          e_flux_vP[k][i][j] += h_south[k][i][j] * vv[k][i][j]  * 0.5*(pi[i][j]+pi[i][jm1]) * avg_fac_e*dt;
        }
      }
    }
  }
  
  // If we haven't converged in the required number
  // of iterations, return 0 to indicate this
  if (iters == maxiters)
  {
    return 0;
  }
  else
  {
    return iters;
  }
}












/**
 *
 * optimizeSOR
 *
 * Optimizes the relaxation parameter to be used for the surface pressure calculation. Takes the current 
 * iteration's variables (uu,vv,hh) and buffers (uu_buf,vv_buf) for temporary storage, along with a guess
 * for the surface pressure (pi). Note that uu_buf, vv_buf and pi will be modified within this function.
 * The optimized rp is returned.
 *
 */
real optimizeSOR (real *** uu, real *** vv, real *** hh, real ** pi, real *** uu_buf, real *** vv_buf, real ** pi_buf, real dt, real rp)
{
  real rp_min = 0;
  real rp_max = 0;
  real rp_acc = 0;
  uint n_rp = 0;
  uint iters_rp = 0;
  uint n_miniters_rp = 0;
  real miniters_rp = 0;
  real drp = 0;
  real rp_opt_range = rp_opt_max - rp_opt_min;

  if (debug)
  {
    printf("BEGINNING OPTIMIZATION\n");
    fflush(stdout);
  }
  
  // This is a crude iteration scheme to optimize the relaxation parameter. We
  // repeatedly call surfPressure over a range of values of rp, successively
  // narrowing our search window until we converge.
  rp_min = rp_opt_min;
  rp_max = rp_opt_max;
  rp_acc = rp_max - rp_min;
  while (rp_acc > rp_acc_max)
  {
    // rp increment
    drp = (rp_max-rp_min)/N_rp;
    
    // Iterate through values of rp to find the one that yields the fewest iterations
    // required for convergence
    miniters_rp = maxiters;
    for (n_rp = 1; n_rp < N_rp; n_rp ++)
    {
      // Relaxation parameter to test
      rp = rp_min + n_rp*drp;
      
      // Configure the pressure solve input identically for each iteration
      memcpy(*(*uu_buf),*(*uu),N*sizeof(real));
      memcpy(*(*vv_buf),*(*vv),N*sizeof(real));
      memcpy(*pi_buf,*pi,Nx*Ny*sizeof(real));
      
      if (debug)
      {
        printf("Trying rp=%lf\n",rp);
        fflush(stdout);
      }
      
      // Determine how many iterations were required to solve for the pressure
      iters_rp = surfPressure(uu_buf,vv_buf,hh,dt,pi_buf,rp,false);
      if (iters_rp == 0)
      {
        iters_rp = maxiters;
      }
      
      // Find the value of rp that yields the minimum number of iterations
      if (iters_rp < miniters_rp)
      {
        miniters_rp = iters_rp;
        n_miniters_rp = n_rp;
      }
    }
    
    // Define new rp range over which to search, centered on the minimum of rp
    rp_max = rp_min + (n_miniters_rp+1)*drp;
    rp_min = rp_min + (n_miniters_rp-1)*drp;
    
    // Calculate accuracy with which rp is currently constrained
    rp_acc = 0.5*(rp_max-rp_min);
    
    // This code adapts the user-specified range of rp to search beyond the limits if the optimal rp
    // lies at the edges of the range
    if ((rp_acc <= rp_acc_max))
    {
      
      // Optimal rp is at the bottom of the range
      if ((rp_min-rp_opt_min <= rp_acc_max) && (rp_opt_min-1.0 > rp_acc_max))
      {
        if (debug)
        {
          printf("DECREASING rp_opt_min\n");
          fflush(stdout);
        }
        
        rp_opt_min -= fmin(rp_opt_range/2,rp_opt_min-1.0);
        rp_opt_max = rp_opt_min + rp_opt_range;
        rp_min = rp_opt_min;
        rp_max = rp_opt_max;
        rp_acc = rp_opt_range;
      }
      
      // Optimal rp is at the top of the range
      else if ((rp_opt_max-rp_max <= rp_acc_max) && (2.0-rp_opt_max > rp_acc_max))
      {
        if (debug)
        {
          printf("INCREASING rp_opt_max\n");
          fflush(stdout);
        }
        
        rp_opt_max += fmin(rp_opt_range/2,2.0-rp_opt_max);
        rp_opt_min = rp_opt_max - rp_opt_range;
        rp_min = rp_opt_min;
        rp_max = rp_opt_max;
        rp_acc = rp_opt_range;
      }
      
    }
    
  }
  
  // Finally, determine the optimal rp
  rp = 0.5*(rp_min+rp_max);
  
  if (debug)
  {
    printf("OPTIMIZATION COMPLETE: new rp=%lf\n",rp);
    fflush(stdout);
  }
  
  return rp;
}
