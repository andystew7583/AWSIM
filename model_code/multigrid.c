/**
 * multigrid.c
 *
 * MultiGrid pressure solver implementation for AWSIM.
 *
 */
#include <math.h>
#include <time.h>

#include "multigrid.h"

typedef enum transfer_MG
{
  TRANSFER_POWER2_MG,
  TRANSFER_GENERAL_MG
}
transfer_MG;

// Data storage structure for MultiGrid scheme
typedef struct data_MG
{
  real ** pi; // Nx x Ny matrix containing current solution
  real ** pi_temp; // Nx x Ny temporary storage matrix
  real ** pi_rhs; // Nx x Ny matrix containing right-hand side of weighted Poisson equation

  uint Nx; // Grid size in x (first dimension)
  uint Ny; // Grid size in y (second dimension)
  real dx; // Grid spacing in x
  real dy; // Grid spacing in y

  real ** Hc; // Water column thickness on cell centers
  real ** Hw; // Water column thickness on cell western edges
  real ** Hs; // Water column thickness on cell southern edges

  real ** Ow; // Nx x Ny matrix containing east-west operators on western edges of grid cells
  real ** Os; // Nx x Ny matrix containing north-south operators on southern edges of grid cells
  real ** Osum; // Nx x Ny matrix containing operator sum around edges of each grid cell
  real ** _Osum; // Nx x Ny matrix containing reciprocal of operator sum around edges of each grid cell

  real ** wmm; // Nx x Ny interpolation weight matrices. Supply weights to be used in summing
  real ** wmp; // coarse-grid elements to the southwest (wmm), northwest (wmp), northeast (wpp)
  real ** wpm; // and southeast (wpm) of each fine-grid element.
  real ** wpp;

  transfer_MG transfer;

  real * im1_vec; // Vectors to store indices of adjacent grid points,
  real * ip1_vec; // pre-computed for efficiency
  real * jm1_vec;
  real * jp1_vec;
}
data_MG;

static data_MG * mg_grids = NULL;
static real * F_MG = NULL;
static uint Ngrids = 0;
static bool debug = false;
static bool use_fullMG = false;
static real pi_tol = 0;
static uint maxiters = 0;
static real omega_WJ = 0;
static bool use_wall_ew = false;
static bool use_wall_ns = false;

static void restrict_MG ( data_MG * fine_grid,
                          real **   pi_f,
                          data_MG * coarse_grid,
                          real **   pi_c );

/**
 * alloc_MG_matrix
 *
 * Allocates an Nx x Ny matrix for the MultiGrid solver.
 *
 * mat - Pointer to matrix to allocate.
 * Nx - Grid size in x (first dimension).
 * Ny - Grid size in y (second dimension).
 *
 * Returns true if allocation succeeds, otherwise returns false.
 *
 */
static bool alloc_MG_matrix (real *** mat, uint Nx, uint Ny)
{
  *mat = matalloc(Nx,Ny);
  if (*mat == NULL)
  {
    fprintf(stderr,"ERROR: Unable to allocate memory\r\n");
    return false;
  }

  return true;
}

/**
 * alloc_MG_vector
 *
 * Allocates a vector for the MultiGrid solver.
 *
 * vec - Pointer to vector to allocate.
 * N - Vector length.
 *
 * Returns true if allocation succeeds, otherwise returns false.
 *
 */
static bool alloc_MG_vector (real ** vec, uint N)
{
  *vec = vecalloc(N);
  if (*vec == NULL)
  {
    fprintf(stderr,"ERROR: Unable to allocate memory\r\n");
    return false;
  }

  return true;
}

/**
 * coarsen_size_MG
 *
 * Calculates the size of the next coarser grid in one direction.
 *
 * N - Fine grid size in one direction.
 *
 * Returns the coarse grid size. Even grids are coarsened exactly by a factor of
 * two, while odd grids retain the extra point on the coarser grid.
 *
 */
static uint coarsen_size_MG (uint N)
{
  return (N > 1) ? (N+1)/2 : 1;
}

/**
 * power2_transfer_MG
 *
 * Checks whether a fine/coarse grid pair is an exact 2:1 coarsening.
 *
 * fine_grid - Fine MultiGrid level.
 * coarse_grid - Coarse MultiGrid level.
 *
 * Returns true if the fast power-of-two transfer operators may be used.
 *
 */
static bool power2_transfer_MG (data_MG * fine_grid, data_MG * coarse_grid)
{
  return (fine_grid->Nx == 2*coarse_grid->Nx) && (fine_grid->Ny == 2*coarse_grid->Ny);
}

/**
 * calcFaceThickness_MG
 *
 * Calculates water column thickness on cell faces for a MultiGrid level.
 *
 * mg_grid - MultiGrid level containing cell-center and face thickness matrices.
 *
 */
static void calcFaceThickness_MG (data_MG * mg_grid)
{
  uint i,j,im1,jm1;

  for (i = 0; i < mg_grid->Nx; i ++)
  {
    im1 = (i+mg_grid->Nx-1) % mg_grid->Nx;

    for (j = 0; j < mg_grid->Ny; j ++)
    {
      jm1 = (j+mg_grid->Ny-1) % mg_grid->Ny;

      mg_grid->Hw[i][j] = 0.5 * (mg_grid->Hc[i][j]+mg_grid->Hc[im1][j]);
      mg_grid->Hs[i][j] = 0.5 * (mg_grid->Hc[i][j]+mg_grid->Hc[i][jm1]);
    }
  }
}

/**
 * init_MG
 *
 * Initializes data and parameters required by the MultiGrid solver.
 *
 * Nx - Grid size in x (first dimension).
 * Ny - Grid size in y (second dimension).
 * Lx - Domain size in x.
 * Ly - Domain size in y.
 * Hc - Nx x Ny matrix containing water column thickness on cell centers.
 * useWallEW - Whether to apply wall boundary conditions in x.
 * useWallNS - Whether to apply wall boundary conditions in y.
 * useFullMG - Whether to use full MultiGrid iteration.
 * tol - Error tolerance for the pressure solve.
 * maxIters - Maximum number of MultiGrid iterations.
 * omega - Weighted Jacobi relaxation parameter.
 * debugOutput - Whether to write debug output to stdout.
 *
 * Returns true if initialization succeeds, otherwise returns false.
 *
 */
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
               bool    debugOutput )
{
  uint Nx_m,Ny_m,F_len;
  int m;
  uint i,j,ip1,jp1;

  use_fullMG = useFullMG;
  pi_tol = tol;
  maxiters = maxIters;
  omega_WJ = omega;
  debug = debugOutput;
  use_wall_ew = useWallEW;
  use_wall_ns = useWallNS;

  // Count grid levels by coarsening until the hierarchy reaches a single cell.
  Ngrids = 1;
  Nx_m = Nx;
  Ny_m = Ny;
  while ((Nx_m > 1) || (Ny_m > 1))
  {
    Nx_m = coarsen_size_MG(Nx_m);
    Ny_m = coarsen_size_MG(Ny_m);
    Ngrids ++;
  }

  // Allocate workspace for the exact solve used on the coarsest 1-D grid.
  F_len = (Nx > Ny) ? Nx : Ny;
  if (!alloc_MG_vector(&F_MG,F_len+1))
  {
    return false;
  }

  mg_grids = (data_MG *) malloc(Ngrids*sizeof(data_MG));
  if (mg_grids == NULL)
  {
    fprintf(stderr,"ERROR: Unable to allocate memory\r\n");
    return false;
  }

  // Construct grid sizes from finest to coarsest, stored so index Ngrids-1 is finest.
  Nx_m = Nx;
  Ny_m = Ny;
  for (m = Ngrids-1; m >=0; m --)
  {
    mg_grids[m].Nx = Nx_m;
    mg_grids[m].Ny = Ny_m;
    mg_grids[m].dx = Lx/mg_grids[m].Nx;
    mg_grids[m].dy = Ly/mg_grids[m].Ny;
    mg_grids[m].transfer = TRANSFER_GENERAL_MG;

    Nx_m = coarsen_size_MG(Nx_m);
    Ny_m = coarsen_size_MG(Ny_m);
  }

  // Mark each fine/coarse pair for the fastest valid transfer operator.
  for (m = Ngrids-1; m > 0; m --)
  {
    mg_grids[m].transfer = power2_transfer_MG(&mg_grids[m],&mg_grids[m-1]) ? TRANSFER_POWER2_MG : TRANSFER_GENERAL_MG;
  }

  if (debug)
  {
    printf("Using %s MultiGrid hierarchy\n", (ispow2u(Nx) && ispow2u(Ny)) ? "fast power-of-two" : "general");
    fflush(stdout);
  }

  for (m = Ngrids-1; m >=0; m --)
  {
    // Allocate all solution, operator, and transfer arrays for this grid level.
    if ( !alloc_MG_matrix(&mg_grids[m].pi_temp,mg_grids[m].Nx,mg_grids[m].Ny)
      || !alloc_MG_matrix(&mg_grids[m].Hc,mg_grids[m].Nx,mg_grids[m].Ny)
      || !alloc_MG_matrix(&mg_grids[m].Hw,mg_grids[m].Nx,mg_grids[m].Ny)
      || !alloc_MG_matrix(&mg_grids[m].Hs,mg_grids[m].Nx,mg_grids[m].Ny)
      || !alloc_MG_matrix(&mg_grids[m].Ow,mg_grids[m].Nx,mg_grids[m].Ny)
      || !alloc_MG_matrix(&mg_grids[m].Os,mg_grids[m].Nx,mg_grids[m].Ny)
      || !alloc_MG_matrix(&mg_grids[m].Osum,mg_grids[m].Nx,mg_grids[m].Ny)
      || !alloc_MG_matrix(&mg_grids[m]._Osum,mg_grids[m].Nx,mg_grids[m].Ny)
      || !alloc_MG_matrix(&mg_grids[m].wmm,mg_grids[m].Nx,mg_grids[m].Ny)
      || !alloc_MG_matrix(&mg_grids[m].wpm,mg_grids[m].Nx,mg_grids[m].Ny)
      || !alloc_MG_matrix(&mg_grids[m].wmp,mg_grids[m].Nx,mg_grids[m].Ny)
      || !alloc_MG_matrix(&mg_grids[m].wpp,mg_grids[m].Nx,mg_grids[m].Ny)
      || !alloc_MG_vector(&mg_grids[m].im1_vec,mg_grids[m].Nx)
      || !alloc_MG_vector(&mg_grids[m].ip1_vec,mg_grids[m].Nx)
      || !alloc_MG_vector(&mg_grids[m].jm1_vec,mg_grids[m].Ny)
      || !alloc_MG_vector(&mg_grids[m].jp1_vec,mg_grids[m].Ny) )
    {
      return false;
    }

    if (m < Ngrids-1)
    {
      if ( !alloc_MG_matrix(&mg_grids[m].pi,mg_grids[m].Nx,mg_grids[m].Ny)
        || !alloc_MG_matrix(&mg_grids[m].pi_rhs,mg_grids[m].Nx,mg_grids[m].Ny) )
      {
        return false;
      }
    }
    else
    {
      mg_grids[m].pi = NULL;
      mg_grids[m].pi_rhs = NULL;
    }

    // Restrict water column thickness onto all coarser grids.
    if (m == Ngrids-1)
    {
      memcpy(*(mg_grids[m].Hc),*Hc,Nx*Ny*sizeof(real));
    }
    else
    {
      restrict_MG(&mg_grids[m+1],mg_grids[m+1].Hc,&mg_grids[m],mg_grids[m].Hc);
    }

    calcFaceThickness_MG(&mg_grids[m]);

    // Construct the weighted Poisson operator on southern cell faces.
    for (i = 0; i < mg_grids[m].Nx; i ++)
    {
      for (j = 0; j < mg_grids[m].Ny; j ++)
      {
        mg_grids[m].Os[i][j] = mg_grids[m].Hs[i][j] / SQUARE(mg_grids[m].dy);
      }

      if (useWallNS)
      {
        mg_grids[m].Os[i][0] = 0;
      }
    }

    // Construct the weighted Poisson operator on western cell faces.
    for (j = 0; j < mg_grids[m].Ny; j ++)
    {
      for (i = 0; i < mg_grids[m].Nx; i ++)
      {
        mg_grids[m].Ow[i][j] = mg_grids[m].Hw[i][j] / SQUARE(mg_grids[m].dx);
      }

      if (useWallEW)
      {
        mg_grids[m].Ow[0][j] = 0;
      }
    }

    // Pre-compute the operator sum and reciprocal for Jacobi relaxation.
    for (i = 0; i < mg_grids[m].Nx; i ++)
    {
      ip1 = (i + mg_grids[m].Nx + 1) % mg_grids[m].Nx;

      for (j = 0; j < mg_grids[m].Ny; j ++)
      {
        jp1 = (j + mg_grids[m].Ny + 1) % mg_grids[m].Ny;

        mg_grids[m].Osum[i][j] = mg_grids[m].Ow[i][j] + mg_grids[m].Ow[ip1][j] + mg_grids[m].Os[i][j] + mg_grids[m].Os[i][jp1];
        mg_grids[m]._Osum[i][j] = 1 / mg_grids[m].Osum[i][j];
      }
    }

    // Pre-compute periodic neighboring indices used by the operator stencils.
    for (i = 0; i < mg_grids[m].Nx; i ++)
    {
      mg_grids[m].im1_vec[i] = (i+mg_grids[m].Nx-1) % mg_grids[m].Nx;
      mg_grids[m].ip1_vec[i] = (i+mg_grids[m].Nx+1) % mg_grids[m].Nx;
    }
    for (j = 0; j < mg_grids[m].Ny; j ++)
    {
      mg_grids[m].jm1_vec[j] = (j+mg_grids[m].Ny-1) % mg_grids[m].Ny;
      mg_grids[m].jp1_vec[j] = (j+mg_grids[m].Ny+1) % mg_grids[m].Ny;
    }

    // Pre-compute interpolation weights for the exact 2:1 transfer operator.
    for (i = 0; i < mg_grids[m].Nx; i ++)
    {
      for (j = 0; j < mg_grids[m].Ny; j ++)
      {
        mg_grids[m].wmm[i][j] = ((i%2==1) ? 0.75 : 0.25) * ((j%2==1) ? 0.75 : 0.25);
        mg_grids[m].wmp[i][j] = ((i%2==1) ? 0.75 : 0.25) * ((j%2==1) ? 0.25 : 0.75);
        mg_grids[m].wpm[i][j] = ((i%2==1) ? 0.25 : 0.75) * ((j%2==1) ? 0.75 : 0.25);
        mg_grids[m].wpp[i][j] = ((i%2==1) ? 0.25 : 0.75) * ((j%2==1) ? 0.25 : 0.75);
      }
    }

    // Modify interpolation weights next to northern/southern walls.
    if (useWallNS)
    {
      for (i = 0; i < mg_grids[m].Nx; i ++)
      {
        mg_grids[m].wmm[i][0] = 0;
        mg_grids[m].wpm[i][0] = 0;
        mg_grids[m].wmp[i][0] = ((i%2==1) ? 0.75 : 0.25);
        mg_grids[m].wpp[i][0] = ((i%2==1) ? 0.25 : 0.75);
        mg_grids[m].wmp[i][mg_grids[m].Ny-1] = 0;
        mg_grids[m].wpp[i][mg_grids[m].Ny-1] = 0;
        mg_grids[m].wmm[i][mg_grids[m].Ny-1] = ((i%2==1) ? 0.75 : 0.25);
        mg_grids[m].wpm[i][mg_grids[m].Ny-1] = ((i%2==1) ? 0.25 : 0.75);
      }
    }

    // Modify interpolation weights next to eastern/western walls.
    if (useWallEW)
    {
      for (j = 0; j < mg_grids[m].Ny; j ++)
      {
        mg_grids[m].wmm[0][j] = 0;
        mg_grids[m].wmp[0][j] = 0;
        mg_grids[m].wpm[0][j] = ((j%2==1) ? 0.75 : 0.25);
        mg_grids[m].wpp[0][j] = ((j%2==1) ? 0.25 : 0.75);
        mg_grids[m].wpm[mg_grids[m].Nx-1][j] = 0;
        mg_grids[m].wpp[mg_grids[m].Nx-1][j] = 0;
        mg_grids[m].wmm[mg_grids[m].Nx-1][j] = ((j%2==1) ? 0.75 : 0.25);
        mg_grids[m].wmp[mg_grids[m].Nx-1][j] = ((j%2==1) ? 0.25 : 0.75);
      }
    }

    // Modify corner interpolation weights when both boundary directions are walled.
    if (useWallNS && useWallEW)
    {
      mg_grids[m].wmm[0][0] = 0;
      mg_grids[m].wmp[0][0] = 0;
      mg_grids[m].wpm[0][0] = 0;
      mg_grids[m].wpp[0][0] = 1;

      mg_grids[m].wmm[0][mg_grids[m].Ny-1] = 0;
      mg_grids[m].wmp[0][mg_grids[m].Ny-1] = 0;
      mg_grids[m].wpm[0][mg_grids[m].Ny-1] = 1;
      mg_grids[m].wpp[0][mg_grids[m].Ny-1] = 0;

      mg_grids[m].wmm[mg_grids[m].Nx-1][0] = 0;
      mg_grids[m].wmp[mg_grids[m].Nx-1][0] = 1;
      mg_grids[m].wpm[mg_grids[m].Nx-1][0] = 0;
      mg_grids[m].wpp[mg_grids[m].Nx-1][0] = 0;

      mg_grids[m].wmm[mg_grids[m].Nx-1][mg_grids[m].Ny-1] = 1;
      mg_grids[m].wmp[mg_grids[m].Nx-1][mg_grids[m].Ny-1] = 0;
      mg_grids[m].wpm[mg_grids[m].Nx-1][mg_grids[m].Ny-1] = 0;
      mg_grids[m].wpp[mg_grids[m].Nx-1][mg_grids[m].Ny-1] = 0;
    }
  }

  return true;
}

/**
 * exactSolve_MG
 *
 * Constructs an exact solution to the Poisson equation for a vector.
 *
 * Nx - Grid size in x (first dimension). Must be 1 if Ny>1.
 * Ny - Grid size in y (second dimension). Must be 1 if Nx>1.
 * pi - Nx x Ny matrix to store solution
 * pi_rhs - Nx x Ny matrix containing right-hand side of weighted Poisson equation
 * Ow - Nx x Ny matrix containing east-west operators on western edges of grid cells
 * Os - Nx x Ny matrix containing north-south operators on southern edges of grid cells
 *
 */
static void exactSolve_MG (  uint      Nx,
                             uint      Ny,
                             real **   pi,
                             real **   pi_rhs,
                             real **   Ow,
                             real **   Os  )
{
  // Looping variables
  uint i,j;
  
  // Nothing to do in this case
  if (Nx == 1 && Ny == 1)
  {
    pi[0][0] = 0;
    return;
  }
  
  // Grid construction ensures that either Nx==1 and/or Ny==1 on smallest grid
  if (Nx == 1)
  {
    // Integrate RHS once to get "fluxes" across cell southern edges
    F_MG[0] = 0;
    for (j = 1; j <= Ny; j ++)
    {
      F_MG[j] = F_MG[j-1] + pi_rhs[0][j-1];
    }
    
    // Integrate "flux" to get pressure in each grid cell
    pi[0][0] = 0;
    for (j = 1; j < Ny; j ++)
    {
      pi[0][j] = pi[0][j-1] + F_MG[j]/Os[0][j];
    }
  }
  else
  {
    // Integrate RHS once to get "fluxes" across cell southern edges
    F_MG[0] = 0;
    for (i = 1; i <= Nx; i ++)
    {
      F_MG[i] = F_MG[i-1] + pi_rhs[i-1][0];
    }
    
    // Integrate "flux" to get pressure in each grid cell
    pi[0][0] = 0;
    for (i = 1; i < Nx; i ++)
    {
      pi[i][0] = pi[i-1][0] + F_MG[i]/Ow[i][0];
    }
  }
}









/**
 * relax_MG
 *
 * Iterates the solution of the weighted Poisson equation once using a weighted Jacobi scheme.
 *
 * Nx - Grid size in x (first dimension)
 * Ny - Grid size in y (second dimension)
 * pi - Nx x Ny matrix containing current solution. Will be modified to store updated solution.
 * pi_prev - Nx x Ny matrix. Will be modified to store current solution.
 * pi_rhs - Nx x Ny matrix containing right-hand side of weighted Poisson equation
 * Ow - Nx x Ny matrix containing east-west operators on western edges of grid cells
 * Os - Nx x Ny matrix containing north-south operators on southern edges of grid cells
 * _Osum - Nx x Ny matrix containing reciprocal of operator sum around edges of each grid cell
 *
 */
static void relax_MG (   uint      Nx,
                         uint      Ny,
                         real **   pi,
                         real **   pi_prev,
                         real **   pi_rhs,
                         real **   Ow,
                         real **   Os,
                         real **   _Osum )
{
  uint i,j,im1,ip1,jm1,jp1;
  
  // Copy current solution to pi_prev matrix
  memcpy(*pi_prev,*pi,Nx*Ny*sizeof(real));
  
  
  // Loop through all indices and perform a single weighted Jacobi iteration, storing the result in pi
  for (i = 0; i < Nx; i ++)
  {
    im1 = (i+Nx-1) % Nx;
    ip1 = (i+Nx+1) % Nx;
    
    for (j = 0; j < Ny; j ++)
    {
      jm1 = (j+Ny-1) % Ny;
      jp1 = (j+Ny+1) % Ny;
      
      // N.B. This code is periodic in y, but the operator Os is set such that the wall BCs are included
      pi[i][j] = omega_WJ * _Osum[i][j] *  ( Os[i][jp1]*pi[i][jp1] + Os[i][j]*pi[i][jm1] + Ow[ip1][j]*pi[ip1][j] + Ow[i][j]*pi[im1][j] - pi_rhs[i][j] )
               + (1-omega_WJ) * pi[i][j];
      
    }
  }
}















/**
 * residual_MG
 *
 * Computes the residual (i.e. error) in the solution to the weighted Poisson equation.
 *
 * Nx - Grid size in x (first dimension)
 * Ny - Grid size in y (second dimension)
 * pi - Nx x Ny matrix containing current solution.
 * res - Nx x Ny matrix. Will be modified to store the residual.
 * pi_rhs - Nx x Ny matrix containing right-hand side of weighted Poisson equation
 * Ow - Nx x Ny matrix containing east-west operators on western edges of grid cells
 * Os - Nx x Ny matrix containing north-south operators on southern edges of grid cells
 * Osum - Nx x Ny matrix containing operator sum around edges of each grid cell
 *
 */
static void residual_MG (  uint      Nx,
                           uint      Ny,
                           real **   pi,
                           real **   res,
                           real **   pi_rhs,
                           real **   Ow,
                           real **   Os,
                           real **   Osum )
{
  // For looping
  uint i,j,im1,ip1,jm1,jp1;
  
  
  // Loop through all indices and compute residual between weighted Laplacian operator and right-hand side
  for (i = 0; i < Nx; i ++)
  {
    im1 = (i+Nx-1) % Nx;
    ip1 = (i+Nx+1) % Nx;
    
    for (j = 0; j < Ny; j ++)
    {
      jm1 = (j+Ny-1) % Ny;
      jp1 = (j+Ny+1) % Ny;
      
      // N.B. This code is periodic in y, but the operator Os is set such that the wall BCs are included
      res[i][j] =  Os[i][jp1]*pi[i][jp1] + Os[i][j]*pi[i][jm1] + Ow[ip1][j]*pi[ip1][j] + Ow[i][j]*pi[im1][j] - Osum[i][j]*pi[i][j] - pi_rhs[i][j];
    }
  }
}













/**
 * restrict_power2_MG
 *
 * Restriction operator for exact 2:1 grid pairs.
 *
 * fine_grid - Fine MultiGrid level.
 * pi_f - Fine-resolution matrix containing input data.
 * pi_c - Coarse-resolution matrix to store output data.
 *
 */
static void restrict_power2_MG ( data_MG * fine_grid,
                                 real **   pi_f,
                                 real **   pi_c )
{
  // For looping
  uint i,j;
  
  
  // Restriction operator is simple average of the four fine-grid elements that lie
  // within each coarse-grid cell
  for (i = 0; i < fine_grid->Nx/2; i ++)
  {
    for (j = 0; j < fine_grid->Ny/2; j ++)
    {
      pi_c[i][j] = 0.25 * (pi_f[2*i][2*j] + pi_f[2*i+1][2*j] + pi_f[2*i][2*j+1] + pi_f[2*i+1][2*j+1]);
    }
  }
}

/**
 * restrict_general_MG
 *
 * Restriction operator for grid pairs that are not exact 2:1 coarsenings.
 *
 * fine_grid - Fine MultiGrid level.
 * pi_f - Fine-resolution matrix containing input data.
 * coarse_grid - Coarse MultiGrid level.
 * pi_c - Coarse-resolution matrix to store output data.
 *
 */
static void restrict_general_MG ( data_MG * fine_grid,
                                  real **   pi_f,
                                  data_MG * coarse_grid,
                                  real **   pi_c )
{
  uint ic,jc,ifine,jfine;
  real x0,x1,y0,y1,xf0,xf1,yf0,yf1;
  real wx,wy,w,wsum;

  memset(*pi_c,0,coarse_grid->Nx*coarse_grid->Ny*sizeof(real));

  // Treat each coarse cell as a box in fine-grid index space.
  for (ic = 0; ic < coarse_grid->Nx; ic ++)
  {
    x0 = ((real) ic) * fine_grid->Nx / coarse_grid->Nx;
    x1 = ((real) (ic+1)) * fine_grid->Nx / coarse_grid->Nx;

    for (jc = 0; jc < coarse_grid->Ny; jc ++)
    {
      y0 = ((real) jc) * fine_grid->Ny / coarse_grid->Ny;
      y1 = ((real) (jc+1)) * fine_grid->Ny / coarse_grid->Ny;
      wsum = 0;

      // Sum all fine cells with nonzero area overlap with this coarse cell.
      for (ifine = (uint) floor(x0); ifine < (uint) ceil(x1); ifine ++)
      {
        xf0 = (real) ifine;
        xf1 = (real) (ifine+1);
        wx = fmin(x1,xf1) - fmax(x0,xf0);

        if (wx <= 0)
        {
          continue;
        }

        for (jfine = (uint) floor(y0); jfine < (uint) ceil(y1); jfine ++)
        {
          yf0 = (real) jfine;
          yf1 = (real) (jfine+1);
          wy = fmin(y1,yf1) - fmax(y0,yf0);

          if (wy <= 0)
          {
            continue;
          }

          w = wx * wy;
          pi_c[ic][jc] += w * pi_f[ifine][jfine];
          wsum += w;
        }
      }

      // Normalize by actual overlap area to preserve constants exactly.
      if (wsum > 0)
      {
        pi_c[ic][jc] /= wsum;
      }
    }
  }
}

/**
 * restrict_MG
 *
 * Restriction operator for MultiGrid solver. Takes a matrix and restricts it to
 * a coarser grid.
 *
 * fine_grid - Fine MultiGrid level.
 * pi_f - Fine-resolution matrix containing input data.
 * coarse_grid - Coarse MultiGrid level.
 * pi_c - Coarse-resolution matrix to store output data.
 *
 */
static void restrict_MG ( data_MG * fine_grid,
                          real **   pi_f,
                          data_MG * coarse_grid,
                          real **   pi_c )
{
  if (fine_grid->transfer == TRANSFER_POWER2_MG)
  {
    restrict_power2_MG(fine_grid,pi_f,pi_c);
  }
  else
  {
    restrict_general_MG(fine_grid,pi_f,coarse_grid,pi_c);
  }
}













/**
 * interp_index_MG
 *
 * Maps an interpolation index to a valid grid index.
 *
 * idx - Raw index supplied by the interpolation stencil.
 * N - Grid size in one direction.
 * use_wall - Whether this direction uses wall boundary conditions.
 *
 * Returns a valid array index, using clamping at walls and periodic wrapping
 * otherwise.
 *
 */
static uint interp_index_MG (int idx, uint N, bool use_wall)
{
  if (N <= 1)
  {
    return 0;
  }

  if (use_wall)
  {
    if (idx < 0)
    {
      return 0;
    }
    if (idx >= (int) N)
    {
      return N-1;
    }
    return (uint) idx;
  }

  return (uint) ((idx + (int) N) % (int) N);
}

/**
 * interpolate_power2_MG
 *
 * Interpolation operator for exact 2:1 grid pairs.
 *
 * fine_grid - Fine MultiGrid level.
 * pi_f - Fine-resolution matrix to store output data.
 * pi_c - Coarse-resolution matrix containing input data.
 *
 */
static void interpolate_power2_MG ( data_MG * fine_grid,
                                    real **   pi_f,
                                    real **   pi_c )
{
  // For looping
  uint i,j,im,ip,jm,jp;

  
  for (i = 0; i < fine_grid->Nx; i ++)
  {
    // Identify coarse-grid i-points to the east and west of this fine-grid i-point
    ip = (((i+1)-(i+1)%2)/2 + fine_grid->Nx/2) % (fine_grid->Nx/2);
    im = (ip-1 + fine_grid->Nx/2) % (fine_grid->Nx/2);
    
    for (j = 0; j < fine_grid->Ny; j ++)
    {
      // Identify coarse-grid j-points to the north and south of this fine-grid j-point
      jp = (((j+1)-(j+1)%2)/2 + fine_grid->Ny/2) % (fine_grid->Ny/2);
      jm = (jp-1 + fine_grid->Ny/2) % (fine_grid->Ny/2);
      
      // Interpolate using pre-calculated weights
      pi_f[i][j] = fine_grid->wmm[i][j]*pi_c[im][jm] + fine_grid->wpm[i][j]*pi_c[ip][jm] + fine_grid->wmp[i][j]*pi_c[im][jp] + fine_grid->wpp[i][j]*pi_c[ip][jp];
    }
  }
}

/**
 * interpolate_general_MG
 *
 * Interpolation operator for grid pairs that are not exact 2:1 coarsenings.
 *
 * fine_grid - Fine MultiGrid level.
 * pi_f - Fine-resolution matrix to store output data.
 * coarse_grid - Coarse MultiGrid level.
 * pi_c - Coarse-resolution matrix containing input data.
 *
 */
static void interpolate_general_MG ( data_MG * fine_grid,
                                     real **   pi_f,
                                     data_MG * coarse_grid,
                                     real **   pi_c )
{
  uint i,j,im,ip,jm,jp;
  int im_raw,jm_raw;
  real x,y,tx,ty;


  for (i = 0; i < fine_grid->Nx; i ++)
  {
    // Map the fine-grid cell center to coarse-grid cell-center coordinates.
    x = (((real) i) + 0.5) * coarse_grid->Nx / fine_grid->Nx - 0.5;

    // Clamp interpolation to the nearest coarse point at walls.
    if (use_wall_ew && (x <= 0))
    {
      im_raw = 0;
      tx = 0;
    }
    else if (use_wall_ew && (x >= coarse_grid->Nx-1))
    {
      im_raw = coarse_grid->Nx-1;
      tx = 0;
    }
    else
    {
      im_raw = (int) floor(x);
      tx = x - im_raw;
    }

    // Convert raw interpolation indices to valid array indices.
    im = interp_index_MG(im_raw,coarse_grid->Nx,use_wall_ew);
    ip = interp_index_MG(im_raw+1,coarse_grid->Nx,use_wall_ew);

    for (j = 0; j < fine_grid->Ny; j ++)
    {
      // Map the fine-grid cell center to coarse-grid cell-center coordinates.
      y = (((real) j) + 0.5) * coarse_grid->Ny / fine_grid->Ny - 0.5;

      // Clamp interpolation to the nearest coarse point at walls.
      if (use_wall_ns && (y <= 0))
      {
        jm_raw = 0;
        ty = 0;
      }
      else if (use_wall_ns && (y >= coarse_grid->Ny-1))
      {
        jm_raw = coarse_grid->Ny-1;
        ty = 0;
      }
      else
      {
        jm_raw = (int) floor(y);
        ty = y - jm_raw;
      }

      jm = interp_index_MG(jm_raw,coarse_grid->Ny,use_wall_ns);
      jp = interp_index_MG(jm_raw+1,coarse_grid->Ny,use_wall_ns);

      // Bilinear interpolation from the four neighboring coarse-grid values.
      pi_f[i][j] = (1-tx)*(1-ty)*pi_c[im][jm]
                 + tx*(1-ty)*pi_c[ip][jm]
                 + (1-tx)*ty*pi_c[im][jp]
                 + tx*ty*pi_c[ip][jp];
    }
  }
}

/**
 * interpolate_MG
 *
 * Interpolation operator for MultiGrid solver. Takes a matrix and interpolates
 * it to a finer grid.
 *
 * fine_grid - Fine MultiGrid level.
 * pi_f - Fine-resolution matrix to store output data.
 * coarse_grid - Coarse MultiGrid level.
 * pi_c - Coarse-resolution matrix containing input data.
 *
 */
static void interpolate_MG ( data_MG * fine_grid,
                             real **   pi_f,
                             data_MG * coarse_grid,
                             real **   pi_c )
{
  if (fine_grid->transfer == TRANSFER_POWER2_MG)
  {
    interpolate_power2_MG(fine_grid,pi_f,pi_c);
  }
  else
  {
    interpolate_general_MG(fine_grid,pi_f,coarse_grid,pi_c);
  }
}















/**
 * correct_MG
 *
 * Correction operator for MultiGrid solver. Subtracts correction from the current solution.
 *
 * Nx - Fine grid size in x (first dimension)
 * Ny - Fine grid size in y (second dimension)
 * pi - Nx x Ny matrix containing current solution
 * cor - Nx x Ny matrix containing correction
 *
 */
static void correct_MG ( uint      Nx,
                         uint      Ny,
                         real **   pi,
                         real **   cor )
{
  // For looping
  uint i,j;
  

  // Apply the coarse-grid correction to the current fine-grid solution.
  for (i = 0; i < Nx; i ++)
  {
    for (j = 0; j < Ny; j ++)
    {
      pi[i][j] -= cor[i][j];
    }
  }
  
}












/**
 * vcycle_MG
 *
 * Performs a MultiGrid 'V' cycle.
 *
 */
static void vcycle_MG (data_MG * mg_grids, uint n)
{
  // Pointer to finest grid in vector of grids
  data_MG * mg_grid = mg_grids + n;
  
  // Smallest grid; solve exactly
  if (n == 0)
  {
    // Perform initial weighted Jacobi iteration
    exactSolve_MG(  mg_grid->Nx,
                    mg_grid->Ny,
                    mg_grid->pi,
                    mg_grid->pi_rhs,
                    mg_grid->Ow,
                    mg_grid->Os     );
    
    return;
  }

  // Perform initial weighted Jacobi iteration
  relax_MG( mg_grid->Nx,
              mg_grid->Ny,
              mg_grid->pi,
              mg_grid->pi_temp,
              mg_grid->pi_rhs,
              mg_grid->Ow,
              mg_grid->Os,
              mg_grid->_Osum   );

  // Compute residual of weighted Poisson equation and store in pi_temp
  residual_MG(  mg_grid->Nx,
                mg_grid->Ny,
                mg_grid->pi,
                mg_grid->pi_temp,
                mg_grid->pi_rhs,
                mg_grid->Ow,
                mg_grid->Os,
                mg_grid->Osum   );

  // Restrict residual to coarser grid and set it as the right-hand side for the coarser grid solution
  restrict_MG( mg_grid,
               mg_grid->pi_temp,
               &mg_grids[n-1],
               mg_grids[n-1].pi_rhs  );

  // Prior for coarser grid solution is zero
  memset(*(mg_grids[n-1].pi),0,(mg_grids[n-1].Nx)*(mg_grids[n-1].Ny)*sizeof(real));

  // Step down to solve on coarser grid
  vcycle_MG(mg_grids,n-1);

  // Interpolate correction back to this grid and store in pi_temp
  interpolate_MG( mg_grid,
                  mg_grid->pi_temp,
                  &mg_grids[n-1],
                  mg_grids[n-1].pi );

  // Subtract correction
  correct_MG( mg_grid->Nx,
              mg_grid->Ny,
              mg_grid->pi,
              mg_grid->pi_temp );

  // Perform final weighted Jacobi iteration
  relax_MG(   mg_grid->Nx,
              mg_grid->Ny,
              mg_grid->pi,
              mg_grid->pi_temp,
              mg_grid->pi_rhs,
              mg_grid->Ow,
              mg_grid->Os,
              mg_grid->_Osum   );
  
}











/**
 * max_residual
 *
 * Computes the maximum residual (i.e. error) in the pressure.
 *
 * Nx - Grid size in x (first dimension)
 * Ny - Grid size in y (second dimension)
 * pi - Nx x Ny matrix containing current solution.
 
 * pi_rhs - Nx x Ny matrix containing right-hand side of weighted Poisson equation
 * Ow - Nx x Ny matrix containing east-west operators on western edges of grid cells
 * Os - Nx x Ny matrix containing north-south operators on southern edges of grid cells
 * _Osum - Nx x Ny matrix containing reciprocal of operator sum around edges of each grid cell
 *
 */
static real max_residual ( uint      Nx,
                           uint      Ny,
                           real **   pi,
                           real **   pi_rhs,
                           real **   Ow,
                           real **   Os,
                           real **   _Osum )
{
  // For looping
  uint i,j,im1,ip1,jm1,jp1;
  real maxres = 0;
  real res = 0;
  
  
  // Loop through all indices and compute residual between weighted Laplacian operator and right-hand side
  for (i = 0; i < Nx; i ++)
  {
    im1 = (i+Nx-1) % Nx;
    ip1 = (i+Nx+1) % Nx;
    
    for (j = 0; j < Ny; j ++)
    {
      jm1 = (j+Ny-1) % Ny;
      jp1 = (j+Ny+1) % Ny;
      
      // N.B. This code is periodic in y, but the operator Os is set such that the wall BCs are included
      res = pi[i][j] - _Osum[i][j] * ( Os[i][jp1]*pi[i][jp1] + Os[i][j]*pi[i][jm1] + Ow[ip1][j]*pi[ip1][j] + Ow[i][j]*pi[im1][j] - pi_rhs[i][j] );
      maxres = fmax(maxres,fabs(res));
    }
  }
  
  return maxres;
  
}









/**
 * iterate_MG
 *
 * Repeatedly uses MultiGrid V-cycles to solve the weighted Poisson equation on a given grid.
 *
 */
static uint iterate_MG (data_MG * mg_grids, uint n)
{
  // To record iterations
  int iters = 0;
  
  // To calculate convergence
  real maxres = pi_tol+1;
  
  // Pointer to MultiGrid grid that we should iterate
  data_MG * mg_grid = mg_grids + n;

  // Perform MultiGrid iteration
  iters = 0;
  while ((maxres > pi_tol) && (iters < maxiters))
  {
    // Start multigrid V-cycle
    vcycle_MG(mg_grids,n);

    // Calculate current residual
    maxres = max_residual(  mg_grids[n].Nx,
                            mg_grids[n].Ny,
                            mg_grids[n].pi,
                            mg_grids[n].pi_rhs,
                            mg_grids[n].Ow,
                            mg_grids[n].Os,
                            mg_grids[n]._Osum   );
    
    // Keep count of iterations
    iters ++;
  }
  
  // Debug output
  if (debug)
  {
    printf("Grid size: %u x %u\n",mg_grid->Nx,mg_grid->Ny);
    printf("Iterations: %u\n",iters);
    printf("Error: %e\n",maxres);
    fflush(stdout);
  }
  
  return iters;
}













/**
 * full_MG
 *
 * Performs a full MultiGrid cycle, consisting of a series of V-cycles of increasing length.
 *
 * Returns number of iterations required in the deepest V-cycle, 
 * or maxiters if convergence was not achieved.
 *
 */
static uint full_MG (data_MG * mg_grids)
{
  // To record number of iterations
  uint iters = 0;
  
  // For looping
  uint n;

  // First, restrict RHS to all coarser grids
  for (n = Ngrids-1; n > 0; n --)
  {
    restrict_MG( &mg_grids[n],
                 mg_grids[n].pi_rhs,
                 &mg_grids[n-1],
                 mg_grids[n-1].pi_rhs );
  }

  // Solve on the coarsest grid
  vcycle_MG(mg_grids,0);

  // Now perform successive V-cycles, starting from finer and finer grids
  for (n = 1; n < Ngrids; n ++)
  {
    // Interpolate correction back to this grid and store in pi_temp
    interpolate_MG( &mg_grids[n],
                    mg_grids[n].pi,
                    &mg_grids[n-1],
                    mg_grids[n-1].pi );
    
    // Perform the V-cycle
    iters = iterate_MG(mg_grids,n);
  }
  
  return iters;
}



















/**
 * solve_MG
 *
 * Inverts the thickness-weighted Poisson equation iteratively using a multi-grid method.
 *
 * The Nx x Ny matrix pi is used as a prior for the iterative procedure, and is modified to
 * store the updated value of pi when this function returns.
 *
 * The Nx x Ny matrix pi_rhs stores the right-hand side of the Poisson equation.
 *
 * Returns the number of iterations required for convergence, or 'maxiters' if convergence was
 * not achieved within the stipulated number of iterations.
 *
 */
uint solve_MG (real ** pi, real ** pi_rhs)
{
  // To record iterations
  int iters = 0;
  
  // For timing
  clock_t start;
  clock_t end;
  
  // Pointer to finest grid in vector of grids
  data_MG * mg_grid = mg_grids + Ngrids - 1;
  
  // Point finest-grid solution and rhs matrices to pi and pi_rhs
  mg_grid->pi = pi;
  mg_grid->pi_rhs = pi_rhs;
  
  // Initialize timer
  start = clock();

  if (use_fullMG)
  {
    iters = full_MG(mg_grids);
  }
  else
  {
    iters = iterate_MG(mg_grids,Ngrids-1);
  }
  
  // Stop timer
  end = clock();
  
  // Debug output
  if (debug)
  {
    printf("Time: %lu\n",end-start);
    fflush(stdout);
  }
  
  return iters;
}
