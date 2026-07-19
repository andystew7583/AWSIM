/**
 * sor.c
 *
 * Successive over-relaxation pressure solver implementation for AWSIM.
 *
 */
#include <math.h>
#include <time.h>

#include "sor.h"

static uint sor_Nx = 0;
static uint sor_Ny = 0;
static real pi_tol = 0;
static uint maxiters = 0;
static bool debug = false;

static real ** Ow = NULL; // "Operators" participating in weighted Laplacian
static real ** Os = NULL;
static real ** Osum = NULL;
static real ** _Osum = NULL;
static real * im1_vec = NULL; // Vectors to store indices of adjacent grid points,
static real * ip1_vec = NULL; // pre-computed for efficiency
static real * jm1_vec = NULL;
static real * jp1_vec = NULL;

static bool alloc_SOR_matrix (real *** mat, uint Nx, uint Ny)
{
  *mat = matalloc(Nx,Ny);
  if (*mat == NULL)
  {
    fprintf(stderr,"ERROR: Unable to allocate memory\r\n");
    return false;
  }

  return true;
}

static bool alloc_SOR_vector (real ** vec, uint N)
{
  *vec = vecalloc(N);
  if (*vec == NULL)
  {
    fprintf(stderr,"ERROR: Unable to allocate memory\r\n");
    return false;
  }

  return true;
}

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
                bool    debugOutput )
{
  uint i,j,ip1,jp1;
  real dxsq = SQUARE(dx);
  real dysq = SQUARE(dy);

  sor_Nx = Nx;
  sor_Ny = Ny;
  pi_tol = tol;
  maxiters = maxIters;
  debug = debugOutput;

  if ( !alloc_SOR_matrix(&Ow,Nx,Ny)
    || !alloc_SOR_matrix(&Os,Nx,Ny)
    || !alloc_SOR_matrix(&Osum,Nx,Ny)
    || !alloc_SOR_matrix(&_Osum,Nx,Ny)
    || !alloc_SOR_vector(&im1_vec,Nx)
    || !alloc_SOR_vector(&ip1_vec,Nx)
    || !alloc_SOR_vector(&jm1_vec,Ny)
    || !alloc_SOR_vector(&jp1_vec,Ny) )
  {
    return false;
  }

  // Define north-south operators for the pressure solve
  for (i = 0; i < Nx; i ++)
  {
    for (j = 0; j < Ny; j ++)
    {
      Os[i][j] = Hs[i][j] / dysq;
    }

    // No north-south periodic operator if there is a wall
    if (useWallNS)
    {
      Os[i][0] = 0;
    }
  }

  // Define east-west operators for the pressure solve
  for (j = 0; j < Ny; j ++)
  {
    for (i = 0; i < Nx; i ++)
    {
      Ow[i][j] = Hw[i][j] / dxsq;
    }

    // No east-west periodic operator if there is a wall
    if (useWallEW)
    {
      Ow[0][j] = 0;
    }
  }

  // Define sum of operators around each cell for computational efficiency
  for (i = 0; i < Nx; i ++)
  {
    ip1 = (i + Nx + 1) % Nx;

    for (j = 0; j < Ny; j ++)
    {
      jp1 = (j + Ny + 1) % Ny;

      Osum[i][j] = Ow[i][j] + Ow[ip1][j] + Os[i][j] + Os[i][jp1];
      _Osum[i][j] = 1 / Osum[i][j];
    }
  }

  // Indexing for adjacent gridpoints - to be used in pressure solve for efficiency
  for (i = 0; i < Nx; i ++)
  {
    im1_vec[i] = (i+Nx-1) % Nx;
    ip1_vec[i] = (i+Nx+1) % Nx;
  }
  for (j = 0; j < Ny; j ++)
  {
    jm1_vec[j] = (j+Ny-1) % Ny;
    jp1_vec[j] = (j+Ny+1) % Ny;
  }

  return true;
}

/**
 * solve_SOR
 *
 * Inverts the thickness-weighted Poisson equation iteratively using successive over-relaxation.
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
uint solve_SOR (real ** pi, real ** pi_rhs, real rp)
{
  // To record iterations
  int iters = 0;

  // To calculate convergence
  real maxdiff;
  real diff = 0;
  real pi_prev = 0;

  // For timing
  clock_t start;
  clock_t end;

  // Looping variables
  int i,j,im1,ip1,jp1,jm1;

  // Initialize timer
  start = clock();

  // Perform SOR iteration
  maxdiff = pi_tol + 1;
  iters = 0;
  while ((maxdiff > pi_tol) && (iters < maxiters))
  {
    maxdiff = 0;


    for (i = 0; i < sor_Nx; i ++)
    {
      im1 = im1_vec[i];
      ip1 = ip1_vec[i];

      for (j = 0; j < sor_Ny; j ++)
      {
        jm1 = jm1_vec[j];
        jp1 = jp1_vec[j];

        // Store current grid value of pi
        pi_prev = pi[i][j];

        // N.B. This code is periodic in y, but the operator Os is set such that the wall BCs are included
        pi[i][j] = (1-rp)*pi[i][j]
                 + rp * _Osum[i][j]
                      *  ( Os[i][jp1]*pi[i][jp1] + Os[i][j]*pi[i][jm1] + Ow[ip1][j]*pi[ip1][j] + Ow[i][j]*pi[im1][j] - pi_rhs[i][j] );

        // Calculate the absolute difference between iterations
        diff = fabs(pi[i][j]-pi_prev);
        maxdiff = fmax(diff,maxdiff);
      }
    }

    iters ++;
  }

  // Stop timer
  end = clock();
  if (debug)
  {
    printf("Time: %lu\n",end-start);
    printf("Iterations: %u\n",iters);
    printf("Error: %e\n",maxdiff);
    fflush(stdout);
  }

  return iters;
}
