/**
  * AWSIM.c
  *
  * Andrew Stewart
  *
  * Main execution file for the AndreW Stewart Isopycnal Model (AWSIM).
  *
  */
#include <time.h>
#include <math.h>

#include "rk.h"
#include "ab.h"
#include "defs.h"

#ifdef ALLOW_FFTW
#include <fftw3.h>
#endif

// Must match number of input parameters defined via "setParam" below
#define NPARAMS 91

// Avoids memory errors associated with usual definition of bool
typedef int mybool;

// Uncomment to print debug messages
const mybool debug = true;

// Work arrays for the 'tderiv' function
real * vars_w = NULL;     // Current iteration data
real *** uu_w = NULL;     // Matrices that store current iteration data,
real *** vv_w = NULL;     // including space for ghost points
real *** hh_w = NULL;
real *** bb_w = NULL;
real *** eta_w = NULL;
real *** uu_d = NULL;     // Matrices that store current iteration data,
real *** vv_d = NULL;     // with no ghost points
real *** hh_d = NULL;
real *** bb_d = NULL;
real *** dt_uu_w = NULL;   // Time derivatives of current iteration
real *** dt_vv_w = NULL;
real *** dt_hh_w = NULL;
real *** dt_bb_w = NULL;
real ** hh_q = NULL;      // Height values on q-gridpoints
real ** zeta = NULL;      // Relative vorticity
real ** qq = NULL;        // Potential vorticity
real *** h_west = NULL;    // Layer thickness on western cell face
real *** h_south = NULL;   // Layer thickness on southern cell face
real ** d2hx = NULL; // Second-order derivatives of h
real ** d2hy = NULL;
real ** huu = NULL;       // Mass-weighted velocities
real ** hvv = NULL;
real ** alpha_w = NULL;   // Arakawa parameters
real ** beta_w = NULL;
real ** gamma_w = NULL;
real ** delta_w = NULL;
real ** epsilon_w = NULL;
real ** phi_w = NULL;
real ** MM_B = NULL;       // Montgomery potential contribution to Bernoulli Potential
real ** KE_B = NULL;       // Kinetic energy contribution to Bernoulli potential
real ** pp = NULL;         // Pressure in each model layer
real ** A2_q = NULL;       // Laplacian viscosity on cell corners
real ** A2_h = NULL;       // Laplacian viscosity on cell centers
real ** A4sqrt_q = NULL;   // Biharmonic viscosity on cell corners
real ** A4sqrt_h = NULL;   // Biharmonic viscosity on cell centers
real ** DD_T = NULL;       // Horizontal tension
real ** DD_S = NULL;       // Horizontal shearing strain
real ** DD4_T = NULL;      // Biharmonic horizontal tension
real ** DD4_S = NULL;      // Biharmonic shearing strain
real ** UU4 = NULL;        // Intermediate "velocity" fields for biharmonic viscosity operator
real ** VV4 = NULL;
real ** Omega_z_w = NULL; // Vertical component of rotation vector
real ** hhs_w = NULL;     // Surface elevation
real ** hhb_w = NULL;     // Bottom elevation
real ** dh_dx = NULL;     // Derivatives of layer thickness
real ** dh_dy = NULL;

// Required for surface/bottom momentum forcing
real ** usq_surf = NULL;
real ** vsq_surf = NULL;
real ** usq_bot = NULL;
real ** vsq_bot = NULL;
real ** uabs_surf = NULL;
real ** uabs_bot = NULL;
real *** hFsurf_west = NULL;
real *** hFsurf_south= NULL;
real *** hFbot_west = NULL;
real *** hFbot_south= NULL;
real *** eta_west = NULL;
real *** eta_south = NULL;
real ** hhs_west = NULL;
real ** hhs_south = NULL;
real ** hhb_west = NULL;
real ** hhb_south = NULL;

// Required for tracer evolution
real ** udb_dx = NULL;
real ** vdb_dy = NULL;
real ** db_dx = NULL;
real ** db_dy = NULL;
real ** d2bx = NULL;
real ** d2by = NULL;
real ** hub = NULL;
real ** hvb = NULL;
real ** hz = NULL;
real ** BB4 = NULL;     // Intermediate "tracer" field for biharmonic diffusion operator

// Used by surfPressure
real ** pi_rhs = NULL;
real ** Hc = NULL; // Water column thickness on cell centers and edges
real ** Hw = NULL;
real ** Hs = NULL;
real ** Ow = NULL; // "Operators" participating in weighted Laplacian
real ** Os = NULL;
real ** Osum = NULL;
real ** _Osum = NULL;
real * im1_vec = NULL; // Vectors to store indices of adjacent grid points,
real * ip1_vec = NULL; // pre-computed for efficiency
real * jm1_vec = NULL;
real * jp1_vec = NULL;

// For calculation of E and Z in calcEZ
real ** Zdens = NULL;
real *** KEdens = NULL;
real *** PEdens = NULL;

// Relaxation matrices
real *** uRelax = NULL;
real *** vRelax = NULL;
real *** hRelax = NULL;
real *** eRelax = NULL;
real *** bRelax = NULL;
real *** uTime = NULL;
real *** vTime = NULL;
real ** hTime = NULL;
real *** eTime = NULL;
real *** bTime = NULL;

// Diabatic velocity
real wDiaPeriod = 0;                     // Period over which diabatic velocity is applied
uint wDiaNrecs = 1;                      // Number of temporal records supplied in diabatic velocity files
real *** wdia = NULL;
real *** wdia_u = NULL;
real *** wdia_v = NULL;
real **** wdia_ff = NULL;

// Grid size - number of grid boxes in the domain
uint Nx = 0;       // Number of interior gridpoints in the x-direction
uint Ny = 0;       // Number of interior gridpoints in the y-direction
uint Nlay = 0;     // Number of isopycnal layers
uint N = 0;        // Number of gridpoints for each interior variable (u,v,h)
uint Ntotal = 0;   // Total number of gridpoints passed to 'tderiv'
uint N_g = 0;   // As above, but including ghost points
uint Ntotal_g = 0;
const uint Ng = 3; // Number of ghost points at y-edges of domain
real dt = 0;            // Time step
real dx = 0;       // Grid spacing
real dy = 0;       // Grid spacing
real dxsq = 0;     // Grid spacing squared
real dysq = 0;     // Grid spacing squared

// Input parameters
mybool useRL = true;                    // Whether to use the rigid lid approximation
mybool useTrad = true;                  // Whether to use the traditional approximation
mybool useWallEW = false;               // Whether to have a wall the eastern/western boundary
mybool useWallNS = false;               // Whether to have a wall at the northern/southern boundary
mybool useTracer = false;               // Whether to include a tracer variable in each layer
mybool useBuoyancy = false;             // Whether to use the tracer as a buoyancy variable
real * gg = NULL;                       // Gravity/reduced gravity at each layer interface
real * geff = NULL;                     // Effective gravity in each layer (basically g times density of that layer)
real ** hhs = NULL;                     // Surface topography elevation
real ** hhb = NULL;                     // Bottom topography elevation
real ** Omega_x = NULL;                 // Horizontal rotation component in the x-direction
real ** Omega_y = NULL;                 // Horizontal rotation component in the y-direction
real ** Omega_z = NULL;                 // Vertical rotation component
real h0 = 0.0;                          // Minimum fluid height (Salmon layer thickness)
real hsml = 0;                          // Surface "mixed layer" thickness
real hbbl = 0;                          // Bottom "boundary layer" thickness
real hmin_surf = 0;                     // Minimum layer thickness for surface momentum forcing
real hmin_bot = 0;                      // Minimum layer thickness for bottom momentum forcing
real A2const = 0.0;                     // Constant Laplacian viscosity coefficient
real A4const = 0.0;                     // Constant biharmonic viscosity coefficient
real A2smag = 0.0;                      // Smagorinsky Laplacian viscosity coefficient
real A4smag = 0.0;                      // Smagorinsky biharmonic viscosity coefficient
real useA2 = false;                     // Laplacian viscosity flag
real useA4 = false;                     // Biharmonic viscosity flag
real K2 = 0.0;                          // Laplacian tracer diffusion coefficient
real K4 = 0.0;                          // Biharmonic tracer diffusion coefficient
real *** taux = NULL;                   // x-component of wind stress
real *** tauy = NULL;                   // y-component of wind stress
real tauPeriod = 0;                     // Period over which wind forcing is applied
uint tauNrecs = 1;                      // Number of temporal records supplied in wind forcing files
mybool useWind = false;                 // Will be set true if wind stress is input
mybool useRelax = false;                // Will be set true if a relaxation time scale file is specified
mybool useWDia = false;                 // Will be set true if diabatic velocity is active
real rDrag = 0;                         // Linear bottom drag coefficient
real CdBot = 0;                         // Quadratic bottom drag coefficient
real rSurf = 0;                         // Linear surface drag coefficient
real CdSurf = 0;                        // Quadratic surface drag coefficient
real ** uLid = NULL;                    // u-velocity of rigid lid - used in the calculation of surface drag
real ** vLid = NULL;                    // v-velocity of rigid lid - used in the calculation of surface drag
real ** uLid_g = NULL;                  // u-velocity of rigid lid including ghost points
real ** vLid_g = NULL;                  // v-velocity of rigid lid including ghost points
real ** Fbaro_x = NULL;                 // x-component of imposed barotropic forcing
real ** Fbaro_y = NULL;                 // y-component of imposed barotropic forcing
mybool useFbaro = false;                // Will be set true if barotropic forcing is input

// Pressure solve parameters
mybool use_MG = false;                  // Set true to use MultiGrid rather than SOR
mybool use_fullMG = false;              // Set true to use full MultiGrid algorithm
real pi_tol = 1e-5;                     // SOR iteration error tolerance
uint maxiters = 10000;                  // Maximum number of iterations

// SOR-specific parameters
real rp_acc_max = 0.001;                // Accuracy to use in optimizing relaxation parameter
real rp_opt_max = 2.0;                  // Max value of rp over which to optimize
real rp_opt_min = 1.0;                  // Min value of rp over which to optimize
uint N_rp = 10;                         // Determines discretization used in optimizing rp
uint rp_opt_freq = 1000;                // Number of time steps between optimizations of rp
uint n_first_optim = 10;                // Number of iterations to wait before first optimization
real * pi_prev = NULL;                  // Pressure in previous SOR iteration

// MultiGrid-specific parameters
uint Npx = 0;
uint Npy = 0;
uint Ngrids = 0;
real *** pi_MG = NULL;
real omega_WJ = 2.0/3.0;
real * F_MG = NULL;
uint F_len = 0;

// Thickness advection parameters
real KT00_sigma = 1.4;

// Random forcing parameters
real RF_F0_rot = 0;                     // RMS amplitude of rotational component of random forcing vector
real RF_F0_div = 0;                     // RMS amplitude of divergent component of random forcing vector
real RF_tau = 0;                        // Random forcing autocorrelation time scale
mybool useRandomForcing = false;        // Whether random forcing is to be used
mybool useThicWeightedRF = false;       // Whether random forcing is applied to the momentum equation
                                        // (i.e. d(hu)/dt = F + ...) or the velocity equation (i.e. du/dt = F ...)
mybool useBarotropicRF = true;          // Whether the same random forcing field should be used for all \
                                        // isopycnal layers, or evolved separately for each layer
real *** RFu = NULL;                    // Random forcing tendencies in real space
real *** RFv = NULL;
real *** RFmask_fft = NULL;             // Mask for amplitudes of random forcing components in spectral space
real *** RFmask_rot = NULL;             // Masks for rotational and divergent components of random forcing in real space
real *** RFmask_div = NULL;
#ifdef ALLOW_FFTW
fftw_plan fft_backward = NULL;          // FFTW plans for inverse 2D FFT
fftw_complex * RFrot = NULL;          // Rotational and divergent components of random forcing tendencies in real space
fftw_complex * RFdiv = NULL;
fftw_complex * RFrot_fft = NULL;      // Rotational and divergent components of random forcing tendencies in spectral space
fftw_complex * RFdiv_fft = NULL;
#endif

// Default parameter values
const real hhs0 = 0;
const real hhb0 = -1000;
const real g0 = 9.81;
const real Omega_x0 = 0;
const real Omega_y0 = 0;
const real Omega_z0 = 0;

// Parameters for averaging u-momentum equation
real avg_fac_hu = 1;       // Fraction of time step used for averaging
uint n_avg_hu = 0;         // Keeps track of number of averages computed
real dt_avg_hu = 0;        // Save time step for averaged output
uint n_prev_avg_hu = 0;    // Model time step number corresponding to previous average output
real t_next_avg_hu = 0;    // Next save time for averaged output
real avg_len_hu = 0;       // True length of time average
real *** hu_tend_q = NULL;
real *** hu_tend_gradM = NULL;
real *** hu_tend_gradKE = NULL;
real *** hu_tend_dhdt = NULL;
real *** hu_tend_A2 = NULL;
real *** hu_tend_A4 = NULL;
real *** hu_tend_rDrag = NULL;
real *** hu_tend_rSurf = NULL;
real *** hu_tend_CdBot = NULL;
real *** hu_tend_CdSurf = NULL;
real *** hu_tend_wind = NULL;
real *** hu_tend_buoy = NULL;
real *** hu_tend_relax = NULL;
real *** hu_tend_wdia = NULL;
real *** hu_tend_rand = NULL;
real *** hu_tend_Fbaro = NULL;

// Parameters for averaging v-momentum equation
real avg_fac_hv = 1;       // Fraction of time step used for averaging
uint n_avg_hv = 0;         // Keeps track of number of averages computed
real dt_avg_hv = 0;        // Save time step for averaged output
uint n_prev_avg_hv = 0;    // Model time step number corresponding to previous average output
real t_next_avg_hv = 0;    // Next save time for averaged output
real avg_len_hv = 0;       // True length of time average
real *** hv_tend_q = NULL;
real *** hv_tend_gradM = NULL;
real *** hv_tend_gradKE = NULL;
real *** hv_tend_dhdt = NULL;
real *** hv_tend_A2 = NULL;
real *** hv_tend_A4 = NULL;
real *** hv_tend_rDrag = NULL;
real *** hv_tend_rSurf = NULL;
real *** hv_tend_CdBot = NULL;
real *** hv_tend_CdSurf = NULL;
real *** hv_tend_wind = NULL;
real *** hv_tend_buoy = NULL;
real *** hv_tend_relax = NULL;
real *** hv_tend_wdia = NULL;
real *** hv_tend_rand = NULL;
real *** hv_tend_Fbaro = NULL;

// Parameters for averaging h-equation
real avg_fac_h = 1;        // Fraction of time step used for averaging
uint n_avg_h = 0;         // Keeps track of number of averages computed
real dt_avg_h = 0;        // Save time step for averaged output
uint n_prev_avg_h = 0;    // Model time step number corresponding to previous average output
real t_next_avg_h = 0;    // Next save time for averaged output
real avg_len_h = 0;       // True length of time average
real *** h_tend_adv = NULL;
real *** h_tend_relax = NULL;

// Parameters for averaging energy diagnostics
real avg_fac_e = 1;        // Fraction of time step used for averaging
uint n_avg_e = 0;         // Keeps track of number of averages computed
real dt_avg_e = 0;        // Save time step for averaged output
uint n_prev_avg_e = 0;    // Model time step number corresponding to previous average output
real t_next_avg_e = 0;    // Next save time for averaged output
real avg_len_e = 0;       // True length of time average
real *** e_flux_uP = NULL;
real *** e_flux_uKE = NULL;
real *** e_flux_uPE = NULL;
real *** e_flux_vP = NULL;
real *** e_flux_vKE = NULL;
real *** e_flux_vPE = NULL;
real *** e_tend_adv = NULL;
real *** e_tend_gradM = NULL;
real *** e_tend_wind = NULL;
real *** e_tend_rDrag = NULL;
real *** e_tend_rSurf = NULL;
real *** e_tend_CdBot = NULL;
real *** e_tend_CdSurf = NULL;
real *** e_tend_A2 = NULL;
real *** e_tend_A4 = NULL;
real *** e_tend_wdiaPE = NULL;
real *** e_tend_wdiaKE = NULL;
real *** e_tend_Fbaro = NULL;
real *** e_tend_buoy = NULL;
real *** e_tend_rand = NULL;
real *** e_tend_relax = NULL;

// Time-stepping method to use
uint timeSteppingScheme = TIMESTEPPING_AB3;

// Momentum discretization method to use
uint momentumScheme = MOMENTUM_AL81;

// Thickness discretization method to use
uint thicknessScheme = THICKNESS_AL81;

// Tracer discretization method to use
uint tracerScheme = TRACER_AL81;

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
  
  real * im1_vec; // Vectors to store indices of adjacent grid points,
  real * ip1_vec; // pre-computed for efficiency
  real * jm1_vec;
  real * jp1_vec;
}
data_MG;

// To store MultiGrid solver grids
data_MG * mg_grids = NULL;





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
void exactSolve_MG (  uint      Nx,
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
void relax_MG (   uint      Nx,
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
  
#pragma parallel
  
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
void residual_MG (  uint      Nx,
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
  
#pragma parallel
  
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
 * restrict_MG
 *
 * Restriction operator for MultiGrid solver. Takes a matrix and restricts it to a coarser
 * grid.
 *
 * Nx - Fine grid size in x (first dimension)
 * Ny - Fine grid size in y (second dimension)
 * pi_f - Nx x Ny fine-resolution matrix containing input data
 * pi_c - Nx/2 x Ny/2 coarse-resolution matrix to store output data
 *
 * Note: Currently only designed for coarsening to a grid with exactly half as many grid cells
 * in each direction.
 *
 */
void restrict_MG (  uint      Nx,
                    uint      Ny,
                    real **   pi_f,
                    real **   pi_c )
{
  // For looping
  uint i,j;
  
#pragma parallel
  
  // Restriction operator is simple average of the four fine-grid elements that lie
  // within each coarse-grid cell
  for (i = 0; i < Nx/2; i ++)
  {
    for (j = 0; j < Ny/2; j ++)
    {
      pi_c[i][j] = 0.25 * (pi_f[2*i][2*j] + pi_f[2*i+1][2*j] + pi_f[2*i][2*j+1] + pi_f[2*i+1][2*j+1]);
    }
  }
}













/**
 * interpolate_MG
 *
 * Interpolation operator for MultiGrid solver. Takes a matrix and restricts it to a coarser
 * grid.
 *
 * Nx - Fine grid size in x (first dimension)
 * Ny - Fine grid size in y (second dimension)
 * pi_f - Nx x Ny fine-resolution matrix to store output data
 * pi_c - Nx/2 x Ny/2 coarse-resolution matrix containing input data
 *
 * wmm,wmp,wpm,wpp - Nx x Ny interpolation weight matrices. Supply weights to be used in summing
 *                   coarse-grid elements to the southwest (wmm), northwest (wmp), northeast (wpp)
 *                   and southeast (wpm) of each fine-grid element.
 *
 * Note: Currently only designed for interpolating to a grid with exactly twice as many grid cells
 * in each direction.
 *
 */
void interpolate_MG ( uint      Nx,
                      uint      Ny,
                      real **   pi_f,
                      real **   pi_c,
                      real **   wmm,
                      real **   wmp,
                      real **   wpm,
                      real **   wpp   )
{
  // For looping
  uint i,j,im,ip,jm,jp;

#pragma parallel
  
  for (i = 0; i < Nx; i ++)
  {
    // Identify coarse-grid i-points to the east and west of this fine-grid i-point
    ip = (((i+1)-(i+1)%2)/2 + Nx/2) % (Nx/2);
    im = (ip-1 + Nx/2) % (Nx/2);
    
    for (j = 0; j < Ny; j ++)
    {
      // Identify coarse-grid j-points to the north and south of this fine-grid j-point
      jp = (((j+1)-(j+1)%2)/2 + Ny/2) % (Ny/2);
      jm = (jp-1 + Ny/2) % (Ny/2);
      
      // Interpolate using pre-calculated weights
      pi_f[i][j] = wmm[i][j]*pi_c[im][jm] + wpm[i][j]*pi_c[ip][jm] + wmp[i][j]*pi_c[im][jp] + wpp[i][j]*pi_c[ip][jp];
    }
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
 * Note: Currently only designed for coarsening to a grid with exactly half as many grid cells
 * in each direction.
 *
 */
void correct_MG ( uint      Nx,
                  uint      Ny,
                  real **   pi,
                  real **   cor )
{
  // For looping
  uint i,j;
  
#pragma parallel

  // Restriction operator is simple average of the four fine-grid elements that lie
  // within each coarse-grid cell
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
void vcycle_MG (data_MG * mg_grids, uint n)
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
  restrict_MG(  mg_grid->Nx,
                mg_grid->Ny,
                mg_grid->pi_temp,
                mg_grids[n-1].pi_rhs  );

  // Prior for coarser grid solution is zero
  memset(*(mg_grids[n-1].pi),0,(mg_grids[n-1].Nx)*(mg_grids[n-1].Ny)*sizeof(real));

  // Step down to solve on coarser grid
  vcycle_MG(mg_grids,n-1);

  // Interpolate correction back to this grid and store in pi_temp
  interpolate_MG( mg_grid->Nx,
                  mg_grid->Ny,
                  mg_grid->pi_temp,
                  mg_grids[n-1].pi,
                  mg_grid->wmm,
                  mg_grid->wmp,
                  mg_grid->wpm,
                  mg_grid->wpp   );

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
real max_residual ( uint      Nx,
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
  
#pragma parallel
  
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
uint iterate_MG (data_MG * mg_grids, uint n)
{
  // To record iterations
  int iters = 0;
  
  // To calculate convergence
  real maxres = pi_tol+1;
  real diff = 0;
  
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
uint full_MG (data_MG * mg_grids)
{
  // To record number of iterations
  uint iters;
  
  // For looping
  uint n;

  // First, restrict RHS to all coarser grids
  for (n = Ngrids-1; n > 0; n --)
  {
    restrict_MG(  mg_grids[n].Nx,
                  mg_grids[n].Ny,
                  mg_grids[n].pi_rhs,
                  mg_grids[n-1].pi_rhs  );
  }

  // Solve on the coarsest grid
  vcycle_MG(mg_grids,0);

  // Now perform successive V-cycles, starting from finer and finer grids
  for (n = 1; n < Ngrids; n ++)
  {
    // Interpolate correction back to this grid and store in pi_temp
    interpolate_MG( mg_grids[n].Nx,
                    mg_grids[n].Ny,
                    mg_grids[n].pi,
                    mg_grids[n-1].pi,
                    mg_grids[n].wmm,
                    mg_grids[n].wmp,
                    mg_grids[n].wpm,
                    mg_grids[n].wpp   );
    
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
  real temp;
  real diff = 0;
  
  // For timing
  clock_t start;
  clock_t end;
  
  // Looping variables
  int i,j,k,im1,ip1,jp1,jm1;
  
  // Initialize timer
  start = clock();
  
  // Perform SOR iteration
  maxdiff = pi_tol + 1;
  iters = 0;
  while ((maxdiff > pi_tol) && (iters < maxiters))
  {
    maxdiff = 0;
    
#pragma parallel
    
    for (i = 0; i < Nx; i ++)
    {
      im1 = im1_vec[i];
      ip1 = ip1_vec[i];
      
      for (j = 0; j < Ny; j ++)
      {
        jm1 = jm1_vec[j];
        jp1 = jp1_vec[j];
        
        // Store current grid value of pi
        *pi_prev = pi[i][j];
        
        // N.B. This code is periodic in y, but the operator Os is set such that the wall BCs are included
        pi[i][j] = (1-rp)*pi[i][j]
                 + rp * _Osum[i][j]
                      *  ( Os[i][jp1]*pi[i][jp1] + Os[i][j]*pi[i][jm1] + Ow[ip1][j]*pi[ip1][j] + Ow[i][j]*pi[im1][j] - pi_rhs[i][j] );
        
        // Calculate the absolute difference between iterations
        diff = fabs(pi[i][j]-*pi_prev);
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













/**
 *
 * calcFaceThickness
 *
 * Calculates layer thickness on western and southern cell faces. The calculation
 * depends on which definition of the interpolated layer thickness (or equivalently
 * the kinetic energy in the Hamiltonian; Salmon (2004)) is used. The flag use_ghost
 * tells this function whether to treat hh as an Nx x Ny matrix containing only points
 * that lie within the computational domain, or to treat hh as an Nx+2*Ng x Ny+2*Ng
 * matrix with periodic/wall boundary conditions implicitly accounted for via
 * suitably-prescribed ghost points.
 *
 * Velocities may be specified on cell faces. These will be used to determine the layer 
 * thicknesses if third-order upwinding is being used. If uu or vv are NULL pointers 
 * then the layer thicknesses will be determined via centered differencing.
 *
 * Note that if there are walls on the boundaries then the thickness will still be
 * calculated on those boundaries by "wrapping" the cell center thickness periodically
 * across the boundaries. This is fine because h on walls doesn't influence the evolution
 * of the model.
 *
 */
void calcFaceThickness (real ** hh, real ** h_west, real ** h_south, mybool use_ghost, uint Nx, uint Ny, real ** uu, real ** vv)
{
  int i,j,im1,ip1,jm1,jp1,imin,imax,jmin,jmax,n;
  real hh_xm, hh_xp, hh_ym, hh_yp; // Intermediate h-estimates for the Kurganov-Tadmor scheme
  
  // Min/max i/j indices depend on whether arrays include ghost points
  imin = use_ghost ? 1 : 0;
  jmin = use_ghost ? 1 : 0;
  imax = use_ghost ? Nx+2*Ng-1 : Nx;
  jmax = use_ghost ? Ny+2*Ng-1 : Ny;
  
  // Flags whether velocities have been set, meaning that third-order upwinding should
  // be used. If false then thicknesses will be computed via centered averages, even
  // if the discretization method specifies third-order upwinding.
  mybool vel_flag = (uu != NULL) && (vv != NULL);
  

  // Arakawa and Lamb (1981) centered scheme. Also the default if no velocities are provided as input parameters.
  if ((thicknessScheme == THICKNESS_AL81) || (!vel_flag))
  {
    
#pragma parallel
    
    // Compute thickness on cell faces
    for (i = imin; i < imax; i ++)
    {
      im1 = use_ghost ? i-1 : (i+Nx-1) % Nx;
      
      for (j = jmin; j < jmax; j ++)
      {
        jm1 = use_ghost ? j-1 : (j+Ny-1) % Ny;
        
        h_west[i][j] = 0.5 * (hh[i][j]+hh[im1][j]);
        h_south[i][j] = 0.5 * (hh[i][j]+hh[i][jm1]);
      }
    }
    
  } // end if (thicknessScheme == THICKNESS_AL81 || (!vel_flag))
  
  
  // Hollingsworth and Kallberg (1983) extended central differencing formulation
  if (thicknessScheme == THICKNESS_HK83)
  {
    
#pragma parallel
    
    // Compute thickness on cell faces
    for (i = imin; i < imax; i ++)
    {
      im1 = use_ghost ? i-1 : (i+Nx-1) % Nx;
      ip1 = use_ghost ? i+1 : (i+Nx+1) % Nx;
      
      for (j = jmin; j < jmax; j ++)
      {
        jm1 = use_ghost ? j-1 : (j+Ny-1) % Ny;
        jp1 = use_ghost ? j+1 : (j+Ny+1) % Ny;
        

        // i end-points are only calculated if ghost points aren't used. If ghost points are used, the near-boundary
        // points will be calculated correctly via the interior formula combined with the ghost points.
        if (useWallEW && (i==0))
        {
          h_south[i][j] = (1.0/3) * ( 1.25 * (hh[i][j]+hh[i][jm1]) + 0.25 * (hh[ip1][j]+hh[ip1][jm1]) );
        }
        else if (useWallEW && (i==Nx-1))
        {
          h_south[i][j] = (1.0/3) * ( 1.25 * (hh[i][j]+hh[i][jm1]) + 0.25 * (hh[im1][j]+hh[im1][jm1]) );
        }
        else
        {
          h_south[i][j] = (1.0/3) * ( hh[i][j]+hh[i][jm1] + 0.25 * (hh[ip1][j]+hh[im1][j]+hh[ip1][jm1]+hh[im1][jm1]) );
        }
     
        // j end-points are only calculated if ghost points aren't used. If ghost points are used, the near-boundary
        // points will be calculated correctly via the interior formula combined with the ghost points.
        if (useWallNS && (j == 0))
        {
          h_west[i][j] = use_ghost ? 0 : (1.0/3) * ( 1.25*(hh[i][j]+hh[im1][j]) + 0.25*(hh[i][jp1]+hh[im1][jp1]) );
        }
        else if (useWallNS && (j == Ny-1))
        {
          h_west[i][j] = use_ghost ? 0 : (1.0/3) * ( 1.25*(hh[i][j]+hh[im1][j]) + 0.25*(hh[i][jm1]+hh[im1][jm1]) );
        }
        else
        {
          h_west[i][j] = (1.0/3) * ( hh[i][j]+hh[im1][j] + 0.25*(hh[i][jp1]+hh[im1][jp1]+hh[i][jm1]+hh[im1][jm1]) );
        }
      }
    }
    
  } // end if (thicknessScheme == THICKNESS_HK83)
  
  
  // Third-order upwinding
  if (vel_flag && (thicknessScheme == THICKNESS_UP3))
  {
    
#pragma parallel
    
    // Compute second derivative for third-order upwinding
    for (i = imin; i < imax; i ++)
    {
      im1 = use_ghost ? i-1 : (i+Nx-1) % Nx;
      ip1 = use_ghost ? i+1 : (i+Nx+1) % Nx;
      
      for (j = jmin; j < jmax; j ++)
      {
        jm1 = use_ghost ? j-1 : (j+Ny-1) % Ny;
        jp1 = use_ghost ? j+1 : (j+Ny+1) % Ny;
        
        d2hx[i][j] = (hh[ip1][j]-2*hh[i][j]+hh[im1][j]);
        d2hy[i][j] = (hh[i][jp1]-2*hh[i][j]+hh[i][jm1]);
      }
    }
    
#pragma parallel
    
    // Compute thickness on cell faces
    for (i = imin; i < imax; i ++)
    {
      im1 = use_ghost ? i-1 : (i+Nx-1) % Nx;
      ip1 = use_ghost ? i+1 : (i+Nx+1) % Nx;
      
      for (j = jmin; j < jmax; j ++)
      {
        jm1 = use_ghost ? j-1 : (j+Ny-1) % Ny;
        jp1 = use_ghost ? j+1 : (j+Ny+1) % Ny;
        

        h_south[i][j] = 0.5 * (hh[i][j]+hh[i][jm1]);
        h_south[i][j] -= (1/6.0) * ( (vv[i][j] > 0) ? d2hy[i][jm1] : d2hy[i][j] );
        h_west[i][j] = 0.5 * (hh[i][j]+hh[im1][j]);
        h_west[i][j] -= (1/6.0) * ( (uu[i][j] > 0) ? d2hx[im1][j] : d2hx[i][j] );
      }
    }
    
  } // end if (thicknessScheme == THICKNESS_UP3)
  
  
  if (vel_flag && (thicknessScheme == THICKNESS_KT00))
  {
    
#pragma parallel
  
    // Determine limited slopes at cell centres
    for (i = imin; i < imax; i ++)
    {
      im1 = use_ghost ? i-1 : (i+Nx-1) % Nx;
      ip1 = use_ghost ? i+1 : (i+Nx+1) % Nx;
      
      for (j = jmin; j < jmax; j ++)
      {
        jm1 = use_ghost ? j-1 : (j+Ny-1) % Ny;
        jp1 = use_ghost ? j+1 : (j+Ny+1) % Ny;
        
        dh_dx[i][j] = minmod( KT00_sigma * (hh[ip1][j]-hh[i][j]),
                               0.5 * (hh[ip1][j]-hh[im1][j]),
                               KT00_sigma * (hh[i][j]-hh[im1][j]) );
        dh_dy[i][j] = minmod( KT00_sigma * (hh[i][jp1]-hh[i][j]),
                               0.5 * (hh[i][jp1]-hh[i][jm1]),
                               KT00_sigma * (hh[i][j]-hh[i][jm1]) );
      }
    }
    
#pragma parallel
    
    // Interpolate hh to cell faces
    for (i = imin; i < imax; i ++)
    {
      im1 = use_ghost ? i-1 : (i+Nx-1) % Nx;
      ip1 = use_ghost ? i+1 : (i+Nx+1) % Nx;
      
      for (j = jmin; j < jmax; j ++)
      {
        jm1 = use_ghost ? j-1 : (j+Ny-1) % Ny;
        jp1 = use_ghost ? j+1 : (j+Ny+1) % Ny;

        // Interpolate h to cell west faces
        hh_xm = hh[im1][j] + 0.5*dh_dx[im1][j];
        hh_xp = hh[i][j] - 0.5*dh_dx[i][j];
        h_west[i][j] = (uu[i][j] > 0) ? hh_xm : hh_xp;

        // Interpolate h to cell south faces
        hh_ym = hh[i][jm1] + 0.5*dh_dy[i][jm1];
        hh_yp = hh[i][j] + 0.5*dh_dy[i][j];
        h_south[i][j] = (vv[i][j] > 0) ? hh_ym : hh_yp;
      }
    }

  } // end if (thicknessScheme == THICKNESS_KT00)

  
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
  int i,j,k,im1,imin,ip1,jp1,jm1,jmin;
  
  // Set right-hand side of Poisson equation
  memset(*pi_rhs,0,Nx*Ny*sizeof(real));
  
  for (k = 0; k < Nlay; k ++)
  {

#pragma parallel
    
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

#pragma parallel
    
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
    
#pragma parallel
    
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
    
#pragma parallel
    
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
          e_flux_uP[k][i][j] += h_west[k][i][j] * uu[k][i][j] * 0.5*(pi[i][j]+pp[im1][j]) * avg_fac_e*dt;
          e_flux_vP[k][i][j] += h_south[k][i][j] * vv[k][i][j]  * 0.5*(pi[i][j]+pp[i][jm1]) * avg_fac_e*dt;
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
  uint rp_iters = 0;
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




/**
 *
 * setGhostPoints
 *
 * Extends the input matrices uu, vv, hh, and bb (if tracers are active) by 2*Ng points in the x and y directions
 * and sets those points in order to enforce boundary conditions. The results are stored in
 * uu_g, vv_g, hh_g, and bb_g (if tracers are active) respectively.
 *
 */
void setGhostPoints (real *** uu, real *** vv, real *** hh, real *** bb,
                     real *** uu_g, real *** vv_g, real *** hh_g, real *** bb_g)
{
  // Looping variables
  int i,j,k;
  
  // Create work arrays and set ghost points in each layer
  for (k = 0; k < Nlay; k ++)
  {
    
#pragma parallel
    
    // Copy to work arrays with room for ghost points
    for (i = 0; i < Nx; i ++)
    {
      for (j = 0; j < Ny; j ++)
      {
        uu_g[k][i+Ng][j+Ng] = uu[k][i][j];
        vv_g[k][i+Ng][j+Ng] = vv[k][i][j];
        hh_g[k][i+Ng][j+Ng] = hh[k][i][j];
        if (useTracer)
        {
          bb_g[k][i+Ng][j+Ng] = bb[k][i][j];
        }
      }
    }
    
#pragma parallel
    
    // Extend grid in the y-direction.
    // Loop over just i-indices within the computational domain.
    for (i = Ng; i < Ng+Nx; i ++)
    {
      
      // With N/S walls
      if (useWallNS)
      {
        // Explicitly set northern boundary velocities to zero
        vv_g[k][i][Ng] = 0;
        vv_g[k][i][Ng+Ny] = 0;
        
        // Set ghost points
        for (j = 0; j < Ng; j ++)
        {
          // Velocities and layer thicknesses are just reflected at the wall.
          // The layer thickness ghost points are used for calculating h on q-points.
          // The velocity ghost points are used to average u^2 in the Bernoulli potential (when HK83 discretization is used),
          // and are consistent with harmonic and biharmonic no-stress boundary conditions.
          uu_g[k][i][j] = uu_g[k][i][2*Ng-j-1];
          vv_g[k][i][j] = -vv_g[k][i][2*Ng-j];
          hh_g[k][i][j] = hh_g[k][i][2*Ng-j-1];
          uu_g[k][i][Ny+Ng+j] = uu_g[k][i][Ny+Ng-j-1];
          vv_g[k][i][Ny+Ng+j] = -vv_g[k][i][Ny+Ng-j];
          hh_g[k][i][Ny+Ng+j] = hh_g[k][i][Ny+Ng-j-1];
          if (useTracer)
          {
            bb_g[k][i][j] = bb_g[k][i][2*Ng-j-1];
            bb_g[k][i][Ny+Ng+j] = bb_g[k][i][Ny+Ng-j-1];
          }
        }
      } // With N/S walls
      
      // Without N/S walls
      else
      {
        // Set ghost points
        for (j = 0; j < Ng; j ++)
        {
          // Velocities and layer thicknesses are just reflected at the wall.
          // The layer thickness ghost points are used for calculating h on q-points.
          // The velocity ghost points are used to average u^2 in the Bernoulli potential (when HK83 discretization is used),
          // and are consistent with harmonic and biharmonic no-stress boundary conditions.
          uu_g[k][i][j] = uu_g[k][i][Ny+j];
          vv_g[k][i][j] = vv_g[k][i][Ny+j];
          hh_g[k][i][j] = hh_g[k][i][Ny+j];
          uu_g[k][i][Ny+Ng+j] = uu_g[k][i][Ng+j];
          vv_g[k][i][Ny+Ng+j] = vv_g[k][i][Ng+j];
          hh_g[k][i][Ny+Ng+j] = hh_g[k][i][Ng+j];
          if (useTracer)
          {
            bb_g[k][i][j] = bb_g[k][i][Ny+j];
            bb_g[k][i][Ny+Ng+j] = bb_g[k][i][Ng+j];
          }
        }
      } // Without N/S walls
      
    }
    
#pragma parallel
    
    // Extend grid in the x-direction
    // Loop over ALL j-points, including ghost points.
    for (j = 0; j < Ny+2*Ng; j ++)
    {
      
      // With E/W walls
      if (useWallEW)
      {
        // Set velocities to zero explicitly on the walls
        uu_g[k][Ng][j] = 0;
        uu_g[k][Nx+Ng][j] = 0;
        
        // Set ghost points
        for (i = 0; i < Ng; i ++)
        {
          uu_g[k][i][j] = -uu_g[k][2*Ng-i][j];
          vv_g[k][i][j] = vv_g[k][2*Ng-i-1][j];
          hh_g[k][i][j] = hh_g[k][2*Ng-i-1][j];
          uu_g[k][Nx+Ng+i][j] = -uu_g[k][Nx+Ng-i][j];
          vv_g[k][Nx+Ng+i][j] = vv_g[k][Nx+Ng-i-1][j];
          hh_g[k][Nx+Ng+i][j] = hh_g[k][Nx+Ng-i-1][j];
          if (useTracer)
          {
            bb_g[k][i][j] = bb_g[k][2*Ng-i-1][j];
            bb_g[k][Nx+Ng+i][j] = bb_g[k][Nx+Ng-i-1][j];
          }
        }
      } // With E/W walls
      
      // Without E/W walls
      else
      {
        // Set ghost points
        for (i = 0; i < Ng; i ++)
        {
          uu_g[k][i][j] = uu_g[k][Nx+i][j];
          vv_g[k][i][j] = vv_g[k][Nx+i][j];
          hh_g[k][i][j] = hh_g[k][Nx+i][j];
          uu_g[k][Nx+Ng+i][j] = uu_g[k][Ng+i][j];
          vv_g[k][Nx+Ng+i][j] = vv_g[k][Ng+i][j];
          hh_g[k][Nx+Ng+i][j] = hh_g[k][Ng+i][j];
          if (useTracer)
          {
            bb_g[k][i][j] = bb_g[k][Nx+i][j];
            bb_g[k][Nx+Ng+i][j] = bb_g[k][Ng+i][j];
          }
        }
      } // Without E/W walls

    }
  }
}








/**
 *
 * calcEta
 *
 * Calculates the upper surface elevation for each layer. The result is stored in eta.
 *
 */
void calcEta (real ** hhb, real *** hh, real *** eta, uint Nlay, uint Nx, uint Ny)
{
  // Looping variables
  int i,j,k;
  
  for (k = Nlay; k >= 0; k --)
  {
#pragma parallel
    for (i = 0; i < Nx; i ++)
    {
      for (j = 0; j < Ny; j ++)
      {
        if (k == Nlay)
        {
          eta[k][i][j] = hhb[i][j];
        }
        else
        {
          eta[k][i][j] = eta[k+1][i][j] + hh[k][i][j];
        }
      }
    }
  }
}








/**
 * tderiv
 *
 * Calculates the time derivatives at all u, v and h gridpoints using the
 * Arakawa-based formulation.
 *
 */
void tderiv (const real t, const real * data, real * dt_data, const uint numvars)
{
  int i, j, k, i0, ip1, im1, ip2, im2, jp1, jm1, j0, jp2, jm2;
  real * pdata = NULL;
  
  // For vorticity points at northern/southern boundaries
  int jn = Ny + Ng, js = Ng;
  real Fnw, Fsw, Fne, Fse, G_n, G_s;
  
  // For storing tendencies
  real rhs_h, rhs_u, rhs_v, rhs_b;
  
  // "effective" surface mixed layer and bottom boundary layer thicknesses,
  // accounting for minimimum layer thicknesses
  real heff_west, heff_south, z_sml_west, z_sml_south, z_bbl_west, z_bbl_south, hF_west_sum, hF_south_sum;
  
  // Surface and bottom velocities - used to compute quadratic drag
  real u_surf, v_surf, u_bot, v_bot;
  
  // Used to compute Smagorinsky viscosities
  real DD_q, DD_h, A2smag_fac, A4smag_fac, A4sqrt_const;
  
  // Used to compute tracer fluxes on cell faces
  real bb_xm, bb_xp, bb_ym, bb_yp;
  
  // Used for diabatic velocity calculations
  real wdia_hRelax;
  real uu_p, uu_m, du_p, du_m;
  real vv_p, vv_m, dv_p, dv_m;
  real bb_p, bb_m, db_p, db_m;
  
  // Used for wind stress interpolation
  real taux_interp, tauy_interp;
  real nflt;
  int nm1, np1;
  
  // Used for energy diagnostics
  real pres_0, pres_im1, pres_jm1;
  
  
  
  
  
  // Construct matrices pointing to 'data' input vector
  pdata = (real *) data;
  vec2mat3(pdata,&uu_d,Nlay,Nx,Ny);
  vec2mat3(pdata+=N,&vv_d,Nlay,Nx,Ny);
  vec2mat3(pdata+=N,&hh_d,Nlay,Nx,Ny);
  if (useTracer)
  {
    vec2mat3(pdata+=N,&bb_d,Nlay,Nx,Ny);
  }
  
  // Construct matrices pointing to 'dt_data' input vector
  pdata = dt_data;
  vec2mat3(pdata,&dt_uu_w,Nlay,Nx,Ny);
  vec2mat3(pdata+=N,&dt_vv_w,Nlay,Nx,Ny);
  vec2mat3(pdata+=N,&dt_hh_w,Nlay,Nx,Ny);
  if (useTracer)
  {
    vec2mat3(pdata+=N,&dt_bb_w,Nlay,Nx,Ny);
  }
  
  // Create work arrays and set ghost points in each layer
  setGhostPoints (uu_d,vv_d,hh_d,bb_d,uu_w,vv_w,hh_w,bb_w);
  
  // Calculate layer interface heights - needed for pressure calculation
  calcEta(hhb_w,hh_w,eta_w,Nlay,Nx+2*Ng,Ny+2*Ng);
  
  // Compute layer thicknesses on cell faces
  for (k = 0; k < Nlay; k ++)
  {
    calcFaceThickness(hh_w[k],h_west[k],h_south[k],true,Nx,Ny,uu_w[k],vv_w[k]);
  }
  
  
  //////////////////////////////////
  ///// SML/BBL THICKNESS CODE /////
  //////////////////////////////////
  
  
  // If a non-zero SML/BBL is specified then we need to compute layer interface heights on
  // cell faces
  if ((hsml > 0) || (hbbl > 0))
  {
    calcEta(hhb_west,h_west,eta_west,Nlay,Nx+2*Ng,Ny+2*Ng);
    calcEta(hhb_south,h_south,eta_south,Nlay,Nx+2*Ng,Ny+2*Ng);
  }
  
  // If a non-zero SML is specified then calculate "forcing thickness" of each layer at the surface
  if (hsml > 0)
  {

#pragma parallel

    // Loop over all points in the domain, including ghost points
    for (i = 1; i < Nx+2*Ng-1; i ++)
    {
      for (j = 1; j < Ny+2*Ng-1; j ++)
      {

        z_sml_west = hhs_west[i][j] - hsml;
        z_sml_south = hhs_south[i][j] - hsml;
        hF_west_sum = 0;
        hF_south_sum = 0;
        for (k = 0; k < Nlay; k ++)
        {
          // "forcing" layer thicknesses - fraction of the layer thickness that lies within the
          // "surface mixed layer", and so will receive momentum forcing.
          if (hF_west_sum < 1.0)
          {
            z_sml_west -= fmin(hmin_surf,h_west[k][i][j]);
            heff_west = fmax(eta_west[k][i][j],z_sml_west) - fmax(eta_west[k+1][i][j],z_sml_west);
            heff_west = fmax(heff_west-hmin_surf,0);
            hFsurf_west[k][i][j] =  heff_west / hsml;
            hF_west_sum += hFsurf_west[k][i][j];
          }
          else
          {
            hFsurf_west[k][i][j] = 0;
          }
          if (hF_south_sum < 1.0)
          {
            z_sml_south -= fmin(hmin_surf,h_south[k][i][j]);
            heff_south = fmax(eta_south[k][i][j],z_sml_south) - fmax(eta_south[k+1][i][j],z_sml_south);
            heff_south = fmax(heff_south-hmin_surf,0);
            hFsurf_south[k][i][j] = heff_south / hsml;
            hF_south_sum += hFsurf_south[k][i][j];
          }
          else
          {
            hFsurf_south[k][i][j] = 0;
          }
        }

      }
    }
  }
  
  // If a non-zero SML is specified then calculate "forcing thickness" of each layer at the surface
  if (hbbl > 0)
  {
      
#pragma parallel

    // Loop over all points in the domain, including ghost points
    for (i = 1; i < Nx+2*Ng-1; i ++)
    {
      for (j = 1; j < Ny+2*Ng-1; j ++)
      {
        
        z_bbl_west = hhb_west[i][j]+hbbl;
        z_bbl_south = hhb_south[i][j]+hbbl;
        hF_west_sum = 0;
        hF_south_sum = 0;
        for (k = Nlay-1; k >= 0; k --)
        {
          // "forcing" layer thicknesses - fraction of the layer thickness that lies within the
          // "bottom boundary layer", and so will receive momentum forcing.
          if (hF_west_sum < 1.0)
          {
            z_bbl_west += fmin(hmin_bot,h_west[k][i][j]);
            heff_west = fmin(eta_west[k][i][j],z_bbl_west) - fmin(eta_west[k+1][i][j],z_bbl_west);
            heff_west = fmax(heff_west-hmin_bot,0);
            hFbot_west[k][i][j] =  heff_west / hbbl;
            hF_west_sum += hFsurf_west[k][i][j];
          }
          else
          {
            hFbot_west[k][i][j] = 0;
          }
          if (hF_south_sum < 1.0)
          {
            z_bbl_south += fmin(hmin_bot,h_south[k][i][j]);
            heff_south = fmin(eta_south[k][i][j],z_bbl_south) - fmin(eta_south[k+1][i][j],z_bbl_south);
            heff_south = fmax(heff_south-hmin_bot,0);
            hFbot_south[k][i][j] = heff_south / hbbl;
            hF_south_sum += hFsurf_south[k][i][j];
          }
          else
          {
            hFbot_south[k][i][j] = 0;
          }
        }
        
      }
    }

  }
  

  ///////////////////////////////
  ///// QUADRATIC DRAG CODE /////
  ///////////////////////////////
  

  /// Calculate absolute velocity at the surface by averaging over SML
  if (CdSurf > 0)
  {
    
#pragma parallel

    // Determine squared "surface" velocities on u/v points
    for (i = 0; i < Nx+2*Ng; i ++)
    {
      for (j = 0; j < Ny+2*Ng; j ++)
      {
        // If we're using an SML then average velocity vertically over SML depth
        if (hsml > 0)
        {
          u_surf = 0;
          v_surf = 0;
          for (k = 0; k < Nlay; k ++)
          {
            // "surface" velocity for computing drag - just the velocity averaged vertically over
            // the top hsml of the water column
            u_surf += uu_w[k][i][j] * hFsurf_west[k][i][j];
            v_surf += vv_w[k][i][j] * hFsurf_south[k][i][j];
          }
        }
        // Otherwise just use velocity from top layer
        else
        {
          u_surf = uu_w[0][i][j];
          v_surf = vv_w[0][i][j];
        }
        
        // Store squared surface velocity components, relative to lid velocity components
        usq_surf[i][j] = (u_surf-uLid_g[i][j])*(u_surf-uLid_g[i][j]);
        vsq_surf[i][j] = (v_surf-vLid_g[i][j])*(v_surf-vLid_g[i][j]);
      }
    }
    
    // Calculate absolute "surface" velocity
    for (i = 0; i < Nx+2*Ng-1; i ++)
    {
      i0 = i;
      ip1 = i+1;
      
      for (j = 0; j < Ny+2*Ng-1; j ++)
      {
        j0 = j;
        jp1 = j+1;
        
        uabs_surf[i][j] = sqrt( 0.5*(usq_surf[i0][j0]+usq_surf[ip1][j0]) + 0.5*(vsq_surf[i0][j0]+vsq_surf[i0][jp1]) );
      }
    }
    
  } // end if (CdSurf > 0)
  
    /// Calculate absolute velocity at the surface by averaging over BBL
  if (CdBot > 0)
  {
    
#pragma parallel
    
    // Determine squared "bottom" velocities on u/v points
    for (i = 0; i < Nx+2*Ng; i ++)
    {
      for (j = 0; j < Ny+2*Ng; j ++)
      {
        if (hbbl > 0)
        {
          u_bot = 0;
          v_bot = 0;
          for (k = 0; k < Nlay; k ++)
          {
            // "bottom" velocity for computing drag - just the velocity averaged vertically over the
            // bottom hbbl of the water column
            u_bot += uu_w[k][i][j] * hFbot_west[k][i][j];
            v_bot += vv_w[k][i][j] * hFbot_south[k][i][j];
          }
        }
        else
        {
          u_bot = uu_w[Nlay-1][i][j];
          v_bot = vv_w[Nlay-1][i][j];
        }
          
        // Store squared bottom velocities
        usq_bot[i][j] = u_bot*u_bot;
        vsq_bot[i][j] = v_bot*v_bot;
      }
    }
    
    // Calculate absolute "bottom" velocity
    for (i = 0; i < Nx+2*Ng-1; i ++)
    {
      i0 = i;
      ip1 = i+1;
      
      for (j = 0; j < Ny+2*Ng-1; j ++)
      {
        j0 = j;
        jp1 = j+1;
        
        uabs_bot[i][j] = sqrt( 0.5*(usq_bot[i0][j0]+usq_bot[ip1][j0]) + 0.5*(vsq_bot[i0][j0]+vsq_bot[i0][jp1]) );
      }
    }    
    
  } // End if (CdBot > 0)
  
  
  
  //////////////////////////////////
  ///// DIABATIC VELOCITY CODE /////
  //////////////////////////////////

  if (useWDia)
  {
    
    // Set up indices for linear interpolation (in time) of diapycnal velocities
    if ((wDiaNrecs==1) || (wDiaPeriod<=0))
    {
      nflt = 0;
      np1 = 1;
      nm1 = 0;
    }
    else
    {
      nflt = (fmod(t,wDiaPeriod)/wDiaPeriod) * wDiaNrecs; // Floating point "index" in time
      nm1 = (int) floor(nflt); // Indices of adjacent wind stress time points
      np1 = (int) ceil(nflt);
      if (np1 == nm1)
      {
        np1 = (nm1 + 1); // Catch cases in which nflt is an exact integer
      }
    }
    
    // Calculate diapycnal velocities at  centers of cell upper faces
    for (i = 0; i < Nx; i ++)
    {
      // Offsets account for vectors/matrices that include ghost points.
      i0 = i + Ng;
      
      for (j = 0; j < Ny; j ++)
      {
        // Offsets account for vectors/matrices that include ghost points.
        j0 = j + Ng;

        // Initialize velocity at sea floor
        // N.B. This code also does linear interpolation in time when wdia_ff varies in time
        wdia[Nlay][i][j] = (np1-nflt)*wdia_ff[nm1][Nlay][i][j] + (nflt-nm1)*wdia_ff[np1%wDiaNrecs][Nlay][i][j];
        
        // Keeps track of diabatic velocity due to thickness relaxation,
        // which must be calculated cumulatively
        wdia_hRelax = 0;
        
        // Assign diabatic velocity across each layer interface
        for (k = Nlay-1; k >= 0; k --)
        {
          // Account for fixed fluxes (will be zero if no fixed flux is specified)
          // N.B. This code also does linear interpolation in time when wdia_ff varies in time
          wdia[k][i][j] = (np1-nflt)*wdia_ff[nm1][k][i][j] + (nflt-nm1)*wdia_ff[np1%wDiaNrecs][k][i][j];
          
          // Account for thickness relaxation
          if (useRelax && (hTime[i][j] > 0))
          {
            wdia_hRelax += (hh_w[k][i0][j0] - hRelax[k][i][j]) / hTime[i][j];
            wdia[k][i][j] += wdia_hRelax;
          }
          
          // Account for interface relaxation
          if (useRelax && (eTime[k][i][j] > 0))
          {
            wdia[k][i][j] += (eta_w[k][i0][j0] - eRelax[k][i][j]) / eTime[k][i][j];
          }
        }
      }
    }
      
    
    // Calculate diapycnal velocities at u and v points
    for (i = 0; i < Nx; i ++)
    {
      // At domain boundaries we just average across the boundary, because
      // boundary-normal component of velocity is guaranteed to be zero
      im1 = (i+Nx-1) % Nx;
      
      for (j = 0; j < Ny; j ++)
      {
        jm1 = (j+Ny-1) % Ny;
        
        for (k = 0; k < Nlay+1; k ++)
        {
          wdia_u[k][i][j] = 0.5 * (wdia[k][i][j] + wdia[k][im1][j]);
          wdia_v[k][i][j] = 0.5 * (wdia[k][i][j] + wdia[k][i][jm1]);
        }
      }
    }
  
  }
    
  
  
  
  /////////////////////////////////
  ///// MAIN LOOP OVER LAYERS /////
  /////////////////////////////////


  for (k = 0; k < Nlay; k ++)
  {
    
#pragma parallel
    
    // Generate q, u* and v*
    for (i = 1; i < Nx+2*Ng-1; i ++)
    {
      i0 = i;
      im1 = i-1;
      ip1 = i+1;
      
      for (j = 1; j < Ny+2*Ng-1; j ++)
      {
        j0 = j;
        jp1 = j+1;
        jm1 = j-1;
        
        // Interior approximation to the relative vorticity
        // Note that this will enforce zeta=0 on the boundaries via the ghost points
        zeta[i0][j0] = (uu_w[k][i0][jm1]-uu_w[k][i0][j0]) / dy + (vv_w[k][i0][j0]-vv_w[k][im1][j0]) / dx;
        
        // Approximation to h at the q-gridpoints
        hh_q[i0][j0] = 0.25 * (hh_w[k][i0][j0]+hh_w[k][im1][jm1]+hh_w[k][im1][j0]+hh_w[k][i0][jm1]);
        
        // Potential vorticity
        qq[i0][j0] = (2*Omega_z_w[i0][j0] + zeta[i0][j0]) / hh_q[i0][j0];
        
        // Zonal/meridional mass fluxes at u-/v-gridpoints
        huu[i0][j0] = uu_w[k][i0][j0]*h_west[k][i0][j0];
        hvv[i0][j0] = vv_w[k][i0][j0]*h_south[k][i0][j0];
      }
    }
    
#pragma parallel
    
    // Generate PV coefficients and potential function
    for (i = 1; i < Nx+2*Ng-1; i ++)
    {
      ip1 = i+1;
      im1 = i-1;
      
      for (j = 1; j < Ny+2*Ng-1; j ++)
      {
        jp1 = j+1;
        jm1 = j-1;
        
        // Calculate pressure. Here we perform this calculation cumulatively, using the pressure set
        // by the previous layer (k-1) and adding to it to get the pressure in the current layer (k).
        if (k == 0)
        {
          // In rigid lid case we'll add the surface pressure later
          if (useRL)
          {
            pp[i][j] = 0;
          }
          else
          {
            pp[i][j] = gg[0]*eta_w[0][i][j];
          }
        }
        else
        {
          pp[i][j] += gg[k]*eta_w[k][i][j];
        }
        
        // Potential function
        switch (momentumScheme)
        {
          case MOMENTUM_HK83:
          {
            KE_B[i][j] = ( SQUARE(uu_w[k][ip1][j]) + SQUARE(uu_w[k][i][j])
                           + 0.25 * ( SQUARE(uu_w[k][ip1][jp1]) + SQUARE(uu_w[k][i][jp1])
                                     + SQUARE(uu_w[k][ip1][jm1]) + SQUARE(uu_w[k][i][jm1]) )
                           + SQUARE(vv_w[k][i][jp1]) + SQUARE(vv_w[k][i][j])
                           + 0.25 * ( SQUARE(vv_w[k][ip1][jp1]) + SQUARE(vv_w[k][ip1][j])
                                     + SQUARE(vv_w[k][im1][jp1]) + SQUARE(vv_w[k][im1][j]) )
                          ) / 6;
            break;
          }
            // This version doesn't use ghost h-points, but is subject to an
            // internal symmetric computational instability (Hollingsworth and Kallberg 1983)
          case MOMENTUM_AL81:
          case MOMENTUM_S75e:
          {
            KE_B[i][j] = ( SQUARE(uu_w[k][ip1][j]) + SQUARE(uu_w[k][i][j]) ) / 4
                        + ( SQUARE(vv_w[k][i][jp1]) + SQUARE(vv_w[k][i][j]) ) / 4;
            break;
          }
          default:
          {
            fprintf(stderr,"ERROR: Unknown spatial discretization method\n");
            break;
          }
        }
        
        // Montgomery potential
        MM_B[i][j] = pp[i][j];
        
        // Add Salmon (2002) term
        MM_B[i][j] -= geff[k] * POW4(h0) / POW3(hh_w[k][i][j]) / 3;
        
        // PV coefficients
        if ((momentumScheme == MOMENTUM_AL81) || (momentumScheme == MOMENTUM_HK83))
        {
          alpha_w[i][j] = (2*qq[ip1][jp1]+qq[i][jp1]+2*qq[i][j]+qq[ip1][j]) / 24;
          beta_w[i][j] = (qq[i][jp1]+2*qq[im1][jp1]+qq[im1][j]+2*qq[i][j]) / 24;
          gamma_w[i][j] = (2*qq[i][jp1]+qq[im1][jp1]+2*qq[im1][j]+qq[i][j]) / 24;
          delta_w[i][j] = (qq[ip1][jp1]+2*qq[i][jp1]+qq[i][j]+2*qq[ip1][j]) / 24;
          epsilon_w[i][j] = (qq[ip1][jp1]+qq[i][jp1]-qq[i][j]-qq[ip1][j]) / 24;
          phi_w[i][j] = (-qq[ip1][jp1]+qq[i][jp1]+qq[i][j]-qq[ip1][j]) / 24;
        }
      }
      
    }
    
    // Horizontal tension and shearing strain - required for Laplacian/biharmonic viscosity tensor,
    // and for definition of Smagornisky viscosities
    if (useA2 || useA4)
    {
      for (i = 1; i < Nx+2*Ng-1; i ++)
      {
        ip1 = i+1;
        im1 = i-1;
        
        for (j = 1; j < Ny+2*Ng-1; j ++)
        {
          jp1 = j+1;
          jm1 = j-1;
      
          DD_T[i][j] = (uu_w[k][ip1][j]-uu_w[k][i][j])/dx - (vv_w[k][i][jp1]-vv_w[k][i][j])/dy;
          DD_S[i][j] = (uu_w[k][i][j]-uu_w[k][i][jm1])/dy + (vv_w[k][i][j]-vv_w[k][im1][j])/dx;
        }
      }
    }

    //  Calculate Laplacian and biharmonic viscosities
    if ((A2smag > 0) || (A4smag > 0))
    {
      // These constants are multiplied by the deformation rate to calculate the Smagorinsky viscosities
      A2smag_fac = SQUARE(A2smag*fmax(dx,dy)/_PI);
      A4smag_fac = SQUARE(A4smag/_PI) * POW4(fmax(dx,dy)) / 8;
      
      for (i = 1; i < Nx+2*Ng-1; i ++)
      {
        ip1 = i+1;
        im1 = i-1;
        
        for (j = 1; j < Ny+2*Ng-1; j ++)
        {
          jp1 = j+1;
          jm1 = j-1;
          
          // Cell corners
          DD_q = sqrt( SQUARE(0.25 * (DD_T[i][j] + DD_T[im1][j] + DD_T[i][jm1] + DD_T[im1][jm1])) + SQUARE(DD_S[i][j]) );
          A2_q[i][j] = fmax(A2const,A2smag_fac*DD_q);
          A4sqrt_q[i][j] = sqrt(fmax(A4const,A4smag_fac*DD_q));
          
          // Cell centers
          DD_h = sqrt( SQUARE(DD_T[i][j]) + SQUARE(0.25 * (DD_S[i][j] + DD_S[ip1][j] + DD_S[i][jp1] + DD_S[ip1][jp1])) );
          A2_h[i][j] = fmax(A2const,A2smag_fac*DD_h);
          A4sqrt_h[i][j] = sqrt(fmax(A4const,A4smag_fac*DD_h));
        }
      }
    }
    // No Smagorinsky viscosities, so just set viscosities to constant values everywhere
    else
    {
      if (useA2)
      {
        for (i = 1; i < Nx+2*Ng-1; i ++)
        {
          for (j = 1; j < Ny+2*Ng-1; j ++)
          {
            A2_h[i][j] = A2const;
            A2_q[i][j] = A2const;
          }
        }
      }
      if (useA4)
      {
        A4sqrt_const = sqrt(A4const);
        for (i = 1; i < Nx+2*Ng-1; i ++)
        {
          for (j = 1; j < Ny+2*Ng-1; j ++)
          {
            A4sqrt_h[i][j] = A4sqrt_const;
            A4sqrt_q[i][j] = A4sqrt_const;
          }
        }
      }
    }
    
    // Calculate intermediate "velocities" (approximately the Laplacians of u and v) if we're using biharmonic viscosity
    if (useA4)
    {
      for (i = 1; i < Nx+2*Ng-1; i ++)
      {
        i0 = i;
        ip1 = i+1;
        im1 = i-1;
        
        for (j = 1; j < Ny+2*Ng-1; j ++)
        {
          j0 = j;
          jp1 = j+1;
          jm1 = j-1;
          
          UU4[i0][j0] = (
                           ( A4sqrt_h[i0][j0]*hh_w[k][i0][j0]*DD_T[i0][j0] - A4sqrt_h[im1][j0]*hh_w[k][im1][j0]*DD_T[im1][j0] ) / dx
                         + ( A4sqrt_q[i0][jp1]*hh_q[i0][jp1]*DD_S[i0][jp1] - A4sqrt_q[i0][j0]*hh_q[i0][j0]*DD_S[i0][j0] ) / dy
                        )
                      / h_west[k][i0][j0];
          VV4[i0][j0] = (
                          ( A4sqrt_q[ip1][j0]*hh_q[ip1][j0]*DD_S[ip1][j0] - A4sqrt_q[i0][j0]*hh_q[i0][j0]*DD_S[i0][j0] ) / dx
                        - ( A4sqrt_h[i0][j0]*hh_w[k][i0][j0]*DD_T[i0][j0] - A4sqrt_h[i0][jm1]*hh_w[k][i0][jm1]*DD_T[i0][jm1] ) / dy
                        )
                      / h_south[k][i0][j0];
          
          // Old version with no thickness weighting, as in Griffies and Hallberg (2000)
          /*
          UU4[i0][j0] = (
                         ( A4sqrt_h[i0][j0]*DD_T[i0][j0] - A4sqrt_h[im1][j0]*DD_T[im1][j0] ) / dx
                         + ( A4sqrt_q[i0][jp1]*DD_S[i0][jp1] - A4sqrt_q[i0][j0]*DD_S[i0][j0] ) / dy
                         );
          VV4[i0][j0] = (
                         ( A4sqrt_q[ip1][j0]*DD_S[ip1][j0] - A4sqrt_q[i0][j0]*DD_S[i0][j0] ) / dx
                         - ( A4sqrt_h[i0][j0]*DD_T[i0][j0] - A4sqrt_h[i0][jm1]*DD_T[i0][jm1] ) / dy
                         );
          */
        }
      }
    }

    // Biharmonic horizontal tension and shearing strain - required to construct biharmonic viscous stress tensor
    if (useA4)
    {
      for (i = 1; i < Nx+2*Ng-1; i ++)
      {
        ip1 = i+1;
        im1 = i-1;
        
        for (j = 1; j < Ny+2*Ng-1; j ++)
        {
          jp1 = j+1;
          jm1 = j-1;
          
          DD4_T[i][j] = (UU4[ip1][j]-UU4[i][j])/dx - (VV4[i][jp1]-VV4[i][j])/dy;
          DD4_S[i][j] = (UU4[i][j]-UU4[i][jm1])/dy + (VV4[i][j]-VV4[im1][j])/dx;
        }
      }
    }
    
    // Advective tracer tendency terms
    if (useTracer)
    {
      
#pragma parallel
      
      for (i = 1; i < Nx+2*Ng-1; i ++)
      {
        im1 = i-1;
        ip1 = i+1;
        
        for (j = 1; j < Ny+2*Ng-1; j ++)
        {
          jm1 = j-1;
          jp1 = j+1;
          
          
          // Centered differencing
          if (tracerScheme == TRACER_AL81)
          {
            // Tracer gradients and products with velocities
            db_dx[i][j] = (bb_w[k][i][j] - bb_w[k][im1][j]) / dx;
            db_dy[i][j] = (bb_w[k][i][j] - bb_w[k][i][jm1]) / dy;
            udb_dx[i][j] = uu_w[k][i][j] * db_dx[i][j];
            vdb_dy[i][j] = vv_w[k][i][j] * db_dy[i][j];
          }
          
          // Product of layer thickness with layer mid-depth
          if (useBuoyancy)
          {
            hz[i][j] = hh_w[k][i][j] * 0.5 * (eta_w[k][i][j] + eta_w[k+1][i][j]);
          }
          
          // If using biharmonic diffusion then calculate intermediate "tracer" quantity required for biharmonic operator
          if (K4 > 0)
          {
            BB4[i][j] = (
                          ( h_west[k][ip1][j]*(bb_w[k][ip1][j]-bb_w[k][i][j]) - h_west[k][i][j]*(bb_w[k][i][j]-bb_w[k][im1][j]) ) / dxsq
                        + ( h_south[k][i][jp1]*(bb_w[k][i][jp1]-bb_w[k][i][j]) - h_south[k][i][j]*(bb_w[k][i][j]-bb_w[k][i][jm1]) ) / dysq
                        )
                      / hh_w[k][i][j];

          }
          
        }
      }
      
      // Third-order upwinding
      if (tracerScheme == TRACER_UP3)
      {
        
#pragma parallel
        
        // Compute second derivative for third-order upwinding
        for (i = 1; i < Nx+2*Ng-1; i ++)
        {
          im1 = i-1;
          ip1 = i+1;
          
          for (j = 1; j < Ny+2*Ng-1; j ++)
          {
            jm1 = j-1;
            jp1 = j+1;
            
            d2bx[i][j] = (bb_w[k][ip1][j]-2*bb_w[k][i][j]+bb_w[k][im1][j]);
            d2by[i][j] = (bb_w[k][i][jp1]-2*bb_w[k][i][j]+bb_w[k][i][jm1]);
          }
        }
        
#pragma parallel
        
        // Compute total tracer on cell faces
        for (i = 1; i < Nx+2*Ng-1; i ++)
        {
          im1 = i-1;
          ip1 = i+1;
          
          for (j = 1; j < Ny+2*Ng-1; j ++)
          {
            jm1 = j-1;
            jp1 = j+1;
            
            hvb[i][j] = 0.5 * (bb_w[k][i][j]+bb_w[k][i][jm1]);
            hvb[i][j] -= (1/6.0) * ( (vv_w[k][i][j] > 0) ? d2by[i][jm1] : d2by[i][j] );
            hvb[i][j] *= hvv[i][j];
            hub[i][j] = 0.5 * (bb_w[k][i][j]+bb_w[k][im1][j]);
            hub[i][j] -= (1/6.0) * ( (uu_w[k][i][j] > 0) ? d2bx[im1][j] : d2bx[i][j] );
            hub[i][j] *= huu[i][j];
          }
        }
        
      } // end if (tracerScheme == THICKNESS_UP3)
      
      
      if (tracerScheme == TRACER_KT00)
      {
        
#pragma parallel
        
        // Determine limited slopes at cell centres
        for (i = 1; i < Nx+2*Ng-1; i ++)
        {
          im1 = i-1;
          ip1 = i+1;
          
          for (j = 1; j < Ny+2*Ng-1; j ++)
          {
            jm1 = j-1;
            jp1 = j+1;
            
            db_dx[i][j] = minmod( KT00_sigma * (bb_w[k][ip1][j]-bb_w[k][i][j]),
                                 0.5 * (bb_w[k][ip1][j]-bb_w[k][im1][j]),
                                 KT00_sigma * (bb_w[k][i][j]-bb_w[k][im1][j]) );
            db_dy[i][j] = minmod( KT00_sigma * (bb_w[k][i][jp1]-bb_w[k][i][j]),
                                 0.5 * (bb_w[k][i][jp1]-bb_w[k][i][jm1]),
                                 KT00_sigma * (bb_w[k][i][j]-bb_w[k][i][jm1]) );
          }
        }
        
#pragma parallel
        
        // Interpolate hb to cell faces
        for (i = 1; i < Nx+2*Ng-1; i ++)
        {
          im1 = i-1;
          ip1 = i+1;
          
          for (j = 1; j < Ny+2*Ng-1; j ++)
          {
            jm1 = j-1;
            jp1 = j+1;
            
            // Interpolate h to cell west faces
            bb_xm = bb_w[k][im1][j] + 0.5*db_dx[im1][j];
            bb_xp = bb_w[k][i][j] - 0.5*db_dx[i][j];
            hub[i][j] = (uu_w[k][i][j] > 0) ? bb_xm : bb_xp;
            hub[i][j] *= huu[i][j];
            
            // Interpolate h to cell south faces
            bb_ym = bb_w[k][i][jm1] + 0.5*db_dy[i][jm1];
            bb_yp = bb_w[k][i][j] + 0.5*db_dy[i][j];
            hvb[i][j] = (vv_w[k][i][j] > 0) ? bb_ym : bb_yp;
            hvb[i][j] *= hvv[i][j];
          }
        }
        
      } // end if (tracerScheme == TRACER_KT00)

      
    } // end if (useTracer)

#pragma parallel
    
    // Calculate time derivatives
    for (i = 0; i < Nx; i ++)
    {
      // Offsets account for the fact that all vectors/matrices except dt_... matrices include ghost points.
      i0 = i + Ng;
      ip1 = i + Ng + 1;
      ip2 = i + Ng + 2;
      im1 = i + Ng - 1;
      im2 = i + Ng - 2;
      
      for (j = 0; j < Ny; j ++)
      {
        // Offsets account for the fact that all vectors/matrices except dt_... matrices include ghost points.
        j0 = j + Ng;
        jp1 = j + Ng + 1;
        jm1 = j + Ng - 1;
        jp2 = j + Ng + 2;
        jm2 = j + Ng - 2;
        
        
        
        ///////////////////////////////
        /// Interpolate wind stress ///
        ///////////////////////////////
        
        // Interpolate user-supplied wind stress in time
        if (useWind)
        {
          if ((tauNrecs==1) || (tauPeriod<=0))
          {
            taux_interp = taux[0][i][j];
            tauy_interp = tauy[0][i][j];
          }
          else
          {
            nflt = (fmod(t,tauPeriod)/tauPeriod) * tauNrecs; // Floating point "index" in time
            nm1 = (int) floor(nflt); // Indices of adjacent wind stress time points
            np1 = (int) ceil(nflt);            
            if (np1 == nm1)
            {
              np1 = (nm1 + 1); // Catch cases in which nflt is an exact integer
            }
            // Interpolate linearly
            taux_interp = (np1-nflt)*taux[nm1][i][j] + (nflt-nm1)*taux[np1%tauNrecs][i][j];
            tauy_interp = (np1-nflt)*tauy[nm1][i][j] + (nflt-nm1)*tauy[np1%tauNrecs][i][j];
          }
        }

        
        
        
        //////////////////////////
        /// Thickness tendency ///
        //////////////////////////
        
        // Initialize
        dt_hh_w[k][i][j] = 0;
        
        
        // Advection
        rhs_h = (huu[i0][j0]-huu[ip1][j0]) / dx + (hvv[i0][j0]-hvv[i0][jp1]) / dy;
        dt_hh_w[k][i][j] += rhs_h;
        // Do time-averaging if required
        if (dt_avg_h > 0)
        {
          h_tend_adv[k][i][j] += avg_fac_h*rhs_h*dt;
        }
        // Need to add h tendency to hu and hv tendencies
        // NOTE: This formulation assumes centered averaging of h to u and v points
        if (dt_avg_hu > 0)
        {
          hu_tend_dhdt[k][i][j] += 0.5 * uu_w[k][i0][j0]*avg_fac_hu*rhs_h*dt;
          hu_tend_dhdt[k][(i+Nx+1) % Nx][j] += 0.5 * uu_w[k][ip1][j0]*avg_fac_hu*rhs_h*dt;
        }
        if (dt_avg_hv > 0)
        {
          hv_tend_dhdt[k][i][j] += 0.5 * vv_w[k][i0][j0]*avg_fac_hv*rhs_h*dt;
          hv_tend_dhdt[k][i][(j+Ny+1) % Ny] += 0.5 * vv_w[k][i0][jp1]*avg_fac_hv*rhs_h*dt;
        }
        // Also contributes to the KE tendency
        // N.B. This formulation ignores grid point locations
        if (dt_avg_e > 0)
        {
          e_tend_adv[k][i][j] += 0.5*(SQUARE(uu_w[k][i0][j0])+SQUARE(vv_w[k][i0][j0]))*rhs_h*avg_fac_e*dt;
        }
        
        // Add diabatic velocity terms. Negative values mean no relaxation.
        if (useWDia)
        {
          rhs_h = - (wdia[k][i][j]-wdia[k+1][i][j]);
          dt_hh_w[k][i][j] += rhs_h;

          // Do time-averaging if required
          if (dt_avg_h > 0)
          {
            h_tend_relax[k][i][j] += avg_fac_h*rhs_h*dt;
          }
          
          // Need to add h tendency to hu and hv tendencies
          // NOTE: This formulation assumes centered averaging of h to u and v points
          if (dt_avg_hu > 0)
          {
            hu_tend_wdia[k][i][j] += 0.5 * uu_w[k][i0][j0]*avg_fac_hu*rhs_h*dt;
            hu_tend_wdia[k][(i+Nx+1) % Nx][j] += 0.5 * uu_w[k][ip1][j0]*avg_fac_hu*rhs_h*dt;
          }
          if (dt_avg_hv > 0)
          {
            hv_tend_wdia[k][i][j] += 0.5 * vv_w[k][i0][j0]*avg_fac_hv*rhs_h*dt;
            hv_tend_wdia[k][i][(j+Ny+1) % Ny] += 0.5 * vv_w[k][i0][jp1]*avg_fac_hv*rhs_h*dt;
          }
          
          // Calculate energy tendencies due to diapycnal fluxes
          if (dt_avg_e > 0)
          {
            // N.B. This formulation ignores grid point locations
            e_tend_wdiaKE[k][i][j] += 0.5*(SQUARE(uu_w[k][i0][j0])+SQUARE(vv_w[k][i0][j0]))*rhs_h*avg_fac_e*dt;
            
            e_tend_wdiaPE[k][i][j] += - geff[k]*(eta_w[k][i0][j0]*wdia[k][i0][j0]-eta_w[k+1][i0][j0]*wdia[k+1][i0][j0])*avg_fac_e*dt
                                  - (-geff[k]*POW4(h0)/POW3(hh_w[k][i][j])/3) * (-rhs_h) * avg_fac_e*dt;
          }
        }

        ///////////////////////////////
        /// Zonal momentum tendency ///
        ///////////////////////////////

        // Initialize
        dt_uu_w[k][i][j] = 0;
        
        // Gradient of Montgomery potential
        rhs_u = - (MM_B[i0][j0]-MM_B[im1][j0]) / dx;
        dt_uu_w[k][i][j] += rhs_u;
        if (dt_avg_hu > 0)
        {
          hu_tend_gradM[k][i][j] += h_west[k][i0][j0] * avg_fac_hu*rhs_u*dt;
        }
        if (dt_avg_e > 0)
        {
          e_tend_gradM[k][i][j] += rhs_u*h_west[k][i0][j0]*uu_w[k][i0][j0]*avg_fac_e*dt;
        }
        
        // Gradient of kinetic energy
        rhs_u = - (KE_B[i0][j0]-KE_B[im1][j0]) / dx;
        dt_uu_w[k][i][j] += rhs_u;
        if (dt_avg_hu > 0)
        {
          hu_tend_gradKE[k][i][j] += h_west[k][i0][j0] * avg_fac_hu*rhs_u*dt;
        }
        if (dt_avg_e > 0)
        {
          e_tend_adv[k][i][j] += rhs_u*h_west[k][i0][j0]*uu_w[k][i0][j0]*avg_fac_e*dt;
        }
        
        // Add Coriolis/vorticity advection terms
        switch (momentumScheme)
        {
          case MOMENTUM_AL81:
          case MOMENTUM_HK83:
          {
            rhs_u             = alpha_w[i0][j0]*hvv[i0][jp1]
                              + beta_w[i0][j0]*hvv[im1][jp1]
                              + gamma_w[i0][j0]*hvv[im1][j0]
                              + delta_w[i0][j0]*hvv[i0][j0]
                              - epsilon_w[i0][j0]*huu[ip1][j0]
                              + epsilon_w[im1][j0]*huu[im1][j0];
            break;
          }
          case MOMENTUM_S75e:
          {
            rhs_u             = 0.25 * qq[i0][j0] * (hvv[i0][j0] + hvv[im1][j0])
                              + 0.25 * qq[i0][jp1] * (hvv[i0][jp1] + hvv[im1][jp1]);
            break;
          }
          default:
          {
            fprintf(stderr,"ERROR: Unknown spatial discretization method\n");
            break;
          }
        }
        dt_uu_w[k][i][j] += rhs_u;
        if (dt_avg_hu > 0)
        {
          hu_tend_q[k][i][j] += h_west[k][i0][j0] * avg_fac_hu*rhs_u*dt;
        }
        if (dt_avg_e > 0)
        {
          e_tend_adv[k][i][j] += rhs_u*h_west[k][i0][j0]*uu_w[k][i0][j0]*avg_fac_e*dt;
        }
        
        // Add diabatic advection
        if (useWDia)
        {
          uu_p = (k == 0) ? 0 : uu_w[k-1][i0][j0];
          uu_m = (k == Nlay-1) ? 0 : uu_w[k+1][i0][j0];
          du_p = (wdia_u[k][i][j] > 0) ? 0 : uu_p-uu_w[k][i0][j0];
          du_m = (wdia_u[k+1][i][j] > 0) ? uu_w[k][i0][j0]-uu_m : 0;
          rhs_u = - (wdia_u[k][i][j] * du_p + wdia_u[k+1][i][j] * du_m);
          dt_uu_w[k][i][j] += rhs_u / h_west[k][i0][j0];
          if (dt_avg_hu > 0)
          {
            hu_tend_wdia[k][i][j] += avg_fac_hu*rhs_u*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_wdiaKE[k][i][j] += rhs_u*uu_w[k][i0][j0]*avg_fac_e*dt;
          }
        }
        
        // Add Laplacian viscosity
        if (useA2)
        {
          rhs_u = (
                      ( A2_h[i0][j0]*hh_w[k][i0][j0]*DD_T[i0][j0] - A2_h[im1][j0]*hh_w[k][im1][j0]*DD_T[im1][j0] ) / dx
                    + ( A2_q[i0][jp1]*hh_q[i0][jp1]*DD_S[i0][jp1] - A2_q[i0][j0]*hh_q[i0][j0]*DD_S[i0][j0] ) / dy
                  );
          dt_uu_w[k][i][j] += rhs_u / h_west[k][i0][j0];
          if (dt_avg_hu > 0)
          {
            hu_tend_A2[k][i][j] += avg_fac_hu*rhs_u*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_A2[k][i][j] += rhs_u*uu_w[k][i0][j0]*avg_fac_e*dt;
          }
        }

        // Add biharmonic viscosity
        if (useA4)
        {
          rhs_u = - (
                        ( A4sqrt_h[i0][j0]*hh_w[k][i0][j0]*DD4_T[i0][j0] - A4sqrt_h[im1][j0]*hh_w[k][im1][j0]*DD4_T[im1][j0] ) / dx
                      + ( A4sqrt_q[i0][jp1]*hh_q[i0][jp1]*DD4_S[i0][jp1] - A4sqrt_q[i0][j0]*hh_q[i0][j0]*DD4_S[i0][j0] ) / dy
                    );
          dt_uu_w[k][i][j] += rhs_u / h_west[k][i0][j0];
          if (dt_avg_hu > 0)
          {
            hu_tend_A4[k][i][j] += avg_fac_hu*rhs_u*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_A4[k][i][j] += rhs_u*uu_w[k][i0][j0]*avg_fac_e*dt;
          }
        }
        
        // Add wind stress
        if (useWind)
        {
          // Compute tendency
          rhs_u = taux_interp * hFsurf_west[k][i0][j0];
          dt_uu_w[k][i][j] += rhs_u / h_west[k][i0][j0];
          if (dt_avg_hu > 0)
          {
            hu_tend_wind[k][i][j] += avg_fac_hu*rhs_u*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_wind[k][i][j] += rhs_u*uu_w[k][i0][j0]*avg_fac_e*dt;
          }
        }
        
        // Add barotropic forcing
        if (useFbaro)
        {
          // Compute tendency
          rhs_u = Fbaro_x[i][j];
          dt_uu_w[k][i][j] += rhs_u;
          if (dt_avg_hu > 0)
          {
            hu_tend_Fbaro[k][i][j] += h_west[k][i0][j0] * avg_fac_hu*rhs_u*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_Fbaro[k][i][j] += rhs_u*h_west[k][i0][j0]*uu_w[k][i0][j0]*avg_fac_e*dt;
          }
        }
        
        // Add linear bottom drag
        if (rDrag > 0)
        {
          rhs_u = - rDrag * uu_w[k][i0][j0] * hFbot_west[k][i0][j0];
          dt_uu_w[k][i][j] += rhs_u / h_west[k][i0][j0];
          if (dt_avg_hu > 0)
          {
            hu_tend_rDrag[k][i][j] += avg_fac_hu*rhs_u*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_rDrag[k][i][j] += rhs_u*uu_w[k][i0][j0]*avg_fac_e*dt;
          }
        }
        
        // Add linear surface drag
        if (rSurf > 0)
        {
          rhs_u = - rSurf * (uu_w[k][i0][j0] - uLid[i][j]) * hFsurf_west[k][i0][j0];
          dt_uu_w[k][i][j] += rhs_u / h_west[k][i0][j0];
          if (dt_avg_hu > 0)
          {
            hu_tend_rSurf[k][i][j] += avg_fac_hu*rhs_u*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_rSurf[k][i][j] += rhs_u*uu_w[k][i0][j0]*avg_fac_e*dt;
          }
        }
        
        // Add quadratic bottom drag
        if (CdBot > 0)
        {
          rhs_u = - CdBot * 0.5*(uabs_bot[i0][j0]+uabs_bot[im1][j0]) * uu_w[k][i0][j0] * hFbot_west[k][i0][j0];
          dt_uu_w[k][i][j] += rhs_u / h_west[k][i0][j0];
          if (dt_avg_hu > 0)
          {
            hu_tend_CdBot[k][i][j] += avg_fac_hu*rhs_u*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_CdBot[k][i][j] += rhs_u*uu_w[k][i0][j0]*avg_fac_e*dt;
          }
        }
        
        // Add quadratic surface drag
        if (CdSurf > 0)
        {
          rhs_u = - CdSurf * 0.5*(uabs_surf[i0][j0]+uabs_surf[im1][j0]) * (uu_w[k][i0][j0] - uLid[i][j]) * hFsurf_west[k][i0][j0];
          dt_uu_w[k][i][j] += rhs_u / h_west[k][i0][j0];
          if (dt_avg_hu > 0)
          {
            hu_tend_CdSurf[k][i][j] += avg_fac_hu*rhs_u*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_CdSurf[k][i][j] += rhs_u*uu_w[k][i0][j0]*avg_fac_e*dt;
          }
        }
        
        // Acceleration due to lateral buoyancy gradient
        if (useTracer && useBuoyancy)
        {
          rhs_u = 0.5*(hz[i0][j0]+hz[im1][j0]) * db_dx[i0][j0];
          dt_uu_w[k][i][j] += rhs_u / h_west[k][i0][j0];
          if (dt_avg_hu > 0)
          {
            hu_tend_buoy[k][i][j] += avg_fac_hu*rhs_u*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_buoy[k][i][j] += rhs_u*uu_w[k][i0][j0]*avg_fac_e*dt;
          }
        }
        
        // Acceleration due to random forcing
        if (useRandomForcing)
        {
          // Compute tendency
          rhs_u = RFu[k][i][j];
          
          // Apply tendency depending on whether it is applied to the velocity (i.e. du/dt = F + ...)
          // or to the momentum (i.e. d(hu)/dt = F + ...
          if (useThicWeightedRF)
          {
            dt_uu_w[k][i][j] += rhs_u / h_west[k][i0][j0];
            if (dt_avg_hu > 0)
            {
              hu_tend_rand[k][i][j] += avg_fac_hu*rhs_u*dt;
            }
            if (dt_avg_e > 0)
            {
              e_tend_rand[k][i][j] += rhs_u*uu_w[k][i0][j0]*avg_fac_e*dt;
            }
          }
          else
          {
            dt_uu_w[k][i][j] += rhs_u;
            if (dt_avg_hu > 0)
            {
              hu_tend_rand[k][i][j] += h_west[k][i0][j0] * avg_fac_hu*rhs_u*dt;
            }
            if (dt_avg_e > 0)
            {
              e_tend_rand[k][i][j] += rhs_u*h_west[k][i0][j0]*uu_w[k][i0][j0]*avg_fac_e*dt;
            }
          }
        }
        
        // MUST BE LAST CONTRIBUTION TO du/dt
        // Add relaxation terms. Negative values mean no relaxation.
        if (useRelax && (uTime[k][i][j] > 0))
        {
          rhs_u = - (uu_w[k][i0][j0] - uRelax[k][i][j]) / uTime[k][i][j];
          dt_uu_w[k][i][j] += rhs_u;
          
          if (dt_avg_hu > 0)
          {
            hu_tend_relax[k][i][j] += h_west[k][i0][j0] * avg_fac_hu*rhs_u*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_relax[k][i][j] += rhs_u*h_west[k][i0][j0]*uu_w[k][i0][j0]*avg_fac_e*dt;
          }
        }
        
      

        ////////////////////////////////////
        /// Meridional momentum tendency ///
        ////////////////////////////////////
        
        // Initialize
        dt_vv_w[k][i][j] = 0;
        
        // Gradient of Montgomery potential
        rhs_v = - (MM_B[i0][j0]-MM_B[i0][jm1]) / dy;
        dt_vv_w[k][i][j] += rhs_v;
        if (dt_avg_hv > 0)
        {
          hv_tend_gradM[k][i][j] += h_south[k][i0][j0] * avg_fac_hv*rhs_v*dt;
        }
        if (dt_avg_e > 0)
        {
          e_tend_gradM[k][i][j] += rhs_v*h_south[k][i0][j0]*vv_w[k][i0][j0]*avg_fac_e*dt;
        }
        
        // Gradient of kinetic energy
        rhs_v = - (KE_B[i0][j0]-KE_B[i0][jm1]) / dy;
        dt_vv_w[k][i][j] += rhs_v;
        if (dt_avg_hv > 0)
        {
          hv_tend_gradKE[k][i][j] += h_south[k][i0][j0] * avg_fac_hv*rhs_v*dt;
        }
        if (dt_avg_e > 0)
        {
          e_tend_adv[k][i][j] += rhs_v*h_south[k][i0][j0]*vv_w[k][i0][j0]*avg_fac_e*dt;
        }
        
        // Add Coriolis/vorticity advection terms
        switch (momentumScheme)
        {
          case MOMENTUM_AL81:
          case MOMENTUM_HK83:
          {
            rhs_v             = phi_w[i0][jm1]*hvv[i0][jm1]
                              - phi_w[i0][j0]*hvv[i0][jp1]
                              - gamma_w[ip1][j0]*huu[ip1][j0]
                              - delta_w[i0][j0]*huu[i0][j0]
                              - alpha_w[i0][jm1]*huu[i0][jm1]
                              - beta_w[ip1][jm1]*huu[ip1][jm1];
            break;
          }
          case MOMENTUM_S75e:
          {
            rhs_v         = - ( 0.25 * qq[i0][j0] * (huu[i0][j0] + huu[i0][jm1])
                              + 0.25 * qq[ip1][j0] * (huu[ip1][j0] + huu[ip1][jm1]) );
            break;
          }
          default:
          {
            fprintf(stderr,"ERROR: Unknown spatial discretization method\n");
            break;
          }
        }
        dt_vv_w[k][i][j] += rhs_v;
        if (dt_avg_hv > 0)
        {
          hv_tend_q[k][i][j] += h_south[k][i0][j0] * avg_fac_hv*rhs_v*dt;
        }
        if (dt_avg_e > 0)
        {
          e_tend_adv[k][i][j] += rhs_v*h_south[k][i0][j0]*vv_w[k][i0][j0]*avg_fac_e*dt;
        }
        
        // Add diabatic advection
        if (useWDia)
        {
          vv_p = (k == 0) ? 0 : vv_w[k-1][i0][j0];
          vv_m = (k == Nlay-1) ? 0 : vv_w[k+1][i0][j0];
          dv_p = (wdia_v[k][i][j] > 0) ? 0 : vv_p-vv_w[k][i0][j0];
          dv_m = (wdia_v[k+1][i][j] > 0) ? vv_w[k][i0][j0]-vv_m : 0;
          rhs_v = - (wdia_v[k][i][j] * dv_p + wdia_v[k+1][i][j] * dv_m);
          dt_vv_w[k][i][j] += rhs_v / h_south[k][i0][j0];
          if (dt_avg_hv > 0)
          {
            hv_tend_wdia[k][i][j] += avg_fac_hv*rhs_v*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_wdiaKE[k][i][j] += rhs_v*vv_w[k][i0][j0]*avg_fac_e*dt;
          }
        }
        
        // Add Laplacian viscosity
        if (useA2)
        {
          rhs_v = (
                       ( A2_q[ip1][j0]*hh_q[ip1][j0]*DD_S[ip1][j0] - A2_q[i0][j0]*hh_q[i0][j0]*DD_S[i0][j0] ) / dx
                      - ( A2_h[i0][j0]*hh_w[k][i0][j0]*DD_T[i0][j0] - A2_h[i0][jm1]*hh_w[k][i0][jm1]*DD_T[i0][jm1] ) / dy
                   );
          dt_vv_w[k][i][j] += rhs_v / h_south[k][i0][j0];
          if (dt_avg_hv > 0)
          {
            hv_tend_A2[k][i][j] += avg_fac_hv*rhs_v*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_A2[k][i][j] += rhs_v*vv_w[k][i0][j0]*avg_fac_e*dt;
          }
        }
        
        // Add biharmonic viscosity
        if (useA4)
        {
          rhs_v = - (
                        ( A4sqrt_q[ip1][j0]*hh_q[ip1][j0]*DD4_S[ip1][j0] - A4sqrt_q[i0][j0]*hh_q[i0][j0]*DD4_S[i0][j0] ) / dx
                      - ( A4sqrt_h[i0][j0]*hh_w[k][i0][j0]*DD4_T[i0][j0] - A4sqrt_h[i0][jm1]*hh_w[k][i0][jm1]*DD4_T[i0][jm1] ) / dy
                    );
          dt_vv_w[k][i][j] += rhs_v / h_south[k][i0][j0];
          if (dt_avg_hv > 0)
          {
            hv_tend_A4[k][i][j] += avg_fac_hv*rhs_v*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_A4[k][i][j] += rhs_v*vv_w[k][i0][j0]*avg_fac_e*dt;
          }
        }
        
        // Add wind stress
        if (useWind)
        {
          rhs_v = tauy_interp * hFsurf_south[k][i0][j0];
          dt_vv_w[k][i][j] += rhs_v / h_south[k][i0][j0];
          if (dt_avg_hv > 0)
          {
            hv_tend_wind[k][i][j] += avg_fac_hv*rhs_v*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_wind[k][i][j] += rhs_v*vv_w[k][i0][j0]*avg_fac_e*dt;
          }
        }
        
        // Add barotropic forcing
        if (useFbaro)
        {
          rhs_v = Fbaro_y[i][j];
          dt_vv_w[k][i][j] += rhs_v;
          if (dt_avg_hv > 0)
          {
            hv_tend_Fbaro[k][i][j] += h_south[k][i0][j0] * avg_fac_hv*rhs_v*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_Fbaro[k][i][j] += rhs_v*h_south[k][i0][j0]*vv_w[k][i0][j0]*avg_fac_e*dt;
          }
        }
        
        // Add linear bottom drag
        if (rDrag > 0)
        {
          rhs_v = - rDrag * vv_w[k][i0][j0] * hFbot_south[k][i0][j0];
          dt_vv_w[k][i][j] += rhs_v / h_south[k][i0][j0];
          if (dt_avg_hv > 0)
          {
            hv_tend_rDrag[k][i][j] += avg_fac_hv*rhs_v*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_rDrag[k][i][j] += rhs_v*vv_w[k][i0][j0]*avg_fac_e*dt;
          }
        }
        
        // Add linear surface drag
        if (rSurf > 0)
        {
          rhs_v = - rSurf * (vv_w[k][i0][j0] - vLid[i][j]) * hFsurf_south[k][i0][j0];
          dt_vv_w[k][i][j] += rhs_v / h_south[k][i0][j0];
          if (dt_avg_hv > 0)
          {
            hv_tend_rSurf[k][i][j] += avg_fac_hv*rhs_v*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_rSurf[k][i][j] += rhs_v*vv_w[k][i0][j0]*avg_fac_e*dt;
          }
        }
        
        // Add quadratic bottom drag
        if (CdBot > 0)
        {
          rhs_v = - CdBot * 0.5*(uabs_bot[i0][j0]+uabs_bot[i0][jm1]) * vv_w[k][i0][j0] * hFbot_south[k][i0][j0];
          dt_vv_w[k][i][j] += rhs_v / h_south[k][i0][j0];
          if (dt_avg_hv > 0)
          {
            hv_tend_CdBot[k][i][j] += avg_fac_hv*rhs_v*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_CdBot[k][i][j] += rhs_v*vv_w[k][i0][j0]*avg_fac_e*dt;
          }
        }
        
        // Add quadratic surface drag
        if (CdSurf > 0)
        {
          rhs_v = - CdSurf * 0.5*(uabs_surf[i0][j0]+uabs_surf[i0][jm1]) * (vv_w[k][i0][j0] - vLid[i][j]) * hFsurf_south[k][i0][j0];
          dt_vv_w[k][i][j] += rhs_v / h_south[k][i0][j0];
          if (dt_avg_hv > 0)
          {
            hv_tend_CdSurf[k][i][j] += avg_fac_hv*rhs_v*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_CdSurf[k][i][j] += rhs_v*vv_w[k][i0][j0]*avg_fac_e*dt;
          }
        }
        
        // Acceleration due to lateral buoyancy gradient
        if (useTracer && useBuoyancy)
        {
          rhs_v = 0.5*(hz[i0][j0]+hz[i0][jm1]) * db_dy[i0][j0];
          dt_vv_w[k][i][j] += rhs_v / h_south[k][i0][j0];
          if (dt_avg_hv > 0)
          {
            hv_tend_buoy[k][i][j] += avg_fac_hv*rhs_v*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_buoy[k][i][j] += rhs_v*vv_w[k][i0][j0]*avg_fac_e*dt;
          }
        }
        
        // Acceleration due to random forcing
        if (useRandomForcing)
        {
          // Compute tendency
          rhs_v = RFv[k][i][j];
          
          // Apply tendency depending on whether it is applied to the velocity (i.e. dv/dt = F + ...)
          // or to the momentum (i.e. d(hv)/dt = F + ...
          if (useThicWeightedRF)
          {
            dt_vv_w[k][i][j] += rhs_v / h_south[k][i0][j0];
            if (dt_avg_hv > 0)
            {
              hv_tend_rand[k][i][j] += avg_fac_hv*rhs_v*dt;
            }
            if (dt_avg_e > 0)
            {
              e_tend_rand[k][i][j] += rhs_v*vv_w[k][i0][j0]*avg_fac_e*dt;
            }
          }
          else
          {
            dt_vv_w[k][i][j] += rhs_v;
            if (dt_avg_hv > 0)
            {
              hv_tend_rand[k][i][j] += h_south[k][i0][j0] * avg_fac_hv*rhs_v*dt;
            }
            if (dt_avg_e > 0)
            {
              e_tend_rand[k][i][j] += rhs_v*h_south[k][i0][j0]*vv_w[k][i0][j0]*avg_fac_e*dt;
            }
          }
        }
        
        // MUST BE LAST CONTRIBUTION TO dv/dt
        // Add relaxation terms. Negative values mean no relaxation.
        if (useRelax && (vTime[k][i][j] > 0))
        {
          rhs_v -= (vv_w[k][i0][j0] - vRelax[k][i][j]) / vTime[k][i][j];
          dt_vv_w[k][i][j] += rhs_v;
          
          if (dt_avg_hu > 0)
          {
            hv_tend_relax[k][i][j] += h_south[k][i0][j0] * avg_fac_hv*rhs_v*dt;
          }
          if (dt_avg_e > 0)
          {
            e_tend_relax[k][i][j] += rhs_v*h_south[k][i0][j0]*vv_w[k][i0][j0]*avg_fac_e*dt;
          }
        }
        
        
        
        
        
        // Additional energy budget terms
        if (dt_avg_e > 0)
        {
          // Pressure flux
          pres_0 = pp[i0][j0] - geff[k]*0.5*(eta_w[k][i0][j0]+eta_w[k+1][i0][j0]);
          pres_im1 = pp[im1][j0] - geff[k]*0.5*(eta_w[k][im1][j0]+eta_w[k+1][im1][j0]);
          pres_jm1 = pp[i0][jm1] - geff[k]*0.5*(eta_w[k][i0][jm1]+eta_w[k+1][i0][jm1]);
          e_flux_uP[k][i][j] += h_west[k][i0][j0] * uu_w[k][i0][j0] * 0.5*(pres_0+pres_im1) * avg_fac_e*dt;
          e_flux_vP[k][i][j] += h_south[k][i0][j0] * vv_w[k][i0][j0]  * 0.5*(pres_0+pres_jm1) * avg_fac_e*dt;
          
          // PE flux
          e_flux_uPE[k][i][j] += h_west[k][i0][j0] * uu_w[k][i0][j0] * 0.5*((MM_B[i0][j0]-pres_0)+(MM_B[im1][j0]-pres_im1)) * avg_fac_e*dt;
          e_flux_vPE[k][i][j] += h_south[k][i0][j0] * vv_w[k][i0][j0]  * 0.5*((MM_B[i0][j0]-pres_0)+(MM_B[i0][jm1]-pres_jm1)) * avg_fac_e*dt;
          //e_flux_uPE[k][i][j] += h_west[k][i0][j0] * uu_w[k][i0][j0] * 0.5*(MM_B[i0][j0]+MM_B[im1][j0]) * avg_fac_e*dt;
          //e_flux_vPE[k][i][j] += h_south[k][i0][j0] * vv_w[k][i0][j0]  * 0.5*(MM_B[i0][j0]+MM_B[i0][jm1]) * avg_fac_e*dt;
          
          // KE flux
          e_flux_uKE[k][i][j] += h_west[k][i0][j0] * uu_w[k][i0][j0] * 0.5*(KE_B[i0][j0]+KE_B[im1][j0]) * avg_fac_e*dt;
          e_flux_vKE[k][i][j] += h_south[k][i0][j0] * vv_w[k][i0][j0]  * 0.5*(KE_B[i0][j0]+KE_B[i0][jm1]) * avg_fac_e*dt;
        }
        
        
        
        
        
        
        
        
        
        ///////////////////////
        /// Tracer tendency ///
        ///////////////////////
        
        if (useTracer)
        {
          // Initialize
          dt_bb_w[k][i][j] = 0;
          
          // Advection
          if (tracerScheme == TRACER_AL81)
          {
            dt_bb_w[k][i][j] = dt_bb_w[k][i][j] - 0.5 * (udb_dx[ip1][j0] + udb_dx[i0][j0] + vdb_dy[i0][jp1] + vdb_dy[i0][j0]);
          }
          else
          {
            // Calculate advective terms via flux form for total tracer, i,e. as ( div (hub) - b div(hu) ) / h
            dt_bb_w[k][i][j] = dt_bb_w[k][i][j] - ( (hub[ip1][j0]-hub[i0][j0]) / dx + (hvb[i0][jp1]-hvb[i0][j0]) / dy ) / hh_w[k][i0][j0];
            dt_bb_w[k][i][j] = dt_bb_w[k][i][j] + bb_w[k][i0][j0] * ( (huu[ip1][j0]-huu[i0][j0]) / dx + (hvv[i0][jp1]-hvv[i0][j0]) / dy ) / hh_w[k][i0][j0];
            
            // N.B. the above is equivalent to:
            // dt_bb_w[k][i][j] = dt_bb_w[k][i][j] - (
            //                                           huu[ip1][j0]*(b_west[ip1][j0]-bb_w[k][i0][j0])/dx
            //                                        +  huu[i0][j0]*(bb_w[k][i0][j0]-b_west[i0][j0])/dx
            //                                        +  hvv[i0][jp1]*(b_south[ip1][j0]-bb_w[k][i0][j0])/dy
            //                                        +  hvv[i0][j0]*(bb_w[k][i0][j0]-b_south[i0][j0])/dy
            //                                       )
            //                                       / hh_w[k][i0][j0];
          }
          
          // Add diabatic advection
          if (useWDia)
          {
            bb_p = (k == 0) ? 0 : bb_w[k-1][i0][j0];
            bb_m = (k == Nlay-1) ? 0 : bb_w[k+1][i0][j0];
            db_p = (wdia[k][i][j] > 0) ? 0 : bb_p-bb_w[k][i0][j0];
            db_m = (wdia[k+1][i][j] > 0) ? bb_w[k][i0][j0]-bb_m : 0;
            dt_bb_w[k][i][j] = dt_bb_w[k][i][j] - (wdia[k][i][j] * db_p + wdia[k+1][i][j] * db_m) / hh_w[k][i0][j0];
          }
          
          // Add Laplacian diffusion
          if (K2 > 0)
          {
            dt_bb_w[k][i][j] = dt_bb_w[k][i][j] + K2 *
            (
                ( h_west[k][ip1][j0]*(bb_w[k][ip1][j0]-bb_w[k][i0][j0]) - h_west[k][i0][j0]*(bb_w[k][i0][j0]-bb_w[k][im1][j0]) ) / dxsq
              + ( h_south[k][i0][jp1]*(bb_w[k][i0][jp1]-bb_w[k][i0][j0]) - h_south[k][i0][j0]*(bb_w[k][i0][j0]-bb_w[k][i0][jm1]) ) / dysq
            )
            / hh_w[k][i0][j0];
          }
          
          // Add biharmonic diffusion
          if (K4 > 0)
          {
            dt_bb_w[k][i][j] = dt_bb_w[k][i][j] - K4 *
            (
                ( h_west[k][ip1][j0]*(BB4[ip1][j0]-BB4[i0][j0]) - h_west[k][i0][j0]*(BB4[i0][j0]-BB4[im1][j0]) ) / dxsq
              + ( h_south[k][i0][jp1]*(BB4[i0][jp1]-BB4[i0][j0]) - h_south[k][i0][j0]*(BB4[i0][j0]-BB4[i0][jm1]) ) / dysq
            )
            / hh_w[k][i0][j0];
          }
          
          // Add relaxation terms. Negative values mean no relaxation.
          if (useRelax && (bTime[k][i][j] > 0))
          {
            dt_bb_w[k][i][j] -= (bb_w[k][i0][j0] - bRelax[k][i][j]) / bTime[k][i][j];
          }
          
        } // end if (useTracer)
        
        
      } // end j-loop
      
    } // end i-loop
    
    
    
    // Cross-boundary velocities fixed and equal to zero.
    // Only necessary for southern/western walls, as only
    // these gridpoints are stored in the time-stepped matrices.
    if (useWallNS)
    {
      for (i = 0; i < Nx; i ++)
      {
        dt_vv_w[k][i][0] = 0;
      }
    }
    if (useWallEW)
    {
      for (j = 0; j < Ny; j ++)
      {
        dt_uu_w[k][0][j] = 0;
      }
    }
    
  } // end k-loop

} // end tderiv










/**
 *
 * calcEZ
 *
 * Calculates total energy over all layers and total potential enstrophy in each layer.
 *
 */
void calcEZ (const int t, real *** uu, real *** vv, real *** hh, real *** bb, real * KE, real * PE, real * E, real * Z)
{
  // Looping variables
  int i,j,k,ip1,im1,jm1,jp1;
  
  // To store intermediate quantities
  real h_q = 0;
  real pv = 0;
  real gtot = 0;
  
  // Repurpose the eta_w matrix to calculate eta here.
  // N.B. Here we treat eta as an Nlay x Nx x Ny matrix
  // instead of an Nlay * (Nx+2*Ng) x (Ny+2*Ng) matrix.
  real *** eta = eta_w;
  
  // Calculate the potential enstrophy in each layer
  for (k = 0; k < Nlay; k ++)
  {
    
#pragma parallel
    
    // NOTE: We only set vorticities on an Nx x Ny grid, though the PV grid may be
    // Nx x Ny+1 or Nx+1 x Ny or Nx+1 x Ny+1 depending on the periodicity of the domain.
    // However, if the domain is non-periodic in either direction then the additional
    // vorticities would just be zero anyway due to free-slip BCs.
    for (i = 0; i < Nx; i ++)
    {
      im1 = (i+Nx-1) % Nx;
      
      for (j = 0; j < Ny; j ++)
      {
        jm1 = (j+Ny-1) % Ny;
        
        if (useWallNS && (j == 0))
        {
          h_q = (hh[k][im1][j]+hh[k][i][j]) / 2;
          zeta[i][j] = 0;
        }
        else if (useWallEW && (i == 0))
        {
          h_q = (hh[k][i][jm1]+hh[k][i][j]) / 2;
          zeta[i][j] = 0;
        }
        else
        {
          h_q = (hh[k][i][j]+hh[k][im1][jm1]+hh[k][im1][j]+hh[k][i][jm1]) / 4;
          zeta[i][j] = (uu[k][i][jm1]-uu[k][i][j]) / dy + (vv[k][i][j]-vv[k][im1][j]) / dx;
        }
        pv = (2*Omega_z[i][j] + zeta[i][j]) / h_q;
        Zdens[i][j] = dx*dy * 0.5 * h_q * SQUARE(pv);
      }
    }
    
    // Sum potential enstrophy densities within this layer
    Z[k] = kahanSum (*Zdens,Nx*Ny);
    
  } // End k-loop
  
  // Calculate layer interface heights - needed for the energy
  calcEta(hhb,hh,eta,Nlay,Nx,Ny);
  
  // Calculate the energy
  for (k = 0; k < Nlay; k ++)
  {
    
#pragma parallel
    
    // In calculating the total energy we always assume periodicity in both directions. KE is defined on each h-point,
    // and has contributions from u^2_west, u^2_east, v^2_north and v^2_south on that grid cell. If a wall is present
    // on one edge of the cell then the normal velocity is guaranteed to be zero, and so makes no contribution to the KE.
    for (i = 0; i < Nx; i ++)
    {
      ip1 = (i+Nx+1) % Nx;
      im1 = (i+Nx-1) % Nx;
      
      for (j = 0; j < Ny; j ++)
      {
        
        jp1 = (j+Ny+1) % Ny;
        jm1 = (j+Ny-1) % Ny;
        
        switch (thicknessScheme)
        {
          case THICKNESS_AL81:
          case THICKNESS_UP3:
          case THICKNESS_KT00:
          {
            KEdens[k][i][j] = 0.125 * dx*dy * ( (hh[k][i][j]+hh[k][ip1][j]) * SQUARE(uu[k][ip1][j])
                                 + (hh[k][im1][j]+hh[k][i][j]) * SQUARE(uu[k][i][j])
                                 + (hh[k][i][j]+hh[k][i][jp1]) * SQUARE(vv[k][i][jp1])
                                 + (hh[k][i][jm1]+hh[k][i][j]) * SQUARE(vv[k][i][j]) );
            break;
          }
            
          // The HK83 stencil extends beyond the walls of the domain, so we require additional checks.
          case THICKNESS_HK83:
          {
            
            KEdens[k][i][j] = 0;
            
            if (useWallNS && j==0)
            {
              KEdens[k][i][j] +=
                  (hh[k][i][j]+hh[k][ip1][j]) * ( 1.25*SQUARE(uu[k][ip1][j]) + 0.25*SQUARE(uu[k][ip1][jp1])  )
                + (hh[k][im1][j]+hh[k][i][j]) * ( 1.25*SQUARE(uu[k][i][j])  + 0.25*SQUARE(uu[k][i][jp1]) );
            }
            else if (useWallNS && j==Ny-1)
            {
              KEdens[k][i][j] +=
                  (hh[k][i][j]+hh[k][ip1][j]) * ( 1.25*SQUARE(uu[k][ip1][j]) + 0.25*SQUARE(uu[k][ip1][jm1]) )
                + (hh[k][im1][j]+hh[k][i][j]) * ( 1.25*SQUARE(uu[k][i][j]) + 0.25*SQUARE(uu[k][i][jm1]) );
            }
            else
            {
              KEdens[k][i][j] +=
                  (hh[k][i][j]+hh[k][ip1][j]) * ( SQUARE(uu[k][ip1][j]) + 0.25*SQUARE(uu[k][ip1][jp1]) + 0.25*SQUARE(uu[k][ip1][jm1]) )
                + (hh[k][im1][j]+hh[k][i][j]) * ( SQUARE(uu[k][i][j])  + 0.25*SQUARE(uu[k][i][jp1]) + 0.25*SQUARE(uu[k][i][jm1]) );
            }
            
            
            if (useWallEW && i==0)
            {
              KEdens[k][i][j] +=
                  (hh[k][i][j]+hh[k][i][jp1]) * ( 1.25*SQUARE(vv[k][i][jp1]) + 0.25*SQUARE(vv[k][ip1][jp1]) )
                + (hh[k][i][jm1]+hh[k][i][j]) * ( 1.25*SQUARE(vv[k][i][j]) + 0.25*SQUARE(vv[k][ip1][j]) );
            }
            else if (useWallEW && i==Nx-1)
            {
              KEdens[k][i][j] +=
                  (hh[k][i][j]+hh[k][i][jp1]) * ( 1.25*SQUARE(vv[k][i][jp1]) + 0.25*SQUARE(vv[k][im1][jp1]) )
                + (hh[k][i][jm1]+hh[k][i][j]) * ( 1.25*SQUARE(vv[k][i][j]) + 0.25*SQUARE(vv[k][im1][j]) );
            }
            else
            {
              KEdens[k][i][j] +=
                  (hh[k][i][j]+hh[k][i][jp1]) * ( SQUARE(vv[k][i][jp1]) + 0.25*SQUARE(vv[k][ip1][jp1]) + 0.25*SQUARE(vv[k][im1][jp1]) )
                + (hh[k][i][jm1]+hh[k][i][j]) * ( SQUARE(vv[k][i][j]) + 0.25*SQUARE(vv[k][ip1][j]) + 0.25*SQUARE(vv[k][im1][j]) );
            }
            
            KEdens[k][i][j] *= (dx*dy/12.0);
            
            break;
          }
            
          default:
          {
            fprintf(stderr,"ERROR: Unknown spatial discretization method\n");
            break;
          }
        }
        
        // If we have an active buoyancy variable then modify the effective gravity accordingly
        if (useTracer && useBuoyancy)
        {
          gtot = geff[k] - bb[k][i][j];
        }
        else
        {
          gtot = geff[k];
        }
        
        // GPE
        PEdens[k][i][j] = dx*dy * gtot * hh[k][i][j] * (eta[k][i][j] - 0.5*hh[k][i][j]);
        
        // Add Salmon term
        PEdens[k][i][j] += dx*dy * gtot * POW4(h0) / POW2(hh[k][i][j]) / 6;
        
      }
    }
  }
  
  // Sum all energy densities to get total energy
  *KE = kahanSum(*(*KEdens),Nlay*Nx*Ny);
  *PE = kahanSum(*(*PEdens),Nlay*Nx*Ny);
  *E = *KE + *PE;
}












/**
 *
 * correctThickness
 *
 * Enforces layer thicknesses that sum to the total water column thickness.
 *
 */
void correctThickness (real *** hh)
{
  int i, j, k;
  real h_err, h_tot;
  
  for (i = 0; i < Nx; i ++)
  {
    for (j = 0; j < Ny; j ++)
    {
      // Calculate sum of layer thicknesses
      h_tot = 0;
      for (k = 0; k < Nlay; k ++)
      {
        h_tot += hh[k][i][j];
      }
      
      // Error vs true water column thickness
      h_err = h_tot - Hc[i][j];
      
      // Apply weighted correction
      for (k = 0; k < Nlay; k ++)
      {
        hh[k][i][j] *= (1 - h_err / h_tot);
      }
    }
  }
}













/**
 *
 * constructOutputName
 *
 * Constructs an output file path from the destination directory (outdir), 
 * variable id (varid, defined at the top of this file), isopycnal layer (k),
 * and output index (n). The resulting file path is stored in outfile. If k<0
 * then no isopycnal layer index is specified in the output file.
 *
 */
void constructOutputName (char * outdir, int varid, int k, uint n, char * outfile)
{
  char nstr[MAX_PARAMETER_FILENAME_LENGTH];
  char kstr[MAX_PARAMETER_FILENAME_LENGTH];
  
  // Create a string for the current iteration number
  sprintf(nstr,"%d",n);
  
  // Create a string for the current isopycnal layer
  if (k >= 0)
  {
    sprintf(kstr,"%d",k);
  }
  else
  {
    sprintf(kstr,"");
  }
  
  // Initialize output file path
  strcpy(outfile,outdir);
  strcat(outfile,"/");
  
  // Append variable name
  switch (varid)
  {
      
    // Instantaneous output
    case VARID_U:
    {
      strcat(outfile,OUTN_U);
      break;
    }
    case VARID_V:
    {
      strcat(outfile,OUTN_V);
      break;
    }
    case VARID_H:
    {
      strcat(outfile,OUTN_H);
      break;
    }
    case VARID_B:
    {
      strcat(outfile,OUTN_B);
      break;
    }
    case VARID_PI:
    {
      strcat(outfile,OUTN_PI);
      break;
    }
    case VARID_W:
    {
      strcat(outfile,OUTN_W);
      break;
    }
      
    // Time-averaged quantities
    case VARID_U_AVG:
    {
      strcat(outfile,OUTN_U_AVG);
      break;
    }
    case VARID_V_AVG:
    {
      strcat(outfile,OUTN_V_AVG);
      break;
    }
    case VARID_H_AVG:
    {
      strcat(outfile,OUTN_H_AVG);
      break;
    }
    case VARID_B_AVG:
    {
      strcat(outfile,OUTN_B_AVG);
      break;
    }
    case VARID_M_AVG:
    {
      strcat(outfile,OUTN_M_AVG);
      break;
    }
    case VARID_PI_AVG:
    {
      strcat(outfile,OUTN_PI_AVG);
      break;
    }
    case VARID_W_AVG:
    {
      strcat(outfile,OUTN_W_AVG);
      break;
    }
    case VARID_HU_AVG:
    {
      strcat(outfile,OUTN_HU_AVG);
      break;
    }
    case VARID_HV_AVG:
    {
      strcat(outfile,OUTN_HV_AVG);
      break;
    }
    case VARID_HUU_AVG:
    {
      strcat(outfile,OUTN_HUU_AVG);
      break;
    }
    case VARID_HVV_AVG:
    {
      strcat(outfile,OUTN_HVV_AVG);
      break;
    }
    case VARID_HUV_AVG:
    {
      strcat(outfile,OUTN_HUV_AVG);
      break;
    }
     
    // U-momentum output
    case VARID_UMOM_Q:
    {
      strcat(outfile,OUTN_UMOM_Q);
      break;
    }
    case VARID_UMOM_GRADM:
    {
      strcat(outfile,OUTN_UMOM_GRADM);
      break;
    }
    case VARID_UMOM_GRADKE:
    {
      strcat(outfile,OUTN_UMOM_GRADKE);
      break;
    }
    case VARID_UMOM_DHDT:
    {
      strcat(outfile,OUTN_UMOM_DHDT);
      break;
    }
    case VARID_UMOM_A2:
    {
      strcat(outfile,OUTN_UMOM_A2);
      break;
    }
    case VARID_UMOM_A4:
    {
      strcat(outfile,OUTN_UMOM_A4);
      break;
    }
    case VARID_UMOM_RDRAG:
    {
      strcat(outfile,OUTN_UMOM_RDRAG);
      break;
    }
    case VARID_UMOM_RSURF:
    {
      strcat(outfile,OUTN_UMOM_RSURF);
      break;
    }
    case VARID_UMOM_CDBOT:
    {
      strcat(outfile,OUTN_UMOM_CDBOT);
      break;
    }
    case VARID_UMOM_CDSURF:
    {
      strcat(outfile,OUTN_UMOM_CDSURF);
      break;
    }
    case VARID_UMOM_WIND:
    {
      strcat(outfile,OUTN_UMOM_WIND);
      break;
    }
    case VARID_UMOM_BUOY:
    {
      strcat(outfile,OUTN_UMOM_BUOY);
      break;
    }
    case VARID_UMOM_RELAX:
    {
      strcat(outfile,OUTN_UMOM_RELAX);
      break;
    }
    case VARID_UMOM_WDIA:
    {
      strcat(outfile,OUTN_UMOM_WDIA);
      break;
    }
    case VARID_UMOM_RAND:
    {
      strcat(outfile,OUTN_UMOM_RAND);
      break;
    }
    case VARID_UMOM_FBARO:
    {
      strcat(outfile,OUTN_UMOM_FBARO);
      break;
    }
      
    // V-momentum output
    case VARID_VMOM_Q:
    {
      strcat(outfile,OUTN_VMOM_Q);
      break;
    }
    case VARID_VMOM_GRADM:
    {
      strcat(outfile,OUTN_VMOM_GRADM);
      break;
    }
    case VARID_VMOM_GRADKE:
    {
      strcat(outfile,OUTN_VMOM_GRADKE);
      break;
    }
    case VARID_VMOM_DHDT:
    {
      strcat(outfile,OUTN_VMOM_DHDT);
      break;
    }
    case VARID_VMOM_A2:
    {
      strcat(outfile,OUTN_VMOM_A2);
      break;
    }
    case VARID_VMOM_A4:
    {
      strcat(outfile,OUTN_VMOM_A4);
      break;
    }
    case VARID_VMOM_RDRAG:
    {
      strcat(outfile,OUTN_VMOM_RDRAG);
      break;
    }
    case VARID_VMOM_RSURF:
    {
      strcat(outfile,OUTN_VMOM_RSURF);
      break;
    }
    case VARID_VMOM_CDBOT:
    {
      strcat(outfile,OUTN_VMOM_CDBOT);
      break;
    }
    case VARID_VMOM_CDSURF:
    {
      strcat(outfile,OUTN_VMOM_CDSURF);
      break;
    }
    case VARID_VMOM_WIND:
    {
      strcat(outfile,OUTN_VMOM_WIND);
      break;
    }
    case VARID_VMOM_BUOY:
    {
      strcat(outfile,OUTN_VMOM_BUOY);
      break;
    }
    case VARID_VMOM_RELAX:
    {
      strcat(outfile,OUTN_VMOM_RELAX);
      break;
    }
    case VARID_VMOM_WDIA:
    {
      strcat(outfile,OUTN_VMOM_WDIA);
      break;
    }
    case VARID_VMOM_RAND:
    {
      strcat(outfile,OUTN_VMOM_RAND);
      break;
    }
    case VARID_VMOM_FBARO:
    {
      strcat(outfile,OUTN_VMOM_FBARO);
      break;
    }
      
    // Thickness equation output
    case VARID_THIC_ADV:
    {
      strcat(outfile,OUTN_THIC_ADV);
      break;
    }
    case VARID_THIC_RELAX:
    {
      strcat(outfile,OUTN_THIC_RELAX);
      break;
    }
      
    // Energy equation diagnostics
    case VARID_ENERGY_UPFLUX:
    {
      strcat(outfile,OUTN_ENERGY_UPFLUX);
      break;
    }
    case VARID_ENERGY_UKEFLUX:
    {
      strcat(outfile,OUTN_ENERGY_UKEFLUX);
      break;
    }
    case VARID_ENERGY_UPEFLUX:
    {
      strcat(outfile,OUTN_ENERGY_UPEFLUX);
      break;
    }
    case VARID_ENERGY_VPFLUX:
    {
      strcat(outfile,OUTN_ENERGY_VPFLUX);
      break;
    }
    case VARID_ENERGY_VKEFLUX:
    {
      strcat(outfile,OUTN_ENERGY_VKEFLUX);
      break;
    }
    case VARID_ENERGY_VPEFLUX:
    {
      strcat(outfile,OUTN_ENERGY_VPEFLUX);
      break;
    }
    case VARID_ENERGY_ADV:
    {
      strcat(outfile,OUTN_ENERGY_ADV);
      break;
    }
    case VARID_ENERGY_GRADM:
    {
      strcat(outfile,OUTN_ENERGY_GRADM);
      break;
    }
    case VARID_ENERGY_WIND:
    {
      strcat(outfile,OUTN_ENERGY_WIND);
      break;
    }
    case VARID_ENERGY_RDRAG:
    {
      strcat(outfile,OUTN_ENERGY_RDRAG);
      break;
    }
    case VARID_ENERGY_RSURF:
    {
      strcat(outfile,OUTN_ENERGY_RSURF);
      break;
    }
    case VARID_ENERGY_CDBOT:
    {
      strcat(outfile,OUTN_ENERGY_CDBOT);
      break;
    }
    case VARID_ENERGY_CDSURF:
    {
      strcat(outfile,OUTN_ENERGY_CDSURF);
      break;
    }
    case VARID_ENERGY_A2:
    {
      strcat(outfile,OUTN_ENERGY_A2);
      break;
    }
    case VARID_ENERGY_A4:
    {
      strcat(outfile,OUTN_ENERGY_A4);
      break;
    }
    case VARID_ENERGY_WDIAPE:
    {
      strcat(outfile,OUTN_ENERGY_WDIAPE);
      break;
    }
    case VARID_ENERGY_WDIAKE:
    {
      strcat(outfile,OUTN_ENERGY_WDIAKE);
      break;
    }
    case VARID_ENERGY_FBARO:
    {
      strcat(outfile,OUTN_ENERGY_FBARO);
      break;
    }
    case VARID_ENERGY_BUOY:
    {
      strcat(outfile,OUTN_ENERGY_BUOY);
      break;
    }
    case VARID_ENERGY_RAND:
    {
      strcat(outfile,OUTN_ENERGY_RAND);
      break;
    }
    case VARID_ENERGY_RELAX:
    {
      strcat(outfile,OUTN_ENERGY_RELAX);
      break;
    }
      
    default:
    {
      fprintf(stderr,"ERROR: Unknown variable ID (%d) specified in constructOutputName\n",varid);
      break;
    }
  }
  
  // Append layer and output indices
  strcat(outfile,kstr);
  strcat(outfile,"_n=");
  strcat(outfile,nstr);
  strcat(outfile,".dat");
}














/**
 *
 * writeOutputFile
 *
 * Convenience function to write a 2D matrix to a specified output file.
 *
 */
mybool writeOutputFile (char * outfile, real ** mat, uint m, uint n)
{
  FILE * outfd = NULL;
  
  outfd = fopen(outfile,"w");
  if (outfd == NULL)
  {
    fprintf(stderr,"ERROR: Could not open output file: %s\n",outfile);
    return false;
  }
  printMatrix(outfd,mat,m,n);
  fclose(outfd);
  return true;
}











/**
 *
 * writeModelState
 *
 * Writes model variables to output files. Returns false if there was a write error,
 * or true if the write was successful.
 *
 */
mybool writeModelState (const int t, const int n, real *** uu, real *** vv, real *** hh, real *** bb, real ** pi, real *** wdia, char * outdir)
{
  uint k = 0;
  char outfile[MAX_PARAMETER_FILENAME_LENGTH];
  
  // Loop over layers
  for (k = 0; k < Nlay; k ++)
  {
    // Write iteration data for u
    constructOutputName(outdir,VARID_U,k,n,outfile);
    if (!writeOutputFile(outfile,uu[k],Nx,Ny)) return false;
    
    // Write iteration data for v
    constructOutputName(outdir,VARID_V,k,n,outfile);
    if (!writeOutputFile(outfile,vv[k],Nx,Ny)) return false;
    
    // Write iteration data for h
    constructOutputName(outdir,VARID_H,k,n,outfile);
    if (!writeOutputFile(outfile,hh[k],Nx,Ny)) return false;
    
    // Write iteration data for b
    if (useTracer)
    {
      constructOutputName(outdir,VARID_B,k,n,outfile);
      if (!writeOutputFile(outfile,bb[k],Nx,Ny)) return false;
    }
    
    // Write iteration data for wdia
    constructOutputName(outdir,VARID_W,k,n,outfile);
    if (!writeOutputFile(outfile,wdia[k],Nx,Ny)) return false;
  }
  
  // Write iteration data for wdia at the sea floor (should be zero in most circumstances,
  // but it's conceivable that a flux through the sea floor would be desireable for some applications)
  k = Nlay;
  constructOutputName(outdir,VARID_W,k,n,outfile);
  if (!writeOutputFile(outfile,wdia[k],Nx,Ny)) return false;

  // Write iteration data for surface pressure
  if (useRL)
  {
    constructOutputName (outdir,VARID_PI,-1,n,outfile);
    if (!writeOutputFile (outfile,pi,Nx,Ny)) return false;
  }

  return true;
}












/**
 *
 * writeAverageState
 *
 * Writes time-averaged model variables to output files. Returns false if there was a write error,
 * or true if the write was successful.
 *
 */
mybool writeAverageState (const int t, const int n, real *** uu, real *** vv, real *** hh, real *** MM, real *** bb, real ** pi, real *** wdia, real *** hu, real *** hv, real *** huu, real *** hvv, real *** huv, char * outdir)
{
  uint k = 0;
  char outfile[MAX_PARAMETER_FILENAME_LENGTH];
  
  // Loop over layers
  for (k = 0; k < Nlay; k ++)
  {
    // Write average data for u
    constructOutputName(outdir,VARID_U_AVG,k,n,outfile);
    if (!writeOutputFile(outfile,uu[k],Nx,Ny)) return false;
    
    // Write average data for v
    constructOutputName(outdir,VARID_V_AVG,k,n,outfile);
    if (!writeOutputFile(outfile,vv[k],Nx,Ny)) return false;
    
    // Write average data for h
    constructOutputName(outdir,VARID_H_AVG,k,n,outfile);
    if (!writeOutputFile(outfile,hh[k],Nx,Ny)) return false;
    
    // Write average data for M
    constructOutputName(outdir,VARID_M_AVG,k,n,outfile);
    if (!writeOutputFile(outfile,MM[k],Nx,Ny)) return false;
    
    // Write average data for wdia
    constructOutputName(outdir,VARID_W_AVG,k,n,outfile);
    if (!writeOutputFile(outfile,wdia[k],Nx,Ny)) return false;
    
    // Write average data for b
    if (useTracer)
    {
      constructOutputName(outdir,VARID_B_AVG,k,n,outfile);
      if (!writeOutputFile(outfile,bb[k],Nx,Ny)) return false;
    }
    
    // Write average data for hu
    constructOutputName(outdir,VARID_HU_AVG,k,n,outfile);
    if (!writeOutputFile(outfile,hu[k],Nx,Ny)) return false;
    
    // Write average data for hv
    constructOutputName(outdir,VARID_HV_AVG,k,n,outfile);
    if (!writeOutputFile(outfile,hv[k],Nx,Ny)) return false;
    
    // Write average data for huu
    constructOutputName(outdir,VARID_HUU_AVG,k,n,outfile);
    if (!writeOutputFile(outfile,huu[k],Nx,Ny)) return false;
    
    // Write average data for hvv
    constructOutputName(outdir,VARID_HVV_AVG,k,n,outfile);
    if (!writeOutputFile(outfile,hvv[k],Nx,Ny)) return false;
    
    // Write average data for huv
    constructOutputName(outdir,VARID_HUV_AVG,k,n,outfile);
    if (!writeOutputFile(outfile,huv[k],Nx,Ny)) return false;
  }
  
  // Write average data for wdia at the sea floor
  k = Nlay;
  constructOutputName(outdir,VARID_W_AVG,k,n,outfile);
  if (!writeOutputFile(outfile,wdia[k],Nx,Ny)) return false;
  
  // Write iteration data for surface pressure
  if (useRL)
  {
    constructOutputName (outdir,VARID_PI_AVG,-1,n,outfile);
    if (!writeOutputFile (outfile,pi,Nx,Ny)) return false;
  }
  
  return true;
}














/**
 *
 * writeUMomentumAverages
 *
 * Calculates averages of terms in U-momentum equation, writes averages to output
 * files, and resets averaging buffers and parameters.
 *
 */
mybool writeUMomentumAverages (uint n, char * outdir)
{
  uint i,j,k;
  char outfile[MAX_PARAMETER_FILENAME_LENGTH];
  
#pragma parallel
  
  // Divide by averaging step length to compute averages
  for (i = 0; i < Nx; i ++)
  {
    for (j = 0; j < Ny; j ++)
    {
      for (k = 0; k < Nlay; k ++)
      {
        hu_tend_q[k][i][j] /= avg_len_hu;
        hu_tend_gradM[k][i][j] /= avg_len_hu;
        hu_tend_gradKE[k][i][j] /= avg_len_hu;
        hu_tend_dhdt[k][i][j] /= avg_len_hu;
        hu_tend_A2[k][i][j] /= avg_len_hu;
        hu_tend_A4[k][i][j] /= avg_len_hu;
        hu_tend_rDrag[k][i][j] /= avg_len_hu;
        hu_tend_rSurf[k][i][j] /= avg_len_hu;
        hu_tend_CdBot[k][i][j] /= avg_len_hu;
        hu_tend_CdSurf[k][i][j] /= avg_len_hu;
        hu_tend_wind[k][i][j] /= avg_len_hu;
        hu_tend_buoy[k][i][j] /= avg_len_hu;
        hu_tend_relax[k][i][j] /= avg_len_hu;
        hu_tend_wdia[k][i][j] /= avg_len_hu;
        hu_tend_rand[k][i][j] /= avg_len_hu;
        hu_tend_Fbaro[k][i][j] /= avg_len_hu;
      }
    }
  }
  
  // Save averages to output files
  for (k = 0; k < Nlay; k ++)
  {
    constructOutputName(outdir,VARID_UMOM_Q,k,n,outfile);
    if (!writeOutputFile(outfile,hu_tend_q[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_UMOM_GRADM,k,n,outfile);
    if (!writeOutputFile(outfile,hu_tend_gradM[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_UMOM_GRADKE,k,n,outfile);
    if (!writeOutputFile(outfile,hu_tend_gradKE[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_UMOM_DHDT,k,n,outfile);
    if (!writeOutputFile(outfile,hu_tend_dhdt[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_UMOM_A2,k,n,outfile);
    if (!writeOutputFile(outfile,hu_tend_A2[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_UMOM_A4,k,n,outfile);
    if (!writeOutputFile(outfile,hu_tend_A4[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_UMOM_RDRAG,k,n,outfile);
    if (!writeOutputFile(outfile,hu_tend_rDrag[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_UMOM_RSURF,k,n,outfile);
    if (!writeOutputFile(outfile,hu_tend_rSurf[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_UMOM_CDBOT,k,n,outfile);
    if (!writeOutputFile(outfile,hu_tend_CdBot[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_UMOM_CDSURF,k,n,outfile);
    if (!writeOutputFile(outfile,hu_tend_CdSurf[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_UMOM_WIND,k,n,outfile);
    if (!writeOutputFile(outfile,hu_tend_wind[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_UMOM_BUOY,k,n,outfile);
    if (!writeOutputFile(outfile,hu_tend_buoy[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_UMOM_RELAX,k,n,outfile);
    if (!writeOutputFile(outfile,hu_tend_relax[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_UMOM_WDIA,k,n,outfile);
    if (!writeOutputFile(outfile,hu_tend_wdia[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_UMOM_RAND,k,n,outfile);
    if (!writeOutputFile(outfile,hu_tend_rand[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_UMOM_FBARO,k,n,outfile);
    if (!writeOutputFile(outfile,hu_tend_Fbaro[k],Nx,Ny)) return false;
  }
  
#pragma parallel
  
  // Reset averaging buffers
  for (i = 0; i < Nx; i ++)
  {
    for (j = 0; j < Ny; j ++)
    {
      for (k = 0; k < Nlay; k ++)
      {
        hu_tend_q[k][i][j] = 0;
        hu_tend_gradM[k][i][j] = 0;
        hu_tend_gradKE[k][i][j] = 0;
        hu_tend_dhdt[k][i][j] = 0;
        hu_tend_A2[k][i][j] = 0;
        hu_tend_A4[k][i][j] = 0;
        hu_tend_rDrag[k][i][j] = 0;
        hu_tend_rSurf[k][i][j] = 0;
        hu_tend_CdBot[k][i][j] = 0;
        hu_tend_CdSurf[k][i][j] = 0;
        hu_tend_wind[k][i][j] = 0;
        hu_tend_buoy[k][i][j] = 0;
        hu_tend_relax[k][i][j] = 0;
        hu_tend_wdia[k][i][j] = 0;
        hu_tend_rand[k][i][j] = 0;
        hu_tend_Fbaro[k][i][j] = 0;
      }
    }
  }
  
  return true;
}









/**
 *
 * writeVMomentumAverages
 *
 * Calculates averages of terms in V-momentum equation, writes averages to output
 * files, and resets averaging buffers and parameters.
 *
 */
mybool writeVMomentumAverages (uint n, char * outdir)
{
  uint i,j,k;
  char outfile[MAX_PARAMETER_FILENAME_LENGTH];
  
#pragma parallel
  
  // Divide by averaging step length to compute averages
  for (i = 0; i < Nx; i ++)
  {
    for (j = 0; j < Ny; j ++)
    {
      for (k = 0; k < Nlay; k ++)
      {
        hv_tend_q[k][i][j] /= avg_len_hv;
        hv_tend_gradM[k][i][j] /= avg_len_hv;
        hv_tend_gradKE[k][i][j] /= avg_len_hv;
        hv_tend_dhdt[k][i][j] /= avg_len_hv;
        hv_tend_A2[k][i][j] /= avg_len_hv;
        hv_tend_A4[k][i][j] /= avg_len_hv;
        hv_tend_rDrag[k][i][j] /= avg_len_hv;
        hv_tend_rSurf[k][i][j] /= avg_len_hv;
        hv_tend_CdBot[k][i][j] /= avg_len_hv;
        hv_tend_CdSurf[k][i][j] /= avg_len_hv;
        hv_tend_wind[k][i][j] /= avg_len_hv;
        hv_tend_buoy[k][i][j] /= avg_len_hv;
        hv_tend_relax[k][i][j] /= avg_len_hv;
        hv_tend_wdia[k][i][j] /= avg_len_hv;
        hv_tend_rand[k][i][j] /= avg_len_hv;
        hv_tend_Fbaro[k][i][j] /= avg_len_hv;
      }
    }
  }
  
  // Save averages to output files
  for (k = 0; k < Nlay; k ++)
  {
    constructOutputName(outdir,VARID_VMOM_Q,k,n,outfile);
    if (!writeOutputFile(outfile,hv_tend_q[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_VMOM_GRADM,k,n,outfile);
    if (!writeOutputFile(outfile,hv_tend_gradM[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_VMOM_GRADKE,k,n,outfile);
    if (!writeOutputFile(outfile,hv_tend_gradKE[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_VMOM_DHDT,k,n,outfile);
    if (!writeOutputFile(outfile,hv_tend_dhdt[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_VMOM_A2,k,n,outfile);
    if (!writeOutputFile(outfile,hv_tend_A2[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_VMOM_A4,k,n,outfile);
    if (!writeOutputFile(outfile,hv_tend_A4[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_VMOM_RDRAG,k,n,outfile);
    if (!writeOutputFile(outfile,hv_tend_rDrag[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_VMOM_RSURF,k,n,outfile);
    if (!writeOutputFile(outfile,hv_tend_rSurf[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_VMOM_CDBOT,k,n,outfile);
    if (!writeOutputFile(outfile,hv_tend_CdBot[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_VMOM_CDSURF,k,n,outfile);
    if (!writeOutputFile(outfile,hv_tend_CdSurf[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_VMOM_WIND,k,n,outfile);
    if (!writeOutputFile(outfile,hv_tend_wind[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_VMOM_BUOY,k,n,outfile);
    if (!writeOutputFile(outfile,hv_tend_buoy[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_VMOM_RELAX,k,n,outfile);
    if (!writeOutputFile(outfile,hv_tend_relax[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_VMOM_WDIA,k,n,outfile);
    if (!writeOutputFile(outfile,hv_tend_wdia[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_VMOM_RAND,k,n,outfile);
    if (!writeOutputFile(outfile,hv_tend_rand[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_VMOM_FBARO,k,n,outfile);
    if (!writeOutputFile(outfile,hv_tend_Fbaro[k],Nx,Ny)) return false;
  }
  
#pragma parallel
  
  // Reset averaging buffers
  for (i = 0; i < Nx; i ++)
  {
    for (j = 0; j < Ny; j ++)
    {
      for (k = 0; k < Nlay; k ++)
      {
        hv_tend_q[k][i][j] = 0;
        hv_tend_gradM[k][i][j] = 0;
        hv_tend_gradKE[k][i][j] = 0;
        hv_tend_dhdt[k][i][j] = 0;
        hv_tend_A2[k][i][j] = 0;
        hv_tend_A4[k][i][j] = 0;
        hv_tend_rDrag[k][i][j] = 0;
        hv_tend_rSurf[k][i][j] = 0;
        hv_tend_CdBot[k][i][j] = 0;
        hv_tend_CdSurf[k][i][j] = 0;
        hv_tend_wind[k][i][j] = 0;
        hv_tend_buoy[k][i][j] = 0;
        hv_tend_relax[k][i][j] = 0;
        hv_tend_wdia[k][i][j] = 0;
        hv_tend_rand[k][i][j] = 0;
        hv_tend_Fbaro[k][i][j] = 0;
      }
    }
  }
  
  return true;
}










/**
 *
 * writeThicknessAverages
 *
 * Calculates averages of terms in the thickness equation, writes averages to output
 * files, and resets averaging buffers and parameters.
 *
 */
mybool writeThicknessAverages (uint n, char * outdir)
{
  uint i,j,k;
  char outfile[MAX_PARAMETER_FILENAME_LENGTH];
  
#pragma parallel
  
  // Divide by averaging step length to compute averages
  for (i = 0; i < Nx; i ++)
  {
    for (j = 0; j < Ny; j ++)
    {
      for (k = 0; k < Nlay; k ++)
      {
        h_tend_adv[k][i][j] /= avg_len_h;
        h_tend_relax[k][i][j] /= avg_len_h;
      }
    }
  }
  
  // Save averages to output files
  for (k = 0; k < Nlay; k ++)
  {
    constructOutputName(outdir,VARID_THIC_ADV,k,n,outfile);
    if (!writeOutputFile(outfile,h_tend_adv[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_THIC_RELAX,k,n,outfile);
    if (!writeOutputFile(outfile,h_tend_relax[k],Nx,Ny)) return false;
  }
  
#pragma parallel
  
  // Reset averaging buffers
  for (i = 0; i < Nx; i ++)
  {
    for (j = 0; j < Ny; j ++)
    {
      for (k = 0; k < Nlay; k ++)
      {
        h_tend_adv[k][i][j] = 0;
        h_tend_relax[k][i][j] = 0;
      }
    }
  }
  
  return true;
}







/**
 *
 * writeEnergyAverages
 *
 * Calculates averages of terms in the energy equation, writes averages to output
 * files, and resets averaging buffers and parameters.
 *
 */
mybool writeEnergyAverages (uint n, char * outdir)
{
  uint i,j,k;
  char outfile[MAX_PARAMETER_FILENAME_LENGTH];
  
#pragma parallel
  
  // Divide by averaging step length to compute averages
  for (i = 0; i < Nx; i ++)
  {
    for (j = 0; j < Ny; j ++)
    {
      for (k = 0; k < Nlay; k ++)
      {
        e_flux_uP[k][i][j] /= avg_len_e;
        e_flux_uKE[k][i][j] /= avg_len_e;
        e_flux_uPE[k][i][j] /= avg_len_e;
        e_flux_vP[k][i][j] /= avg_len_e;
        e_flux_vKE[k][i][j] /= avg_len_e;
        e_flux_vPE[k][i][j] /= avg_len_e;
        e_tend_adv[k][i][j] /= avg_len_e;
        e_tend_gradM[k][i][j] /= avg_len_e;
        e_tend_wind[k][i][j] /= avg_len_e;
        e_tend_rDrag[k][i][j] /= avg_len_e;
        e_tend_rSurf[k][i][j] /= avg_len_e;
        e_tend_CdBot[k][i][j] /= avg_len_e;
        e_tend_CdSurf[k][i][j] /= avg_len_e;
        e_tend_A2[k][i][j] /= avg_len_e;
        e_tend_A4[k][i][j] /= avg_len_e;
        e_tend_wdiaPE[k][i][j] /= avg_len_e;
        e_tend_wdiaKE[k][i][j] /= avg_len_e;
        e_tend_Fbaro[k][i][j] /= avg_len_e;
        e_tend_buoy[k][i][j] /= avg_len_e;
        e_tend_rand[k][i][j] /= avg_len_e;
        e_tend_relax[k][i][j] /= avg_len_e;
      }
    }
  }
  
  // Save averages to output files
  for (k = 0; k < Nlay; k ++)
  {
    constructOutputName(outdir,VARID_ENERGY_UPFLUX,k,n,outfile);
    if (!writeOutputFile(outfile,e_flux_uP[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_ENERGY_UKEFLUX,k,n,outfile);
    if (!writeOutputFile(outfile,e_flux_uKE[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_ENERGY_UPEFLUX,k,n,outfile);
    if (!writeOutputFile(outfile,e_flux_uPE[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_ENERGY_VPFLUX,k,n,outfile);
    if (!writeOutputFile(outfile,e_flux_vP[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_ENERGY_VKEFLUX,k,n,outfile);
    if (!writeOutputFile(outfile,e_flux_vKE[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_ENERGY_VPEFLUX,k,n,outfile);
    if (!writeOutputFile(outfile,e_flux_vPE[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_ENERGY_ADV,k,n,outfile);
    if (!writeOutputFile(outfile,e_tend_adv[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_ENERGY_GRADM,k,n,outfile);
    if (!writeOutputFile(outfile,e_tend_gradM[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_ENERGY_WIND,k,n,outfile);
    if (!writeOutputFile(outfile,e_tend_wind[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_ENERGY_RDRAG,k,n,outfile);
    if (!writeOutputFile(outfile,e_tend_rDrag[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_ENERGY_RSURF,k,n,outfile);
    if (!writeOutputFile(outfile,e_tend_rSurf[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_ENERGY_CDBOT,k,n,outfile);
    if (!writeOutputFile(outfile,e_tend_CdBot[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_ENERGY_CDSURF,k,n,outfile);
    if (!writeOutputFile(outfile,e_tend_CdSurf[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_ENERGY_A2,k,n,outfile);
    if (!writeOutputFile(outfile,e_tend_A2[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_ENERGY_A4,k,n,outfile);
    if (!writeOutputFile(outfile,e_tend_A4[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_ENERGY_WDIAPE,k,n,outfile);
    if (!writeOutputFile(outfile,e_tend_wdiaPE[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_ENERGY_WDIAKE,k,n,outfile);
    if (!writeOutputFile(outfile,e_tend_wdiaKE[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_ENERGY_FBARO,k,n,outfile);
    if (!writeOutputFile(outfile,e_tend_Fbaro[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_ENERGY_BUOY,k,n,outfile);
    if (!writeOutputFile(outfile,e_tend_buoy[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_ENERGY_RAND,k,n,outfile);
    if (!writeOutputFile(outfile,e_tend_rand[k],Nx,Ny)) return false;
    constructOutputName(outdir,VARID_ENERGY_RELAX,k,n,outfile);
    if (!writeOutputFile(outfile,e_tend_relax[k],Nx,Ny)) return false;
  }
  
#pragma parallel
  
  // Reset averaging buffers
  for (i = 0; i < Nx; i ++)
  {
    for (j = 0; j < Ny; j ++)
    {
      for (k = 0; k < Nlay; k ++)
      {
        e_flux_uP[k][i][j] = 0;
        e_flux_uKE[k][i][j] = 0;
        e_flux_uPE[k][i][j] = 0;
        e_flux_vP[k][i][j] = 0;
        e_flux_vKE[k][i][j] = 0;
        e_flux_vPE[k][i][j] = 0;
        e_tend_adv[k][i][j] = 0;
        e_tend_gradM[k][i][j] = 0;
        e_tend_wind[k][i][j] = 0;
        e_tend_rDrag[k][i][j] = 0;
        e_tend_rSurf[k][i][j] = 0;
        e_tend_CdBot[k][i][j] = 0;
        e_tend_CdSurf[k][i][j] = 0;
        e_tend_A2[k][i][j] = 0;
        e_tend_A4[k][i][j] = 0;
        e_tend_wdiaPE[k][i][j] = 0;
        e_tend_wdiaKE[k][i][j] = 0;
        e_tend_Fbaro[k][i][j] = 0;
        e_tend_buoy[k][i][j] = 0;
        e_tend_rand[k][i][j] = 0;
        e_tend_relax[k][i][j] = 0;
      }
    }
  }
  
  return true;
}












/**
 *  printUsage
 *
 *  Prints an information message describing correct usage of the
 *  program in terms of required parameters.
 *
 */
void printUsage()
{
    printf
    (
     "USAGE: AWSIM <infile> <outdir>\n"
     "  \n"
     "  where <infile> is the name of the input parameter file, and\n"
     "  <outdir> is the directory into which output files should be\n"
     "  written.\n"
     "  \n"
     "  The input file must specify a series of input parameters\n"
     "  in name-value form, separated by whitespace, i.e.\n"
     "  <name> <value> [<name> <value> [...]]\n"
     "  \n"
     "  All matrix parameters should be formatted in column-major\n"
     "  order.\n"
     "  \n"
     "  Nlay                Number of shallow water layers.\n"
     "  Nx                  Number of x-gridpoints.\n"
     "  Ny                  Number of y-gridpoints.\n"
     "  Nt                  Number of time-values.\n"
     "                      Exactly two of {Nt, tmax, dt} must be specified.\n"
     "  Lx                  Zonal domain size. If set to zero then x grid spacing dx will be set\n"
     "                      equal to y grid spacing dy, and Lx=Nx*dx. Optional, default is 0.\n"
     "  Ly                  Meridional domain size. Must be > 0.\n"
     "  tmin                Time at which integration should start.\n"
     "                      Optional. If unset or negative then set to\n"
     "                      the pickup time of the restart, if specified. If tmin\n"
     "                      is set and >=0  then it overrides the pickup time.\n"
     "  tmax                Time at which integration should stop.\n"
     "                      Must be >0.\n"
     "                      Exactly two of {Nt, tmax, dt} must be specified.\n"
     "  dt                  Time step size.\n"
     "                      Exactly two of {Nt, tmax, dt} must be specified.\n"
     "  savefrequency       Data storage frequency, in units of time.\n"
     "                      Must be >0. Optional - default is the time step.\n"
     "  savefreqAvg         Averaged output frequency, in units of time. For values <=0 no\n"
     "                      averaged products will be calculated. Optional - default is 0.\n"
     "  savefreqUMom        Frequency of averaged output of u-momentum equation terms,\n"
     "                      in units of time. For values <=0 no averaged products will be\n"
     "                      calculated. Optional - default is 0.\n"
     "  savefreqVMom        Frequency of averaged output of v-momentum equation terms,\n"
     "                      in units of time. For values <=0 no averaged products will be\n"
     "                      calculated. Optional - default is 0.\n"
     "  savefreqThic        Frequency of averaged output of thickness equation terms,\n"
     "                      in units of time. For values <=0 no averaged products will be\n"
     "                      calculated. Optional - default is 0.\n"
     "  savefreqEnergy      Frequency of averaged diagnostics of energy budget terms,\n"
     "                      in units of time. For values <=0 no averaged products will be\n"
     "                      calculated. Optional - default is 0. N.B. The energy\n"
     "                      diagnostics neglect barotropic forcing, horizontal buoyancy\n"
     "                      effects and random forcing.\n"
     "  savefreqEZ          Storage frequency for energy and potential\n"
     "                      enstrophy, in units of time.\n"
     "                      Negative or zero values disable this diagnostic.\n"
     "                      Optional - default is 0.\n"
     "  restart             Set to 0 to initialize the simulation from t=0\n"
     "                      using input from the xInitFile parameters. Set to\n"
     "                      1 to pick up from a previous run using output files\n"
     "                      specified by startIdx as input files for this run.\n"
     "                      Optional - default is 0.\n"
     "  startIdx            Specifies the index of the output files to use as\n"
     "                      input files for this run. The output index counts\n"
     "                      the number of output files written in any given\n"
     "                      simulation. If the tmin parameter is not specified\n"
     "                      then the simulation start time is calculated as\n"
     "                      startIdx*savefrequency. Optional - default is 0.\n"
     "  timeSteppingScheme  Selects time stepping scheme. See SWdefs.h for\n"
     "                      numerical identifiers. Default is TIMESTEPPING_AB3.\n"
     "  momentumScheme      Selects momentum advection scheme. See SWdefs.h for\n"
     "                      numerical identifiers. Default is MOMENTUM_AL81.\n"
     "  thicknessScheme     Selects thickness advection scheme. See SWdefs.h for\n"
     "                      numerical identifiers. Default is THICKNESS_AL81.\n"
     "  tracerScheme        Selects tracer advection scheme. See SWdefs.h for\n"
     "                      numerical identifiers. Default is TRACER_AL81.\n"
     "  useRL               Set to 0 (default) to use a free surface, or any\n"
     "                      other value to use a rigid lid.\n"
     "  useTrad             Set to 0 (default) to use the traditional approximation,\n"
     "                      or any other value to include non-traditional terms.\n"
     "  useWallEW           Set to 0 (default) to remove western and eastern walls,\n"
     "                      or any other value to include them.\n"
     "  useWallNS           Set to 0 (default) to remove northern and southern walls,\n"
     "                      or any other value to include them.\n"
     "  useTracer           Set !=0 to allow passive/active tracer advection.\n"
     "                      Optional - default is 0 (no tracer).\n"
     "  useBuoyancy         Set !=0 to make the tracer field active as a buoyancy variable.\n"
     "                      Only activated if useTracer!=0, if useRL!=0 and if Nlay==1 for\n"
     "                      ease of implementation. Optional - default is 0 (passive tracer).\n"
     "  h0                  Minimum layer thickness. Must be >=0.\n"
     "                      Optional - default is 0.0.\n"
     "  hsml                Surface \"mixed layer\" thickness. Momentum forcing at the\n"
     "                      surface will be distributed over a layer of thickness hsml\n"
     "                      if more than one isopycnal layer lies within hsml of the\n"
     "                      surface. Must be >=0.\n"
     "                      Optional - default is 0.0, in which case all momentum\n"
     "                      input will be applied to the uppermost layer.\n"
     "  hbbl                Bottom \"boundary layer\" thickness. Momentum forcing at the\n"
     "                      ocean bed will be distributed over a layer of thickness hbbl\n"
     "                      if more than one isopycnal layer lies within hbbl of the\n"
     "                      ocean bed. Must be >=0.\n"
     "                      Optional - default is 0.0, in which case all momentum\n"
     "                      input will be applied to the lowermost layer.\n"
     "  hmin_surf           Minimum layer thickness for momentum forcing. Only the fraction\n"
     "                      of the layer thickness that exceeds hmin will receive surface\n"
     "                      momentum forcing. Must be >=0. Optional - default is 0.0\n"
     "  hmin_bot            Minimum layer thickness for momentum forcing. Only the fraction\n"
     "                      of the layer thickness that exceeds hmin will receive\n"
     "                      bottom momentum forcing. Must be >=0. Optional - default is 0.0\n"
     "  A2                  Constant Laplacian viscosity. Must be >= 0.\n"
     "                      Optional - default is 0.0.\n"
     "  A4                  Constant Biharmonic viscosity. Must be >= 0.\n"
     "                      Optional - default is 0.0.\n"
     "  A2smag              Smagorinsky Laplacian viscosity. Must be >= 0. Recommended range\n"
     "                      is 3--4. Optional - default is 0.0.\n"
     "  A4smag              Smagorinsky Biharmonic viscosity. Must be >= 0. Recommended range\n"
     "                      is 3--4. Optional - default is 0.0.\n"
     "  K2                  Constant Laplacian tracer diffusion coefficient. Must be >= 0.\n"
     "                      Optional - default is 0.0.\n"
     "  K4                  Constant Biharmonic tracer diffusion coefficient. Must be >= 0.\n"
     "                      Optional - default is 0.0.\n"
     "  uInitFile           String file name containing initial u-velocities.\n"
     "                      Default is u_k=0 everywhere for all k.\n"
     "                      Size: Nlay x Nx x Ny.\n"
     "  vInitFile           String file name containing initial v-velocities.\n"
     "                      Default is v_k=0 everywhere for all k.\n"
     "                      Size: Nlay x Nx x Ny.\n"
     "  hInitFile           String file name containing initial layer thicknesses.\n"
     "                      Default is h_k=1000/Nlay everywhere for all k.\n"
     "                      Size: Nlay x Nx x Ny.\n"
     "  bInitFile           String file name containing initial tracer concentrations.\n"
     "                      Default is b=0 everywhere for all k.\n"
     "                      Size: Nlay x Nx x Ny.\n"
     "  hsFile              String file name containing surface elevation.\n"
     "                      Size: Nx x Ny. Default: hs=0 everywhere. \n"
     "                      Not used if useRL=false.\n"
     "  hbFile              String file name containing bottom elevation.\n"
     "                      Size: Nx x Ny. Default: hb=-1000 everywhere.\n"
     "  OmegaxFile          String file name containing x component of rotation.\n"
     "                      Size: Nx x Ny. Default: Omega_x=0 everywhere.\n"
     "  OmegayFile          String file name containing y component of rotation.\n"
     "                      Size: Nx x Ny. Default: Omega_y=0 everywhere.\n"
     "  OmegazFile          String file name containing vertical component of rotation.\n"
     "                      Size: Nx+1 x Ny+1. Default: Omega_z=0 everywhere.\n"
     "  gFile               String file name containing Nlay x 1 vector of reduced\n"
     "                      gravities at each layer's upper interface. Default is\n"
     "                      g'=9.81 everywhere. \n"
     "  tauxFile            String file name containing x-component of surface wind.\n"
     "                      stress normalized by density.\n"
     "                      Size: Nx x Ny x tauNrecs. Default: taux=0 everywhere.\n"
     "  tauyFile            String file name containing y-component of surface wind.\n"
     "                      stress normalized by density.\n"
     "                      Size: Nx x Ny x tauNrecs. Default: tauy=0 everywhere.\n"
     "  tauPeriod           Period over which the imposed wind stress varies.\n"
     "                      If <=0 then only the first records in tauxFile and tauyFile\n"
     "                      will be used, i.e. wind stress is steady. Default is 0.\n"
     "  tauNrecs            Number of temporal records supplied in tauxFile and tauyFile.\n"
     "                      Default is 1. Must be > 0.\n"
     "  wDiaFile            String file name containing fixed component of diabatic\n"
     "                      velocity across layer interfaces. Note that velocities will.\n"
     "                      be modified if thickness and/or interface restoring is\n"
     "                      applied. Size: Nx x Ny x Nlay+1. Default: 0 everywhere.\n"
     "  wDiaPeriod          Period over which the imposed diapycnal velocity varies.\n"
     "                      If <=0 then only the first record in wDiaFile\n"
     "                      will be used, i.e. diapycnal velocity is steady. Default is 0.\n"
     "  wDiaNrecs           Number of temporal records supplied in wDiaFile.\n"
     "                      Default is 1. Must be > 0.\n"
     "  linDragCoeff        Linear bottom drag coefficient (m/s). Must be >= 0.\n"
     "                      Optional - default is 0.0.\n"
     "  quadDragCoeff       Quadratic bottom drag coefficient (dimensionless).\n"
     "                      Must be >= 0. Optional - default is 0.0.\n"
     "  linDragSurf         Linear surface drag coefficient (m/s). Must be >= 0.\n"
     "                      Optional - default is 0.0.\n"
     "  quadDragSurf        Quadratic surface drag coefficient (dimensionless).\n"
     "                      Must be >= 0. Optional - default is 0.0.\n"
     "  uLidFile            String file name containing x-component of rigid lid\n"
     "                      velocity (m/s), relative to which the surface drag will be\n"
     "                      calculated. Optional - default is 0.0.\n"
     "                      Size: Nx x Ny.\n"
     "  vLidFile            String file name containing v-component of rigid lid\n"
     "                      velocity (m/s), relative to which the surface drag will be\n"
     "                      calculated. Optional - default is 0.0.\n"
     "                      Size: Nx x Ny.\n"
     "  FbaroXFile          String file name containing barotropic acceleration (m/s^2)\n"
     "                      in x that will be applied throughout the water column.\n"
     "                      Default is zero everywhere. Size: Nx x Ny.\n"
     "  FbaroYFile          String file name containing barotropic acceleration (m/s^2)\n"
     "                      in y that will be applied throughout the water column.\n"
     "                      Default is zero everywhere. Size: Nx x Ny.\n"
     "  uRelaxFile          String file name containing relaxation values for\n"
     "                      u-velocities. Default is 0 everywhere for all k.\n"
     "                      Size: Nlay x Nx x Ny.\n"
     "  vRelaxFile          String file name containing relaxation values for\n"
     "                      v-velocities. Default is 0 everywhere for all k.\n"
     "                      Size: Nlay x Nx x Ny.\n"
     "  hRelaxFile          String file name containing relaxation values for\n"
     "                      layer thicknesses. Default is 1000/Nlay everywhere.\n"
     "                      Size: Nlay x Nx x Ny.\n"
     "  eRelaxFile          String file name containing relaxation values for\n"
     "                      layer upper-surface elevations. Default layer thicknesses of\n"
     "                      1000/Nlay everywhere. Size: Nlay x Nx x Ny.\n"
     "  bRelaxFile          String file name containing relaxation values for\n"
     "                      tracer concentration. Default is 0 everywhere.\n"
     "                      Size: Nlay x Nx x Ny.\n"
     "  uTimeFile           String file name containing relaxation time scales\n"
     "                      for u-velocities. Zero or negative values \n"
     "                      imply no relaxation. Optional, default is -1\n"
     "                      (no relaxation) everywhere. Size: Nlay x Nx x Ny.\n"
     "  vTimeFile           String file name containing relaxation time scales\n"
     "                      for v-velocities. Zero or  negative values \n"
     "                      imply no relaxation. Optional, default is -1\n"
     "                      (no relaxation) everywhere. Size: Nlay x Nx x Ny.\n"
     "  hTimeFile           String file name containing relaxation time scales\n"
     "                      for layer thicknesses. Zero or negative values\n"
     "                      imply no relaxation. Optional, default is -1\n"
     "                      (no relaxation) everywhere. Size: Nx x Ny.\n"
     "  eTimeFile           String file name containing relaxation time scales\n"
     "                      for layer upper surfaces. Zero or negative values\n"
     "                      imply no relaxation. Optional, default is -1\n"
     "                      (no relaxation) everywhere. Size: Nlay x Nx x Ny.\n"
     "  bTimeFile           String file name containing relaxation time scales\n"
     "                      for tracer concentrations. Zero or negative values\n"
     "                      imply no relaxation.  Optional, default is -1\n"
     "                      (no relaxation) everywhere. Size: Nlay x Nx x Ny.\n"
     "  \n"
     "  Pressure solve parameters:\n"
     "  \n"
     "  use_MG              Set to 0 (default) to use the point-SOR algorithm,\n"
     "                      or any other value to use the MultiGrid algorithm.\n"
     "                      Note: To use MultiGrid, Nx and Ny must be powers of 2.\n"
     "  use_fullMG          Set to 0 (default) to use vanilla MultiGrid, or any\n"
     "                      other value to use the Full MultiGrid algorithm.\n"
     "                      Full MultiGrid is theoretically faster, but in\n"
     "                      practice vanilla MultiGrid can perform better.\n"
     "  tol                 L_infinity accuracy with which to converge surface\n"
     "                      iterative solution for surface pressure.\n"
     "                      Must be > 0. Default is 10^-8 m^2/s^2.\n"
     "  omega_MG            Weighting parameter for Multigrid Jacobi iteration.\n"
     "                      Default is 2/3.\n"
     "  maxiters            Maxmimum allowable iterations in any given pressure\n"
     "                      solve. Must be > 0. Default is 10,000.\n"
     "  SOR_rp_max          Upper bound of search range to use when optimizing\n"
     "                      SOR relaxation parameter.\n"
     "                      Must lie in [0,2]. Default is 2.0.\n"
     "  SOR_rp_min          Lower bound of search range to use when optimizing\n"
     "                      SOR relaxation parameter.\n"
     "                      Must lie in [0,SOR_rp_max). Default is 1.0.\n"
     "  SOR_rp_acc          Accuracy with which to determine optimal\n"
     "                      SOR relaxation parameter.\n"
     "                      Must be > 0. Default is 0.001.\n"
     "  SOR_opt_freq        Frequency (as a number of iterations) with which to\n"
     "                      re-optimize the SOR relaxation parameter.\n"
     "                      Must be > 0. Default is 1000.\n"
     "  \n"
     "  Random forcing parameters:\n"
     "  \n"
     "  useRandomForcing    Set to 1 to use random forcing, or 0 (default)\n"
     "                      to deactivate.\n"
     "  useThicWeightedRF   Set to 1 to apply random forcing to the momentum\n"
     "                      equation tendency, or 0 (default) to apply to the\n"
     "                      velocity equation tendency.\n"
     "  useBarotropicRF     Set to 1 (default) to apply the same random forcing\n"
     "                      field in all isopycnal layers, or 0 to calculate a\n"
     "                      separate random forcing field for each layer.\n"
     "  RF_F0_rot           RMS amplitude of the rotational component of random forcing.\n"
     "                      Has units of m^2/s^2 if useThickWeightedRF=true, or units of\n"
     "                      m/s^2 if useThickWeightedRF=false.\n"
     "  RF_F0_div           RMS amplitude of the divergent component of random forcing.\n"
     "                      Has units of m^2/s^2 if useThickWeightedRF=true, or units of\n"
     "                      m/s^2 if useThickWeightedRF=false.\n"
     "  RF_tau              Autocorrelation time scale for the random forcing\n"
     "                      function. Must be > 0.\n"
     "  RF_fftMaskFile      String file name containing a mask for the amplitudes\n"
     "                      of the forcing function in spectral space, for each\n"
     "                      isopycnal layer. Units are arbitrary because this mask\n\n"
     "                      is normalized before multiplying the forcing function.\n"
     "                      Default is 1 everywhere. Size: Nlay x Nx x Ny.\n"
     "  RF_rotMaskFile      String file name containing a mask for the rotational\n"
     "                      component of the forcing function in real space, for each\n"
     "                      isopycnal layer. Default is 1 everywhere.\n"
     "                      Size: Nlay x Nx x Ny.\n"
     "  RF_divMaskFile      String file name containing a mask for the divergent\n"
     "                      component of the forcing function in real space, for each\n"
     "                      isopycnal layer. Default is 1 everywhere.\n"
     "                      Size: Nlay x Nx x Ny.\n"
   );
}







                               

/**
 * main
 *
 * Program entry point. Initialises all of the required memory,
 * reads in input parameters, performs iterations of the required
 * numerical method, and finally cleans up.
 *
 */
int main (int argc, char ** argv)
{
  
  // Time parameters
  real tmin = -1;         // Set to -1 initially - indicates that this parameter has not been set
  real tmax = 0;          // End time
  real t = 0;             // Time counter
  int Nt = 0;             // Number of steps
  real dt_s = 0;          // Data save time step
  real dt_EZ = 0;         // E/P save time step
  real t_next = 0;        // Next save time
  real t_next_EZ = 0;     // E/P next save time
  int tParamCnt = 0;      // Time parameter counter
  
  // Iteration parameters
  uint n = 0;
  uint n_saves = 0;
  uint n_EZ = 0;
  mybool restart = false;
  uint n0 = 0;
  
  // Time-averaging parameters
  real avg_fac = 0;       // Fraction of time step used for averaging
  uint n_avg = 0;         // Keeps track of number of averages computed
  real dt_avg = 0;        // Save time step for averaged output
  real t_next_avg = 0;    // Next save time for averaged output
  
  // Physical size of the domain
  real Lx = 0;
  real Ly = 0;
  
  // File descriptors for input parameter file and output data file
  char * infname = argv[1];
  FILE * infile = NULL;
  char * outdir = argv[2];
  char EZfilename[MAX_PARAMETER_FILENAME_LENGTH];
  FILE * EZfile = NULL;
  time_t now; // To store current time
  char tfilename[MAX_PARAMETER_FILENAME_LENGTH];
  FILE * tfile = NULL;
  char inbuf[256];
  mybool readerr = false;
  char fmt_str[80];
  
  // Looping variables
  int i = 0;
  int j = 0;
  int k = 0;
  int m = 0;
  int im1 = 0;
  int ip1 = 0;
  int jm1 = 0;
  int jp1 = 0;

  // Parameters involved in optimization of the relaxation parameters
  real rp = 0;
  
  // Parameters required for random forcing code
  real RF_phase = 0;            // Random phase
  real RF_Tratio = 0;           // Ratio of time step to autocorrelation timescale
  real RF_mask_norm = 0;        // For normalizing mask in wavenumber space
  int RF_idx = 0;               // For indexing fftw vectors
  int RF_idx0 = 0;
  int RF_idx_im1 = 0;
  int RF_idx_ip1 = 0;
  int RF_idx_jm1 = 0;
  int RF_idx_jp1 = 0;
  real *** RFamp_rot_fft = NULL;    // Spectral amplitudes of forcing
  real *** RFamp_div_fft = NULL;    // Spectral amplitudes of forcing
  real ** RF_kk_recip = NULL;   // Absolute wavenumber
  
  // Storage arrays for velocity, height and PV
  real * vars = NULL;
  real *** uu = NULL;
  real *** vv = NULL;
  real *** hh = NULL;
  real *** bb = NULL;
  real ** pi = NULL;
  
  // Buffers for data to be written to output files
  real * vars_buf = NULL;
  real *** uu_buf = NULL;
  real *** vv_buf = NULL;
  real *** hh_buf = NULL;
  real *** bb_buf = NULL;
  real ** pi_buf = NULL;
  
  // Storage for updated variables
  real * vars_out = NULL;
  real *** uu_out = NULL;
  real *** vv_out = NULL;
  real *** hh_out = NULL;
  real *** bb_out = NULL;
  
  // For storing/calculating time-averaged quantities
  real *** uu_avg = NULL;
  real *** vv_avg = NULL;
  real *** hh_avg = NULL;
  real *** MM_avg = NULL;
  real *** bb_avg = NULL;
  real ** pi_avg = NULL;
  real *** wdia_avg = NULL;
  real *** hu_avg = NULL;
  real *** hv_avg = NULL;
  real *** huu_avg = NULL;
  real *** hvv_avg = NULL;
  real *** huv_avg = NULL;
  real *** udt = NULL;
  real *** vdt = NULL;
  real *** hdt = NULL;
  real *** Mdt = NULL;
  real *** bdt = NULL;
  real ** pidt = NULL;
  real *** wdt = NULL;
  real *** hudt = NULL;
  real *** hvdt = NULL;
  real *** huudt = NULL;
  real *** hvvdt = NULL;
  real *** huvdt = NULL;
  
  // Work arrays and storage arrays for time derivatives
  real * k13 = NULL;
  real * k24 = NULL;
  real * dt_vars = NULL;
  real * dt_vars_1 = NULL;
  real * dt_vars_2 = NULL;
  real * pvars = NULL;
  
  // Used for calculation of energy and potential enstrophy
  real KE = 0;
  real PE = 0;
  real E = 0;
  real * Z = 0;
  
  // Stores data required for parsing input parameters
  paramdata params[NPARAMS];
  int paramcntr = 0;
  
  // Filename holders for input parameter arrays
  char uInitFile[MAX_PARAMETER_FILENAME_LENGTH];
  char vInitFile[MAX_PARAMETER_FILENAME_LENGTH];
  char hInitFile[MAX_PARAMETER_FILENAME_LENGTH];
  char bInitFile[MAX_PARAMETER_FILENAME_LENGTH];
  char hsFile[MAX_PARAMETER_FILENAME_LENGTH];
  char hbFile[MAX_PARAMETER_FILENAME_LENGTH];
  char OmegaxFile[MAX_PARAMETER_FILENAME_LENGTH];
  char OmegayFile[MAX_PARAMETER_FILENAME_LENGTH];
  char OmegazFile[MAX_PARAMETER_FILENAME_LENGTH];
  char gFile[MAX_PARAMETER_FILENAME_LENGTH];
  char tauxFile[MAX_PARAMETER_FILENAME_LENGTH];
  char tauyFile[MAX_PARAMETER_FILENAME_LENGTH];
  char wDiaFile[MAX_PARAMETER_FILENAME_LENGTH];
  char uLidFile[MAX_PARAMETER_FILENAME_LENGTH];
  char vLidFile[MAX_PARAMETER_FILENAME_LENGTH];
  char FbaroXFile[MAX_PARAMETER_FILENAME_LENGTH];
  char FbaroYFile[MAX_PARAMETER_FILENAME_LENGTH];
  char uRelaxFile[MAX_PARAMETER_FILENAME_LENGTH];
  char vRelaxFile[MAX_PARAMETER_FILENAME_LENGTH];
  char hRelaxFile[MAX_PARAMETER_FILENAME_LENGTH];
  char eRelaxFile[MAX_PARAMETER_FILENAME_LENGTH];
  char bRelaxFile[MAX_PARAMETER_FILENAME_LENGTH];
  char uTimeFile[MAX_PARAMETER_FILENAME_LENGTH];
  char vTimeFile[MAX_PARAMETER_FILENAME_LENGTH];
  char hTimeFile[MAX_PARAMETER_FILENAME_LENGTH];
  char eTimeFile[MAX_PARAMETER_FILENAME_LENGTH];
  char bTimeFile[MAX_PARAMETER_FILENAME_LENGTH];
  char RF_fftMaskFile[MAX_PARAMETER_FILENAME_LENGTH];
  char RF_rotMaskFile[MAX_PARAMETER_FILENAME_LENGTH];
  char RF_divMaskFile[MAX_PARAMETER_FILENAME_LENGTH];
  
  // Default file name parameters - zero-length strings
  uInitFile[0] = '\0';
  vInitFile[0] = '\0';
  hInitFile[0] = '\0';
  bInitFile[0] = '\0';
  hsFile[0] = '\0';
  hbFile[0] = '\0';
  OmegaxFile[0] = '\0';
  OmegayFile[0] = '\0';
  OmegazFile[0] = '\0';
  gFile[0] = '\0';
  tauxFile[0] = '\0';
  tauyFile[0] = '\0';
  wDiaFile[0] = '\0';
  uLidFile[0] = '\0';
  vLidFile[0] = '\0';
  FbaroXFile[0] = '\0';
  FbaroYFile[0] = '\0';
  uRelaxFile[0] = '\0';
  vRelaxFile[0] = '\0';
  hRelaxFile[0] = '\0';
  eRelaxFile[0] = '\0';
  bRelaxFile[0] = '\0';
  uTimeFile[0] = '\0';
  vTimeFile[0] = '\0';
  hTimeFile[0] = '\0';
  eTimeFile[0] = '\0';
  bTimeFile[0] = '\0';
  RF_fftMaskFile[0] = '\0';
  RF_rotMaskFile[0] = '\0';
  RF_divMaskFile[0] = '\0';

  // Define input parameter data
  setParam(params,paramcntr++,"Nlay","%u",&Nlay,false);
  setParam(params,paramcntr++,"Nx","%u",&Nx,false);
  setParam(params,paramcntr++,"Ny","%u",&Ny,false);
  setParam(params,paramcntr++,"Nt","%u",&Nt,true);
  setParam(params,paramcntr++,"tmin",FLT_FMT,&tmin,true);
  setParam(params,paramcntr++,"tmax",FLT_FMT,&tmax,false);
  setParam(params,paramcntr++,"Lx",FLT_FMT,&Lx,true);
  setParam(params,paramcntr++,"Ly",FLT_FMT,&Ly,false);
  setParam(params,paramcntr++,"dt",FLT_FMT,&dt,true);
  setParam(params,paramcntr++,"savefrequency",FLT_FMT,&dt_s,true);
  setParam(params,paramcntr++,"savefreqAvg",FLT_FMT,&dt_avg,true);
  setParam(params,paramcntr++,"savefreqUMom",FLT_FMT,&dt_avg_hu,true);
  setParam(params,paramcntr++,"savefreqVMom",FLT_FMT,&dt_avg_hv,true);
  setParam(params,paramcntr++,"savefreqThic",FLT_FMT,&dt_avg_h,true);
  setParam(params,paramcntr++,"savefreqEnergy",FLT_FMT,&dt_avg_e,true);
  setParam(params,paramcntr++,"savefreqEZ",FLT_FMT,&dt_EZ,true);
  setParam(params,paramcntr++,"restart","%d",&restart,true);
  setParam(params,paramcntr++,"startIdx","%u",&n0,true);
  setParam(params,paramcntr++,"timeSteppingScheme","%u",&timeSteppingScheme,true);
  setParam(params,paramcntr++,"momentumScheme","%u",&momentumScheme,true);
  setParam(params,paramcntr++,"thicknessScheme","%u",&thicknessScheme,true);
  setParam(params,paramcntr++,"tracerScheme","%u",&tracerScheme,true);
  setParam(params,paramcntr++,"useRL","%d",&useRL,true);
  setParam(params,paramcntr++,"useTrad","%d",&useTrad,true);
  setParam(params,paramcntr++,"useWallEW","%d",&useWallEW,true);
  setParam(params,paramcntr++,"useWallNS","%d",&useWallNS,true);
  setParam(params,paramcntr++,"useTracer","%d",&useTracer,true);
  setParam(params,paramcntr++,"useBuoyancy","%d",&useBuoyancy,true);
  setParam(params,paramcntr++,"h0",FLT_FMT,&h0,true);
  setParam(params,paramcntr++,"hsml",FLT_FMT,&hsml,true);
  setParam(params,paramcntr++,"hbbl",FLT_FMT,&hbbl,true);
  setParam(params,paramcntr++,"hmin_surf",FLT_FMT,&hmin_surf,true);
  setParam(params,paramcntr++,"hmin_bot",FLT_FMT,&hmin_bot,true);
  setParam(params,paramcntr++,"A2",FLT_FMT,&A2const,true);
  setParam(params,paramcntr++,"A4",FLT_FMT,&A4const,true);
  setParam(params,paramcntr++,"A2smag",FLT_FMT,&A2smag,true);
  setParam(params,paramcntr++,"A4smag",FLT_FMT,&A4smag,true);
  setParam(params,paramcntr++,"K2",FLT_FMT,&K2,true);
  setParam(params,paramcntr++,"K4",FLT_FMT,&K4,true);
  setParam(params,paramcntr++,"uInitFile","%s",uInitFile,true);
  setParam(params,paramcntr++,"vInitFile","%s",vInitFile,true);
  setParam(params,paramcntr++,"hInitFile","%s",hInitFile,true);
  setParam(params,paramcntr++,"bInitFile","%s",bInitFile,true);
  setParam(params,paramcntr++,"hsFile","%s",hsFile,true);
  setParam(params,paramcntr++,"hbFile","%s",hbFile,true);
  setParam(params,paramcntr++,"OmegaxFile","%s",OmegaxFile,true);
  setParam(params,paramcntr++,"OmegayFile","%s",OmegayFile,true);
  setParam(params,paramcntr++,"OmegazFile","%s",OmegazFile,true);
  setParam(params,paramcntr++,"gFile","%s",gFile,true);
  setParam(params,paramcntr++,"tauxFile","%s",tauxFile,true);
  setParam(params,paramcntr++,"tauyFile","%s",tauyFile,true);
  setParam(params,paramcntr++,"tauPeriod",FLT_FMT,&tauPeriod,true);
  setParam(params,paramcntr++,"tauNrecs","%u$",&tauNrecs,true);
  setParam(params,paramcntr++,"wDiaFile","%s",wDiaFile,true);
  setParam(params,paramcntr++,"wDiaPeriod",FLT_FMT,&wDiaPeriod,true);
  setParam(params,paramcntr++,"wDiaNrecs","%u$",&wDiaNrecs,true);
  setParam(params,paramcntr++,"linDragCoeff",FLT_FMT,&rDrag,true);
  setParam(params,paramcntr++,"quadDragCoeff",FLT_FMT,&CdBot,true);
  setParam(params,paramcntr++,"linDragSurf",FLT_FMT,&rSurf,true);
  setParam(params,paramcntr++,"quadDragSurf",FLT_FMT,&CdSurf,true);
  setParam(params,paramcntr++,"uLidFile","%s",&uLidFile,true);
  setParam(params,paramcntr++,"vLidFile","%s",&vLidFile,true);
  setParam(params,paramcntr++,"FbaroXFile","%s",&FbaroXFile,true);
  setParam(params,paramcntr++,"FbaroYFile","%s",&FbaroYFile,true);
  setParam(params,paramcntr++,"uRelaxFile","%s",&uRelaxFile,true);
  setParam(params,paramcntr++,"vRelaxFile","%s",&vRelaxFile,true);
  setParam(params,paramcntr++,"hRelaxFile","%s",&hRelaxFile,true);
  setParam(params,paramcntr++,"eRelaxFile","%s",&eRelaxFile,true);
  setParam(params,paramcntr++,"bRelaxFile","%s",&bRelaxFile,true);
  setParam(params,paramcntr++,"uTimeFile","%s",&uTimeFile,true);
  setParam(params,paramcntr++,"vTimeFile","%s",&vTimeFile,true);
  setParam(params,paramcntr++,"hTimeFile","%s",&hTimeFile,true);
  setParam(params,paramcntr++,"eTimeFile","%s",&eTimeFile,true);
  setParam(params,paramcntr++,"bTimeFile","%s",&bTimeFile,true);
  
  // Pressure solve parameters
  setParam(params,paramcntr++,"use_MG","%d",&use_MG,true);
  setParam(params,paramcntr++,"use_fullMG","%d",&use_fullMG,true);
  setParam(params,paramcntr++,"tol",FLT_FMT,&pi_tol,true);
  setParam(params,paramcntr++,"omega_MG",FLT_FMT,&omega_WJ,true);
  setParam(params,paramcntr++,"maxiters","%u",&maxiters,true);
  setParam(params,paramcntr++,"SOR_rp_max",FLT_FMT,&rp_opt_max,true);
  setParam(params,paramcntr++,"SOR_rp_min",FLT_FMT,&rp_opt_min,true);
  setParam(params,paramcntr++,"SOR_rp_acc",FLT_FMT,&rp_acc_max,true);
  setParam(params,paramcntr++,"SOR_opt_freq","%u",&rp_opt_freq,true);
  
  // Random forcing parameters
  setParam(params,paramcntr++,"useRandomForcing","%d",&useRandomForcing,true);
  setParam(params,paramcntr++,"useThicWeightedRF","%d",&useThicWeightedRF,true);
  setParam(params,paramcntr++,"useBarotropicRF","%d",&useBarotropicRF,true);
  setParam(params,paramcntr++,"RF_F0_rot",FLT_FMT,&RF_F0_rot,true);
  setParam(params,paramcntr++,"RF_F0_div",FLT_FMT,&RF_F0_div,true);
  setParam(params,paramcntr++,"RF_tau",FLT_FMT,&RF_tau,true);
  setParam(params,paramcntr++,"RF_fftMaskFile","%s",&RF_fftMaskFile,true);
  setParam(params,paramcntr++,"RF_rotMaskFile","%s",&RF_rotMaskFile,true);
  setParam(params,paramcntr++,"RF_divMaskFile","%s",&RF_divMaskFile,true);
  
  // Check that a file name has been specified
  if (argc < 3)
  {
    fprintf(stderr,"ERROR: No input parameter file name supplied\n");
    printUsage();
    return 0;
  }
  
  // Calculate elapsed time and write to a file
  strcpy(tfilename,outdir);
  strcat(tfilename,"/");
  strcat(tfilename,TFILE);
  tfile = fopen(tfilename,"w");
  if (tfile != NULL)
  {
    time(&now);
    fprintf(tfile,"Program started at %s\n", ctime(&now));
    fflush(tfile);
  }
  
  
  //////////////////////////////////////////
  ///// BEGIN READING INPUT PARAMETERS /////
  //////////////////////////////////////////
  
  // Attempt to read input parameter data. Errors will be printed
  // to stderr, and will result in 'false' being returned.
  if (!readParams(infname,params,NPARAMS,stderr))
  {
    printUsage();
    return 0;
  }

  // Check that required parameters take legal values
  if ((Nlay == 0) || (Nx == 0) || (Ny == 0) || (dt < 0) || (dt_s < 0) ||
      (tmax < 0.0) || (Lx < 0.0) || (Ly <= 0.0) || (h0 < 0.0) || (hsml < 0.0) || (hbbl < 0.0) || (hmin_surf < 0.0) || (hmin_bot < 0.0) ||
      (A2const < 0.0) || (A4const < 0.0) || (A2smag < 0.0) || (A4smag < 0.0) || (K2 < 0.0) || (K4 < 0.0) ||
      (rDrag < 0.0) || (rSurf < 0.0) || (CdBot < 0.0) || (CdSurf < 0.0) ||
      (pi_tol <= 0) || (rp_opt_max > 2.0) || (rp_opt_min < 0.0) ||
      (rp_opt_min >= rp_opt_max) || (rp_acc_max <= 0) || (rp_opt_freq <= 0) ||
      (maxiters <= 0) || (tauNrecs <= 0) || (wDiaNrecs <= 0))
  {
    fprintf(stderr,"ERROR: Invalid input parameter values\n");
    printUsage();
    return 0;
  }
  
  // Check number of time parameters
  tParamCnt = ((Nt > 0) ? 1 : 0) + ((tmax > 0) ? 1 : 0) + ((dt > 0) ? 1 : 0);
  if (tParamCnt != 2)
  {
    fprintf(stderr,"Must specify exactly two of {Nt, tmax, dt}\n");
    printUsage();
    return 0;
  }
  
  // MultiGrid only implemented for power-of-2 grid sizes
  if (useRL && use_MG && (!ispow2u(Nx) || !ispow2u(Ny)))
  {
    fprintf(stderr,"ERROR: To use MultiGrid methods, grid dimensions must be powers of 2.\n");
    printUsage();
    return 0;
  }
  
  /// Buoyancy only implemented for 1-layer case with rigid lid
  if (useTracer && useBuoyancy && ((!useRL) || (Nlay>1)))
  {
    fprintf(stderr,"ERROR: Active (buoyancy) tracer only available in one-layer setup with rigid lid.\n");
    printUsage();
    return 0;
  }
  
  /// Random forcing autocorrelation time scale must be positive
  if (useRandomForcing && (RF_tau<=0))
  {
    fprintf(stderr,"ERROR: Random forcing autocorrelation time scale must be positive.\n");
    printUsage();
    return 0;
  }
  
  // Flag to determine whether any wind forcing is prescribed
  useWind = (strlen(tauxFile) > 0)  || (strlen(tauyFile) > 0);
  
  // Flag to determine whether any barotropic forcing is prescribed
  useFbaro = (strlen(FbaroXFile) > 0)  || (strlen(FbaroYFile) > 0);
  
  // Flag to determine whether relaxation is in use
  useRelax = (strlen(uTimeFile) > 0) || (strlen(vTimeFile) > 0) || (strlen(hTimeFile) > 0)  || (strlen(eTimeFile) > 0) || (useTracer && (strlen(bTimeFile) > 0));
  
  // Flag to determine whether diabatic velocity is in use
  useWDia = useRelax || (strlen(wDiaFile) > 0);

  // Gridpoint counters
  N = Nlay*Nx*Ny; // Number of gridpoints for each interior variable (u,v,h)
  Ntotal = useTracer ? 4*N : 3*N; // Total number of gridpoints passed to 'tderiv'
  N_g = Nlay*(Nx+2*Ng)*(Ny+2*Ng); // As above, but including ghost points
  Ntotal_g = useTracer ? 4*N_g : 3*N_g;
  
  // Grid dimensions
  dy = Ly/Ny;
  if (Lx == 0)
  {
    dx = dy;
    Lx = Nx*dx;
  }
  else
  {
    dx = Lx/Nx;
  }
  dxsq = dx*dx;
  dysq = dy*dy;

  
  // Set the integration start time
  // Note: If dt_s is not specified then this will just be zero
  if (tmin < 0)
  {
    if (restart)
    {
      tmin = n0*dt_s;
      if (tmin >= tmax)
      {
        fprintf(stderr,"Integration start time (=startIdx*savefrequency) exceeds integration end time (tmax)");
        printUsage();
        return 0;
      }
    }
    else
    {
      tmin = 0;
    }
  }
  else
  {
    if (tmin >= tmax)
    {
      fprintf(stderr,"Integration start time (tmin) exceeds integration end time (tmax)");
      printUsage();
      return 0;
    }
  }
  t = tmin;
  
  // Determine the unknown time parameter
  if (Nt == 0)
  {
    Nt = (uint) ceil((tmax-tmin)/dt);
  }
  if (tmax == 0)
  {
    tmax = tmin + dt*Nt;
  }
  if (dt == 0)
  {
    dt = (tmax-tmin)/Nt;
  }

  // Determine output frequency and first output time
  if (dt_s == 0)
  {
    dt_s = dt;
  }
  t_next = tmin + dt_s;
  
  // If restarting then initialize output count accordingly
  if (restart)
  {
    n_saves = n0;
  }
  else
  {
    n_saves = 0;
  }
  
  // Initialize next output time for EZ diagnostics and EZ counter
  if (dt_EZ > 0)
  {
    t_next_EZ = tmin + dt_EZ;
    n_EZ = round(tmin/dt_EZ) + 1;
  }
  
  // Initialize next output time for averaged diagnostics
  if (dt_avg > 0)
  {
    t_next_avg = tmin + dt_avg;
    n_avg = round(tmin/dt_avg) + 1;
  }
  
  // Initialize next output time for U-momentum diagnostics
  if (dt_avg_hu > 0)
  {
    t_next_avg_hu = tmin + dt_avg_hu;
    n_prev_avg_hu = 0;
    n_avg_hu = round(tmin/dt_avg_hu) + 1;
  }
  
  // Initialize next output time for V-momentum diagnostics
  if (dt_avg_hv > 0)
  {
    t_next_avg_hv = tmin + dt_avg_hv;
    n_prev_avg_hv = 0;
    n_avg_hv = round(tmin/dt_avg_hv) + 1;
  }
  
  // Initialize next output time for thickness diagnostics
  if (dt_avg_h > 0)
  {
    t_next_avg_h = tmin + dt_avg_h;
    n_prev_avg_h = 0;
    n_avg_h = round(tmin/dt_avg_h) + 1;
  }
  
  // Initialize next output time for V-momentum diagnostics
  if (dt_avg_e > 0)
  {
    t_next_avg_e = tmin + dt_avg_e;
    n_prev_avg_e = 0;
    n_avg_e = round(tmin/dt_avg_e) + 1;
  }

  //  Multigrid-specific parameters
  if (useRL && use_MG)
  {
    // Number of powers of 2 constituting each grid dimension
    Npx = log2u(Nx);
    Npy = log2u(Ny);
    
    // Ngrids -> Number of grids in the multigrid scheme
    // F_len -> Length of vector that will constitute the very smallest grid
    if (Npx > Npy)
    {
      Ngrids = Npy + 1;
      F_len = pow2u(Npx-Npy);
    }
    else
    {
      Ngrids = Npx + 1;
      F_len = pow2u(Npy-Npx);
    }
  }
  
  // Viscosity flags
  useA2 = (A2const > 0) || (A2smag > 0);
  useA4 = (A4const > 0) || (A4smag > 0);
  

  
  ////////////////////////////////////////
  ///// END READING INPUT PARAMETERS /////
  ////////////////////////////////////////
  


  ///////////////////////////////////
  ///// BEGIN MEMORY ALLOCATION /////
  ///////////////////////////////////
  
  VECALLOC(vars,Ntotal);
  VECALLOC(vars_out,Ntotal);
  VECALLOC(vars_buf,Ntotal);
  VECALLOC(k13,Ntotal);
  VECALLOC(k24,Ntotal);
  VECALLOC(dt_vars,Ntotal);
  VECALLOC(dt_vars_1,Ntotal);
  VECALLOC(dt_vars_2,Ntotal);
  MATALLOC(hhs,Nx,Ny);
  MATALLOC(hhb,Nx,Ny);
  MATALLOC(Omega_x,Nx,Ny);
  MATALLOC(Omega_y,Nx,Ny);
  MATALLOC(Omega_z,Nx+1,Ny+1);
  VECALLOC(gg,Nlay);
  VECALLOC(geff,Nlay);
  VECALLOC(Z,Nlay);
  MATALLOC3(taux,tauNrecs,Nx,Ny);
  MATALLOC3(tauy,tauNrecs,Nx,Ny);
  MATALLOC(uLid,Nx,Ny);
  MATALLOC(vLid,Nx,Ny);
  MATALLOC(uLid_g,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(vLid_g,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(Fbaro_x,Nx,Ny);
  MATALLOC(Fbaro_y,Nx,Ny);
  
  // Divide up the 'vars' array between the dependent variables
  pvars = vars;
  vec2mat3(pvars,&uu,Nlay,Nx,Ny);
  vec2mat3(pvars+=N,&vv,Nlay,Nx,Ny);
  vec2mat3(pvars+=N,&hh,Nlay,Nx,Ny);
  if (useTracer)
  {
    vec2mat3(pvars+=N,&bb,Nlay,Nx,Ny);
  }
  
  // Buffers for variables to be written to output files
  pvars = vars_buf;
  vec2mat3(pvars,&uu_buf,Nlay,Nx,Ny);
  vec2mat3(pvars+=N,&vv_buf,Nlay,Nx,Ny);
  vec2mat3(pvars+=N,&hh_buf,Nlay,Nx,Ny);
  if (useTracer)
  {
    vec2mat3(pvars+=N,&bb_buf,Nlay,Nx,Ny);
  }
  
  // Buffers for output data
  pvars = vars_out;
  vec2mat3(pvars,&uu_out,Nlay,Nx,Ny);
  vec2mat3(pvars+=N,&vv_out,Nlay,Nx,Ny);
  vec2mat3(pvars+=N,&hh_out,Nlay,Nx,Ny);
  if (useTracer)
  {
    vec2mat3(pvars+=N,&bb_out,Nlay,Nx,Ny);
  }
  
  // Work arrays for model variables, including space for ghost points
  VECALLOC(vars_w,Ntotal_g);
  pvars = vars_w;
  vec2mat3(pvars,&uu_w,Nlay,Nx+2*Ng,Ny+2*Ng);
  vec2mat3(pvars+=N_g,&vv_w,Nlay,Nx+2*Ng,Ny+2*Ng);
  vec2mat3(pvars+=N_g,&hh_w,Nlay,Nx+2*Ng,Ny+2*Ng);
  if (useTracer)
  {
    vec2mat3(pvars+=N_g,&bb_w,Nlay,Nx+2*Ng,Ny+2*Ng);
  }
  
  // Work arrays for model variables, NOT including space for ghost points
  pvars = vars_buf;
  vec2mat3(pvars,&uu_d,Nlay,Nx,Ny);
  vec2mat3(pvars+=N,&vv_d,Nlay,Nx,Ny);
  vec2mat3(pvars+=N,&hh_d,Nlay,Nx,Ny);
  if (useTracer)
  {
    vec2mat3(pvars+=N,&bb_d,Nlay,Nx,Ny);
  }
  
  // Work arrays for time derivatives of model variables. Here they are just
  // assigned to point to 'vars_buf' in order to initialize all of the pointers.
  // When 'tderiv' is called, these pointers will be reassigned.
  pvars = vars_buf;
  vec2mat3(pvars,&dt_uu_w,Nlay,Nx,Ny);
  vec2mat3(pvars+=N,&dt_vv_w,Nlay,Nx,Ny);
  vec2mat3(pvars+=N,&dt_hh_w,Nlay,Nx,Ny);
  if (useTracer)
  {
    vec2mat3(pvars+=N,&dt_bb_w,Nlay,Nx,Ny);
  }
  
  // For calculating online averaged of model variables and their products.
  if (dt_avg > 0)
  {
    MATALLOC3(uu_avg,Nlay,Nx,Ny);
    MATALLOC3(vv_avg,Nlay,Nx,Ny);
    MATALLOC3(hh_avg,Nlay,Nx,Ny);
    MATALLOC3(MM_avg,Nlay,Nx,Ny);
    if (useTracer)
    {
      MATALLOC3(bb_avg,Nlay,Nx,Ny);
    }
    MATALLOC(pi_avg,Nx,Ny);
    MATALLOC3(wdia_avg,Nlay+1,Nx,Ny);
    MATALLOC3(hu_avg,Nlay,Nx,Ny);
    MATALLOC3(huu_avg,Nlay,Nx,Ny);
    MATALLOC3(hv_avg,Nlay,Nx,Ny);
    MATALLOC3(hvv_avg,Nlay,Nx,Ny);
    MATALLOC3(huv_avg,Nlay,Nx,Ny);
    MATALLOC3(udt,Nlay,Nx,Ny);
    MATALLOC3(vdt,Nlay,Nx,Ny);
    MATALLOC3(hdt,Nlay,Nx,Ny);
    MATALLOC3(Mdt,Nlay,Nx,Ny);
    if (useTracer)
    {
      MATALLOC3(bdt,Nlay,Nx,Ny);
    }
    MATALLOC(pidt,Nx,Ny);
    MATALLOC3(wdt,Nlay+1,Nx,Ny);
    MATALLOC3(hudt,Nlay,Nx,Ny);
    MATALLOC3(huudt,Nlay,Nx,Ny);
    MATALLOC3(hvdt,Nlay,Nx,Ny);
    MATALLOC3(hvvdt,Nlay,Nx,Ny);
    MATALLOC3(huvdt,Nlay,Nx,Ny);
  }
  
  // For averaging u-momentum equation
  if (dt_avg_hu > 0)
  {
    MATALLOC3(hu_tend_q,Nlay,Nx,Ny);
    MATALLOC3(hu_tend_gradM,Nlay,Nx,Ny);
    MATALLOC3(hu_tend_gradKE,Nlay,Nx,Ny);
    MATALLOC3(hu_tend_dhdt,Nlay,Nx,Ny);
    MATALLOC3(hu_tend_A2,Nlay,Nx,Ny);
    MATALLOC3(hu_tend_A4,Nlay,Nx,Ny);
    MATALLOC3(hu_tend_rDrag,Nlay,Nx,Ny);
    MATALLOC3(hu_tend_rSurf,Nlay,Nx,Ny);
    MATALLOC3(hu_tend_CdBot,Nlay,Nx,Ny);
    MATALLOC3(hu_tend_CdSurf,Nlay,Nx,Ny);
    MATALLOC3(hu_tend_wind,Nlay,Nx,Ny);
    MATALLOC3(hu_tend_buoy,Nlay,Nx,Ny);
    MATALLOC3(hu_tend_relax,Nlay,Nx,Ny);
    MATALLOC3(hu_tend_wdia,Nlay,Nx,Ny);
    MATALLOC3(hu_tend_rand,Nlay,Nx,Ny);
    MATALLOC3(hu_tend_Fbaro,Nlay,Nx,Ny);
  }
  
  // For averaging v-momentum equation
  if (dt_avg_hv > 0)
  {
    MATALLOC3(hv_tend_q,Nlay,Nx,Ny);
    MATALLOC3(hv_tend_gradM,Nlay,Nx,Ny);
    MATALLOC3(hv_tend_gradKE,Nlay,Nx,Ny);
    MATALLOC3(hv_tend_dhdt,Nlay,Nx,Ny);
    MATALLOC3(hv_tend_A2,Nlay,Nx,Ny);
    MATALLOC3(hv_tend_A4,Nlay,Nx,Ny);
    MATALLOC3(hv_tend_rDrag,Nlay,Nx,Ny);
    MATALLOC3(hv_tend_rSurf,Nlay,Nx,Ny);
    MATALLOC3(hv_tend_CdBot,Nlay,Nx,Ny);
    MATALLOC3(hv_tend_CdSurf,Nlay,Nx,Ny);
    MATALLOC3(hv_tend_wind,Nlay,Nx,Ny);
    MATALLOC3(hv_tend_buoy,Nlay,Nx,Ny);
    MATALLOC3(hv_tend_relax,Nlay,Nx,Ny);
    MATALLOC3(hv_tend_wdia,Nlay,Nx,Ny);
    MATALLOC3(hv_tend_rand,Nlay,Nx,Ny);
    MATALLOC3(hv_tend_Fbaro,Nlay,Nx,Ny);
  }
  
  // For averaging thickness equation
  if (dt_avg_h > 0)
  {
    MATALLOC3(h_tend_adv,Nlay,Nx,Ny);
    MATALLOC3(h_tend_relax,Nlay,Nx,Ny);
  }
  
  // For averaging energy budget diagnostics
  if (dt_avg_e > 0)
  {
    MATALLOC3(e_flux_uP,Nlay,Nx,Ny);
    MATALLOC3(e_flux_uKE,Nlay,Nx,Ny);
    MATALLOC3(e_flux_uPE,Nlay,Nx,Ny);
    MATALLOC3(e_flux_vP,Nlay,Nx,Ny);
    MATALLOC3(e_flux_vKE,Nlay,Nx,Ny);
    MATALLOC3(e_flux_vPE,Nlay,Nx,Ny);
    MATALLOC3(e_tend_adv,Nlay,Nx,Ny);
    MATALLOC3(e_tend_gradM,Nlay,Nx,Ny);
    MATALLOC3(e_tend_wind,Nlay,Nx,Ny);
    MATALLOC3(e_tend_rDrag,Nlay,Nx,Ny);
    MATALLOC3(e_tend_rSurf,Nlay,Nx,Ny);
    MATALLOC3(e_tend_CdBot,Nlay,Nx,Ny);
    MATALLOC3(e_tend_CdSurf,Nlay,Nx,Ny);
    MATALLOC3(e_tend_A2,Nlay,Nx,Ny);
    MATALLOC3(e_tend_A4,Nlay,Nx,Ny);
    MATALLOC3(e_tend_wdiaPE,Nlay,Nx,Ny);
    MATALLOC3(e_tend_wdiaKE,Nlay,Nx,Ny);
    MATALLOC3(e_tend_Fbaro,Nlay,Nx,Ny);
    MATALLOC3(e_tend_buoy,Nlay,Nx,Ny);
    MATALLOC3(e_tend_rand,Nlay,Nx,Ny);
    MATALLOC3(e_tend_relax,Nlay,Nx,Ny);
  }
  
  // Work arrays for time derivative function
  MATALLOC(qq,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(hh_q,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(zeta,Nx+2*Ng,Ny+2*Ng);
  MATALLOC3(h_west,Nlay,Nx+2*Ng,Ny+2*Ng);
  MATALLOC3(h_south,Nlay,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(d2hx,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(d2hy,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(dh_dx,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(dh_dy,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(huu,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(hvv,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(alpha_w,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(beta_w,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(gamma_w,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(delta_w,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(epsilon_w,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(phi_w,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(MM_B,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(KE_B,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(pp,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(Omega_z_w,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(hhs_w,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(hhb_w,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(A2_q,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(A2_h,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(A4sqrt_q,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(A4sqrt_h,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(DD_T,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(DD_S,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(DD4_T,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(DD4_S,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(UU4,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(VV4,Nx+2*Ng,Ny+2*Ng);
  MATALLOC3(eta_w,Nlay+1,Nx+2*Ng,Ny+2*Ng);

  // Arrays required for surface/bottom momentum forcing calculation
  MATALLOC(usq_surf,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(vsq_surf,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(usq_bot,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(vsq_bot,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(uabs_surf,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(uabs_bot,Nx+2*Ng,Ny+2*Ng);
  MATALLOC3(hFsurf_west,Nlay,Nx+2*Ng,Ny+2*Ng);
  MATALLOC3(hFsurf_south,Nlay,Nx+2*Ng,Ny+2*Ng);
  MATALLOC3(hFbot_west,Nlay,Nx+2*Ng,Ny+2*Ng);
  MATALLOC3(hFbot_south,Nlay,Nx+2*Ng,Ny+2*Ng);
  MATALLOC3(eta_west,Nlay+1,Nx+2*Ng,Ny+2*Ng);
  MATALLOC3(eta_south,Nlay+1,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(hhs_west,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(hhs_south,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(hhb_west,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(hhb_south,Nx+2*Ng,Ny+2*Ng);
  
  // For tracer tendency calculation
  MATALLOC(udb_dx,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(vdb_dy,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(db_dx,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(db_dy,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(d2bx,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(d2by,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(hub,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(hvb,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(BB4,Nx+2*Ng,Ny+2*Ng);
  MATALLOC(hz,Nx+2*Ng,Ny+2*Ng);
  
  // For calculating E and Z in calcEZ
  MATALLOC3(KEdens,Nlay,Nx,Ny);
  MATALLOC3(PEdens,Nlay,Nx,Ny);
  MATALLOC(Zdens,Nx,Ny);
  
  // Required for surface pressure calculations
  MATALLOC(pi,Nx,Ny);
  MATALLOC(pi_buf,Nx,Ny);
  MATALLOC(pi_rhs,Nx,Ny);
  MATALLOC(Hc,Nx,Ny);
  MATALLOC(Hw,Nx,Ny);
  MATALLOC(Hs,Nx,Ny);
  MATALLOC(Ow,Nx,Ny);
  MATALLOC(Os,Nx,Ny);
  MATALLOC(Osum,Nx,Ny);
  MATALLOC(_Osum,Nx,Ny);
  VECALLOC(pi_prev,Ny);
  VECALLOC(im1_vec,Nx);
  VECALLOC(ip1_vec,Nx);
  VECALLOC(jm1_vec,Ny);
  VECALLOC(jp1_vec,Ny);
  
  // Random forcing
  if (useRandomForcing)
  {
    MATALLOC3(RFu,Nlay,Nx,Ny);
    MATALLOC3(RFv,Nlay,Nx,Ny);
    MATALLOC3(RFmask_fft,Nlay,Nx,Ny);
    MATALLOC3(RFmask_rot,Nlay,Nx,Ny);
    MATALLOC3(RFmask_div,Nlay,Nx,Ny);
#ifdef ALLOW_FFTW
    RFrot = malloc(Nlay*Nx*Ny*sizeof(fftw_complex));
    RFdiv = malloc(Nlay*Nx*Ny*sizeof(fftw_complex));
    RFrot_fft = malloc(Nlay*Nx*Ny*sizeof(fftw_complex));
    RFdiv_fft = malloc(Nlay*Nx*Ny*sizeof(fftw_complex));
#endif
    MATALLOC3(RFamp_rot_fft,Nlay,Nx,Ny);
    MATALLOC3(RFamp_div_fft,Nlay,Nx,Ny);
    MATALLOC(RF_kk_recip,Nx,Ny);
  }
  
  // Relaxation matrices
  if (useRelax)
  {
    MATALLOC3(uRelax,Nlay,Nx,Ny);
    MATALLOC3(vRelax,Nlay,Nx,Ny);
    MATALLOC3(hRelax,Nlay,Nx,Ny);
    MATALLOC3(eRelax,Nlay,Nx,Ny);
    MATALLOC3(uTime,Nlay,Nx,Ny);
    MATALLOC3(vTime,Nlay,Nx,Ny);
    MATALLOC(hTime,Nx,Ny);
    MATALLOC3(eTime,Nlay,Nx,Ny);
    if (useTracer)
    {
      MATALLOC3(bRelax,Nlay,Nx,Ny);
      MATALLOC3(bTime,Nlay,Nx,Ny);
    }
  }
  
  // Diabatic velocity
  MATALLOC3(wdia,Nlay+1,Nx,Ny);
  if (useWDia)
  {
    MATALLOC3(wdia_u,Nlay+1,Nx,Ny);
    MATALLOC3(wdia_v,Nlay+1,Nx,Ny);
    MATALLOC4(wdia_ff,wDiaNrecs,Nlay+1,Nx,Ny);
  }
  
  // For MultiGrid solver
  VECALLOC(F_MG,F_len+1);
  mg_grids = (data_MG *) malloc(Ngrids*sizeof(data_MG));
  for (m = Ngrids-1; m >=0; m --)
  {
    mg_grids[m].Nx = m==Ngrids-1 ? Nx : mg_grids[m+1].Nx/2;
    mg_grids[m].Ny = m==Ngrids-1 ? Ny : mg_grids[m+1].Ny/2;
    mg_grids[m].dx = Lx/mg_grids[m].Nx;
    mg_grids[m].dy = Ly/mg_grids[m].Ny;
    MATALLOC(mg_grids[m].pi,mg_grids[m].Nx,mg_grids[m].Ny);     // These two memory allocations actually won't be used for m=Ngrids-1 because
    MATALLOC(mg_grids[m].pi_rhs,mg_grids[m].Nx,mg_grids[m].Ny); // we'll instead assign these pointers to the code's main pi/pi_rhs matrices
    MATALLOC(mg_grids[m].pi_temp,mg_grids[m].Nx,mg_grids[m].Ny);
    MATALLOC(mg_grids[m].Hc,mg_grids[m].Nx,mg_grids[m].Ny);
    MATALLOC(mg_grids[m].Hw,mg_grids[m].Nx,mg_grids[m].Ny);
    MATALLOC(mg_grids[m].Hs,mg_grids[m].Nx,mg_grids[m].Ny);
    MATALLOC(mg_grids[m].Ow,mg_grids[m].Nx,mg_grids[m].Ny);
    MATALLOC(mg_grids[m].Os,mg_grids[m].Nx,mg_grids[m].Ny);
    MATALLOC(mg_grids[m].Osum,mg_grids[m].Nx,mg_grids[m].Ny);
    MATALLOC(mg_grids[m]._Osum,mg_grids[m].Nx,mg_grids[m].Ny);
    MATALLOC(mg_grids[m].wmm,mg_grids[m].Nx,mg_grids[m].Ny);
    MATALLOC(mg_grids[m].wpm,mg_grids[m].Nx,mg_grids[m].Ny);
    MATALLOC(mg_grids[m].wmp,mg_grids[m].Nx,mg_grids[m].Ny);
    MATALLOC(mg_grids[m].wpp,mg_grids[m].Nx,mg_grids[m].Ny);
    VECALLOC(mg_grids[m].im1_vec,mg_grids[m].Nx);
    VECALLOC(mg_grids[m].ip1_vec,mg_grids[m].Nx);
    VECALLOC(mg_grids[m].jm1_vec,mg_grids[m].Ny);
    VECALLOC(mg_grids[m].jp1_vec,mg_grids[m].Ny);
  }
  
  /////////////////////////////////
  ///// END MEMORY ALLOCATION /////
  /////////////////////////////////
  
  

  
  /////////////////////////////////////
  ///// BEGIN PARAMETER DEFAULTS  /////
  /////////////////////////////////////
  
#pragma parallel
  
  // Horizontal rotation components
  for (i = 0; i < Nx; i ++)
  {
    for (j = 0; j < Ny; j ++)
    {
      Omega_x[i][j] = Omega_x0;
      Omega_y[i][j] = Omega_y0;
    }
  }
  
#pragma parallel
  
  // Vertical rotation component
  for (i = 0; i < Nx+1; i ++)
  {
    for (j = 0; j < Ny+1; j ++)
    {
      Omega_z[i][j] = Omega_z0;
    }
  }

#pragma parallel
  
  // Default surface topography elevation
  for (i = 0; i < Nx; i ++)
  {
    for (j = 0; j < Ny; j ++)
    {
      hhs[i][j] = hhs0;
    }
  }
  
#pragma parallel
  
  // Default bottom topography elevation
  for (i = 0; i < Nx; i ++)
  {
    for (j = 0; j < Ny; j ++)
    {
      hhb[i][j] = hhb0;
    }
  }
  
  // Default reduced gravities
  for (k = 0; k < Nlay; k ++)
  {
    gg[k] = g0;
  }
  
  // Default initial state
  for (k = 0; k < Nlay; k ++)
  {
    
#pragma parallel
    
    for (i = 0; i < Nx; i ++)
    {
      for (j = 0; j < Ny; j ++)
      {
        uu[k][i][j] = 0;
        vv[k][i][j] = 0;
        hh[k][i][j] = -hhb0/Nlay; // Divide ocean depth into Nlay equal layers
        if (useTracer)
        {
          bb[k][i][j] = 0;
        }
      }
    }
    
  }
  
  // Default relaxation targets
  if (useRelax)
  {
    for (k = 0; k < Nlay; k ++)
    {
      
  #pragma parallel
      
      for (i = 0; i < Nx; i ++)
      {
        for (j = 0; j < Ny; j ++)
        {
          uRelax[k][i][j] = 0;
          vRelax[k][i][j] = 0;
          hRelax[k][i][j] = -hhb0/Nlay; // Divide ocean depth into Nlay equal layers
          eRelax[k][i][j] = -k*hhb0/Nlay;
          if (useTracer)
          {
            bRelax[k][i][j] = 0;
          }
        }
      }
      
    }
  }
  
  // Default fixed diabatic velocities
  if (useWDia)
  {
    for (m = 0; m < wDiaNrecs; m ++)
    {
      for (k = 0; k < Nlay+1; k ++)
      {
      
#pragma parallel
      
        for (i = 0; i < Nx; i ++)
        {
          for (j = 0; j < Ny; j ++)
          {
            wdia_ff[m][k][i][j] = 0;
          }
        }
        
      }
    }
  }
  
  // Default time scale is -1 --- corresponds to no relaxation anywhere
  if (useRelax)
  {
    
#pragma parallel

    for (i = 0; i < Nx; i ++)
    {
      for (j = 0; j < Ny; j ++)
      {
        
        hTime[i][j] = -1;
        
        for (k = 0; k < Nlay; k ++)
        {
          uTime[k][i][j] = -1;
          vTime[k][i][j] = -1;
          eTime[k][i][j] = -1;
          if (useTracer)
          {
            bTime[k][i][j] = -1;
          }
        }
        
      }
      
    }
  }
  
#pragma parallel
  
  // Default lid velocity
  for (i = 0; i < Nx; i ++)
  {
    for (j = 0; j < Ny; j ++)
    {
      uLid[i][j] = 0;
      vLid[i][j] = 0;
    }
  }
  
#pragma parallel
  
  // Default barotropic forcing
  for (i = 0; i < Nx; i ++)
  {
    for (j = 0; j < Ny; j ++)
    {
      Fbaro_x[i][j] = 0;
      Fbaro_y[i][j] = 0;
    }
  }
  
  // Default parameters for random forcing
  if (useRandomForcing)
  {
     for (k = 0; k < Nlay; k ++)
     {

#pragma parallel
       
       for (i = 0; i < Nx; i ++)
       {
         for (j = 0; j < Ny; j ++)
         {
           RFu[k][i][j] = 0; // Initialize with zero forcing
           RFv[k][i][j] = 0;
           RFmask_fft[k][i][j] = 1; // Uniform spectral mask
           RFmask_rot[k][i][j] = 1; // Uniform spatial mask for the rotational component
           RFmask_div[k][i][j] = 1; // Uniform spatial mask for the divergent component
         }
       }
     }
  }

  ///////////////////////////////////
  ///// END PARAMETER DEFAULTS  /////
  ///////////////////////////////////
  
  
  
  ////////////////////////////////////////
  ///// BEGIN READING PARAMETER DATA /////
  ////////////////////////////////////////
    
  // Read input matrices and vectors
  if ( ( (strlen(hsFile) > 0)          &&  !readMatrix(hsFile,hhs,Nx,Ny,stderr) ) ||
       ( (strlen(hbFile) > 0)          &&  !readMatrix(hbFile,hhb,Nx,Ny,stderr) ) ||
       ( (strlen(OmegaxFile) > 0)      &&  !readMatrix(OmegaxFile,Omega_x,Nx,Ny,stderr) ) ||
       ( (strlen(OmegayFile) > 0)      &&  !readMatrix(OmegayFile,Omega_y,Nx,Ny,stderr) ) ||
       ( (strlen(OmegazFile) > 0)      &&  !readMatrix(OmegazFile,Omega_z,Nx+1,Ny+1,stderr) ) ||
       ( (strlen(gFile) > 0)           &&  !readVector(gFile,gg,Nlay,stderr) ) ||
       ( (strlen(tauxFile) > 0)        &&  !readMatrix3(tauxFile,taux,tauNrecs,Nx,Ny,stderr) ) ||
       ( (strlen(tauyFile) > 0)        &&  !readMatrix3(tauyFile,tauy,tauNrecs,Nx,Ny,stderr) ) ||
       ( (strlen(wDiaFile) > 0)        &&  !readMatrix4(wDiaFile,wdia_ff,wDiaNrecs,Nlay+1,Nx,Ny,stderr) ) ||
       ( (strlen(uLidFile) > 0)        &&  !readMatrix(uLidFile,uLid,Nx,Ny,stderr) ) ||
       ( (strlen(vLidFile) > 0)        &&  !readMatrix(vLidFile,vLid,Nx,Ny,stderr) ) ||
       ( (strlen(FbaroXFile) > 0)      &&  !readMatrix(FbaroXFile,Fbaro_x,Nx,Ny,stderr) ) ||
       ( (strlen(FbaroYFile) > 0)      &&  !readMatrix(FbaroYFile,Fbaro_y,Nx,Ny,stderr) ) ||
       ( (strlen(uRelaxFile) > 0)      &&  !readMatrix3(uRelaxFile,uRelax,Nlay,Nx,Ny,stderr) ) ||
       ( (strlen(vRelaxFile) > 0)      &&  !readMatrix3(vRelaxFile,vRelax,Nlay,Nx,Ny,stderr) ) ||
       ( (strlen(hRelaxFile) > 0)      &&  !readMatrix3(hRelaxFile,hRelax,Nlay,Nx,Ny,stderr) ) ||
       ( (strlen(eRelaxFile) > 0)      &&  !readMatrix3(eRelaxFile,eRelax,Nlay,Nx,Ny,stderr) ) ||
       ( useTracer && (strlen(bRelaxFile) > 0)      &&  !readMatrix3(bRelaxFile,bRelax,Nlay,Nx,Ny,stderr) ) ||
       ( (strlen(uTimeFile) > 0)       &&  !readMatrix3(uTimeFile,uTime,Nlay,Nx,Ny,stderr) ) ||
       ( (strlen(vTimeFile) > 0)       &&  !readMatrix3(vTimeFile,vTime,Nlay,Nx,Ny,stderr) ) ||
       ( (strlen(hTimeFile) > 0)       &&  !readMatrix(hTimeFile,hTime,Nx,Ny,stderr) )  ||
       ( (strlen(eTimeFile) > 0)       &&  !readMatrix3(eTimeFile,eTime,Nlay,Nx,Ny,stderr) )  ||
       ( useTracer && (strlen(bTimeFile) > 0)       &&  !readMatrix3(bTimeFile,bTime,Nlay,Nx,Ny,stderr) ) ||
       ( useRandomForcing && (strlen(RF_fftMaskFile) > 0)   &&  !readMatrix3(RF_fftMaskFile,RFmask_fft,Nlay,Nx,Ny,stderr) ) ||
       ( useRandomForcing && (strlen(RF_rotMaskFile) > 0)   &&  !readMatrix3(RF_rotMaskFile,RFmask_rot,Nlay,Nx,Ny,stderr) ) ||
       ( useRandomForcing && (strlen(RF_divMaskFile) > 0)   &&  !readMatrix3(RF_divMaskFile,RFmask_div,Nlay,Nx,Ny,stderr) )   )
  {
    printUsage();
    return 0;
  }
  
  // If a restart is specified then read from the output of a previous run
  if (restart)
  {
    // Check initial state
    for (k = 0; k < Nlay; k ++)
    {
      // Overwrite xInitFile parameters as we won't be using them
      constructOutputName (outdir,VARID_U,k,n0,uInitFile);
      constructOutputName (outdir,VARID_V,k,n0,vInitFile);
      constructOutputName (outdir,VARID_H,k,n0,hInitFile);
      if (useTracer)
      {
        constructOutputName (outdir,VARID_B,k,n0,bInitFile);
      }
      if ( !readMatrix(uInitFile,uu[k],Nx,Ny,stderr) ||
           !readMatrix(vInitFile,vv[k],Nx,Ny,stderr) ||
           !readMatrix(hInitFile,hh[k],Nx,Ny,stderr) ||
           (useTracer && !readMatrix(bInitFile,bb[k],Nx,Ny,stderr)) )
      {
        fprintf(stderr,"Unable to read pickup files.\n");
        printUsage();
        return 0;
      }
    }
  }
  // Read from initialization files
  else
  {
    if ( ( (strlen(uInitFile) > 0)       &&  !readMatrix3(uInitFile,uu,Nlay,Nx,Ny,stderr) ) ||
         ( (strlen(vInitFile) > 0)       &&  !readMatrix3(vInitFile,vv,Nlay,Nx,Ny,stderr) ) ||
         ( (strlen(hInitFile) > 0)       &&  !readMatrix3(hInitFile,hh,Nlay,Nx,Ny,stderr) ) ||
         ( useTracer && (strlen(bInitFile) > 0)       &&  !readMatrix3(bInitFile,bb,Nlay,Nx,Ny,stderr) ) )
    {
      fprintf(stderr,"Unable to read initialization files.\n");
      printUsage();
      return 0;
    }
  }
  
  

  
  // Check initial state
  for (k = 0; k < Nlay; k ++)
  {
    for (i = 0; i < Nx; i ++)
    {
      for (j = 0; j < Ny; j ++)
      {
        if (hh[k][i][j] <= 0)
        {
          fprintf(stderr,"hInitFile may contain only values that are >0.");
          printUsage();
          return 0;
        }
      }
    }
  }
  
  // Rigid-lid only checks
  if (useRL)
  {

#pragma parallel
    
    // Check water column thickness is > 0 if in rigid lid mode
    for (i = 0; i < Nx; i ++)
    {
      for (j = 0; j < Ny; j ++)
      {
        if (hhb[i][j] >= hhs[i][j])
        {
          fprintf(stderr,"Bottom elevation must always be lower than surface elevation.");
          printUsage();
          return 0;
        }
      }
    }
    
  } // End useRL
  
  // Reduced gravities must be positive
  for (k = 0; k < Nlay; k ++)
  {
    if (gg[k] < 0)
    {
      fprintf(stderr,"gFile may contain only values that are >0.");
      printUsage();
      return 0;
    }
  }
  
  // Enforce v=0 at southern boundary if using north/south walls
  if (useWallNS)
  {
    for (k = 0; k < Nlay; k ++)
    {
      for (i = 0; i < Nx; i ++)
      {
        vv[k][i][0] = 0;
      }
    }
  }
  
  // Enforce u=0 at western boundary if using east/west walls
  if (useWallEW)
  {
    for (k = 0; k < Nlay; k ++)
    {
      for (j = 0; j < Ny; j ++)
      {
        uu[k][0][j] = 0;
      }
    }
  }
  
  // Enforce vLid=0 at southern boundary if using north/south walls
  if (useWallNS)
  {
    for (i = 0; i < Nx; i ++)
    {
      vLid[i][0] = 0;
    }
 
  }
  
  // Enforce uLid=0 at western boundary if using east/west walls
  if (useWallEW)
  {
    for (j = 0; j < Ny; j ++)
    {
      uLid[0][j] = 0;
    }
  }
  
  //////////////////////////////////////
  ///// END READING PARAMETER DATA /////
  //////////////////////////////////////



  /////////////////////////////
  ///// BEGIN WORK ARRAYS /////
  /////////////////////////////

#pragma parallel
  
  // Copy to work arrays, leaving room for ghost points
  for (i = 0; i < Nx; i ++)
  {
    for (j = 0; j < Ny; j ++)
    {
      hhb_w[Ng+i][Ng+j] = hhb[i][j];
      Omega_z_w[Ng+i][Ng+j] = Omega_z[i][j];
    }
  }
  
  // If using north/south walls then need Omega_z on the northern wall
  if (useWallNS)
  {
    for (i = 0; i < Nx; i ++)
    {
      Omega_z_w[Ng+i][Ng+Ny] = Omega_z[i][Ny];
    }
  }
  
  // If using east/west walls then need Omega_z on the eastern wall
  if (useWallEW)
  {
    for (j = 0; j < Ny; j ++)
    {
      Omega_z_w[Ng+Nx][Ng+j] = Omega_z[Nx][j];
    }
  }
  
  // Only need northeast corner if we're using north/south and east/west walls
  if (useWallNS && useWallEW)
  {
    Omega_z_w[Nx+Ng][Ng+Ny] = Omega_z[Nx][Ny];
  }

#pragma parallel
  
  // Set y-ghost points
  if (useWallNS)
  {
    // Reflection of bottom topography and nearest-neighbor extrapolation of Omega_z
    for (i = 0; i < Nx+2*Ng; i ++)
    {
      for (j = 0; j < Ng; j ++)
      {
        hhb_w[i][j] = hhb_w[i][2*Ng-j-1];
        hhb_w[i][Ny+Ng+j] = hhb_w[i][Ny+Ng-j-1];
        Omega_z_w[i][j] = Omega_z_w[i][Ng];
        Omega_z_w[i][Ny+Ng+j] = Omega_z_w[i][Ny+Ng];
      }
    }
  }
  else
  {
    // Periodic BC
    for (i = 0; i < Nx+2*Ng; i ++)
    {
      for (j = 0; j < Ng; j ++)
      {
        hhb_w[i][j] = hhb_w[i][Ny+j];
        hhb_w[i][Ny+Ng+j] = hhb_w[i][Ng+j];
        Omega_z_w[i][j] = Omega_z_w[i][Ny+j];
        Omega_z_w[i][Ny+Ng+j] = Omega_z_w[i][Ng+j];
      }
    }
  }
  
#pragma parallel
  
  // Set x-ghost points
  if (useWallEW)
  {
    // Reflection of bottom topography and nearest-neighbor extrapolation of Omega_z
    for (i = 0; i < Ng; i ++)
    {
      for (j = 0; j < Ny+2*Ng; j ++)
      {
        hhb_w[i][j] = hhb_w[2*Ng-i-1][j];
        hhb_w[Nx+Ng+i][j] = hhb_w[Nx+Ng-i-1][j];
        Omega_z_w[i][j] = Omega_z_w[Ng][j];
        Omega_z_w[Nx+Ng+i][j] = Omega_z_w[Nx+Ng][j];
      }
    }
  }
  else
  {
    // Periodic BC
    for (i = 0; i < Ng; i ++)
    {
      for (j = 0; j < Ny+2*Ng; j ++)
      {
        hhb_w[i][j] = hhb_w[Nx+i][j];
        hhb_w[Nx+Ng+i][j] = hhb_w[Ng+i][j];
        Omega_z_w[i][j] = Omega_z_w[Nx+i][j];
        Omega_z_w[Nx+Ng+i][j] = Omega_z_w[Ng+i][j];
      }
    }
  }
  
  // Calculate surface/bottom topography on cell faces
  calcFaceThickness(hhs_w,hhs_west,hhs_south,true,Nx,Ny,NULL,NULL);
  calcFaceThickness(hhb_w,hhb_west,hhb_south,true,Nx,Ny,NULL,NULL);
  
  // Default "forcing" layer thicknesses - fraction of the layer thickness that lies within the
  // "surface mixed layer" and "bottom boundary layer", and so will receive momentum forcing.
  // Default is that only the uppermost (lowermost) layer receives surface (bottom) momentum forcing.
  for (k = 0; k < Nlay; k ++)
  {
    
#pragma parallel
    
    for (i = 0; i < Nx+2*Ng; i ++)
    {
      for (j = 0; j < Ny+2*Ng; j ++)
      {

        if (k == 0)
        {
          hFsurf_west[k][i][j] = 1;
          hFsurf_south[k][i][j] = 1;
        }
        else
        {
          hFsurf_west[k][i][j] = 0;
          hFsurf_south[k][i][j] = 0;
        }
        
        if (k == Nlay-1)
        {
          hFbot_west[k][i][j] = 1;
          hFbot_south[k][i][j] = 1;
        }
        else
        {
          hFbot_west[k][i][j] = 0;
          hFbot_south[k][i][j] = 0;
        }
      }
    }
  }
  
#pragma parallel
  
  // Compute water column thickness on cell centers
  for (i = 0; i < Nx; i ++)
  {
    for (j = 0; j < Ny; j ++)
    {
      Hc[i][j] = hhs[i][j]-hhb[i][j];
    }
  }
  
  // Compute water column thickness on cell faces
  calcFaceThickness(Hc,Hw,Hs,false,Nx,Ny,NULL,NULL);

#pragma parallel
  
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
  
#pragma parallel
  
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
  
#pragma parallel
  
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
  
#pragma parallel
  
  // Extend uLid and vLid to include ghost points
  for (i = 0; i < Nx; i ++)
  {
    for (j = 0; j < Ny; j ++)
    {
      uLid_g[i+Ng][j+Ng] = uLid[i][j];
      vLid_g[i+Ng][j+Ng] = vLid[i][j];
    }
  }
  
#pragma parallel
  
  // Extend uLid and vLid in the y-direction.
  // Loop over just i-indices within the computational domain.
  for (i = Ng; i < Ng+Nx; i ++)
  {
    
    // With N/S walls
    if (useWallNS)
    {
      // Explicitly set northern boundary velocities to zero
      vLid_g[i][Ng] = 0;
      vLid_g[i][Ng+Ny] = 0;
      
      // Set ghost points
      for (j = 0; j < Ng; j ++)
      {
        // Velocities are just reflected at the wall.
        uLid_g[i][j] = uLid_g[i][2*Ng-j-1];
        vLid_g[i][j] = -vLid_g[i][2*Ng-j];
        uLid_g[i][Ny+Ng+j] = uLid_g[i][Ny+Ng-j-1];
        vLid_g[i][Ny+Ng+j] = -vLid_g[i][Ny+Ng-j];
      }
    } // With N/S walls
    
    // Without N/S walls
    else
    {
      // Set ghost points
      for (j = 0; j < Ng; j ++)
      {
        // Velocities set periodically
        uLid_g[i][j] = uLid_g[i][Ny+j];
        vLid_g[i][j] = vLid_g[i][Ny+j];
        uLid_g[i][Ny+Ng+j] = uLid_g[i][Ng+j];
        vLid_g[i][Ny+Ng+j] = vLid_g[i][Ng+j];
      }
    } // Without N/S walls
    
  }
  
#pragma parallel
  
  // Extend uLid and vLid in the x-direction
  // Loop over ALL j-points, including ghost points.
  for (j = 0; j < Ny+2*Ng; j ++)
  {
    
    // With E/W walls
    if (useWallEW)
    {
      // Set velocities to zero explicitly on the walls
      uLid_g[Ng][j] = 0;
      uLid_g[Nx+Ng][j] = 0;
      
      // Set ghost points
      for (i = 0; i < Ng; i ++)
      {
        uLid_g[i][j] = -uLid_g[2*Ng-i][j];
        vLid_g[i][j] = vLid_g[2*Ng-i-1][j];
        uLid_g[Nx+Ng+i][j] = -uLid_g[Nx+Ng-i][j];
        vLid_g[Nx+Ng+i][j] = vLid_g[Nx+Ng-i-1][j];
      }
    } // With E/W walls
    
    // Without E/W walls
    else
    {
      // Set ghost points
      for (i = 0; i < Ng; i ++)
      {
        uLid_g[i][j] = uLid_g[Nx+i][j];
        vLid_g[i][j] = vLid_g[Nx+i][j];
        uLid_g[Nx+Ng+i][j] = uLid_g[Ng+i][j];
        vLid_g[Nx+Ng+i][j] = vLid_g[Ng+i][j];
      }
    } // Without E/W walls
    
  }

  
  // Effective gravity in each layer
  geff[0] = gg[0];
  for (k = 1; k < Nlay; k ++)
  {
    geff[k] = geff[k-1] + gg[k];
  }
  
#pragma parallel
    
  // Initialize diapycnal velocity to zero - it will be set in tderiv
  for (i = 0; i < Nx; i ++)
  {
    for (j = 0; j < Ny; j ++)
    {
      for (k = 0; k < Nlay+1; k ++)
      {
        wdia[k][i][j] = 0;
      }
    }
  }
  
  ///////////////////////////
  ///// END WORK ARRAYS /////
  ///////////////////////////
  
  
  
  ///////////////////////////////////////////
  ///// BEGIN MULTIGRID SOLVER MATRICES /////
  ///////////////////////////////////////////

  if (useRL && use_MG)
  {
    // Loop down through grids (toward coarser resolution) and set up solver at each step
    for (m = Ngrids-1; m >=0; m --)
    {
      // Calculate water column thickness on this grid
      if (m == Ngrids-1)
      {
        memcpy(*(mg_grids[m].Hc),*Hc,Nx*Ny*sizeof(real));
      }
      else
      {
        restrict_MG(mg_grids[m+1].Nx,mg_grids[m+1].Ny,mg_grids[m+1].Hc,mg_grids[m].Hc);
      }
      
      // Compute water column thickness on cell faces
      calcFaceThickness(mg_grids[m].Hc,mg_grids[m].Hw,mg_grids[m].Hs,false,mg_grids[m].Nx,mg_grids[m].Ny,NULL,NULL);
      
#pragma parallel
      
      // Define north-south operators for the pressure solve
      for (i = 0; i < mg_grids[m].Nx; i ++)
      {
        for (j = 0; j < mg_grids[m].Ny; j ++)
        {
          mg_grids[m].Os[i][j] = mg_grids[m].Hs[i][j] / SQUARE(mg_grids[m].dy);
        }
        
        // No north-south periodic operator if there is a wall
        if (useWallNS)
        {
          mg_grids[m].Os[i][0] = 0;
        }
      }
      
#pragma parallel
      
      // Define east-west operators for the pressure solve
      for (j = 0; j < mg_grids[m].Ny; j ++)
      {
        for (i = 0; i < mg_grids[m].Nx; i ++)
        {
          mg_grids[m].Ow[i][j] = mg_grids[m].Hw[i][j] / SQUARE(mg_grids[m].dx);
        }
        
        // No east-west periodic operator if there is a wall
        if (useWallEW)
        {
          mg_grids[m].Ow[0][j] = 0;
        }
      }
      
#pragma parallel
      
      // Define sum of operators around each cell for computational efficiency
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
      
      // Indexing for adjacent gridpoints - to be used in pressure solve for efficiency
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
    
#pragma parallel
      
      // Define weights for interpolation to finer grids
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
  }
  
  /////////////////////////////////////////
  ///// END MULTIGRID SOLVER MATRICES /////
  /////////////////////////////////////////

  
  
  /////////////////////////////////////////////
  ///// BEGIN TIME STEPPING PRELIMINARIES /////
  /////////////////////////////////////////////

  // First we use the pressure solve to ensure that the initial velocity field is non-divergent.
  // This also preconditions pi for future pressure solves. We also optimize the relaxation parameter
  // at this point.
  if (useRL)
  {
    // First, ensure the layer thicknesses sum to the total water column depth
    correctThickness(hh_out);
    
    rp = 0.5*(rp_opt_min+rp_opt_max);
    surfPressure(uu,vv,hh,dt,pi,rp,false); // Needs to be called in order to update uu and vv
  }
  
#ifdef ALLOW_FFTW
  
  // Initialize random forcing variables
  if (useRandomForcing)
  {    
    // First, seed the random number generator
    srand(time(0));
    
    // Create FFTW plan for performing backward FFTs
    fft_backward = fftw_plan_dft_2d(Nx,Ny,RFrot_fft,RFrot,FFTW_BACKWARD,FFTW_PRESERVE_INPUT);
    
    // Ratio of time step to autocorrelation timescale
    RF_Tratio = dt/RF_tau;
    
    // Define absolue wavenumber matrix
    for (i = 0; i < Nx; i ++)
    {
      for (j = 0; j < Ny; j ++)
      {
        // See FFTW documentation for info on layout of components in wavenumber space (i.e. why does this code make sense?)
        RF_kk_recip[i][j] = 1.0 / sqrt( SQUARE(_2PI/Lx * (((int) (i+floor(Nx/2.0)))%Nx-floor(Nx/2.0))) + SQUARE(_2PI/Ly * (((int) (j+floor(Ny/2.0)))%Ny-floor(Ny/2.0))) );
      }
    }
    RF_kk_recip[0][0] = 0; // Avoids infinite mode-0 amplitude
    
    // Loop over layers
    for (k = 0; k < Nlay; k ++)
    {
        
#pragma parallel
        
      // Compute normalization coefficient for spectral mask in this layer
      RF_mask_norm = 0;
      for (i = 0; i < Nx; i ++)
      {
        for (j = 0; j < Ny; j ++)
        {
          RF_mask_norm += RFmask_fft[k][i][j];
        }
      }
      RF_mask_norm = sqrt(0.5*RF_mask_norm); // Factor of 0.5 accounts for the fact that only half of the forcing
                                            // "energy" is preserved when we transform back to real space
      
#pragma parallel
        
      // Evolve rotational and divergent components in spectral space
      for (i = 0; i < Nx; i ++)
      {
        for (j = 0; j < Ny; j ++)
        {
          // Amplitude of white noise forcing term
          RFamp_rot_fft[k][i][j] = RF_F0_rot * RF_kk_recip[i][j] * RFmask_fft[k][i][j] / RF_mask_norm;
          RFamp_div_fft[k][i][j] = RF_F0_div * RF_kk_recip[i][j] * RFmask_fft[k][i][j] / RF_mask_norm;
        }
      }
      
      // If we're using barotropic forcing then we only need to initialize the Fourier
      // components in the uppermost layer
      if ((k == 0) || (!useBarotropicRF))
      {
       
        // Initialize rotational and divergent components in spectral space
        for (i = 0; i < Nx; i ++)
        {
          for (j = 0; j < Ny; j ++)
          {
            RF_idx = k*Nx*Ny + i*Ny + j;
            
            // Generate random phase for rotational component of forcing
            RF_phase = _2PI * (rand() * 1.0 / RAND_MAX);
            // Initialize real part
            RFrot_fft[RF_idx][0] = RFamp_rot_fft[k][i][j] * cos(RF_phase);
            // Initialize complex part
            RFrot_fft[RF_idx][1] = RFamp_rot_fft[k][i][j]  * sin(RF_phase);
            
            // Generate random phase for rotational component of forcing
            RF_phase = _2PI * (rand() * 1.0 / RAND_MAX);
            // Initialize real part
            RFdiv_fft[RF_idx][0] = RFamp_div_fft[k][i][j] * cos(RF_phase);
            // Initialize complex part
            RFdiv_fft[RF_idx][1] = RFamp_div_fft[k][i][j] * sin(RF_phase);
          }
        }
        
        // Transform back to real space
        RF_idx = k*Nx*Ny;
        fftw_execute_dft(fft_backward, RFrot_fft+RF_idx, RFrot+RF_idx);
        fftw_execute_dft(fft_backward, RFdiv_fft+RF_idx, RFdiv+RF_idx);
      }
      else
      {
        // If we're using barotropic random forcing then we just copy the
        // rotational and divergent components in real space from the uppermost layer
        for (i = 0; i < Nx; i ++)
        {
          for (j = 0; j < Ny; j ++)
          {
            RF_idx = k*Nx*Ny + i*Ny + j;
            RF_idx0 = i*Ny + j;
            
            // Only copy the real parts because we don't use the imaginary ones for forcing
            RFrot[RF_idx][0] = RFrot[RF_idx0][0];
            RFdiv[RF_idx][0] = RFdiv[RF_idx0][0];
          }
        }

      } // end if ((k == 0) || (!useBarotropicRF))

    } // end k-loop
    
  } // end if (useRandomForcing)
  
#endif // ALLOW_FFTW

  // Write initial conditions if not restarting the simulation
  if (!writeModelState (t,n_saves,uu,vv,hh,bb,pi,wdia,outdir))
  {
    fprintf(stderr,"Unable to write model initial state");
    printUsage();
    return 0;
  }
  n_saves ++;
  
  // Open the E/Z diagnostic file
  if (dt_EZ > 0)
  {
    strcpy(EZfilename,outdir);
    strcat(EZfilename,"/");
    strcat(EZfilename,EZFILE);
  }
  
  // If we are calculating averaged output then initialize averaging arrays to zero
  if (dt_avg > 0)
  {
#pragma parallel
    for (i = 0; i < Nx; i ++)
    {
      for (j = 0; j < Ny; j ++)
      {
        for (k = 0; k < Nlay; k ++)
        {
          uu_avg[k][i][j] = 0;
          vv_avg[k][i][j] = 0;
          hh_avg[k][i][j] = 0;
          if (useTracer)
          {
            bb_avg[k][i][j] = 0;
          }
          hu_avg[k][i][j] = 0;
          hv_avg[k][i][j] = 0;
          huu_avg[k][i][j] = 0;
          hvv_avg[k][i][j] = 0;
          wdia_avg[k][i][j] = 0;
        }
        wdia_avg[Nlay][i][j] = 0;
        pi_avg[i][j] = 0;
      }
    }
  }
  
  // If we are calculating averaged u-momentum output then initialize averaging arrays to zero
  if (dt_avg_hu > 0)
  {
#pragma parallel
    for (i = 0; i < Nx; i ++)
    {
      for (j = 0; j < Ny; j ++)
      {
        for (k = 0; k < Nlay; k ++)
        {
          hu_tend_q[k][i][j] = 0;
          hu_tend_gradM[k][i][j] = 0;
          hu_tend_gradKE[k][i][j] = 0;
          hu_tend_dhdt[k][i][j] = 0;
          hu_tend_A2[k][i][j] = 0;
          hu_tend_A4[k][i][j] = 0;
          hu_tend_rDrag[k][i][j] = 0;
          hu_tend_rSurf[k][i][j] = 0;
          hu_tend_CdBot[k][i][j] = 0;
          hu_tend_CdSurf[k][i][j] = 0;
          hu_tend_wind[k][i][j] = 0;
          hu_tend_buoy[k][i][j] = 0;
          hu_tend_relax[k][i][j] = 0;
          hu_tend_wdia[k][i][j] = 0;
          hu_tend_rand[k][i][j] = 0;
          hu_tend_Fbaro[k][i][j] = 0;
        }
      }
    }
  }
    
  // If we are calculating averaged u-momentum output then initialize averaging arrays to zero
  if (dt_avg_hv > 0)
  {
#pragma parallel
    for (i = 0; i < Nx; i ++)
    {
      for (j = 0; j < Ny; j ++)
      {
        for (k = 0; k < Nlay; k ++)
        {
          hv_tend_q[k][i][j] = 0;
          hv_tend_gradM[k][i][j] = 0;
          hv_tend_gradKE[k][i][j] = 0;
          hv_tend_dhdt[k][i][j] = 0;
          hv_tend_A2[k][i][j] = 0;
          hv_tend_A4[k][i][j] = 0;
          hv_tend_rDrag[k][i][j] = 0;
          hv_tend_rSurf[k][i][j] = 0;
          hv_tend_CdBot[k][i][j] = 0;
          hv_tend_CdSurf[k][i][j] = 0;
          hv_tend_wind[k][i][j] = 0;
          hv_tend_buoy[k][i][j] = 0;
          hv_tend_wdia[k][i][j] = 0;
          hv_tend_rand[k][i][j] = 0;
          hv_tend_Fbaro[k][i][j] = 0;
        }
      }
    }
  }
  
  // If we are calculating averaged thickness equation output then initialize averaging arrays to zero
  if (dt_avg_h > 0)
  {
#pragma parallel
    for (i = 0; i < Nx; i ++)
    {
      for (j = 0; j < Ny; j ++)
      {
        for (k = 0; k < Nlay; k ++)
        {
          h_tend_adv[k][i][j] = 0;
          h_tend_relax[k][i][j] = 0;
        }
      }
    }
  }
  
  // If we are calculating averaged energy budget diagnostics then initialize averaging arrays to zero
  if (dt_avg_e > 0)
  {
#pragma parallel
    for (i = 0; i < Nx; i ++)
    {
      for (j = 0; j < Ny; j ++)
      {
        for (k = 0; k < Nlay; k ++)
        {
          e_flux_uP[k][i][j] = 0;
          e_flux_uKE[k][i][j] = 0;
          e_flux_uPE[k][i][j] = 0;
          e_flux_vP[k][i][j] = 0;
          e_flux_vKE[k][i][j] = 0;
          e_flux_vPE[k][i][j] = 0;
          e_tend_adv[k][i][j] = 0;
          e_tend_gradM[k][i][j] = 0;
          e_tend_wind[k][i][j] = 0;
          e_tend_rDrag[k][i][j] = 0;
          e_tend_rSurf[k][i][j] = 0;
          e_tend_CdBot[k][i][j] = 0;
          e_tend_CdSurf[k][i][j] = 0;
          e_tend_A2[k][i][j] = 0;
          e_tend_A4[k][i][j] = 0;
          e_tend_wdiaPE[k][i][j] = 0;
          e_tend_wdiaKE[k][i][j] = 0;
          e_tend_Fbaro[k][i][j] = 0;
          e_tend_buoy[k][i][j] = 0;
          e_tend_rand[k][i][j] = 0;
          e_tend_relax[k][i][j] = 0;
        }
      }
    }
  }
  
  // Write first time record to timer file
  if (tfile != NULL)
  {
    fprintf(tfile,EXP_FMT,t);
    fprintf(tfile," ");
    fflush(tfile);
  }
  
  ///////////////////////////////////////////
  ///// END TIME STEPPING PRELIMINARIES /////
  ///////////////////////////////////////////
  
  
  
  ///////////////////////////////
  ///// BEGIN TIME STEPPING /////
  ///////////////////////////////

  // Numerical time-integration loop
  for (n = 1; n <= Nt; n ++)
  {
    // This keeps the time accurate when using single precision
    t = tmin + (n-1)*dt;
    
#ifdef ALLOW_FFTW
    
    // Compute random forcing in real space
    if (useRandomForcing)
    {
      
      // Loop over layers to compute forcing in real space
      for (k = 0; k < Nlay; k ++)
      {
        
#pragma parallel
        
        // Apply mask in real space
        for (i = 0; i < Nx; i ++)
        {
          for (j = 0; j < Ny; j ++)
          {
            RF_idx = k*Nx*Ny + i*Ny + j;
            RFrot[RF_idx][0] *= RFmask_rot[k][i][j];
            RFdiv[RF_idx][0] *= RFmask_div[k][i][j];
          }
        }
        
#pragma parallel
        
        // Calculate tendencies for u and v
        for (i = 0; i < Nx; i ++)
        {
          im1 = (i + Nx - 1) % Nx;
          ip1 = (i + Nx + 1) % Nx;
          
          for (j = 0; j < Ny; j ++)
          {
            jm1 = (j + Ny - 1) % Ny;
            jp1 = (j + Ny + 1) % Ny;
            RF_idx = k*Nx*Ny + i*Ny + j;
            RF_idx_im1 = k*Nx*Ny + im1*Ny + j;
            RF_idx_ip1 = k*Nx*Ny + ip1*Ny + j;
            RF_idx_jm1 = k*Nx*Ny + i*Ny + jm1;
            RF_idx_jp1 = k*Nx*Ny + i*Ny + jp1;
            
            // u-forcing = - d/dy(Frot) + d/dx(Fdiv)
            RFu[k][i][j] = - (RFrot[RF_idx_jp1][0] - RFrot[RF_idx][0]) / dy + (RFdiv[RF_idx][0] - RFdiv[RF_idx_im1][0]) / dx;
            // v-forcing = d/dx(Frot) + d/dy(Fdiv)
            RFv[k][i][j] = (RFrot[RF_idx_ip1][0] - RFrot[RF_idx][0]) / dx + (RFdiv[RF_idx][0] - RFdiv[RF_idx_jm1][0]) / dy;
          }
        }
        
      }
      
    } // end if (useRandomForcing)
    
#endif // ALLOW_FFTW
    
    // Take one time step based on selected scheme
    if (timeSteppingScheme == TIMESTEPPING_AB1)
    {
      ab1 (&t,
           vars,
           vars_out,
           dt_vars,
           dt,
           Ntotal,
           &tderiv             );
    }
    else if (timeSteppingScheme == TIMESTEPPING_AB3)
    {
      // Generate the first couple of iterations
      // with lower-order methods
      if (n == 1)
      {
        ab1 (&t,
             vars,
             vars_out,
             dt_vars,
             dt,
             Ntotal,
             &tderiv             );
      }
      else if (n == 2)
      {
        ab2 (&t,
             vars,
             vars_out,
             dt_vars,
             dt_vars_1,
             dt,
             Ntotal,
             &tderiv             );
        
      }
      // Once we have previous iterations, use third order Adams-Bashforth
      else
      {
        ab3 (   &t,
                vars,
                vars_out,
                dt_vars,
                dt_vars_1,
                dt_vars_2,
                dt,
                Ntotal,
                &tderiv             );
      }
    }
    else
    {
      fprintf(stderr,"ERROR: Unknown time-integration method\n");
      break;
    }
    
    // If we're using a rigid lid then we need to add a correction to
    // the velocity field to account for the surface pressure
    if (useRL)
    {
      
      // First, ensure the layer thicknesses sum to the total water column depth
      correctThickness(hh_out);
      
      // Optimize after the first iteration or whenever we reach the optimization frequency.
      // Optimizing after the first iteration gets us much closer to the optimal rp than
      // using the initial model input, typically.
      if (!use_MG && ((n == n_first_optim) || ((rp_opt_freq > 0) && (n % rp_opt_freq == 0))))
      {
        rp = optimizeSOR(uu_out,vv_out,hh_out,pi,uu_buf,vv_buf,pi_buf,dt,rp);
      }
      
      surfPressure(uu_out,vv_out,hh_out,dt,pi,rp,true);
      
    } // end if (useRL)
    
#ifdef ALLOW_FFTW
    
    // Evolve the random forcing function, if used
    if (useRandomForcing)
    {
      // Evolve forcing layer by layer
      for (k = 0; k < Nlay; k ++)
      {

        // If we're using barotropic forcing then we only need to evolve the Fourier
        // components in the uppermost layer
        if ((k == 0) || (!useBarotropicRF))
        {
          
#pragma parallel
          
          // Evolve rotational and divergent components in spectral space
          for (i = 0; i < Nx; i ++)
          {
            for (j = 0; j < Ny; j ++)
            {
              RF_idx = k*Nx*Ny + i*Ny + j;
              
              // Generate random phase for rotational component of forcing
              RF_phase = _2PI * (rand() * 1.0 / RAND_MAX);
              // Evolve real part
              RFrot_fft[RF_idx][0] = RFrot_fft[RF_idx][0]
                                    - RF_Tratio * RFrot_fft[RF_idx][0] // Exponential decay term
                                    + sqrt(2*RF_Tratio) * RFamp_rot_fft[k][i][j] * cos(RF_phase); // White noise forcing term
              // Evolve complex part
              RFrot_fft[RF_idx][1] = RFrot_fft[RF_idx][1]
                                    - RF_Tratio * RFrot_fft[RF_idx][1] // Exponential decay term
                                    + sqrt(2*RF_Tratio) * RFamp_rot_fft[k][i][j]  * sin(RF_phase); // White noise forcing term
              
              // Generate random phase for rotational component of forcing
              RF_phase = _2PI * (rand() * 1.0 / RAND_MAX);
              // Evolve real part
              RFdiv_fft[RF_idx][0] = RFdiv_fft[RF_idx][0]
                                    - RF_Tratio * RFdiv_fft[RF_idx][0] // Exponential decay term
                                    + sqrt(2*RF_Tratio) * RFamp_div_fft[k][i][j] * cos(RF_phase); // White noise forcing term
              // Evolve complex part
              RFdiv_fft[RF_idx][1] = RFdiv_fft[RF_idx][1]
                                    - RF_Tratio * RFdiv_fft[RF_idx][1] // Exponential decay term
                                    + sqrt(2*RF_Tratio) * RFamp_div_fft[k][i][j] * sin(RF_phase); // White noise forcing term

            }
          }
          
          // Transform back to real space
          RF_idx = k*Nx*Ny;
          fftw_execute_dft(fft_backward, RFrot_fft+RF_idx, RFrot+RF_idx);
          fftw_execute_dft(fft_backward, RFdiv_fft+RF_idx, RFdiv+RF_idx);
        }
        else
        {

#pragma parallel
          
          // If we're using barotropic random forcing then we just copy the
          // rotational and divergent components in real space from the uppermost layer
          for (i = 0; i < Nx; i ++)
          {
            for (j = 0; j < Ny; j ++)
            {
              RF_idx = k*Nx*Ny + i*Ny + j;
              RF_idx0 = i*Ny + j;
              
              // Only copy the real parts because we don't use the imaginary ones for forcing
              RFrot[RF_idx][0] = RFrot[RF_idx0][0];
              RFdiv[RF_idx][0] = RFdiv[RF_idx0][0];
            }
          }
          
        } // end if ((k == 0) || (!useBarotropicRF))
                           
      } // end k-loop
         
    } // end if (useRandomForcing}
    
#endif // ALLOW_FFTW
    
    // If the time step has taken us past a save point (or multiple
    // save points), interpolate and write out the data
    while (t >= t_next)
    {
      // Write out the most recent model state.
      // An older version of the code used interpolation between time steps. In this version we
      // avoid doing this because it breaks conservation properties in the model state.
      // N.B. This uses the most recently-calculated data in the global variable wdia for the diapycnal velocity
      if (!writeModelState (t,n_saves,uu_out,vv_out,hh_out,bb_out,pi,wdia,outdir))
      {
        fprintf(stderr,"Unable to write model state");
        printUsage();
        return 0;
      }
      
      // Increment recorded number of output writes
      n_saves ++;
      
      // Write to time file
      if (tfile != NULL)
      {
        fprintf(tfile,EXP_FMT,t_next);
        fprintf(tfile," ");
        fflush(tfile);
      }
      
      // Update the next save time. We track this using n_saves because if only single
      // precision is used then the save points quickly lose accuracy.
      t_next = n_saves*dt_s;
    }
    
    // If we are time-averaging the model variables then add to the online average
    if (dt_avg > 0)
    {
      
      // Fraction of time step to add to current average
      avg_fac = (t < t_next_avg) ? 1 : 1-(t-t_next_avg)/dt;
      
      // Need layer thicknesses on cell faces to calculate products
      for (k = 0; k < Nlay; k ++)
      {
        calcFaceThickness(hh_out[k],h_west[k],h_south[k],false,Nx,Ny,uu[k],vv[k]);
      }
      
      // Need layer interface heights for pressure calculation
      // N.B. here we use eta_w as an Nlay x Nx x Ny matrix even though it is really
      // an Nlay x (Nx+2*Ng) x (Ny+2*Ng) matrix
      calcEta(hhb,hh_out,eta_w,Nlay,Nx,Ny);

#pragma parallel
      
      // Add to averages
      for (i = 0; i < Nx; i ++)
      {
        im1 = (i+Nx-1) % Nx;
        
        for (j = 0; j < Ny; j ++)
        {
          jm1 = (j+Ny-1) % Ny;
          
          pp[i][j] = 0;
          for (k = 0; k < Nlay; k ++)
          {
            // Calculate Montgomery potential. Here we perform this calculation cumulatively, using the pressure set
            // by the previous layer (k-1) and adding to it to get the pressure in the current layer (k).
            // N.B. Here we use pp as an Nx x Ny matrix, even though it is reall an (Nx+2Ng) x (Ny+2*Ng) matrix
            if (k == 0)
            {
              pp[i][j] = useRL ? pi[i][j] : gg[0]*eta_w[0][i][j];
            }
            else
            {
              pp[i][j] += gg[k]*eta_w[k][i][j];
            }
            
            // Calculate averaging quantities (multiplied by time step size)
            udt[k][i][j] = dt*uu_out[k][i][j];
            vdt[k][i][j] = dt*vv_out[k][i][j];
            hdt[k][i][j] = dt*hh_out[k][i][j];
            wdt[k][i][j] = dt*wdia[k][i][j];
            Mdt[k][i][j] = dt*(pp[i][j] - geff[k] * POW4(h0) / POW3(hh_out[k][i][j]) / 3);
            if (useTracer)
            {
              bdt[k][i][j] = dt*bb_out[k][i][j];
            }
            hudt[k][i][j] = dt*h_west[k][i][j]*uu_out[k][i][j];
            hvdt[k][i][j] = dt*h_south[k][i][j]*vv_out[k][i][j];
            huudt[k][i][j] = dt*h_west[k][i][j]*uu_out[k][i][j]*uu_out[k][i][j];
            hvvdt[k][i][j] = dt*h_south[k][i][j]*vv_out[k][i][j]*vv_out[k][i][j];
            huvdt[k][i][j] = dt * 0.25*(hh_out[k][i][j] + hh_out[k][im1][j] + hh_out[k][i][jm1] + hh_out[k][im1][jm1])
                                * 0.5*(uu_out[k][i][j] + uu_out[k][i][jm1])
                                * 0.5*(vv_out[k][i][j] + vv_out[k][im1][j]);
            
            // Add to running averages
            uu_avg[k][i][j] += avg_fac*udt[k][i][j];
            vv_avg[k][i][j] += avg_fac*vdt[k][i][j];
            hh_avg[k][i][j] += avg_fac*hdt[k][i][j];
            MM_avg[k][i][j] += avg_fac*Mdt[k][i][j];
            wdia_avg[k][i][j] += avg_fac*wdt[k][i][j];
            if (useTracer)
            {
              bb_avg[k][i][j] += avg_fac*bdt[k][i][j];
            }
            hu_avg[k][i][j] += avg_fac*hudt[k][i][j];
            hv_avg[k][i][j] += avg_fac*hvdt[k][i][j];
            huu_avg[k][i][j] += avg_fac*huudt[k][i][j];
            hvv_avg[k][i][j] += avg_fac*hvvdt[k][i][j];
            huv_avg[k][i][j] += avg_fac*huvdt[k][i][j];
          }
          
          // Sea floor diapycnal velocity (typically zero, except in some specific use cases)
          k = Nlay;
          wdia_avg[k][i][j] += avg_fac*wdt[k][i][j];
          
          // Surface pressure
          pidt[i][j] = dt*pi[i][j];
          pi_avg[i][j] += avg_fac*pidt[i][j];
        }
      }
      
      // If we have passed the end of the averaging window then write output files and reset buffers
      if (t >= t_next_avg)
      {
        
#pragma parallel
        
        // Divide by averaging step length to compute averages
        for (i = 0; i < Nx; i ++)
        {
          for (j = 0; j < Ny; j ++)
          {
            for (k = 0; k < Nlay; k ++)
            {
              uu_avg[k][i][j] /= dt_avg;
              vv_avg[k][i][j] /= dt_avg;
              hh_avg[k][i][j] /= dt_avg;
              MM_avg[k][i][j] /= dt_avg;
              wdia_avg[k][i][j] /= dt_avg;
              if (useTracer)
              {
                bb_avg[k][i][j] /= dt_avg;
              }
              hu_avg[k][i][j] /= dt_avg;
              hv_avg[k][i][j] /= dt_avg;
              huu_avg[k][i][j] /= dt_avg;
              hvv_avg[k][i][j] /= dt_avg;
              huv_avg[k][i][j] /= dt_avg;
            }
            
            k = Nlay;
            wdia_avg[k][i][j] /= dt_avg;
            pi_avg[i][j] /= dt_avg;
          }
        }
        
        // Write averaged variables to output files
        if (!writeAverageState (t,n_avg,uu_avg,vv_avg,hh_avg,MM_avg,bb_avg,pi_avg,wdia_avg,hu_avg,hv_avg,huu_avg,hvv_avg,huv_avg,outdir))
        {
          fprintf(stderr,"Unable to write model state");
          printUsage();
          return 0;
        }
        
#pragma parallel
        
        // Reset averaging buffers
        for (i = 0; i < Nx; i ++)
        {
          for (j = 0; j < Ny; j ++)
          {
            for (k = 0; k < Nlay; k ++)
            {
              uu_avg[k][i][j] = (1-avg_fac)*udt[k][i][j];
              vv_avg[k][i][j] = (1-avg_fac)*vdt[k][i][j];
              hh_avg[k][i][j] = (1-avg_fac)*hdt[k][i][j];
              MM_avg[k][i][j] = (1-avg_fac)*Mdt[k][i][j];
              wdia_avg[k][i][j] = (1-avg_fac)*wdt[k][i][j];
              if (useTracer)
              {
                bb_avg[k][i][j] = (1-avg_fac)*bdt[k][i][j];
              }
              hu_avg[k][i][j] = (1-avg_fac)*hudt[k][i][j];
              hv_avg[k][i][j] = (1-avg_fac)*hvdt[k][i][j];
              huu_avg[k][i][j] = (1-avg_fac)*huudt[k][i][j];
              hvv_avg[k][i][j] = (1-avg_fac)*hvvdt[k][i][j];
              huv_avg[k][i][j] = (1-avg_fac)*huvdt[k][i][j];
            }
            
            k=Nlay;
            wdia_avg[k][i][j] = (1-avg_fac)*wdt[k][i][j];
            pi_avg[i][j] = (1-avg_fac)*pidt[i][j];
          }
        }
       
        // Update the next save time. We track this using n_avg because if only single
        // precision is used then the save points quickly lose accuracy.
        n_avg ++;
        t_next_avg = n_avg*dt_avg;
      }
      
    } // end if (dt_avg > 0)
    
    // If we are time-averaging the U-momentum equation and have passed a save checkpoint
    // then write averages to output files and reset averaging buffers
    if ((dt_avg_hu > 0) && (t >= t_next_avg_hu))
    {
      // Calculate true length of averaging period
      avg_len_hu = (n-n_prev_avg_hu)*dt;
      
      // Write to files
      if (!writeUMomentumAverages(n_avg_hu,outdir))
      {
        fprintf(stderr,"Unable to write averaged terms in U momentum equation");
        printUsage();
        return 0;
      }
      
      // Reset parameters for next averaging period
      n_prev_avg_hu = n;
      n_avg_hu ++;
      t_next_avg_hu = n_avg_hu*dt_avg_hu;
    }
    
    // If we are time-averaging the V-momentum equation and have passed a save checkpoint
    // then write averages to output files and reset averaging buffers
    if ((dt_avg_hv > 0) && (t >= t_next_avg_hv))
    {
      // Calculate true length of averaging period
      avg_len_hv = (n-n_prev_avg_hv)*dt;
      
      // Write to files
      if (!writeVMomentumAverages(n_avg_hv,outdir))
      {
        fprintf(stderr,"Unable to write averaged terms in V momentum equation");
        printUsage();
        return 0;
      }
      
      // Reset parameters for next averaging period
      n_prev_avg_hv = n;
      n_avg_hv ++;
      t_next_avg_hv = n_avg_hv*dt_avg_hv;
    }
    
    // If we are time-averaging the thickness equation and have passed a save checkpoint
    // then write averages to output files and reset averaging buffers
    if ((dt_avg_h > 0) && (t >= t_next_avg_h))
    {
      // Calculate true length of averaging period
      avg_len_h = (n-n_prev_avg_h)*dt;
      
      // Write to files
      if (!writeThicknessAverages(n_avg_h,outdir))
      {
        fprintf(stderr,"Unable to write averaged terms in thickness equation");
        printUsage();
        return 0;
      }
      
      // Reset parameters for next averaging period
      n_prev_avg_h = n;
      n_avg_h ++;
      t_next_avg_h = n_avg_h*dt_avg_h;
    }
    
    // If we are time-averaging the energy budget diagnostics and have passed a save checkpoint
    // then write averages to output files and reset averaging buffers
    if ((dt_avg_e > 0) && (t >= t_next_avg_e))
    {
      // Calculate true length of averaging period
      avg_len_e = (n-n_prev_avg_e)*dt;
      
      // Write to files
      if (!writeEnergyAverages(n_avg_e,outdir))
      {
        fprintf(stderr,"Unable to write averaged terms in energy equation");
        printUsage();
        return 0;
      }
      
      // Reset parameters for next averaging period
      n_prev_avg_e = n;
      n_avg_e ++;
      t_next_avg_e = n_avg_e*dt_avg_e;
    }
    
    // If the time step has taken us past an E/Z save point (or multiple
    // save points), write out the latest E and Z. Note that we perform no
    // interpolation here, as doing so would violate conservation of E and Z.
    while ((dt_EZ > 0) && (t >= t_next_EZ))
    {
      // Calculate potential energy and potential enstrophy
      calcEZ(t,uu_out,vv_out,hh_out,bb_out,&KE,&PE,&E,Z);
      
      // Open EZ file to append new data
      if (!restart && (n_EZ == 1))
      {
        EZfile = fopen(EZfilename,"w");
      }
      else
      {
        // N.B. This condition could be reached if the run has been restarted,
        // but the EZ file was deleted. In that case the behavior of fopen is simply
        // to create a new file.
        EZfile = fopen(EZfilename,"a");
      }
      if (EZfile == NULL)
      {
        fprintf(stderr,"ERROR: Could not open E/P file name\n");
        printUsage();
        return 0;
      }
      
      // Write E/Z data to output file
      strcpy(fmt_str,"%.50");
      strcat(fmt_str,EXP_FMT_ARG);
      fprintf(EZfile,EXP_FMT,t_next_EZ);
      fprintf(EZfile," ");
      fprintf(EZfile,fmt_str,KE);
      fprintf(EZfile," ");
      fprintf(EZfile,fmt_str,PE);
      fprintf(EZfile," ");
      fprintf(EZfile,fmt_str,E);
      fprintf(EZfile," ");
      for (k = 0; k < Nlay; k ++)
      {
        fprintf(EZfile,fmt_str,Z[k]);
        fprintf(EZfile," ");
      }
      fprintf(EZfile,"\n");
      fflush(EZfile);
      
      // Close EZ file
      fclose(EZfile);

      // Update the next save time. We track this using n_EZ because if only single
      // precision is used then the save points quickly lose accuracy.
      n_EZ ++;
      t_next_EZ = n_EZ*dt_EZ;
    }
    
    // The final step is to copy the updated variables from vars_out back to vars so
    // that we can begin the next time step.
    memcpy(vars,vars_out,Ntotal*sizeof(real));
  }
  
  // Print completion time
  if (tfile != NULL)
  {
    time(&now);
    fprintf(tfile,"\nProgram completed at %s\n", ctime(&now));
    fclose(tfile);
  }
  
  /////////////////////////////
  ///// END TIME STEPPING /////
  /////////////////////////////
  
  
  
  /////////////////////////
  ///// BEGIN CLEANUP /////
  /////////////////////////
  
  // No cleanup required because the program is exiting anyway, so memory will be cleared and file handles released
  
  ///////////////////////
  ///// END CLEANUP /////
  ///////////////////////
  
  
  
  return 0;
}
