// Header file for aabwdefs.c
#ifndef _SWDEFS_H_
#define _SWDEFS_H_

#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>

#include "nsode.h"

// Define to allow inclusion of FFTW package, required for random forcing
#undef ALLOW_FFTW

// Max name length for files containing input parameter arrays
#define MAX_PARAMETER_FILENAME_LENGTH 256

// Time-integration method identifiers
#define TIMESTEPPING_RK1 0
#define TIMESTEPPING_RK2 1
#define TIMESTEPPING_RK4 2
#define TIMESTEPPING_AB1 3
#define TIMESTEPPING_AB2 4
#define TIMESTEPPING_AB3 5
#define TIMESTEPPING_AB4 6

// Momentum equation discretization scheme identifiers
#define MOMENTUM_AL81 0 // Arakawa and Lamb (1981)
#define MOMENTUM_HK83 1 // Hollingsworth and Kallberg (1983)
#define MOMENTUM_S75e 2 // Sadourny (1975) energy conserving scheme

// Thickness advection scheme identifiers
#define THICKNESS_AL81 0 // Centered differencing. Required to conserve potential enstrophy when combined with Arakawa and Lamb (1981) momentum discretization.
#define THICKNESS_HK83 1// Centered differencing. Required to conserve potential enstrophy when combined with Hollingsworth and Kallberg (1983) momentum discretization.
#define THICKNESS_UP3 2 // Third-order upwinding
#define THICKNESS_KT00 3 // Kurganov and Tadmor (2000) scheme for conservation laws

// Tracer advection scheme identifiers
#define TRACER_AL81 0 // Modified centered differencing. Required to conserve energy if the tracer is serving as a buoyancy variable - see Warneford and Dellar (2014).
#define TRACER_UP3 1 // Third-order upwinding
#define TRACER_KT00 2 // Kurganov and Tadmor (2000) scheme for conservation laws

// Handy constants
#define _PI  3.14159265358979323846
#define _2PI 6.28318530717958647692

// Efficient way to square and cube numbers
#define SQUARE(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

// Alternatives for simple powers
#define POW2(x)  ((x)*(x))
#define POW3(x)  ((x)*(x)*(x))
#define POW4(x)  ((x)*(x)*(x)*(x))
#define POW5(x)  ((x)*(x)*(x)*(x)*(x))
#define POW6(x)  ((x)*(x)*(x)*(x)*(x)*(x))
#define POW7(x)  ((x)*(x)*(x)*(x)*(x)*(x)*(x))
#define POW8(x)  ((x)*(x)*(x)*(x)*(x)*(x)*(x)*(x))
#define POW9(x)  ((x)*(x)*(x)*(x)*(x)*(x)*(x)*(x)*(x))
#define POW10(x) ((x)*(x)*(x)*(x)*(x)*(x)*(x)*(x)*(x)*(x))

#define MATALLOC3(mat3,Nx,Ny,Nz)                                \
    mat3 = matalloc3(Nx,Ny,Nz);                                 \
    if (mat3 == NULL)                                           \
    {                                                           \
      fprintf(stderr,"ERROR: Unable to allocate memory\r\n");   \
      return 0;                                                 \
    }

#define MATALLOC(mat,Nx,Ny)                                     \
    mat = matalloc(Nx,Ny);                                      \
    if (mat == NULL)                                            \
    {                                                           \
        fprintf(stderr,"ERROR: Unable to allocate memory\r\n"); \
        return 0;                                               \
    }

#define VECALLOC(vec,N)                                         \
    vec = vecalloc(N);                                          \
    if (vec == NULL)                                            \
    {                                                           \
        fprintf(stderr,"ERROR: Unable to allocate memory\r\n"); \
        return 0;                                               \
    }

#define CHECKALLOC(ptr,size)                                    \
    ptr = malloc(size);                                         \
    if (ptr == NULL)                                            \
    {                                                           \
        fprintf(stderr,"ERROR: Unable to allocate memory\r\n"); \
        return 0;                                               \
    }

// Output variable IDs
#define VARID_U 0
#define VARID_V 1
#define VARID_H 2
#define VARID_B 3
#define VARID_PI 4

#define VARID_U_AVG 10
#define VARID_V_AVG 11
#define VARID_H_AVG 12
#define VARID_M_AVG 13
#define VARID_B_AVG 14
#define VARID_PI_AVG 15
#define VARID_HU_AVG 16
#define VARID_HV_AVG 17
#define VARID_HUU_AVG 18
#define VARID_HVV_AVG 19
#define VARID_HUV_AVG 20

#define VARID_UMOM_Q 30
#define VARID_UMOM_GRADM 31
#define VARID_UMOM_GRADKE 32
#define VARID_UMOM_DHDT 33
#define VARID_UMOM_A2 34
#define VARID_UMOM_A4 35
#define VARID_UMOM_RDRAG 36
#define VARID_UMOM_RSURF 37
#define VARID_UMOM_CDBOT 38
#define VARID_UMOM_CDSURF 39
#define VARID_UMOM_WIND 40
#define VARID_UMOM_BUOY 41
#define VARID_UMOM_RELAX 42
#define VARID_UMOM_WDIA 43
#define VARID_UMOM_RAND 44

#define VARID_VMOM_Q 50
#define VARID_VMOM_GRADM 51
#define VARID_VMOM_GRADKE 52
#define VARID_VMOM_DHDT 53
#define VARID_VMOM_A2 54
#define VARID_VMOM_A4 55
#define VARID_VMOM_RDRAG 56
#define VARID_VMOM_RSURF 57
#define VARID_VMOM_CDBOT 58
#define VARID_VMOM_CDSURF 59
#define VARID_VMOM_WIND 60
#define VARID_VMOM_BUOY 61
#define VARID_VMOM_RELAX 62
#define VARID_VMOM_WDIA 63
#define VARID_VMOM_RAND 64

#define VARID_THIC_ADV 70
#define VARID_THIC_RELAX 71


// Output filenames
static const char OUTN_U[] = "U";
static const char OUTN_V[] = "V";
static const char OUTN_H[] = "H";
static const char OUTN_B[] = "B";
static const char OUTN_PI[] = "P";

// Averaged output filenames
static const char OUTN_U_AVG[] = "U_avg";
static const char OUTN_V_AVG[] = "V_avg";
static const char OUTN_H_AVG[] = "H_avg";
static const char OUTN_M_AVG[] = "M_avg";
static const char OUTN_B_AVG[] = "B_avg";
static const char OUTN_PI_AVG[] = "P_avg";
static const char OUTN_HU_AVG[] = "HU_avg";
static const char OUTN_HV_AVG[] = "HV_avg";
static const char OUTN_HUU_AVG[] = "HUU_avg";
static const char OUTN_HVV_AVG[] = "HVV_avg";

// U-momentum equation output filenames
static const char OUTN_UMOM_Q[] = "UMom_PVadvection";
static const char OUTN_UMOM_GRADM[] = "UMom_Montgomery";
static const char OUTN_UMOM_GRADKE[] = "UMom_KEgradient";
static const char OUTN_UMOM_DHDT[] = "UMom_dhdt";
static const char OUTN_UMOM_A2[] = "UMom_viscosity";
static const char OUTN_UMOM_A4[] = "UMom_hypervisc";
static const char OUTN_UMOM_RDRAG[] = "UMom_linBotDrag";
static const char OUTN_UMOM_RSURF[] = "UMom_linSurfDrag";
static const char OUTN_UMOM_CDBOT[] = "UMom_quadBotDrag";
static const char OUTN_UMOM_CDSURF[] = "UMom_quadSurfDrag";
static const char OUTN_UMOM_WIND[] = "UMom_windStress";
static const char OUTN_UMOM_BUOY[] = "UMom_buoyForce";
static const char OUTN_UMOM_RELAX[] = "UMom_relaxation";
static const char OUTN_UMOM_WDIA[] = "UMom_diapycnal";
static const char OUTN_UMOM_RAND[] = "UMom_randomForcing";

// V-momentum equation output filenames
static const char OUTN_VMOM_Q[] = "VMom_PVadvection";
static const char OUTN_VMOM_GRADM[] = "VMom_Montgomery";
static const char OUTN_VMOM_GRADKE[] = "VMom_KEgradient";
static const char OUTN_VMOM_DHDT[] = "VMom_dhdt";
static const char OUTN_VMOM_A2[] = "VMom_viscosity";
static const char OUTN_VMOM_A4[] = "VMom_hypervisc";
static const char OUTN_VMOM_RDRAG[] = "VMom_linBotDrag";
static const char OUTN_VMOM_RSURF[] = "VMom_linSurfDrag";
static const char OUTN_VMOM_CDBOT[] = "VMom_quadBotDrag";
static const char OUTN_VMOM_CDSURF[] = "VMom_quadSurfDrag";
static const char OUTN_VMOM_WIND[] = "VMom_windStress";
static const char OUTN_VMOM_BUOY[] = "VMom_buoyForce";
static const char OUTN_VMOM_RELAX[] = "VMom_relaxation";
static const char OUTN_VMOM_WDIA[] = "VMom_diapycnal";
static const char OUTN_VMOM_RAND[] = "VMom_randomForcing";

// Thickness equation output filenames
static const char OUTN_THIC_ADV[] = "Thic_advection";
static const char OUTN_THIC_RELAX[] = "Thic_relaxation";

// Time tracking file and energy/enstrophy diagnostics file
static const char TFILE[] = "time.txt";
static const char EZFILE[] = "EZdiags.dat";

// Input parameter data
typedef struct paramdata
{
    char * name;
    char * type;
    void * pvar;
    bool read;
}
paramdata;

// To define input parameters
void setParam (paramdata * params, uint idx, char * name, char * type, void * pvar, bool read);
bool readParams (char * infname, paramdata * params, uint nparams, FILE * errstrm);

// Data I/O
void printMatrix (FILE * outfile, real ** mat, uint m, uint n);
void printVector (FILE * outfile, real * vec, uint m);
bool readMatrix3 (char * fname, real *** mat3, uint Nr, uint Nc, uint Ns, FILE * errstrm);
bool readMatrix (char * fname, real ** mat, uint m, uint n, FILE * errstrm);
bool readVector (char * fname, real * vec, uint len, FILE * errstrm);

// Wrappers for allocating/freeing arrays
real *** matalloc3 (uint m1, uint m2, uint m3);
void matfree3 (real *** mat);
real ** matalloc (uint m, uint n);
void matfree (real ** mat);
real * vecalloc (uint m);
void vecfree (real * vec);
void vec2mat (real * vec, real *** pmat, uint m, uint n);
void vec2mat3 (real * vec, real **** pmat3, uint Nr, uint Nc, uint Ns);

// Thomas algorithm
void thomas (const real * a, const real * b, real * c, real * d, real * x, unsigned int n);

// Kahan summation
real kahanSum (real * vec, int N);

// Math
bool ispow2u (uint x);
uint log2u (uint x);
uint pow2u (uint x);

// Functions for limiting
real max3 (real v1, real v2, real v3);
real min3 (real v1, real v2, real v3);
real minmod (real v1, real v2, real v3);

#endif
