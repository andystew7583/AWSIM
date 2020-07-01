// Super-header for ODE solver files
#ifndef _NSODE_H_
#define _NSODE_H_


#ifndef _REAL_DEFINED_
#define _REAL_DEFINED_
typedef double real;
#define FLT_FMT "%lf"
#define EXP_FMT "%le"
#define FLT_FMT_ARG "lf"
#define EXP_FMT_ARG "le"
#endif

#ifndef _UINT_DEFINED_
#define _UINT_DEFINED_
typedef unsigned int uint;
#endif

#ifndef _DERIVATIVE_FUNCTION_DEFINED_
#define _DERIVATIVE_FUNCTION_DEFINED_
typedef void	(*DERIVATIVE_FUNCTION) (    const real      t,
                                            const real *    x,
                                            real *          dx_dt,
                                            const uint      numvars);
#endif


#endif
