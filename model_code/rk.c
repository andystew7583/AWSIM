/****************************************************************************************
**
**  File        ::      rk.c
**
**  Author      ::      Andrew Stewart
**
**  Description ::      Provides a series of Runge-Kutta methods of various orders,
**                      for use in numerical integration problems.
**	
****************************************************************************************/
#include "rk.h"

/****************************************************************************************
**
**  Function    ::      rk1
**
**  Purpose     ::      Performs a single iteration of the standard Runge-Kutta first 
**                      order method, incrementing t by the step-size h and approximating
**                      the next set of dependent variable values (x), which are updated
**                      with their new values. 
**
**                      The arrays x and xout MUST each have a length at least
**                      as great as numvars
**
**  Input       ::      PARAMETERS:
**
**                      t
**                          Pointer to the independent variable value.
**                      x
**                          Array of values of dependent variables
**                      xout
**                          Array into which the incremented x-values will be put
**                      h
**                          Step size to take
**                      numvars
**                          Specifies the size of all input arrays
**                      f
**                          Derivative function
**
**  Output      ::      The value of the independent variable (t) will be incremented by h,
**                      and the new values of the dependent variables will be put in
**                      the xout array.
**	
****************************************************************************************/
void rk1 (  real *                  t, 
            real *                  x,
            real *                  xout,
            const real              h, 
            const uint              numvars,
            DERIVATIVE_FUNCTION     f)
{
    // For looping
    uint i;

    // Calculate the derivatives of x at t
    (*f)(*t, x, xout, numvars);
    
    
#pragma parallel
    
    // Increment the dependent variable
    for (i = 0; i < numvars; i ++)
    { 
        xout[i] = xout[i] * h + x[i];
    }

    // Increment t
    *t += h;
}

/****************************************************************************************
**
**  Function    ::      rk2
**
**  Purpose     ::      Performs a single iteration of the standard Runge-Kutta second 
**                      order method, incrementing t by the step-size h and approximating
**                      the next set of dependent variable values (x), which are updated
**                      with their new values. 
**
**                      The arrays x, xout, and k MUST each have a length at least
**                      as great as numvars, and MUST all be distinct arrays that are not
**                      modified as the method proceeds. If these conditions are not met,
**                      the results of this procedure are undefined.
**
**                      k is a temporary storage buffer used
**                      by the Runge-Kutta algorithm. It is passed in as a parameter
**                      so that memory allocation need not happen every time rk2 is
**                      called - the calling code can supply the same buffer
**                      each time, greatly reducing memory allocations.
**
**  Input       ::      PARAMETERS:
**
**                      t
**                          Pointer to the independent variable value
**                      x
**                          Array of values of dependent variables
**                      xout
**                          Array into which the incremented x-values will be put
**                      k
**                          Temporary storage buffer for intermediate values
**                      h
**                          Step size to take
**                      numvars
**                          Specifies the size of all input arrays
**                      f
**                          Derivative function
**
**  Output      ::      The value of the independent variable (t) will be incremented by h,
**                      and the new values of the dependent variables will be put in
**                      the xout array.
**	
****************************************************************************************/
void rk2 (  real *                  t, 
            real *                  x,
            real *                  xout,
            real *                  k,
            const real              h, 
            const uint              numvars,
            DERIVATIVE_FUNCTION     f)
{
    // For looping
    uint i;

    // To prevent multiple calculation
    real h_2 = h / 2;

    // Calculate k1 = h*f(x,t)
    (*f)(*t, x, k, numvars);
    
#pragma parallel
    
    // Prepare the k1 values as parameters for the calculation of k2
    for (i = 0; i < numvars; i ++)
    {
        k[i] = k[i] * h_2 + x[i];
    }

    // Calculate k2
    (*f)(*t + h_2, k, xout, numvars);
    
#pragma parallel
    
    // Calculate the result
    for (i = 0; i < numvars; i ++)
    {
        xout[i] = xout[i] * h + x[i];
    }

    // Increment t
    *t += h;
}

/****************************************************************************************
**
**  Function    ::      rk4
**
**  Purpose     ::      Performs a single iteration of the standard Runge-Kutta fourth 
**                      order method, incrementing t by the step-size h and approximating
**                      the next set of dependent variable values (x), which are updated
**                      with their new values. 
**
**                      The arrays x, xout, k13 and k24 MUST each have a length at least
**                      as great as numvars, and MUST all be distinct arrays that are not
**                      modified as the method proceeds. If these conditions are not met,
**                      the results of this procedure are undefined.
**
**                      k13 and k24 are simply temporary storage buffers used
**                      by the Runge-Kutta algorithm. They are passed in as parameters
**                      so that memory allocation need not happen every time rk4 is
**                      called - the calling code can supply the same set of buffers
**                      each time, greatly reducing memory costs.
**
**  Input       ::      PARAMETERS:
**
**                      t
**                          Pointer to the independent variable value
**                      x
**                          Array of values of dependent variables
**                      xout
**                          Array into which the incremented x-values will be put
**                      k13
**                          Temporary storage buffer for the k1 and k3 values
**                      k24
**                          Temporary storage buffer for the k2 and k4 values
**                      h
**                          Step size to take
**                      numvars
**                          Specifies the size of all input arrays
**                      f
**                          Derivative function
**
**  Output      ::      The value of the independent variable (t) will be incremented by h,
**                      and the new values of the dependent variables will be put in
**                      the xout array.
**	
****************************************************************************************/
void rk4 (  real *                  t, 
            real *                  x,
            real *                  xout,
            real *                  k13,
            real *                  k24,
            const real              h, 
            const uint              numvars,
            DERIVATIVE_FUNCTION     f)
{
    // For looping
    uint    i;

    // Calculate these now to prevent repeated calculation
    real    h_6         = h / 6;
    real    h_3         = h / 3;
    real    h_2         = h / 2;

    // We must calculate k1-k4 for each x[i], but instead of using 4 buffers,
    // we recycle k1 to hold k3 and k2 to hold k4. At each step, we calculate
    // the next set of k-values and add the appropriate quantity to xtemp,
    // which will hold the new values of x at the end.

    // Calculate k1
    (*f)(*t, x, k13, numvars);
    
#pragma parallel

    // Prepare the k1 array as a parameter for generating k2
    for (i = 0; i < numvars; i ++)
    {
        xout[i]	    =       k13[i] * h_6;
        k13[i]      =       k13[i] * h_2 + x[i];
    }

    // Calculate k2
    (*f)(*t + h/2, k13, k24, numvars);
    
#pragma parallel

    // Prepare the k2 array as a parameter for generating k3
    for (i = 0; i < numvars; i ++)
    {
        xout[i]	    +=      k24[i] * h_3;
        k24[i]      =       k24[i] * h_2 + x[i];
    }

    // Calculate k3
    (*f)(*t + h/2, k24, k13, numvars);
    
#pragma parallel

    // Prepare the k3 array as a parameter for generating k4
    for (i = 0; i < numvars; i ++)
    {
        xout[i]	    +=      k13[i] * h_3;
        k13[i]      =       k13[i] * h + x[i];
    }

    // Calculate k4
    (*f)(*t + h, k13, k24, numvars);
    
#pragma parallel

    // Instead of finishing xtemp, just add its current value, plus
    // the modified k4, to the x-values and we have the result
    for (i = 0; i < numvars; i ++)
    {
        xout[i] += x[i] + (k24[i] * h_6);
    }

    // Increment t
    *t += h;
}
