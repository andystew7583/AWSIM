/****************************************************************************************
**
**  File        ::      ab.c
**
**  Author      ::      Andrew Stewart
**
**  Description ::      Provides a series of Adams-Bashforth methods of various orders,
**                      for use in numerical integration problems.
**
****************************************************************************************/
#include "ab.h"


/****************************************************************************************
**
**  Function    ::      ab1
**
**  Purpose     ::      Performs a single iteration of the standard Adams-Bashforth first
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
**                      dxdt
**                          Storage buffer for the time derivative of x
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
void ab1 (  real *                  t,
            real *                  x,
            real *                  xout,
            real *                  dxdt,
            const real              h,
            const uint              numvars,
            DERIVATIVE_FUNCTION     f)
{
    // For looping
    uint i;

    // Calculate the derivatives of x at t
    (*f)(*t, x, dxdt, numvars);
    
#pragma parallel
    
    // Increment the dependent variable
    for (i = 0; i < numvars; i ++)
    {
        xout[i] = dxdt[i] * h + x[i];
    }

    // Increment t
    *t += h;
}


/****************************************************************************************
**
**  Function    ::      ab2
**
**  Purpose     ::      Performs a single iteration of the standard Adams-Bashforth second
**                      order method, incrementing t by the step-size h and approximating
**                      the next set of dependent variable values (x), which are updated
**                      with their new values.
**
**                      The arrays MUST all have a length at least as great as numvars
**                      and must not modified as the method proceeds. They must also all
**                      be distinct arrays, with the exception of x and xout, which may
**                      point to the same array if in-place time-stepping is desired.
**                      If these conditions are not met, the results of this procedure
**                      are undefined.
**
**  Input       ::      PARAMETERS:
**
**                      t
**                          Pointer to the independent variable value
**                      x
**                          Array of values of dependent variables
**                      xout
**                          Array into which the incremented x-values will be put
**                      dxdt
**                          Storage buffer for the time derivative of x
**                      dxdt_1
**                          Array containing the time derivative of x from the
**                          previous iteration
**                      h
**                          Step size to take
**                      numvars
**                          Specifies the size of all input arrays
**                      f
**                          Derivative function
**
**  Output      ::      The value of the independent variable (t) will be incremented by h,
**                      the time derivative of x will be stored in dxdt and dxdt_1, and
**                      the new values of the dependent variables will be put in xout.
**
****************************************************************************************/
void ab2 (  real *                  t,
            real *                  x,
            real *                  xout,
            real *                  dxdt,
            real *                  dxdt_1,
            const real              h,
            const uint              numvars,
            DERIVATIVE_FUNCTION     f)
{
    // For looping
    uint i;

    // To avoid multiple calculation
    real c0 = 3*h/2;
    real c1 = -h/2;

    // Calculate dx/dt = f(x,t)
    (*f)(*t, x, dxdt, numvars);
    
#pragma parallel

    // Calculate the result
    for (i = 0; i < numvars; i ++)
    {
        xout[i] = x[i] + c0*dxdt[i] + c1*dxdt_1[i];
    }

    // The new time derivative becomes the old time derivative
    memcpy(dxdt_1,dxdt,numvars*sizeof(real));

    // Increment t
    *t += h;
}


/****************************************************************************************
**
**  Function    ::      ab3
**
**  Purpose     ::      Performs a single iteration of the standard Adams-Bashforth third
**                      order method, incrementing t by the step-size h and approximating
**                      the next set of dependent variable values (x), which are updated
**                      with their new values.
**
**                      The arrays MUST all have a length at least as great as numvars
**                      and must not modified as the method proceeds. They must also all
**                      be distinct arrays, with the exception of x and xout, which may
**                      point to the same array if in-place time-stepping is desired.
**                      If these conditions are not met, the results of this procedure
**                      are undefined.
**
**  Input       ::      PARAMETERS:
**
**                      t
**                          Pointer to the independent variable value
**                      x
**                          Array of values of dependent variables
**                      xout
**                          Array into which the incremented x-values will be put
**                      dxdt
**                          Storage buffer for the time derivative of x
**                      dxdt_1
**                          Array containing the time derivative of x from the
**                          previous iteration
**                      dxdt_2
**                          Array containing the time derivative of x from the
**                          second-previous iteration
**                      h
**                          Step size to take
**                      numvars
**                          Specifies the size of all input arrays
**                      f
**                          Derivative function
**
**  Output      ::      The value of the independent variable (t) will be incremented by h,
**                      the time derivative of x will be stored in dxdt and dxdt_1, the
**                      content of dxdt_1 will be copied to dxdt_2, and
**                      the new values of the dependent variables will be put in xout.
**
****************************************************************************************/
void ab3 (  real *                  t,
            real *                  x,
            real *                  xout,
            real *                  dxdt,
            real *                  dxdt_1,
            real *                  dxdt_2,
            const real              h,
            const uint              numvars,
            DERIVATIVE_FUNCTION     f)
{
    // For looping
    uint i;

    // To avoid multiple calculation
    real c0 = 23*h/12;
    real c1 = -4*h/3;
    real c2 = 5*h/12;

    // Calculate dx/dt = f(x,t)
    (*f)(*t, x, dxdt, numvars);
    
#pragma parallel

    // Calculate the result
    for (i = 0; i < numvars; i ++)
    {
        xout[i] = x[i] + c0*dxdt[i] + c1*dxdt_1[i] + c2*dxdt_2[i];
    }

    // The new time derivative becomes the old time derivative
    memcpy(dxdt_2,dxdt_1,numvars*sizeof(real));
    memcpy(dxdt_1,dxdt,numvars*sizeof(real));

    // Increment t
    *t += h;
}


/****************************************************************************************
**
**  Function    ::      ab4
**
**  Purpose     ::      Performs a single iteration of the standard Adams-Bashforth fourth
**                      order method, incrementing t by the step-size h and approximating
**                      the next set of dependent variable values (x), which are updated
**                      with their new values.
**
**                      The arrays MUST all have a length at least as great as numvars
**                      and must not modified as the method proceeds. They must also all
**                      be distinct arrays, with the exception of x and xout, which may
**                      point to the same array if in-place time-stepping is desired.
**                      If these conditions are not met, the results of this procedure
**                      are undefined.
**
**  Input       ::      PARAMETERS:
**
**                      t
**                          Pointer to the independent variable value
**                      x
**                          Array of values of dependent variables
**                      xout
**                          Array into which the incremented x-values will be put
**                      dxdt
**                          Storage buffer for the time derivative of x
**                      dxdt_1
**                          Array containing the time derivative of x from the
**                          previous iteration
**                      dxdt_2
**                          Array containing the time derivative of x from the
**                          second-previous iteration
**                      dxdt_3
**                          Array containing the time derivative of x from the
**                          third-previous iteration
**                      h
**                          Step size to take
**                      numvars
**                          Specifies the size of all input arrays
**                      f
**                          Derivative function
**
**  Output      ::      The value of the independent variable (t) will be incremented by h,
**                      the time derivative of x will be stored in dxdt and dxdt_1, the
**                      content of dxdt_1 will be copied to dxdt_2, the content of
**                      dxdt_2 will be copied to dxdt_3, and the new values of the
**                      dependent variables will be put in xout.
**
****************************************************************************************/
void ab4 (  real *                  t,
            real *                  x,
            real *                  xout,
            real *                  dxdt,
            real *                  dxdt_1,
            real *                  dxdt_2,
            real *                  dxdt_3,
            const real              h,
            const uint              numvars,
            DERIVATIVE_FUNCTION     f)
{
    // For looping
    uint i;

    // To avoid multiple calculation
    real c0 = 55*h/24;
    real c1 = -59*h/24;
    real c2 = 37*h/24;
    real c3 = -3*h/8;

    // Calculate dx/dt = f(x,t)
    (*f)(*t, x, dxdt, numvars);
    
#pragma parallel
    
    // Calculate the result
    for (i = 0; i < numvars; i ++)
    {
        xout[i] = x[i] + c0*dxdt[i] + c1*dxdt_1[i] + c2*dxdt_2[i] + c3*dxdt_3[i];
    }

    // The new time derivative becomes the old time derivative
    memcpy(dxdt_3,dxdt_2,numvars*sizeof(real));
    memcpy(dxdt_2,dxdt_1,numvars*sizeof(real));
    memcpy(dxdt_1,dxdt,numvars*sizeof(real));

    // Increment t
    *t += h;
}
