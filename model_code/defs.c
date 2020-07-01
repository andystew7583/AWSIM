/**
 *
 * SWdefs.c
 *
 * Andrew Stewart
 *
 * Contains a series of subroutines that are of practical use for shallow water codes.
 *
 */

#include "defs.h"


/**
 * printMatrix
 *
 * Convenience method to write a matrix of data (mat, dimensions m by n)
 * to a data file (outfile).
 *
 */
void printMatrix (FILE * outfile, real ** mat, uint m, uint n)
{
    uint i,j;

    // Write array of data
    for (j = 0; j < n; j ++)
    {
        for (i = 0; i < m; i ++)
        {
            //fprintf(outfile,EXP_FMT,mat[i][j]);
            //fprintf(outfile," ");
            fwrite(&mat[i][j],sizeof(real),1,outfile);
        }
        //fprintf(outfile,"\r\n");
    }
    //fflush(outfile);
}


/**
 * printVector
 *
 * Convenience method to write a vector of data (vec) of length m
 * to a data file (outfile).
 *
 */
void printVector (FILE * outfile, real * vec, uint m)
{
    uint i;
    
    // Write array of data
    for (i = 0; i < m; i ++)
    {
        fprintf(outfile,EXP_FMT,vec[i]);
        fprintf(outfile," ");
    }
    fprintf(outfile,"\r\n");
    fflush(outfile);
}


/**
 * readVector
 *
 * Reads 'len' real values into the vector 'vec' from the file
 * specified by 'fname'. Any error will result in 'false' being
 * returned, and an error message being printed to the specified
 * FILE/stream 'errstrm'.
 *
 */
bool readVector (char * fname, real * vec, uint len, FILE * errstrm)
{
    // Attempt to open the topography file
    FILE * pfile = NULL;
    int i = 0;    
    bool readok = true;
    
    // Basic error checking
    if (fname == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: NULL file name supplied to 'readVector'\r\n");
        }
        return false;
    }
    if (vec == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: NULL vector supplied to 'readVector'\r\n");
        }
        return false;
    }
    
    // Attempt to open the file
    pfile = fopen(fname,"r");
    
    // Check that the file exists
    if (pfile == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: Could not open %s\r\n",fname);
        }
        return false;
    }
    
    // Read in 'len' real values from the file
    for (i = 0; i < len; i ++)
    {
        //if (fscanf(pfile,FLT_FMT,&vec[i]) == EOF)
        if (fread(&vec[i],sizeof(real),1,pfile) < 1)
        {
            if (errstrm != NULL)
            {
                fprintf(errstrm,"ERROR: Could only read %d of %u values from %s\r\n",i,len,fname);
            }
            readok = false;
            break;
        }
    }
    
    // Close the file and exit
    fclose(pfile);
    return readok;
}


/**
 * readMatrix
 *
 * Reads 'm'x'n' real values into the matrix 'mat' in column-major
 * order (i.e. reads mat[1][1], mat[2][1], mat[3][1], ... , mat[2][1],
 * ... for compatibility with MatLab) from the file specified by
 * 'fname'. Any error will result in 'false' being returned, and an
 * error message being printed to the specified FILE/stream 'errstrm'.
 *
 */
bool readMatrix (char * fname, real ** mat, uint m, uint n, FILE * errstrm)
{
    // Attempt to open the topography file
    FILE * pfile = NULL;
    int i = 0; 
    int j = 0;   
    bool readok = true;
    
    // Basic error checking
    if (fname == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: NULL file name supplied to 'readMatrix'\r\n");
        }
        return false;
    }
    if (mat == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: NULL matrix supplied to 'readMatrix'\r\n");
        }
        return false;
    }
    for (i = 0; i < m; i ++)
    {
        if (mat[i] == NULL)
        {
            if (errstrm != NULL)
            {
                fprintf(errstrm,"ERROR: NULL vector supplied to 'readMatrix'\r\n");
            }
            return false;
        }
    }
    
    // Attempt to open the file
    pfile = fopen(fname,"r");
    
    // Check that the file exists
    if (pfile == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: Could not open %s\r\n",fname);
        }
        return false;
    }
    
    // Read in 'len' real values from the file
    for (j = 0; j < n; j ++)
    {
        for (i = 0; i < m; i ++)
        {
            //if (fscanf(pfile,FLT_FMT,&mat[i][j]) == EOF)
            if (fread(&mat[i][j],sizeof(real),1,pfile) < 1)
            {
                if (errstrm != NULL)
                {
                    fprintf(errstrm,"ERROR: Could only read %d of %u values from %s\r\n",j*m+i,m*n,fname);
                }
                readok = false;
                break;
            }
        }
    }
    
    // Close the file and exit
    fclose(pfile);
    return readok;
}


/**
 * readMatrix3
 *
 * Reads 'Nr'x'Nc'x'Ns' real values into the matrix 'mat3' in column-major
 * order (i.e. reads mat[1][1][1], mat[2][1][1], mat[3][1][1], ... , mat[1][2][1],
 * ... for compatibility with MatLab) from the file specified by
 * 'fname'. Any error will result in 'false' being returned, and an
 * error message being printed to the specified FILE/stream 'errstrm'.
 *
 */
bool readMatrix3 (char * fname, real *** mat3, uint Nr, uint Nc, uint Ns, FILE * errstrm)
{
  FILE * pfile = NULL;
  int i = 0;
  int j = 0;
  int k = 0;
  bool readok = true;
  
  // Basic error checking
  if (fname == NULL)
  {
    if (errstrm != NULL)
    {
      fprintf(errstrm,"ERROR: NULL file name supplied to 'readMatrix3'\r\n");
    }
    return false;
  }
  if (mat3 == NULL)
  {
    if (errstrm != NULL)
    {
      fprintf(errstrm,"ERROR: NULL matrix supplied to 'readMatrix3'\r\n");
    }
    return false;
  }
  for (i = 0; i < Nr; i ++)
  {
    if (mat3[i] == NULL)
    {
      if (errstrm != NULL)
      {
        fprintf(errstrm,"ERROR: NULL matrix supplied to 'readMatrix3'\r\n");
      }
      return false;
    }
    for (j = 0; j < Nc; j ++)
    {
      if (mat3[i][j] == NULL)
      {
        if (errstrm != NULL)
        {
          fprintf(errstrm,"ERROR: NULL vector supplied to 'readMatrix3'\r\n");
        }
        return false;
      }
    }
  }
  
  // Attempt to open the file
  pfile = fopen(fname,"r");
  
  // Check that the file exists
  if (pfile == NULL)
  {
    if (errstrm != NULL)
    {
      fprintf(errstrm,"ERROR: Could not open %s\r\n",fname);
    }
    return false;
  }
  
  // Read data from the file
  for (k = 0; k < Ns; k ++)
  {
    for (j = 0; j < Nc; j ++)
    {
      for (i = 0; i < Nr; i ++)
      {
        //if (fscanf(pfile,"%le",&mat3[i][j][k]) == EOF)
        if (fread(&mat3[i][j][k],sizeof(real),1,pfile) < 1)
        {
          if (errstrm != NULL)
          {
            fprintf(errstrm,"ERROR: Could only read %d of %u values from %s\r\n",k*Nr*Nc+j*Nr+i,Nr*Nc*Ns,fname);
          }
          readok = false;
          break;
        }
      }
    }
  }
  
  // Close the file and exit
  fclose(pfile);
  return readok;
}


/**
 *  matalloc3
 *
 *  Convenience method to allocate a real matrix with dimensions
 *  (Nr,Nc,Ns). Returns a pointer to an array of pointers that point
 *  to an array of pointers that point to the beginning of each row
 *  in the matrix.
 *
 */
real *** matalloc3 (uint Nr, uint Nc, uint Ns)
{
  real * slices = malloc(Nr*Nc*Ns*sizeof(real));
  real ** cols = malloc(Nr*Nc*sizeof(real *));
  real *** rows = malloc(Nr*sizeof(real **));
  
  uint i, j;
  
  if ((slices == NULL) || (cols == NULL) || (rows == NULL))
  {
    return NULL;
  }
  
  memset(slices,0,Nr*Nc*Ns*sizeof(real));
  
  for (i = 0; i < Nr; i ++)
  {
    rows[i] = cols + Nc*i;
    for (j = 0; j < Nc; j ++)
    {
      cols[i*Nc+j] = slices + Ns*j + Nc*Ns*i;
    }
  }
  
  return rows;
}


/**
 *  matfree3
 *
 *  Frees memory used by a matrix allocated with matalloc3.
 *
 */
void matfree3 (real *** mat3)
{
  free(*(*mat3));
  free(*mat3);
  free(mat3);
}


/**
 *  matalloc
 *
 *  Convenience method to allocate a real matrix with dimensions
 *  (m,n). Returns a pointer to an array of pointers that point to the
 *  beginning of each row in the matrix.
 *
 */
real ** matalloc (uint m, uint n)
{
    real * data = malloc(m*n*sizeof(real));
    real ** mat = malloc(m*sizeof(real *));
    uint i;

    if (data == NULL || mat == NULL)
    {
        return NULL;
    }

    memset(data,0,m*n*sizeof(real));
    mat[0] = data;

    for (i = 1; i < m; i ++)
    {
        mat[i] = mat[i-1] + n;
    }

    return mat;
}


/**
 *  matfree
 *
 *  Frees memory used by a matrix allocated with matalloc.
 *
 */
void matfree (real ** mat)
{
    free(*mat);
    free(mat);
}


/**
 *  vecalloc
 *
 *  Convenience method to allocate a real vector of length m.
 *  Returns a pointer to the first element in the allocated vector.
 *
 */
real * vecalloc (uint m)
{
    real * vec = malloc(m*sizeof(real));
    memset(vec,0,m*sizeof(real));
    return vec;
}


/**
 *  vecfree
 *
 *  Frees memory used by a vector allocated with vecalloc.
 *
 */
void vecfree (real * vec)
{
    free(vec);
}


/**
 *  vec2mat
 *
 * Creates a matrix from a vector of length m*n. If the matrix pointed
 * to by pmat (*pmat) is NULL, a vector of pointers is created that must
 * later be freed.
 */
void vec2mat (real * vec, real *** pmat, uint m, uint n)
{
    int i;

    if (*pmat == NULL)
    {
        *pmat = malloc(m*sizeof(real *));
    }

    (*pmat)[0] = vec;

    for (i = 1; i < m; i ++)
    {
        (*pmat)[i] = (*pmat)[i-1] + n;
    }
}


/**
 *  vec2mat3
 *
 * Creates a 3-dimensional matrix from a vector of length Nr*Nc*Ns. If the matrix pointed
 * to by pmat3 (*pmat3) is NULL, a vector of pointers is created that must
 * later be freed.
 */
void vec2mat3 (real * vec, real **** pmat3, uint Nr, uint Nc, uint Ns)
{
  int i,j;
  real *** rows = NULL;
  real ** cols = NULL;
  real * slices = NULL;
  
  // Change nomenclature for clarity
  slices = vec;
  
  // If the matrix doesn't exist then create one
  if (*pmat3 == NULL)
  {
    cols = malloc(Nr*Nc*sizeof(real *));
    rows = malloc(Nr*sizeof(real **));
  }
  else
  {
    rows = *pmat3;
    cols = *rows;
  }
  
  // Can't create a matrix from nothing
  if ((slices == NULL) || (cols == NULL) || (rows == NULL))
  {
    if (pmat3 != NULL)
    {
      *pmat3 = NULL;
    }
    return;
  }
  
  // Assign row and column pointers
  for (i = 0; i < Nr; i ++)
  {
    rows[i] = cols + Nc*i;
    for (j = 0; j < Nc; j ++)
    {
      cols[i*Nc+j] = slices + Ns*j + Nc*Ns*i;
    }
  }
  
  *pmat3 = rows;
}


/**
 * setParam
 *
 * Convenience method to set the values in the paramdata struct at the specified
 * index (idx) in an array (params). The method sets the parameter name (name)
 * parameter type (type, corresponding to the format string required for fscanf),
 * the pointer to the variable that will hold the parameter value (pvar) and
 * whether the parameter should be marked as having been read (read).
 */
void setParam (paramdata * params, uint idx, char * name, char * type, void * pvar, bool read)
{
    params[idx].name = name;
    params[idx].type = type;
    params[idx].pvar = pvar;
    params[idx].read = read;
}


/**
 * readParams
 * 
 * Reads data into the specified array of paramdata structures
 * 'params' of length 'nparams' from the file specified by
 * 'infname'. Any error will result in 'false' being returned, and an
 * error message being printed to the specified FILE/stream 'errstrm'.
 */
bool readParams (char * infname, paramdata * params, uint nparams, FILE * errstrm)
{
    FILE * infile = NULL;
    bool readok = true;
    char inbuf[256];
    int i = 0;
    
    // Basic error checking
    if (infname == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: NULL file name supplied to 'readParams'\r\n");
        }
        return false;
    }
    if (params == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: NULL parameter information supplied to 'readParams'\r\n");
        }
        return false;
    }
    
    // Open the input file and check that it can be read 
    infile = fopen(infname,"r");
    if (infile == NULL)
    {
        if (errstrm != NULL)
        {
            fprintf(errstrm,"ERROR: Could not find input parameter file\r\n");
        }
        return false;
    }
    
    // Scan input file for parameter values
    while (fscanf(infile,"%s",inbuf) != EOF)
    {
        // Check whether the string read matches any parameter names
        for (i = 0; i < nparams; i ++)
        {
            if (strcmp(inbuf,params[i].name) == 0)
            {
                if (fscanf(infile,params[i].type,params[i].pvar) == EOF)
                {
                    if (errstrm != NULL)
                    {
                        fprintf(errstrm,"ERROR: Could not read parameter %s\r\n",inbuf);
                    }
                    readok = false;
                }
                else
                {                    
                    params[i].read = true;
                }
                
                break;
            }
        }
        
        // Parameter name doesn't match any of those in the list
        if (i == nparams)
        {
            if (errstrm != NULL)
            {
                fprintf(errstrm,"ERROR: Parameter not recognised: %s\r\n",inbuf);
            }
            readok = false;
        }
    }
    
    
    // Check that required parameter values have been specified
    for (i = 0; i < nparams; i ++)
    {
        if (!params[i].read)
        {
            if (errstrm != NULL)
            {
                fprintf(errstrm,"ERROR: Required parameter not specified: %s\r\n",params[i].name);
            }
            readok = false;
        }
    }
    
    // We don't need the input file any more
    fclose(infile);
    return readok;
}


/**
 * thomas
 *
 * Performs the thomas algorithm. Here a, b and c are the sub-diagonal, diagonal
 * and super-diagonal in the tridiagonal coefficient matrix, and d is the vector
 * on the right hand side. The solution is filled into x.
 * Warning: will modify c and d!
 */
void thomas (const real * a, const real * b, real * c, real * d, real * x, unsigned int n)
{
  real id;
  int i;
  
  // Modify the coefficients
  c[0] /= b[0];	// Division by zero risk
  d[0] /= b[0];	// Division by zero would imply a singular matrix
  for (i = 1; i < n; i ++)
  {
    id = 1 / (b[i] - c[i-1]*a[i]);     // Division by zero risk
    c[i] *= id;	                       // Last value calculated is redundant
    d[i] = (d[i] - d[i-1]*a[i]) * id;
  }
  
  // Now back substitute
  x[n - 1] = d[n - 1];
  for (i = n-2; i >= 0; i --)
  {
    x[i] = d[i] - c[i] * x[i + 1];
  }
}

/**
 * kahanSum
 *
 * Performs Kahan summation over an input vector 'vec' of length 'N'.
 *
 */
real kahanSum (real * vec, int N)
{
  real sum = 0.0;
  real c = 0.0;                  // A running compensation for lost low-order bits.
  int i;
  real y;
  real t;
  
  for (i = 0; i < N; i ++)
  {
    y = vec[i] - c;     // So far, so good: c is zero.
    t = sum + y;          // Alas, sum is big, y small, so low-order digits of y are lost.
    c = (t - sum) - y;    // (t - sum) recovers the high-order part of y; subtracting y recovers -(low part of y)
    sum = t;              // Algebraically, c should always be zero. Beware overly-aggressive optimizing compilers!
    // Next time around, the lost low part will be added to y in a fresh attempt.
  }
  
  return sum;
}


/**
 * ispow2u
 *
 * Determines whether an integer is a power of 2. Returns true if so, or false otherwise.
 *
 */
bool ispow2u (uint x)
{
  return x && !(x & (x - 1));
}


/**
 * pow2u
 *
 * Calculates an integer p such that p = 2^x. 
 *
 * If the result would be larger than UINT_MAX then 0 is returned.
 *
 */
uint pow2u (uint x)
{
  uint n = 0;
  uint p = 1;
  
  // Bit-shift to find the correct power of 2
  for (n = 0; n < x && n < sizeof(uint)*8; n ++)
  {
    p <<= 1;
  }
  
  return p;
}

/**
 * log2u
 *
 * Calculates an integer n such that x = 2^n. If x is not a power of 2
 * then UINT_MAX is returned.
 *
 */
uint log2u (uint x)
{
  uint n = 0;
  
  // Check that it actually is a power of 2
  if (!ispow2u(x))
  {
    return UINT_MAX;
  }
  
  // Bit-shift to find the correct power of 2
  while (((x & 1) == 0) && (n < sizeof(uint)*8))
  {
    n ++;
    x >>= 1;
  }
  
  return n;
}


/**
 * max3
 *
 * Calculates the maximum of three real values. Not a very pretty implementation,
 * but always employs the minimal number of comparisons (2) to calculate the max.
 *
 */
real max3 (real v1, real v2, real v3)
{
  if (v1 > v2)
  {
    if (v1 > v3)
    {
      return v1;
    }
    else
    {
      return v3;
    }
  }
  else
  {
    if (v2 > v3)
    {
      return v2;
    }
    else
    {
      return v3;
    }
  }
}


/**
 * min3
 *
 * Calculates the minimum of three real values. Not a very pretty implementation,
 * but always employs the minimal number of comparisons (2) to calculate the min.
 *
 */
real min3 (real v1, real v2, real v3)
{
  if (v1 < v2)
  {
    if (v1 < v3)
    {
      return v1;
    }
    else
    {
      return v3;
    }
  }
  else
  {
    if (v2 < v3)
    {
      return v2;
    }
    else
    {
      return v3;
    }
  }
}


/**
 * minmod
 *
 * Implements the minmod function of Kurganov and Tadmor (2000) for three
 * quantities.
 *
 */
real minmod (real v1, real v2, real v3)
{
  real test = v1*v2*v3;
  
  if (v1 < 0 && v2 < 0 && v3 < 0)
  {
    return max3(v1,v2,v3);
  }
  else if (v1 > 0 && v2 > 0 && v3 > 0)
  {
    return min3(v1,v2,v3);
  }
  else
  {
    return 0;
  }
}
