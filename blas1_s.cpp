# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <complex>

using namespace std;

# include "blas0.hpp"
# include "blas1_s.hpp"

//****************************************************************************80

int isamax ( int n, float dx[], int incx )

//****************************************************************************80
//
//  Purpose:
//
//    ISAMAX finds the index of the float vector element of maximum absolute value.
//
//  Discussion:
//
//    WARNING: This index is a 1-based index, not a 0-based index!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, float X[*], the vector to be examined.
//
//    Input, int INCX, the increment between successive entries of SX.
//
//    Output, int ISAMAX, the index of the element of maximum
//    absolute value.
//
{
  float dmax;
  int i;
  int ix;
  int value;

  value = 0;

  if ( n < 1 || incx <= 0 )
  {
    return value;
  }

  value = 1;

  if ( n == 1 )
  {
    return value;
  }

  if ( incx == 1 )
  {
    dmax = fabs ( dx[0] );

    for ( i = 1; i < n; i++ )
    {
      if ( dmax < fabs ( dx[i] ) )
      {
        value = i + 1;
        dmax = fabs ( dx[i] );
      }
    }
  }
  else
  {
    ix = 0;
    dmax = fabs ( dx[0] );
    ix = ix + incx;

    for ( i = 1; i < n; i++ )
    {
      if ( dmax < fabs ( dx[ix] ) )
      {
        value = i + 1;
        dmax = fabs ( dx[ix] );
      }
      ix = ix + incx;
    }
  }

  return value;
}
//****************************************************************************80

float sasum ( int n, float x[], int incx )

//****************************************************************************80
//
//  Purpose:
//
//    SASUM takes the sum of the absolute values of a float vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, float X[*], the vector to be examined.
//
//    Input, int INCX, the increment between successive entries of X.
//    INCX must not be negative.
//
//    Output, float SASUM, the sum of the absolute values of X.
//
{
  int i;
  int j;
  float value;

  value = 0.0;
  j = 0;

  for ( i = 0; i < n; i++ )
  {
    value = value + fabs ( x[j] );
    j = j + incx;
  }

  return value;
}
//****************************************************************************80

void saxpy ( int n, float a, float x[], int incx, float y[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    SAXPY computes float constant times a vector plus a vector.
//
//  Discussion:
//
//    This routine uses unrolled loops for increments equal to one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of elements in X and Y.
//
//    Input, float A, the multiplier of X.
//
//    Input, float X[*], the first vector.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Input/output, float Y[*], the second vector.
//    On output, Y[*] has been replaced by Y[*] + A * X[*].
//
//    Input, int INCY, the increment between successive entries of Y.
//
{
  int i;
  int ix;
  int iy;
  int m;

  if ( n <= 0 )
  {
    return;
  }

  if ( a == 0.0 )
  {
    return;
  }
//
//  Code for unequal increments or equal increments
//  not equal to 1.
//
  if ( incx != 1 || incy != 1 )
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      y[iy] = y[iy] + a * x[ix];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
//
//  Code for both increments equal to 1.
//
  else
  {
    m = n % 4;

    for ( i = 0; i < m; i++ )
    {
      y[i] = y[i] + a * x[i];
    }

    for ( i = m; i < n; i = i + 4 )
    {
      y[i  ] = y[i  ] + a * x[i  ];
      y[i+1] = y[i+1] + a * x[i+1];
      y[i+2] = y[i+2] + a * x[i+2];
      y[i+3] = y[i+3] + a * x[i+3];
    }
  }

  return;
}
//****************************************************************************80

void scopy ( int n, float dx[], int incx, float dy[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    SCOPY copies a float vector X to a vector Y.
//
//  Discussion:
//
//    The routine uses unrolled loops for increments equal to one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of elements in DX and DY.
//
//    Input, float DX[*], the first vector.
//
//    Input, int INCX, the increment between successive entries of DX.
//
//    Output, float DY[*], the second vector.
//
//    Input, int INCY, the increment between successive entries of DY.
//
{
  int i;
  int ix;
  int iy;
  int m;

  if ( n <= 0 )
  {
    return;
  }

  if ( incx == 1 && incy == 1 )
  {
    m = n % 7;

    if ( m != 0 )
    {
      for ( i = 0; i < m; i++ )
      {
        dy[i] = dx[i];
      }
    }

    for ( i = m; i < n; i = i + 7 )
    {
      dy[i] = dx[i];
      dy[i + 1] = dx[i + 1];
      dy[i + 2] = dx[i + 2];
      dy[i + 3] = dx[i + 3];
      dy[i + 4] = dx[i + 4];
      dy[i + 5] = dx[i + 5];
      dy[i + 6] = dx[i + 6];
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( -n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( -n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      dy[iy] = dx[ix];
      ix = ix + incx;
      iy = iy + incy;
    }

  }

  return;
}
//****************************************************************************80

float sdot ( int n, float dx[], int incx, float dy[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    SDOT forms the dot product of two float vectors.
//
//  Discussion:
//
//    This routine uses unrolled loops for increments equal to one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, float DX[*], the first vector.
//
//    Input, int INCX, the increment between successive entries in DX.
//
//    Input, float DY[*], the second vector.
//
//    Input, int INCY, the increment between successive entries in DY.
//
//    Output, float SDOT, the sum of the product of the corresponding
//    entries of DX and DY.
//
{
  float dtemp;
  int i;
  int ix;
  int iy;
  int m;

  dtemp = 0.0;

  if ( n <= 0 )
  {
    return dtemp;
  }
//
//  Code for unequal increments or equal increments
//  not equal to 1.
//
  if ( incx != 1 || incy != 1 )
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      dtemp = dtemp + dx[ix] * dy[iy];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
//
//  Code for both increments equal to 1.
//
  else
  {
    m = n % 5;

    for ( i = 0; i < m; i++ )
    {
      dtemp = dtemp + dx[i] * dy[i];
    }

    for ( i = m; i < n; i = i + 5 )
    {
      dtemp = dtemp + dx[i  ] * dy[i  ]
                    + dx[i+1] * dy[i+1]
                    + dx[i+2] * dy[i+2]
                    + dx[i+3] * dy[i+3]
                    + dx[i+4] * dy[i+4];
    }

  }

  return dtemp;
}
//****************************************************************************80

float snrm2 ( int n, float x[], int incx )

//****************************************************************************80
//
//  Purpose:
//
//    SNRM2 returns the euclidean norm of a float vector.
//
//  Discussion:
//
//     SNRM2 ( X ) = sqrt ( X' * X )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, float X[*], the vector whose norm is to be computed.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Output, float SNRM2, the Euclidean norm of X.
//
{
  float absxi;
  int i;
  int ix;
  float norm;
  float scale;
  float ssq;

  if ( n < 1 || incx < 1 )
  {
    norm = 0.0;
  }
  else if ( n == 1 )
  {
    norm = fabs ( x[0] );
  }
  else
  {
    scale = 0.0;
    ssq = 1.0;
    ix = 0;

    for ( i = 0; i < n; i++ )
    {
      if ( x[ix] != 0.0 )
      {
        absxi = fabs ( x[ix] );
        if ( scale < absxi )
        {
          ssq = 1.0 + ssq * ( scale / absxi ) * ( scale / absxi );
          scale = absxi;
        }
        else
        {
          ssq = ssq + ( absxi / scale ) * ( absxi / scale );
        }
      }
      ix = ix + incx;
    }

    norm = scale * sqrt ( ssq );
  }

  return norm;
}
//****************************************************************************80

void srot ( int n, float x[], int incx, float y[], int incy, float c, float s )

//****************************************************************************80
//
//  Purpose:
//
//    SROT applies a float plane rotation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input/output, float X[*], one of the vectors to be rotated.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Input/output, float Y[*], one of the vectors to be rotated.
//
//    Input, int INCY, the increment between successive elements of Y.
//
//    Input, float C, S, parameters (presumably the cosine and
//    sine of some angle) that define a plane rotation.
//
{
  int i;
  int ix;
  int iy;
  float stemp;

  if ( n <= 0 )
  {
  }
  else if ( incx == 1 && incy == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      stemp = c * x[i] + s * y[i];
      y[i]  = c * y[i] - s * x[i];
      x[i]  = stemp;
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      stemp = c * x[ix] + s * y[iy];
      y[iy] = c * y[iy] - s * x[ix];
      x[ix] = stemp;
      ix = ix + incx;
      iy = iy + incy;
    }

  }

  return;
}
//****************************************************************************80

void srotg ( float *sa, float *sb, float *c, float *s )

//****************************************************************************80
//
//  Purpose:
//
//    SROTG constructs a float Givens plane rotation.
//
//  Discussion:
//
//    Given values A and B, this routine computes
//
//    SIGMA = sign ( A ) if abs ( A ) >  abs ( B )
//          = sign ( B ) if abs ( A ) <= abs ( B );
//
//    R     = SIGMA * ( A * A + B * B );
//
//    C = A / R if R is not 0
//      = 1     if R is 0;
//
//    S = B / R if R is not 0,
//        0     if R is 0.
//
//    The computed numbers then satisfy the equation
//
//    (  C  S ) ( A ) = ( R )
//    ( -S  C ) ( B ) = ( 0 )
//
//    The routine also computes
//
//    Z = S     if abs ( A ) > abs ( B ),
//      = 1 / C if abs ( A ) <= abs ( B ) and C is not 0,
//      = 1     if C is 0.
//
//    The single value Z encodes C and S, and hence the rotation:
//
//    If Z = 1, set C = 0 and S = 1;
//    If abs ( Z ) < 1, set C = sqrt ( 1 - Z * Z ) and S = Z;
//    if abs ( Z ) > 1, set C = 1/ Z and S = sqrt ( 1 - C * C );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input/output, float *SA, *SB,  On input, SA and SB are the values
//    A and B.  On output, SA is overwritten with R, and SB is
//    overwritten with Z.
//
//    Output, float *C, *S, the cosine and sine of the Givens rotation.
//
{
  float r;
  float roe;
  float scale;
  float z;

  if ( fabs ( *sb ) < fabs ( *sa ) )
  {
    roe = *sa;
  }
  else
  {
    roe = *sb;
  }

  scale = fabs ( *sa ) + fabs ( *sb );

  if ( scale == 0.0 )
  {
    *c = 1.0;
    *s = 0.0;
    r = 0.0;
  }
  else
  {
    r = scale * sqrt ( ( *sa / scale ) * ( *sa / scale )
                     + ( *sb / scale ) * ( *sb / scale ) );
    r = r4_sign ( roe ) * r;
    *c = *sa / r;
    *s = *sb / r;
  }

  if ( 0.0 < fabs ( *c ) && fabs ( *c ) <= *s )
  {
    z = 1.0 / *c;
  }
  else
  {
    z = *s;
  }

  *sa = r;
  *sb = z;

  return;
}
//****************************************************************************80

void sscal ( int n, float sa, float x[], int incx )

//****************************************************************************80
//
//  Purpose:
//
//    SSCAL scales a float vector by a constant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, float SA, the multiplier.
//
//    Input/output, float X[*], the vector to be scaled.
//
//    Input, int INCX, the increment between successive entries of X.
//
{
  int i;
  int ix;
  int m;

  if ( n <= 0 )
  {
  }
  else if ( incx == 1 )
  {
    m = n % 5;

    for ( i = 0; i < m; i++ )
    {
      x[i] = sa * x[i];
    }

    for ( i = m; i < n; i = i + 5 )
    {
      x[i]   = sa * x[i];
      x[i+1] = sa * x[i+1];
      x[i+2] = sa * x[i+2];
      x[i+3] = sa * x[i+3];
      x[i+4] = sa * x[i+4];
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    for ( i = 0; i < n; i++ )
    {
      x[ix] = sa * x[ix];
      ix = ix + incx;
    }

  }

  return;
}
//****************************************************************************80

void sswap ( int n, float x[], int incx, float y[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    SSWAP interchanges two float vectors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input/output, float X[*], one of the vectors to swap.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Input/output, float Y[*], one of the vectors to swap.
//
//    Input, int INCY, the increment between successive elements of Y.
//
{
  int i;
  int ix;
  int iy;
  int m;
  float temp;

  if ( n <= 0 )
  {
  }
  else if ( incx == 1 && incy == 1 )
  {
    m = n % 3;

    for ( i = 0; i < m; i++ )
    {
      temp = x[i];
      x[i] = y[i];
      y[i] = temp;
    }

    for ( i = m; i < n; i = i + 3 )
    {
      temp = x[i];
      x[i] = y[i];
      y[i] = temp;

      temp = x[i+1];
      x[i+1] = y[i+1];
      y[i+1] = temp;

      temp = x[i+2];
      x[i+2] = y[i+2];
      y[i+2] = temp;
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      temp = x[ix];
      x[ix] = y[iy];
      y[iy] = temp;
      ix = ix + incx;
      iy = iy + incy;
    }
  }
  return;
}

