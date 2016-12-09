# include <cstdlib>
# include <cmath>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cstring>

using namespace std;

# include "hermite_polynomial.hpp"

//****************************************************************************80

double h_integral ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    H_INTEGRAL evaluates the integral of H(i,x).
//
//  Discussion:
//
//    H(i,x) is the physicist's Hermite polynomial of degree I.
//
//    The integral computed is:
//
//      integral ( -oo < x < +oo ) H(i,x) exp(-x^2) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the integral.  
//    0 <= N.
//
//    Output, double H_INTEGRAL, the value of the integral.
//
{
  const double r8_pi = 3.141592653589793;
  double value;

  if ( ( n % 2 ) == 1 )
  {
    value = 0.0;
  }
  else
  {
    value = r8_factorial2 ( n - 1 ) * sqrt ( r8_pi ) / pow ( 2.0, n / 2 );
  }

  return value;
}
//****************************************************************************80

double *h_polynomial_coefficients ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    H_POLYNOMIAL_COEFFICIENTS: coefficients of H(i,x).
//
//  Discussion:
//
//    H(i,x) is the physicist's Hermite polynomial of degree I.
//
//  First terms:
//
//    N/K     0     1      2      3       4     5      6    7      8    9   10
//
//     0      1
//     1      0     2
//     2     -2     0      4
//     3      0   -12      0      8
//     4     12     0    -48      0      16
//     5      0   120      0   -160       0    32
//     6   -120     0    720      0    -480     0     64
//     7      0 -1680      0   3360       0 -1344      0   128
//     8   1680     0 -13440      0   13440     0  -3584     0    256
//     9      0 30240      0 -80640       0 48384      0 -9216      0 512
//    10 -30240     0 302400      0 -403200     0 161280     0 -23040   0 1024
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//
//    Output, double HERMITE_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the 
//    coefficients of the Hermite polynomials of orders 0 through N.
//
{
  double *c;
  int i;
  int j;

  if ( n < 0 )
  {
    return NULL;
  }

  c = new double[(n+1)*(n+1)];

  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      c[i+j*(n+1)] = 0.0;
    }
  }

  c[0+0*(n+1)] = 1.0;

  if ( n == 0 )
  {
    return c;
  }

  c[1+1*(n+1)] = 2.0;

  for ( i = 2; i <= n; i++ )
  {
    c[i+0*(n+1)] = - 2.0 * ( double ) ( i - 1 ) * c[i-2+0*(n+1)];
    for ( j = 1; j <= i - 2; j++ )
    {
      c[i+j*(n+1)] =   2.0 * c[i-1+(j-1)*(n+1)] 
                     - 2.0 * ( double ) ( i - 1 ) * c[i-2+j*(n+1)];
    }
    c[i+(i-1)*(n+1)] =  2.0 * c[i-1+(i-2)*(n+1)];
    c[i+ i   *(n+1)] =  2.0 * c[i-1+(i-1)*(n+1)];
  }

  return c;
}
//****************************************************************************80

double *h_polynomial_value ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    H_POLYNOMIAL_VALUE evaluates H(i,x).
//
//  Discussion:
//
//    H(i,x) is the physicist's Hermite polynomial of degree I.
//
//  Differential equation:
//
//    Y'' - 2 X Y' + 2 N Y = 0
//
//  First terms:
//
//      1
//      2 X
//      4 X^2     -  2
//      8 X^3     - 12 X
//     16 X^4     - 48 X^2     + 12
//     32 X^5    - 160 X^3    + 120 X
//     64 X^6    - 480 X^4    + 720 X^2    - 120
//    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
//    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
//    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
//   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
//
//  Recursion:
//
//    H(0,X) = 1,
//    H(1,X) = 2*X,
//    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
//
//  Norm:
//
//    Integral ( -oo < X < +oo ) exp ( - X^2 ) * H(N,X)^2 dX
//    = sqrt ( PI ) * 2^N * N!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int M, the number of evaluation points.
//
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//
//    Input, double X[M], the evaluation points.
//
//    Output, double H_POLYNOMIAL_VALUE[M*(N+1)], the values of the first 
//    N+1 Hermite polynomials at the evaluation points.
//
{
  int i;
  int j;
  double *p;

  if ( n < 0 )
  {
    return NULL;
  }

  p = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    p[i+0*m] = 1.0;
  }

  if ( n == 0 )
  {
    return p;
  }

  for ( i = 0; i < m; i++ )
  {
    p[i+1*m] = 2.0 * x[i];
  }

  for ( j = 2; j <= n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      p[i+j*m] = 2.0 * x[i] * p[i+(j-1)*m]
        - 2.0 * ( double ) ( j - 1 ) * p[i+(j-2)*m];
    }
  }
  return p;
}
//****************************************************************************80

void h_polynomial_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    H_POLYNOMIAL_VALUES: tabulated values of H(i,x).
//
//  Discussion:
//
//    H(i,x) is the physicist's Hermite polynomial of degree I.
//
//    In Mathematica, the function can be evaluated by:
//
//      HermiteH[n,x]
//
//  Differential equation:
//
//    Y'' - 2 X Y' + 2 N Y = 0;
//
//  First terms:
//
//      1
//      2 X
//      4 X^2     -  2
//      8 X^3     - 12 X
//     16 X^4     - 48 X^2     + 12
//     32 X^5    - 160 X^3    + 120 X
//     64 X^6    - 480 X^4    + 720 X^2    - 120
//    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
//    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
//    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
//   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
//
//  Recursion:
//
//    H(0,X) = 1,
//    H(1,X) = 2*X,
//    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
//
//  Norm:
//
//    Integral ( -oo < X < +oo ) exp ( - X^2 ) * H(N,X)^2 dX
//    = sqrt ( PI ) * 2^N * N!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the polynomial.
//
//    Output, double &X, the point where the polynomial is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 18

  static double fx_vec[N_MAX] = {
      0.1000000000000000E+01,
      0.1000000000000000E+02,
      0.9800000000000000E+02,
      0.9400000000000000E+03,
      0.8812000000000000E+04,
      0.8060000000000000E+05,
      0.7178800000000000E+06,
      0.6211600000000000E+07,
      0.5206568000000000E+08,
      0.4212712000000000E+09,
      0.3275529760000000E+10,
      0.2432987360000000E+11,
      0.1712370812800000E+12,
      0.0000000000000000E+00,
      0.4100000000000000E+02,
     -0.8000000000000000E+01,
      0.3816000000000000E+04,
      0.3041200000000000E+07 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10, 11,
    12,  5,  5,
     5,  5,  5 };

  static double x_vec[N_MAX] = {
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     0.0E+00,
     0.5E+00,
     1.0E+00,
     3.0E+00,
     1.0E+01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double *h_polynomial_zeros ( int nt )

//****************************************************************************80
//
//  Purpose:
//
//    H_POLYNOMIAL_ZEROS: zeros of H(i,x).
//
//  Discussion:
//
//    H(i,x) is the physicist's Hermite polynomial of degree I.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, the degree of the polynomial.
//
//    Output, double H_POLYNOMIAL_ZEROS[NT], the zeros of the polynomial.
//
{
  double *bj;
  int i;
  const double r8_pi = 3.141592653589793;
  double *wts;
  double *z;

  z = new double[nt];

  for ( i = 0; i < nt; i++ )
  {
    z[i] = 0.0;
  }

  bj = new double[nt];

  for ( i = 0; i < nt; i++ )
  {
    bj[i] = sqrt ( ( double ) ( i + 1 ) / 2.0 );
  }

  wts = new double[nt];
  for ( i = 0; i < nt; i++ )
  {
    wts[i] = 0.0;
  }
  wts[0] = sqrt ( sqrt ( r8_pi ) );

  imtqlx ( nt, z, bj, wts );

  return z;
}
//****************************************************************************80

void h_quadrature_rule ( int nt, double t[], double wts[] )

//****************************************************************************80
//
//  Purpose:
//
//    H_QUADRATURE_RULE: quadrature for H(i,x).
//
//  Discussion:
//
//    H(i,x) is the physicist's Hermite polynomial of degree I.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, the order of the rule.
//
//    Output, double T[NT], WTS[NT], the points and weights 
//    of the rule.
//
{
  double *bj;
  int i;
  const double r8_pi = 3.141592653589793;

  for ( i = 0; i < nt; i++ )
  {
    t[i] = 0.0;
  }

  bj = new double[nt];

  for ( i = 0; i < nt; i++ )
  {
    bj[i] = sqrt ( ( double ) ( i + 1 ) / 2.0 );
  }

  for ( i = 0; i < nt; i++ )
  {
    wts[i] = 0.0;
  }

  wts[0] = sqrt ( sqrt ( r8_pi ) );

  imtqlx ( nt, t, bj, wts );

  for ( i = 0; i < nt; i++ )
  {
    wts[i] = wts[i] * wts[i];
  }
  return;
}
//****************************************************************************80

double he_double_product_integral ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    HE_DOUBLE_PRODUCT_INTEGRAL: integral of He(i,x)*He(j,x)*e^(-x^2/2).
//
//  Discussion:
//
//    He(i,x) represents the probabilist's Hermite polynomial.
//
//    VALUE = integral ( -oo < x < +oo ) He(i,x)*He(j,x) exp(-x^2/2) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical Methods for Stochastic Computations: A Spectral Method Approach,
//    Princeton, 2010,
//    ISBN13: 978-0-691-14212-8,
//    LC: QA274.23.X58.
//
//  Parameters:
//
//    Input, int I, J, the polynomial indices.
//
//    Output, double HE_DOUBLE_PRODUCT_INTEGRAL, the value of the integral.
//
{
  double value;

  if ( i == j )
  {
    value = r8_factorial ( i );
  }
  else
  {
    value = 0.0;
  }
  return value;
}
//****************************************************************************80

double he_integral ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    HE_INTEGRAL evaluates the integral of He(i,x).
//
//  Discussion:
//
//    He(i,x) represents the probabilist's Hermite polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the integral.  
//    0 <= N.
//
//    Output, double HE_INTEGRAL, the value of the integral.
//
{
  const double r8_pi = 3.141592653589793;
  double value;

  if ( ( n % 2 ) == 1 )
  {
    value = 0.0;
  }
  else
  {
    value = r8_factorial2 ( n - 1 ) * sqrt ( 2.0 * r8_pi );
  }

  return value;
}
//****************************************************************************80

double *he_polynomial_coefficients ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    HE_POLYNOMIAL_COEFFICIENTS: coefficients of He(i,x).
//
//  Discussion:
//
//    He(i,x) represents the probabilist's Hermite polynomial.
//
//  First terms:
//
//    N/K     0     1      2      3       4     5      6    7      8    9   10
//
//     0      1
//     1      0     1
//     2     -1     0      1
//     3      0    -3      0      1
//     4      3     0     -6      0       1
//     5      0    15      0    -10       0     1
//     6    -15     0     45      0     -15     0      1
//     7      0  -105      0    105       0   -21      0     1
//     8    105     0   -420      0     210     0    -28     0      1
//     9      0   945      0  -1260       0   378      0   -36      0   1
//    10   -945     0   4725      0   -3150     0    630     0    -45   0    1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//
//    Output, double HE_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients 
//    of the Hermite polynomials of orders 0 through N.
//
{
  double *c;
  int i;
  int j;

  if ( n < 0 )
  {
    return NULL;
  }

  c = new double[(n+1)*(n+1)];

  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      c[i+j*(n+1)] = 0.0;
    }
  }

  c[0+0*(n+1)] = 1.0;

  if ( n == 0 )
  {
    return c;
  }

  c[1+1*(n+1)] = 1.0;

  for ( i = 2; i <= n; i++ )
  {
    c[i+0*(n+1)] =                      - ( double ) ( i - 1 ) * c[i-2+0*(n+1)];
    for ( j = 1; j <= i - 2; j++ )
    {
      c[i+j*(n+1)] = c[i-1+(j-1)*(n+1)] - ( double ) ( i - 1 ) * c[i-2+j*(n+1)];
    }
    c[i+(i-1)*(n+1)] =   c[i-1+(i-2)*(n+1)];
    c[i+ i   *(n+1)] =   c[i-1+(i-1)*(n+1)];
  }

  return c;
}
//****************************************************************************80

double *he_polynomial_value ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    HE_POLYNOMIAL_VALUE evaluates He(i,x).
//
//  Discussion:
//
//    He(i,x) represents the probabilist's Hermite polynomial.
//
//  Differential equation:
//
//    ( exp ( - 0.5 * x^2 ) * He(n,x)' )' + n * exp ( - 0.5 * x^2 ) * He(n,x) = 0
//
//  First terms:
//
//   1
//   X
//   X^2  -  1
//   X^3  -  3 X
//   X^4  -  6 X^2 +   3
//   X^5  - 10 X^3 +  15 X
//   X^6  - 15 X^4 +  45 X^2 -   15
//   X^7  - 21 X^5 + 105 X^3 -  105 X
//   X^8  - 28 X^6 + 210 X^4 -  420 X^2 +  105
//   X^9  - 36 X^7 + 378 X^5 - 1260 X^3 +  945 X
//   X^10 - 45 X^8 + 630 X^6 - 3150 X^4 + 4725 X^2 - 945
//
//  Recursion:
//
//    He(0,X) = 1,
//    He(1,X) = X,
//    He(N,X) = X * He(N-1,X) - (N-1) * He(N-2,X)
//
//  Orthogonality:
//
//    Integral ( -oo < X < +oo ) exp ( - 0.5 * X^2 ) * He(M,X) He(N,X) dX 
//    = sqrt ( 2 * pi ) * N// * delta ( N, M )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
//    NIST Handbook of Mathematical Functions,
//    Cambridge University Press, 2010,
//    ISBN: 978-0521192255,
//    LC: QA331.N57.
//
//  Parameters:
//
//    Input, int M, the number of evaluation points.
//
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//
//    Input, double X[M], the evaluation points.
//
//    Output, double HE_POLYNOMIAL_VALUE[M*(N+1)], the values of the
//    probabilist's Hermite polynomials of index 0 through N.
//
{
  int i;
  int j;
  double *p;

  if ( n < 0 )
  {
    return NULL;
  }

  p = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    p[i+0*m] = 1.0;
  }

  if ( n == 0 )
  {
    return p;
  }

  for ( i = 0; i < m; i++ )
  {
    p[i+1*m] = x[i];
  }

  for ( j = 2; j <= n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      p[i+j*m] = x[i] * p[i+(j-1)*m] - ( double ) ( j - 1 ) * p[i+(j-2)*m];
    }
  }
  return p;
}
//****************************************************************************80

void he_polynomial_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    HE_POLYNOMIAL_VALUES: tabulated values of He(i,x).
//
//  Discussion:
//
//    He(i,x) represents the probabilist's Hermite polynomial.
//
//    In Mathematica, the function can be evaluated by:
//
//      He(n,x) = HermiteH[n,x/Sqrt[2]] / Sqrt [ 2^n ] 
//
//  First terms:
//
//   1
//   X
//   X^2  -  1
//   X^3  -  3 X
//   X^4  -  6 X^2 +   3
//   X^5  - 10 X^3 +  15 X
//   X^6  - 15 X^4 +  45 X^2 -   15
//   X^7  - 21 X^5 + 105 X^3 -  105 X
//   X^8  - 28 X^6 + 210 X^4 -  420 X^2 +  105
//   X^9  - 36 X^7 + 378 X^5 - 1260 X^3 +  945 X
//   X^10 - 45 X^8 + 630 X^6 - 3150 X^4 + 4725 X^2 - 945
//
//  Recursion:
//
//    He(0,X) = 1,
//    He(1,X) = X,
//    He(N,X) = X * He(N-1,X) - (N-1) * He(N-2,X)
//
//  Norm:
//
//    Integral ( -oo < X < +oo ) exp ( - 0.5 * X^2 ) * He(M,X) He(N,X) dX 
//    = sqrt ( 2 * pi ) * N! * delta ( M, N )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the polynomial.
//
//    Output, double &X, the point where the polynomial is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 18

  static double fx_vec[N_MAX] = {
    1.000000000000000E+00, 
    5.000000000000000E+00, 
    24.00000000000000E+00, 
    110.0000000000000E+00,
    478.0000000000000E+00,
    1950.000000000000E+00, 
    7360.000000000000E+00, 
    25100.00000000000E+00, 
    73980.00000000000E+00, 
    169100.0000000000E+00, 
    179680.0000000000E+00, 
   -792600.0000000000E+00, 
   -5939480.000000000E+00, 
    0.000000000000000E+00, 
    6.281250000000000E+00, 
    6.000000000000000E+00, 
    18.00000000000000E+00, 
    90150.00000000000E+00 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10, 11,
    12,  5,  5,
     5,  5,  5 };

  static double x_vec[N_MAX] = {
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     0.0E+00,
     0.5E+00,
     1.0E+00,
     3.0E+00,
     1.0E+01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double *he_polynomial_zeros ( int nt )

//****************************************************************************80
//
//  Purpose:
//
//    HE_POLYNOMIAL_ZEROS: zeros of He(i,x).
//
//  Discussion:
//
//    He(i,x) represents the probabilist's Hermite polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, the degree of the polynomial.
//
//    Output, double HE_POLYNOMIAL_ZEROS[NT], the zeros of the polynomial.
//
{
  double *bj;
  int i;
  const double r8_pi = 3.141592653589793;
  double *wts;
  double *z;

  z = new double[nt];

  for ( i = 0; i < nt; i++ )
  {
    z[i] = 0.0;
  }

  bj = new double[nt];

  for ( i = 0; i < nt; i++ )
  {
    bj[i] = sqrt ( ( double ) ( i + 1 ) / 2.0 );
  }

  wts = new double[nt];

  for ( i = 0; i < nt; i++ )
  {
    wts[i] = 0.0;
  }

  wts[0] = sqrt ( sqrt ( r8_pi ) );

  imtqlx ( nt, z, bj, wts );

  for ( i = 0; i < nt; i++ )
  {
    z[i] = z[i] * sqrt ( 2.0 );
  }

  delete [] bj;
  delete [] wts;

  return z;
}
//****************************************************************************80

void he_quadrature_rule ( int nt, double t[], double wts[] )

//****************************************************************************80
//
//  Purpose:
//
//    HE_QUADRATURE_RULE: quadrature for He(i,x).
//
//  Discussion:
//
//    He(i,x) represents the probabilist's Hermite polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, the order of the rule.
//
//    Output, double T[NT], WTS[NT], the points and weights 
//    of the rule.
//
{
  double *bj;
  int i;
  const double r8_pi = 3.141592653589793;

  for ( i = 0; i < nt; i++ )
  {
    t[i] = 0.0;
  }

  bj = new double[nt];

  for ( i = 0; i < nt; i++ )
  {
    bj[i] = sqrt ( ( double ) ( i + 1 ) / 2.0 );
  }

  for ( i = 0; i < nt; i++ )
  {
    wts[i] = 0.0;
  }
  wts[0] = sqrt ( sqrt ( r8_pi ) );

  imtqlx ( nt, t, bj, wts );

  for ( i = 0; i < nt; i++ )
  {
    t[i] = t[i] * sqrt ( 2.0 );
  }
  for ( i = 0; i < nt; i++ )
  {
    wts[i] = wts[i] * wts[i] * sqrt ( 2.0 );
  }
  delete [] bj;

  return;
}
//****************************************************************************80

double he_triple_product_integral ( int i, int j, int k )

//****************************************************************************80
//
//  Purpose:
//
//    HE_TRIPLE_PRODUCT_INTEGRAL: integral of He(i,x)*He(j,x)*He(k,x)*e^(-x^2/2).
//
//  Discussion:
//
//    He(i,x) represents the probabilist's Hermite polynomial.
//
//    VALUE = integral ( -oo < x < +oo ) He(i,x)*He(j,x)*He(k,x) exp(-x^2/2) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu,
//    Numerical Methods for Stochastic Computations: A Spectral Method Approach,
//    Princeton, 2010,
//    ISBN13: 978-0-691-14212-8,
//    LC: QA274.23.X58.
//
//  Parameters:
//
//    Input, int I, J, K, the polynomial indices.
//
//    Output, double HE_TRIPLE_PRODUCT_INTEGRAL, the value of the integral.
//
{
  int s;
  double value;

  s = ( i + j + k ) / 2;

  if ( s < i || s < j || s < k )
  {
    value = 0.0;
  }
  else if ( ( ( i + j + k ) % 2 ) != 0 )
  {
    value = 0.0;
  }
  else
  {
    value = r8_factorial ( i ) / r8_factorial ( s - i ) 
          * r8_factorial ( j ) / r8_factorial ( s - j ) 
          * r8_factorial ( k ) / r8_factorial ( s - k );
  }

  return value;
}
//****************************************************************************80

double *hen_exponential_product ( int p, double b )

//****************************************************************************80
//
//  Purpose:
//
//    HEN_EXPONENTIAL_PRODUCT: exponential product exp(b*x)*Hen(i,x)*Hen(j,x).
//
//  Discussion:
//
//    Hen(i,x) is the normalized probabilist's Hermite polynomial of degree I.
//
//
//    For polynomial chaos applications, it is of interest to know the
//    value of the integrals of products of exp(B*X) with every possible pair
//    of basis functions.  That is, we'd like to form
//
//      Tij = Integral ( -oo < X < +oo ) 
//        exp(B*X) * Hen(I,X) * Hen(J,X) exp(-0.5*X*X) dx
//
//    We will estimate these integrals using Gauss-Hermite quadrature.
//    Because of the exponential factor exp(B*X), the quadrature will not 
//    be exact.
//
//    However, when B = 0, the quadrature is exact, and moreoever, the
//    table will be the identity matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the maximum degree of the 
//    polyonomial factors.  0 <= P.
//
//    Input, double B, the coefficient of X in the exponential factor.
//
//    Output, double HEN_EXPONENTIAL_PRODUCT[(P+1)*(P+1)], the table of 
//    integrals.  TABLE(I,J) represents the weighted integral of 
//    exp(B*X) * Hen(I,X) * Hen(J,X).
//
{
  double *h_table;
  int i;
  int j;
  int k;
  int order;
  double *table;
  double *w_table;
  double x;
  double *x_table;

  table = new double[(p+1)*(p+1)];
  for ( j = 0; j <= p; j++ )
  {
    for ( i = 0; i <= p; i++ )
    {
      table[i+j*(p+1)] = 0.0;
    }
  }

  order = ( 3 * p + 4 ) / 2;

  x_table = new double[order];
  w_table = new double[order];

  he_quadrature_rule ( order, x_table, w_table );

  for ( k = 0; k < order; k++ )
  {
    x = x_table[k];
    h_table = hen_polynomial_value ( 1, p, x_table+k );
//
//  The following formula is an outer product in H_TABLE.
//
    for ( j = 0; j <= p; j++ )
    {
      for ( i = 0; i <= p; i++ )
      {
        table[i+j*(p+1)] = table[i+j*(p+1)] 
          + w_table[k] * exp ( b * x ) * h_table[i] * h_table[j];
      }
    }
    delete [] h_table;
  }

  delete [] w_table;
  delete [] x_table;

  return table;
}
//****************************************************************************80

double *hen_polynomial_value ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    HEN_POLYNOMIAL_VALUE: evaluates Hen(i,x).
//
//  Discussion:
//
//    Hen(i,x) is the normalized probabilist's Hermite polynomial of degree I.
//
//    These polynomials satisfy the orthonormality condition:
//
//      Integral ( -oo < X < +oo ) exp ( - 0.5 * X^2 ) * Hen(M,X) Hen(N,X) dX 
//      = delta ( N, M )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
//    NIST Handbook of Mathematical Functions,
//    Cambridge University Press, 2010,
//    ISBN: 978-0521192255,
//    LC: QA331.N57.
//
//  Parameters:
//
//    Input, int M, the number of evaluation points.
//
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//
//    Input, double X[M], the evaluation points.
//
//    Output, double HEN_POLYNOMIAL_VALUE[M*(N+1)], the values of the 
//    polynomials of index 0 through N.
//
{
  double fact;
  int i;
  int j;
  double *p;
  const double r8_pi = 3.141592653589793;

  if ( n < 0 )
  {
    return NULL;
  }

  p = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    p[i+0*m] = 1.0;
  }

  if ( n == 0 )
  {
    return p;
  }

  for ( i = 0; i < m; i++ )
  {
    p[i+1*m] = x[i];
  }

  for ( j = 2; j <= n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      p[i+j*m] = x[i] * p[i+(j-1)*m] - ( double ) ( j - 1 ) * p[i+(j-2)*m];
    }
  }
//
//  Normalize.
//
  fact = 1.0;
  for ( j = 0; j <= n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      p[i+j*m] = p[i+j*m] / sqrt ( fact * sqrt ( 2.0 * r8_pi ) );
    }
    fact = fact * ( double ) ( j + 1 );
  }
  return p;
}
//****************************************************************************80

double *hen_power_product ( int p, int e )

//****************************************************************************80
//
//  Purpose:
//
//    HEN_POWER_PRODUCT: power products, x^e*Hen(i,x)*Hen(j,x).
//
//  Discussion:
//
//    Hen(i,x) is the normalized probabilist's Hermite polynomial of degree I.
//
//    For polynomial chaos applications, it is of interest to know the
//    value of the integrals of products of X with every possible pair
//    of basis functions.  That is, we'd like to form
//
//      Tij = Integral ( -oo < X < +oo ) 
//        X^E * Hen(I,X) * Hen(J,X) exp(-0.5*X*X) dx
//
//    We will estimate these integrals using Gauss-Hermite quadrature.
//
//    When E is 0, the computed table should be the identity matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the maximum degree of the polyonomial 
//    factors.  0 <= P.
//
//    Input, int E, the exponent of X in the integrand.
//    0 <= E.
//
//    Output, double HEN_POWER_PRODUCT[(P+1)*(P+1)], the table of integrals.  
//    TABLE(I,J) represents the weighted integral of 
//    X^E * Hen(I,X) * Hen(J,X).
//
{
  double *h_table;
  int i;
  int j;
  int k;
  int order;
  double *table;
  double *w_table;
  double x;
  double *x_table;

  table = new double[(p+1)*(p+1)];
  for ( j = 0; j <= p; j++ )
  {
    for ( i = 0; i <= p; i++ )
    {
      table[i+j*(p+1)] = 0.0;
    }
  }

  order = p + 1 + ( ( e + 1 ) / 2 );

  x_table = new double[order];
  w_table = new double[order];

  he_quadrature_rule ( order, x_table, w_table );

  for ( k = 0; k < order; k++ )
  {
    x = x_table[k];
    h_table = hen_polynomial_value ( 1, p, x_table+k );
//
//  The following formula is an outer product in H_TABLE.
//
    if ( e == 0 )
    {
      for ( j = 0; j <= p; j++ )
      {
        for ( i = 0; i <= p; i++ )
        {
          table[i+j*(p+1)] = table[i+j*(p+1)] 
            + w_table[k] * h_table[i] * h_table[j];
        }
      }
    }
    else
    {
      for ( j = 0; j <= p; j++ )
      {
        for ( i = 0; i <= p; i++ )
        {
          table[i+j*(p+1)] = table[i+j*(p+1)] 
            + w_table[k] * pow ( x, e ) * h_table[i] * h_table[j];
        }
      }
    }
    delete [] h_table;
  }

  delete [] w_table;
  delete [] x_table;

  return table;
}
//****************************************************************************80

double *hf_exponential_product ( int p, double b )

//****************************************************************************80
//
//  Purpose:
//
//    HF_EXPONENTIAL_PRODUCT: exponential products, exp(b*x)*Hf(i,x)*Hf(j,x).
//
//  Discussion:
//
//    Hf(i,x) represents the Hermite function of "degree" I.   
//
//    For polynomial chaos applications, it is of interest to know the
//    value of the integrals of products of exp(B*X) with every possible pair
//    of basis functions.  That is, we'd like to form
//
//      Tij = Integral ( -oo < X < +oo ) exp(B*X) * Hf(I,X) * Hf(J,X) dx
//
//    We will estimate these integrals using Gauss-Hermite quadrature.
//    Because of the exponential factor exp(B*X), the quadrature will not 
//    be exact.
//
//    However, when B = 0, the quadrature is exact, and moreoever, the
//    table will be the identity matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the maximum degree of the polyonomial 
//    factors.  0 <= P.
//
//    Input, double B, the coefficient of X in the exponential factor.
//
//    Output, double HF_EXPONENTIAL_PRODUCT[(P+1)*(P+1)], the table of
//    integrals.  TABLE(I,J) represents the integral of 
//    exp(B*X) * Hf(I,X) * Hf(J,X).
//
{
  double *h_table;
  int i;
  int j;
  int k;
  int order;
  double *table;
  double *w_table;
  double x;
  double *x_table;

  table = new double[(p+1)*(p+1)];
  for ( j = 0; j <= p; j++ )
  {
    for ( i = 0; i <= p; i++ )
    {
      table[i+j*(p+1)] = 0.0;
    }
  }

  order = ( 3 * p + 4 ) / 2;

  x_table = new double[order];
  w_table = new double[order];

  hf_quadrature_rule ( order, x_table, w_table );

  for ( k = 0; k < order; k++ )
  {
    x = x_table[k];
    h_table = hf_function_value ( 1, p, x_table+k );
//
//  The following formula is an outer product in H_TABLE.
//
    for ( j = 0; j <= p; j++ )
    {
      for ( i = 0; i <= p; i++ )
      {
        table[i+j*(p+1)] = table[i+j*(p+1)] 
          + w_table[k] * exp ( b * x ) * h_table[i] * h_table[j];
      }
    }
    delete [] h_table;
  }

  delete [] w_table;
  delete [] x_table;

  return table;
}
//****************************************************************************80

double *hf_function_value ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    HF_FUNCTION_VALUE evaluates Hf(i,x).
//
//  Discussion:
//
//    Hf(i,x) represents the Hermite function of "degree" I.   
//
//    The Hermite function of degree n is related to the physicist's
//    Hermite polynomial H(n,x):
//
//      Hf(n,x) = H(n,x) * exp ( - 0.5 * x^2 ) / sqrt ( 2^n n// sqrt ( pi ) )
//
//    The Hermite functions are orthonormal:
//
//      Integral ( -oo < x < +oo ) Hf(m,x) Hf(n,x) dx = delta ( m, n )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
//    NIST Handbook of Mathematical Functions,
//    Cambridge University Press, 2010,
//    ISBN: 978-0521192255,
//    LC: QA331.N57.
//
//  Parameters:
//
//    Input, int M, the number of evaluation points.
//
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//
//    Input, double X[M], the point at which the polynomials are 
//    to be evaluated.
//
//    Output, double HF_FUNCTION_VALUE[M*(N+1)], the values of the Hermite
//    functions of index 0 through N at the evaluation points.
//
{
  double *f;
  int i;
  int j;
  const double r8_pi = 3.141592653589793;

  f = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    f[i+0*m] = exp ( - 0.5 * x[i] * x[i] ) / sqrt ( sqrt ( r8_pi ) );
  }
  if ( n == 0 )
  {
    return f;
  }

  for ( i = 0; i < m; i++ )
  {
    f[i+1*m] = 2.0 * exp ( - 0.5 * x[i] * x[i] ) * x[i] 
      / sqrt ( 2.0 * sqrt ( r8_pi ) );
  }

  for ( j = 2; j <= n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      f[i+j*m] = ( sqrt ( 2.0 ) * x[i] * f[i+(j-1)*m]
        - sqrt ( ( double ) ( j - 1 ) ) * f[i+(j-2)*m] ) 
        / sqrt ( ( double ) ( j ) );
    }
  }
 
  return f;
}
//****************************************************************************80

void hf_function_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    HF_FUNCTION_VALUES: tabulated values of Hf(i,x).
//
//  Discussion:
//
//    Hf(i,x) represents the Hermite function of "degree" I.   
//
//    In Mathematica, the function can be evaluated by:
//
//      Hf(n,x) = HermiteH[n,x] 
//        * Exp [ -1/2 * x^2] / Sqrt [ 2^n * n! * Sqrt[Pi] ]
//
//    The Hermite functions are orthonormal:
//
//      Integral ( -oo < x < +oo ) Hf(m,x) Hf(n,x) dx = delta ( m, n )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the polynomial.
//
//    Output, double &X, the point where the polynomial is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 23

  static double fx_vec[N_MAX] = {
    0.7511255444649425E+00,  0.0000000000000000E+00, -0.5311259660135985E+00, 
    0.0000000000000000E+00,  0.4599685791773266E+00,  0.0000000000000000E+00, 
    0.4555806720113325E+00,  0.6442883651134752E+00,  0.3221441825567376E+00, 
   -0.2630296236233334E+00, -0.4649750762925110E+00, -0.5881521185179581E-01, 
    0.3905052515434106E+00,  0.2631861423064045E+00, -0.2336911435996523E+00, 
   -0.3582973361472840E+00,  0.6146344487883041E-01,  0.3678312067984882E+00, 
    0.9131969309166278E-01,  0.4385750950032321E+00, -0.2624689527931006E-01, 
    0.5138426125477819E+00,  0.9355563118061758E-01 };

  static int n_vec[N_MAX] = {
    0,  1,  2,  
    3,  4,  5,  
    0,  1,  2,  
    3,  4,  5,  
    6,  7,  8,  
    9, 10, 11,  
   12,  5,  5,  
    5,  5 };

  static double x_vec[N_MAX] = {
    0.0E+00, 0.0E+00, 0.0E+00, 
    0.0E+00, 0.0E+00, 0.0E+00, 
    1.0E+00, 1.0E+00, 1.0E+00, 
    1.0E+00, 1.0E+00, 1.0E+00, 
    1.0E+00, 1.0E+00, 1.0E+00, 
    1.0E+00, 1.0E+00, 1.0E+00, 
    1.0E+00, 0.5E+00, 2.0E+00, 
    3.0E+00, 4.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double *hf_power_product ( int p, int e )

//****************************************************************************80
//
//  Purpose:
//
//    HF_POWER_PRODUCT: power products x^e*Hf(i,x)*Hf(j,x).
//
//  Discussion:
//
//    Hf(i,x) represents the Hermite function of "degree" I.   
//
//    For polynomial chaos applications, it is of interest to know the
//    value of the integrals of products of X with every possible pair
//    of basis functions.  That is, we'd like to form
//
//      Tij = Integral ( -oo < X < +oo ) X^E * Hf(I,X) * Hf(J,X) dx
//
//    We will estimate these integrals using Gauss-Hermite quadrature.
//
//    When E is 0, the computed table should be the identity matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the maximum degree of the polyonomial 
//    factors.  0 <= P.
//
//    Input, int E, the exponent of X in the integrand.
//    0 <= E.
//
//    Output, double HF_POWER_PRODUCT[(P+1)*(P+1)], the table of integrals.  
//    TABLE(I,J) represents the integral of X^E * Hf(I,X) * Hf(J,X).
//
{
  double *h_table;
  int i;
  int j;
  int k;
  int order;
  double *table;
  double *w_table;
  double x;
  double *x_table;

  table = new double[(p+1)*(p+1)];
  for ( j = 0; j <= p; j++ )
  {
    for ( i = 0; i <= p; i++ )
    {
      table[i+j*(p+1)] = 0.0;
    }
  }

  order = p + 1 + ( ( e + 1 ) / 2 );

  x_table = new double[order];
  w_table = new double[order];

  hf_quadrature_rule ( order, x_table, w_table );

  for ( k = 0; k < order; k++ )
  {
    x = x_table[k];
    h_table = hf_function_value ( 1, p, x_table+k );
//
//  The following formula is an outer product in H_TABLE.
//
    if ( e == 0 )
    {
      for ( j = 0; j <= p; j++ )
      {
        for ( i = 0; i <= p; i++ )
        {
          table[i+j*(p+1)] = table[i+j*(p+1)] 
            + w_table[k] * h_table[i] * h_table[j];
        }
      }
    }
    else
    {
      for ( j = 0; j <= p; j++ )
      {
        for ( i = 0; i <= p; i++ )
        {
          table[i+j*(p+1)] = table[i+j*(p+1)] 
            + w_table[k] * pow ( x, e ) * h_table[i] * h_table[j];
        }
      }
    }
    delete [] h_table;
  }

  delete [] w_table;
  delete [] x_table;

  return table;
}
//****************************************************************************80

void hf_quadrature_rule ( int nt, double t[], double wts[] )

//****************************************************************************80
//
//  Purpose:
//
//    HF_QUADRATURE_RULE: quadrature for Hf(i,x).
//
//  Discussion:
//
//    Hf(i,x) represents the Hermite function of "degree" I.   
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, the order of the rule.
//
//    Output, double T[NT], WTS[NT], the points and weights 
//    of the rule.
//
{
  double *bj;
  int i;
  const double r8_pi = 3.141592653589793;

  for ( i = 0; i < nt; i++ )
  {
    t[i] = 0.0;
  }

  bj = new double[nt];

  for ( i = 0; i < nt; i++ )
  {
    bj[i] = sqrt ( ( double ) ( i + 1 ) / 2.0 );
  }

  for ( i = 0; i < nt; i++ )
  {
    wts[i] = 0.0;
  }

  wts[0] = sqrt ( sqrt ( r8_pi ) );

  imtqlx ( nt, t, bj, wts );

  for ( i = 0; i < nt; i++ )
  {
    wts[i] = wts[i] * wts[i] * exp ( t[i] * t[i] );
  }
  return;
}
//****************************************************************************80

double *hn_exponential_product ( int p, double b )

//****************************************************************************80
//
//  Purpose:
//
//    HN_EXPONENTIAL_PRODUCT: exponential products exp(b*x)*Hn(i,x)*Hn(j,x).
//
//  Discussion:
//
//    Hn(i,x) is the normalized physicist's Hermite polynomial of degree I.  
//
//    For polynomial chaos applications, it is of interest to know the
//    value of the integrals of products of exp(B*X) with every possible pair
//    of basis functions.  That is, we'd like to form
//
//      Tij = Integral ( -oo < X < +oo ) 
//        exp(B*X) * Hn(I,X) * Hn(J,X) exp(-X*X) dx
//
//    We will estimate these integrals using Gauss-Hermite quadrature.
//    Because of the exponential factor exp(B*X), the quadrature will not 
//    be exact.
//
//    However, when B = 0, the quadrature is exact, and moreoever, the
//    table will be the identity matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the maximum degree of the polyonomial 
//    factors.  0 <= P.
//
//    Input, double B, the coefficient of X in the exponential factor.
//
//    Output, double HN_EXPONENTIAL_PRODUCT[(P+1)*(P+1)], the table of 
//    integrals.  TABLE(I,J) represents the weighted integral of 
//    exp(B*X) * Hn(I,X) * Hn(J,X).
//
{
  double *h_table;
  int i;
  int j;
  int k;
  int order;
  double *table;
  double *w_table;
  double x;
  double *x_table;

  table = new double[(p+1)*(p+1)];
  for ( j = 0; j <= p; j++ )
  {
    for ( i = 0; i <= p; i++ )
    {
      table[i+j*(p+1)] = 0.0;
    }
  }

  order = ( 3 * p + 4 ) / 2;

  x_table = new double[order];
  w_table = new double[order];

  h_quadrature_rule ( order, x_table, w_table );

  for ( k = 0; k < order; k++ )
  {
    x = x_table[k];
    h_table = hn_polynomial_value ( 1, p, x_table+k );
//
//  The following formula is an outer product in H_TABLE.
//
    for ( j = 0; j <= p; j++ )
    {
      for ( i = 0; i <= p; i++ )
      {
        table[i+j*(p+1)] = table[i+j*(p+1)] 
          + w_table[k] * exp ( b * x ) * h_table[i] * h_table[j];
      }
    }
    delete [] h_table;
  }

  delete [] w_table;
  delete [] x_table;

  return table;
}
//****************************************************************************80

double *hn_polynomial_value ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    HN_POLYNOMIAL_VALUE evaluates Hn(i,x).
//
//  Discussion:
//
//    Hn(i,x) is the normalized physicist's Hermite polynomial of degree I.  
//
//    These polynomials satisfy the orthonormality condition:
//
//      Integral ( -oo < X < +oo ) 
//        exp ( - X^2 ) * Hn(M,X) Hn(N,X) dX = delta ( N, M )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int M, the number of evaluation points.
//
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//
//    Input, double X[M], the evaluation points.
//
//    Output, double HN_POLYNOMIAL_VALUE[M*(N+1)], the values of the first 
//    N+1 Hermite polynomials at the evaluation points.
//
{
  double fact;
  int i;
  int j;
  double *p;
  const double r8_pi = 3.141592653589793;
  double two_power;

  if ( n < 0 )
  {
    return NULL;
  }

  p = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    p[i+0*m] = 1.0;
  }

  if ( n == 0 )
  {
    return p;
  }

  for ( i = 0; i < m; i++ )
  {
    p[i+1*m] = 2.0 * x[i];
  }

  for ( j = 2; j <= n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      p[i+j*m] = 2.0 * x[i] * p[i+(j-1)*m]
        - 2.0 * ( double ) ( j - 1 ) * p[i+(j-2)*m];
    }
  }
//
//  Normalize.
//
  fact = 1.0;
  two_power = 1.0;
  for ( j = 0; j <= n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      p[i+j*m] = p[i+j*m] / sqrt ( fact * two_power * sqrt ( r8_pi ) );
    }
    fact = fact * ( double ) ( j + 1 );
    two_power = two_power * 2.0;
  }
  return p;
}
//****************************************************************************80

double *hn_power_product ( int p, int e )

//****************************************************************************80
//
//  Purpose:
//
//    HN_POWER_PRODUCT: power products x^e*Hn(i,x)*Hn(j,x).
//
//  Discussion:
//
//    Hn(i,x) is the normalized physicist's Hermite polynomial of degree I.  
//
//    For polynomial chaos applications, it is of interest to know the
//    value of the integrals of products of X with every possible pair
//    of basis functions.  That is, we'd like to form
//
//      Tij = Integral ( -oo < X < +oo ) X^E * Hn(I,X) * Hn(J,X) exp(-X*X) dx
//
//    We will estimate these integrals using Gauss-Hermite quadrature.
//
//    When E is 0, the computed table should be the identity matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the maximum degree of the polyonomial 
//    factors.  0 <= P.
//
//    Input, int E, the exponent of X in the integrand.
//    0 <= E.
//
//    Output, double HN_POWER_PRODUCT[(P+1)*(P+1)], the table of 
//    integrals.  TABLE(I,J) represents the weighted integral of 
//    X^E * Hn(I,X) * Hn(J,X).
//
{
  double *h_table;
  int i;
  int j;
  int k;
  int order;
  double *table;
  double *w_table;
  double x;
  double *x_table;

  table = new double[(p+1)*(p+1)];
  for ( j = 0; j <= p; j++ )
  {
    for ( i = 0; i <= p; i++ )
    {
      table[i+j*(p+1)] = 0.0;
    }
  }

  order = p + 1 + ( ( e + 1 ) / 2 );

  x_table = new double[order];
  w_table = new double[order];

  h_quadrature_rule ( order, x_table, w_table );

  for ( k = 0; k < order; k++ )
  {
    x = x_table[k];
    h_table = hn_polynomial_value ( 1, p, x_table+k );
//
//  The following formula is an outer product in H_TABLE.
//
    if ( e == 0 )
    {
      for ( j = 0; j <= p; j++ )
      {
        for ( i = 0; i <= p; i++ )
        {
          table[i+j*(p+1)] = table[i+j*(p+1)] 
            + w_table[k] * h_table[i] * h_table[j];
        }
      }
    }
    else
    {
      for ( j = 0; j <= p; j++ )
      {
        for ( i = 0; i <= p; i++ )
        {
          table[i+j*(p+1)] = table[i+j*(p+1)] 
            + w_table[k] * pow ( x, e ) * h_table[i] * h_table[j];
        }
      }
    }
    delete [] h_table;
  }

  delete [] w_table;
  delete [] x_table;

  return table;
}
