# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

# include "tetrahedron_monte_carlo.hpp"

//****************************************************************************80

double *monomial_value ( int m, int n, int e[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONOMIAL_VALUE evaluates a monomial.
//
//  Discussion:
//
//    This routine evaluates a monomial of the form
//
//      product ( 1 <= i <= m ) x(i)^e(i)
//
//    where the exponents are nonnegative integers.  Note that
//    if the combination 0^0 is encountered, it should be treated
//    as 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int POINT_NUM, the number of points at which the
//    monomial is to be evaluated.
//
//    Input, int E[M], the exponents.
//
//    Input, double X[M*N], the point coordinates.
//
//    Output, double MONOMIAL_VALUE[N], the value of the monomial.
//
{
  int i;
  int j;
  double *v;

  v = new double[n];

  for ( j = 0; j < n; j++ )
  {
    v[j] = 1.0;
  }

  for ( i = 0; i < m; i++ )
  {
    if ( 0 != e[i] )
    {
      for ( j = 0; j < n; j++ )
      {
        v[j] = v[j] * pow ( x[i+j*m], e[i] );
      }
    }
  }

  return v;
}
//****************************************************************************80

double tetrahedron01_monomial_integral ( int e[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON01_MONOMIAL_INTEGRAL: integrals in the unit tetrahedron in 3D.
//
//  Discussion:
//
//    The monomial is F(X,Y,Z) = X^E(1) * Y^E(2) * Z^E(3).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int E[3], the exponents.  
//    Each exponent must be nonnegative.
//
//    Output, double TETRAHEDRON01_MONOMIAL_INTEGRAL, the integral.
//
{
  int i;
  double integral;
  int j;
  int k;
  const int m = 3;

  for ( i = 0; i < m; i++ )
  {
    if ( e[i] < 0 )
    {
      std::cerr << "\n";
      std::cerr << "TETRAHEDRON01_MONOMIAL_INTEGRAL - Fatal error!\n";
      std::cerr << "  All exponents must be nonnegative.\n";
      exit ( 1 );
    }
  }

  k = 0;
  integral = 1.0;

  for ( i = 0; i < m; i++ )
  {
    for ( j = 1; j <= e[i]; j++ )
    {
      k = k + 1;
      integral = integral * double( j ) / double( k );
    }
  }

  for ( i = 0; i < m; i++ )
  {
    k = k + 1;
    integral = integral / double( k );
  }

  return integral;
}
//****************************************************************************80

double *tetrahedron01_sample ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON01_SAMPLE samples the unit tetrahedron in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2015
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Reuven Rubinstein,
//    Monte Carlo Optimization, Simulation, and Sensitivity 
//    of Queueing Networks,
//    Krieger, 1992,
//    ISBN: 0894647644,
//    LC: QA298.R79.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Output, double TETRAHEDRON01_SAMPLE_01[3*N], the points.
//
{
  double *e;
  double e_sum;
  int i;
  int j;
  const int m = 3;
  double *x;

  x = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    e = r8vec_uniform_01_new ( m + 1, seed );

    for ( i = 0; i < m + 1; i++ )
    {
      e[i] = - log ( e[i] );
    }
    e_sum = r8vec_sum ( m + 1, e );

    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = e[i] / e_sum;
    }
    delete [] e;
  }

  return x;
}
//****************************************************************************80

double tetrahedron01_volume ( )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON01_VOLUME returns the volume of the unit tetrahedron in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double TETRAHEDRON01_VOLUME, the volume.
//
{
  double volume;

  volume = 1.0 / 6.0;

  return volume;
}
