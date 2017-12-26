# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

# include "triangle_monte_carlo.hpp"

//****************************************************************************80

double triangle01_area ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE01_AREA returns the area of the unit triangle in 2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 January 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double TRIANGLE01_AREA, the area
//
{
  double area;

  area = 0.5;

  return area;
}
//****************************************************************************80

double triangle01_monomial_integral ( int e[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE01_MONOMIAL_INTEGRAL: monomial integrals in the unit triangle in 2D.
//
//  Discussion:
//
//    The monomial is F(X,Y) = X^E(1) * Y^E(2).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int E[2], the exponents.  
//    Each exponent must be nonnegative.
//
//    Output, double TRIANGLE01_MONOMIAL_INTEGRAL, the integral.
//
{
  int i;
  double integral;
  int j;
  int k;
  const int m = 2;

  for ( i = 0; i < m; i++ )
  {
    if ( e[i] < 0 )
    {
      std::cerr << "\n";
      std::cerr << "TRIANGLE01_MONOMIAL_INTEGRAL - Fatal error!\n";
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

double *triangle01_sample ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE01_SAMPLE samples from the interior of the unit triangle in 2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 January 2014
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
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double TRIANGLE01_SAMPLE[2*N], the points.
//
{
  double *e;
  double e_sum;
  int i;
  int j;
  const int m = 2;
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
