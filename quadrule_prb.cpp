# include <cmath>
# include <cstdlib>
# include <cstring>
# include <iomanip>
# include <iostream>

# include "quadrule.hpp"

int main ( );
void chebyshev_set_test ( );
void chebyshev1_compute_test ( );
void chebyshev1_integral_test ( );
void chebyshev1_set_test ( );
void chebyshev2_compute_test ( );
void chebyshev2_compute_test2 ( );
void chebyshev2_integral_test ( );
void chebyshev2_set_test ( );
void chebyshev3_compute_test ( );
void chebyshev3_integral_test ( );
void chebyshev3_set_test ( );
void clenshaw_curtis_compute_test ( );
void clenshaw_curtis_set_test ( );
void fejer1_compute_test ( );
void fejer1_set_test ( );
void fejer2_compute_test ( );
void fejer2_set_test ( );
void gegenbauer_integral_test ( );
void gegenbauer_ek_compute_test ( );
void gegenbauer_ss_compute_test ( );
void gen_hermite_ek_compute_test ( );
void gen_hermite_integral_test ( );
void gen_laguerre_ek_compute_test ( );
void gen_laguerre_integral_test ( );
void gen_laguerre_ss_compute_test ( );
void hermite_ek_compute_test ( );
void hermite_integral_test ( );
void hermite_set_test ( );
void hermite_ss_compute_test ( );
void hermite_gk16_set_test ( );
void hermite_gk18_set_test ( );
void hermite_gk22_set_test ( );
void hermite_gk24_set_test ( );
void hermite_1_set_test ( );
void hermite_probabilist_set_test ( );
void imtqlx_test ( );
void jacobi_ek_compute_test ( );
void jacobi_integral_test ( );
void jacobi_ss_compute_test ( );
void kronrod_set_test ( );
void laguerre_ek_compute_test ( );
void laguerre_integral_test ( );
void laguerre_set_test ( );
void laguerre_ss_compute_test ( );
void laguerre_1_set_test ( );
void legendre_dr_compute_test ( );
void legendre_ek_compute_test ( );
void legendre_integral_test ( );
void legendre_set_test ( );
void lobatto_compute_test ( );
void lobatto_set_test ( );
void nc_compute_weights_test ( );
void ncc_compute_test ( );
void ncc_set_test ( );
void nco_compute_test ( );
void nco_set_test ( );
void ncoh_compute_test ( );
void ncoh_set_test ( );
void patterson_set_test ( );
void r8_psi_test ( );
void radau_compute_test ( );
void radau_set_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for QUADRULE_PRB.
//
//  Discussion:
//
//    QUADRULE_PRB tests the QUADRULE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2015
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  double alpha;
  int n;
  int order;

  std::cout << "\n";
  std::cout << "QUADRULE_PRB\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the QUADRULE library.\n";

  chebyshev_set_test ( );
  chebyshev1_compute_test ( );
  chebyshev1_integral_test ( );
  chebyshev1_set_test ( );
  chebyshev2_compute_test ( );
  chebyshev2_compute_test2 ( );
  chebyshev2_integral_test ( );
  chebyshev2_set_test ( );
  chebyshev3_compute_test ( );
  chebyshev3_integral_test ( );
  chebyshev3_set_test ( );
  clenshaw_curtis_compute_test ( );
  clenshaw_curtis_set_test ( );
  fejer1_compute_test ( );
  fejer1_set_test ( );
  fejer2_compute_test ( );
  fejer2_set_test ( );
  gegenbauer_integral_test ( );
  gegenbauer_ek_compute_test ( );
  gegenbauer_ss_compute_test ( );
  gen_hermite_ek_compute_test ( );
  gen_hermite_integral_test ( );
  gen_laguerre_ek_compute_test ( );
  gen_laguerre_integral_test ( );
  gen_laguerre_ss_compute_test ( );
  hermite_ek_compute_test ( );
  hermite_integral_test ( );
  hermite_set_test ( );
  hermite_ss_compute_test ( );
  hermite_gk16_set_test ( );
  hermite_gk18_set_test ( );
  hermite_gk22_set_test ( );
  hermite_gk24_set_test ( );
  hermite_1_set_test ( );
  hermite_probabilist_set_test ( );
  imtqlx_test ( );
  jacobi_ek_compute_test ( );
  jacobi_integral_test ( );
  jacobi_ss_compute_test ( );
  kronrod_set_test ( );
  laguerre_ek_compute_test ( );
  laguerre_integral_test ( );
  laguerre_set_test ( );
  laguerre_ss_compute_test ( );
  laguerre_1_set_test ( );
  legendre_dr_compute_test ( );
  legendre_ek_compute_test ( );
  legendre_integral_test ( );
  legendre_set_test ( );
  lobatto_compute_test ( );
  lobatto_set_test ( );
  nc_compute_weights_test ( );
  ncc_compute_test ( );
  ncc_set_test ( );
  nco_compute_test ( );
  nco_set_test ( );
  ncoh_compute_test ( );
  ncoh_set_test ( );
  patterson_set_test ( );
  r8_psi_test ( );
  radau_compute_test ( );
  radau_set_test ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "QUADRULE_PRB\n";
  std::cout << "  Normal end of execution.\n";
  std::cout << "\n";
  timestamp ( );
 
  return 0;
}
//****************************************************************************80

void chebyshev_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV_SET_TEST tests CHEBYSHEV_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "CHEBYSHEV_SET_TEST\n";
  std::cout << "  CHEBYSHEV_SET sets\n";
  std::cout << "  a Chebyshev quadrature rule over [-1,1].\n";

  std::cout << "\n";
  std::cout << "  Index             X                   W\n";
  std::cout << "\n";

  for ( n = 1; n <= 9; n++ )
  {
    if ( n == 8 )
    {
      continue;
    }
    w = new double[n];
    x = new double[n];

    chebyshev_set ( n, x, w );

    std::cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void chebyshev1_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV1_COMPUTE_TEST tests CHEBYSHEV1_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "CHEBYSHEV1_COMPUTE_TEST\n";
  std::cout << "  CHEBYSHEV1_COMPUTE computes\n";
  std::cout << "  a Chebyshev Type 1 quadrature rule over [-1,1].\n";
  std::cout << "\n";
  std::cout << "  Index             X                   W\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    chebyshev1_compute ( n, x, w );

    std::cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void chebyshev1_integral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV1_INTEGRAL_TEST tests CHEBYSHEV1_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int prec;
  double value;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "CHEBYSHEV1_INTEGRAL_TEST\n";
  std::cout << "  CHEBYSHEV1_INTEGRAL evaluates\n";
  std::cout << "  Integral ( -1 < x < +1 ) x^n / sqrt(1-x*x) dx\n";
  std::cout << "\n";
  std::cout << "         N         Value\n";
  std::cout << "\n";

  for ( n = 0; n <= 10; n++ )
  {
    value = chebyshev1_integral ( n );
    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(24) << value << "\n";
  }
 
  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void chebyshev1_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV1_SET_TEST tests CHEBYSHEV1_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "CHEBYSHEV1_SET_TEST\n";
  std::cout << "  CHEBYSHEV1_SET sets\n";
  std::cout << "  a Chebyshev Type 1 quadrature rule over [-1,1].\n";

  std::cout << "\n";
  std::cout << "  Index             X                   W\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    chebyshev1_set ( n, x, w );

    std::cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void chebyshev2_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV2_COMPUTE_TEST tests CHEBYSHEV2_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "CHEBYSHEV2_COMPUTE_TEST\n";
  std::cout << "  CHEBYSHEV2_COMPUTE computes\n";
  std::cout << "  a Chebyshev Type 2 quadrature rule over [-1,1].\n";
  std::cout << "\n";
  std::cout << "  Index             X                   W\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    chebyshev2_compute ( n, x, w );

    std::cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void chebyshev2_compute_test2 ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV2_COMPUTE_TEST uses CHEBYSHEV2_COMPUTE over the semicircle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  double error;
  double exact;
  double *f;
  int i;
  int n = 10;
  int prec;
  double q;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "CHEBYSHEV2_COMPUTE_TEST2\n";
  std::cout << "  Approximate the integral of f(x,y) over the semicircle\n";
  std::cout << "    -1 <= x <= 1, y = sqrt ( 1 - x^2 )\n";
  std::cout << "  using N Chebyshev points.\n";
  std::cout << "  If p(x,y) involves any term of odd degree in y,\n";
  std::cout << "  the estimate will only be approximate.\n";
  std::cout << "\n";
  std::cout << "  Polynomial    N    Integral        Estimate       Error\n";
  std::cout << "\n";

  f = new double[n];
  x = new double[n];
  w = new double[n];

  chebyshev2_compute ( n, x, w );
//
//  f(x,y) = 1
//
  exact = 1.5707963267948966192;
  for ( i = 0; i < n; i++ )
  {
    f[i] = 1.0;
  }
  q = r8vec_dot_product ( n, w, f );
  error = r8_abs ( q - exact );
  std::cout << "  1            "
       << "  " << std::setw(2) << n
       << "  " << std::setw(14) << exact
       << "  " << std::setw(14) << q
       << "  " << std::setw(14) << error << "\n";
//
//  f(x,y) = x
//
  exact = 0.0;
  for ( i = 0; i < n; i++ )
  {
    f[i] = x[i];
  }
  q = r8vec_dot_product ( n, w, f );
  error = r8_abs ( q - exact );
  std::cout << "  x            "
       << "  " << std::setw(2) << n
       << "  " << std::setw(14) << exact
       << "  " << std::setw(14) << q
       << "  " << std::setw(14) << error << "\n";
//
//  f(x,y) = y = sqrt ( 1 - x^2 )
//
  exact = 0.66666666666666666667;
  for ( i = 0; i < n; i++ )
  {
    f[i] = sqrt ( 1.0 - pow ( x[i], 2 ) );
  }
  q = r8vec_dot_product ( n, w, f ) / 2.0;
  error = r8_abs ( q - exact );
  std::cout << "     y         "
       << "  " << std::setw(2) << n
       << "  " << std::setw(14) << exact
       << "  " << std::setw(14) << q
       << "  " << std::setw(14) << error << "\n";
//
//  f(x,y) = x^2
//
  exact = 0.39269908169872415481;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( x[i], 2 );
  }
  q = r8vec_dot_product ( n, w, f );
  error = r8_abs ( q - exact );
  std::cout << "  x^2          "
       << "  " << std::setw(2) << n
       << "  " << std::setw(14) << exact
       << "  " << std::setw(14) << q
       << "  " << std::setw(14) << error << "\n";
//
//  f(x,y) = xy = x * sqrt ( 1 - x^2 )
//
  exact = 0.0;
  for ( i = 0; i < n; i++ )
  {
    f[i] = x[i] * sqrt ( 1.0 - pow ( x[i], 2 ) );
  }
  q = r8vec_dot_product ( n, w, f ) / 2.0;
  error = r8_abs ( q - exact );
  std::cout << "  x  y         "
       << "  " << std::setw(2) << n
       << "  " << std::setw(14) << exact
       << "  " << std::setw(14) << q
       << "  " << std::setw(14) << error << "\n";
//
//  f(x,y) = y^2 -> ( 1 - x^2 )
//
  exact = 0.39269908169872415481;
  for ( i = 0; i < n; i++ )
  {
    f[i] = 1.0 - pow ( x[i], 2 );
  }
  q = r8vec_dot_product ( n, w, f ) / 3.0;
  error = r8_abs ( q - exact );
  std::cout << "     y^2       "
       << "  " << std::setw(2) << n
       << "  " << std::setw(14) << exact
       << "  " << std::setw(14) << q
       << "  " << std::setw(14) << error << "\n";
//
//  f(x,y) = x^3
//
  exact = 0.0;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( x[i], 3 );
  }
  q = r8vec_dot_product ( n, w, f );
  error = r8_abs ( q - exact );
  std::cout << "  x^3          "
       << "  " << std::setw(2) << n
       << "  " << std::setw(14) << exact
       << "  " << std::setw(14) << q
       << "  " << std::setw(14) << error << "\n";
//
//  f(x,y) = x^2 y = x^2 sqrt ( 1 - x^2 )
//
  exact = 0.13333333333333333333;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( x[i], 2 ) * sqrt ( 1.0 - pow ( x[i], 2 ) );
  }
  q = r8vec_dot_product ( n, w, f ) / 2.0;
  error = r8_abs ( q - exact );
  std::cout << "  x^2y         "
       << "  " << std::setw(2) << n
       << "  " << std::setw(14) << exact
       << "  " << std::setw(14) << q
       << "  " << std::setw(14) << error << "\n";
//
//  f(x,y) = x y^2 = x * ( 1 - x^2 )
//
  exact = 0.0;
  for ( i = 0; i < n; i++ )
  {
    f[i] = x[i] * ( 1.0 - pow ( x[i], 2 ) );
  }
  q = r8vec_dot_product ( n, w, f ) / 3.0;
  error = r8_abs ( q - exact );
  std::cout << "  x  y^2       "
       << "  " << std::setw(2) << n
       << "  " << std::setw(14) << exact
       << "  " << std::setw(14) << q
       << "  " << std::setw(14) << error << "\n";
//
//  f(x,y) = y^3
//
  exact = 0.26666666666666666667;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( 1.0 - pow ( x[i], 2 ), 1.5 );
  }
  q = r8vec_dot_product ( n, w, f ) / 4.0;
  error = r8_abs ( q - exact );
  std::cout << "     y^3       "
       << "  " << std::setw(2) << n
       << "  " << std::setw(14) << exact
       << "  " << std::setw(14) << q
       << "  " << std::setw(14) << error << "\n";
//
//  f(x,y) = x^4
//
  exact = 0.19634954084936207740;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( x[i], 4 );
  }
  q = r8vec_dot_product ( n, w, f );
  error = r8_abs ( q - exact );
  std::cout << "  x^4          "
       << "  " << std::setw(2) << n
       << "  " << std::setw(14) << exact
       << "  " << std::setw(14) << q
       << "  " << std::setw(14) << error << "\n";
//
//  f(x,y) = x^2y^2 -> x^2( 1 - x^2 )
//
  exact = 0.065449846949787359135;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( x[i], 2 ) * ( 1.0 - pow ( x[i], 2 ) );
  }
  q = r8vec_dot_product ( n, w, f ) / 3.0;
  error = r8_abs ( q - exact );
  std::cout << "  x^2y^2       "
       << "  " << std::setw(2) << n
       << "  " << std::setw(14) << exact
       << "  " << std::setw(14) << q
       << "  " << std::setw(14) << error << "\n";
//
//  f(x,y) = y^4 -> ( 1 - x^2 )^2
//
  exact = 0.19634954084936207740;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( 1.0 - pow ( x[i], 2 ), 2 );
  }
  q = r8vec_dot_product ( n, w, f ) / 5.0;
  error = r8_abs ( q - exact );
  std::cout << "     y^4       "
       << "  " << std::setw(2) << n
       << "  " << std::setw(14) << exact
       << "  " << std::setw(14) << q
       << "  " << std::setw(14) << error << "\n";
//
//  f(x,y) = x^4y = x^4 sqrt ( 1 - x^2 )
//
  exact = 0.057142857142857142857;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( x[i], 4 ) * sqrt ( 1.0 - pow ( x[i], 2 ) );
  }
  q = r8vec_dot_product ( n, w, f ) / 2.0;
  error = r8_abs ( q - exact );
  std::cout << "  x^4y         "
       << "  " << std::setw(2) << n
       << "  " << std::setw(14) << exact
       << "  " << std::setw(14) << q
       << "  " << std::setw(14) << error << "\n";
//
//  x^2y^3 = x^2 ( 1 - x^2 )^(3/2)
//
  exact = 0.038095238095238095238;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( x[i], 2 ) * pow ( 1.0 - pow ( x[i], 2 ), 1.5 );
  }
  q = r8vec_dot_product ( n, w, f ) / 4.0;
  error = r8_abs ( q - exact );
  std::cout << "  x^2y^3       "
       << "  " << std::setw(2) << n
       << "  " << std::setw(14) << exact
       << "  " << std::setw(14) << q
       << "  " << std::setw(14) << error << "\n";
//
//  f(x,y) = y^5
//
  exact = 0.15238095238095238095;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( 1.0 - pow ( x[i], 2 ), 2.5 );
  }
  q = r8vec_dot_product ( n, w, f ) / 6.0;
  error = r8_abs ( q - exact );
  std::cout << "     y^5       "
       << "  " << std::setw(2) << n
       << "  " << std::setw(14) << exact
       << "  " << std::setw(14) << q
       << "  " << std::setw(14) << error << "\n";
//
//  f(x,y) = x^6
//
  exact = 0.12271846303085129838;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( x[i], 6 );
  }
  q = r8vec_dot_product ( n, w, f );
  error = r8_abs ( q - exact );
  std::cout << "  x^6          "
       << "  " << std::setw(2) << n
       << "  " << std::setw(14) << exact
       << "  " << std::setw(14) << q
       << "  " << std::setw(14) << error << "\n";
//
//  f(x,y) = x^4y^2 -> x^4( 1 - x^2 )
//
  exact = 0.024543692606170259675;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( x[i], 4 ) * ( 1.0 - pow ( x[i], 2 ) );
  }
  q = r8vec_dot_product ( n, w, f ) / 3.0;
  error = r8_abs ( q - exact );
  std::cout << "  x^4y^2       "
       << "  " << std::setw(2) << n
       << "  " << std::setw(14) << exact
       << "  " << std::setw(14) << q
       << "  " << std::setw(14) << error << "\n";
//
//  f(x,y) = x^2y^4 -> x^2( 1 - x^2 )^2
//
  exact = 0.024543692606170259675;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( x[i], 2 ) * pow ( 1.0 - pow ( x[i], 2 ), 2 );
  }
  q = r8vec_dot_product ( n, w, f ) / 5.0;
  error = r8_abs ( q - exact );
  std::cout << "  x^2y^4       "
       << "  " << std::setw(2) << n
       << "  " << std::setw(14) << exact
       << "  " << std::setw(14) << q
       << "  " << std::setw(14) << error << "\n";
//
//  f(x,y) = y^6 -> ( 1 - x^2 )^3
//
  exact = 0.12271846303085129838;
  for ( i = 0; i < n; i++ )
  {
    f[i] = pow ( 1.0 - pow ( x[i], 2 ), 3 );
  }
  q = r8vec_dot_product ( n, w, f ) / 7.0;
  error = r8_abs ( q - exact );
  std::cout << "     y^6       "
       << "  " << std::setw(2) << n
       << "  " << std::setw(14) << exact
       << "  " << std::setw(14) << q
       << "  " << std::setw(14) << error << "\n";

  delete [] f;
  delete [] w;
  delete [] x;

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void chebyshev2_integral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV2_INTEGRAL_TEST tests CHEBYSHEV2_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int prec;
  double value;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "CHEBYSHEV2_INTEGRAL_TEST\n";
  std::cout << "  CHEBYSHEV2_INTEGRAL evaluates\n";
  std::cout << "  Integral ( -1 < x < +1 ) x^n * sqrt(1-x*x) dx\n";
  std::cout << "\n";
  std::cout << "         N         Value\n";
  std::cout << "\n";

  for ( n = 0; n <= 10; n++ )
  {
    value = chebyshev2_integral ( n );
    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(24) << value << "\n";
  }
 
  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void chebyshev2_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV2_SET_TEST tests CHEBYSHEV2_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "CHEBYSHEV2_SET_TEST\n";
  std::cout << "  CHEBYSHEV2_SET sets\n";
  std::cout << "  a Chebyshev Type 2 quadrature rule over [-1,1].\n";

  std::cout << "\n";
  std::cout << "  Index             X                   W\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    chebyshev2_set ( n, x, w );

    std::cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void chebyshev3_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV3_COMPUTE_TEST tests CHEBYSHEV3_COMPUTE
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "CHEBYSHEV3_COMPUTE_TEST\n";
  std::cout << "  CHEBYSHEV3_COMPUTE computes\n";
  std::cout << "  a Chebyshev Type 3 quadrature rule over [-1,1].\n";

  std::cout << "\n";
  std::cout << "  Index             X                   W\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    chebyshev3_compute ( n, x, w );

    std::cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void chebyshev3_integral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV3_INTEGRAL_TEST tests CHEBYSHEV3_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int prec;
  double value;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "CHEBYSHEV3_INTEGRAL_TEST\n";
  std::cout << "  CHEBYSHEV3_INTEGRAL evaluates\n";
  std::cout << "  Integral ( -1 < x < +1 ) x^n / sqrt(1-x*x) dx\n";
  std::cout << "\n";
  std::cout << "         N         Value\n";
  std::cout << "\n";

  for ( n = 0; n <= 10; n++ )
  {
    value = chebyshev3_integral ( n );
    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(24) << value << "\n";
  }
 
  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void chebyshev3_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV3_SET_TEST tests CHEBYSHEV3_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "CHEBYSHEV3_SET_TEST\n";
  std::cout << "  CHEBYSHEV3_SET sets\n";
  std::cout << "  a Chebyshev Type 3 quadrature rule over [-1,1].\n";

  std::cout << "\n";
  std::cout << "  Index             X                   W\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    chebyshev3_set ( n, x, w );

    std::cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void clenshaw_curtis_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CLENSHAW_CURTIS_COMPUTE_TEST tests CLENSHAW_CURTIS_COMPUTE
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2006
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "CLENSHAW_CURTIS_COMPUTE_TEST\n";
  std::cout << "  CLENSHAW_CURTIS_COMPUTE computes\n";
  std::cout << "  a Clenshaw-Curtis quadrature rule over [-1,1].\n";
  std::cout << "\n";
  std::cout << "  Index             X                   W\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    clenshaw_curtis_compute ( n, x, w );

    std::cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void clenshaw_curtis_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CLENSHAW_CURTIS_SET_TEST tests CLENSHAW_CURTIS_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "CLENSHAW_CURTIS_SET_TEST\n";
  std::cout << "  CLENSHAW_CURTIS_SET sets up a Clenshaw-Curtis rule;\n";
  std::cout << "  a Clenshaw-Curtis quadrature rule over [-1,1]\n";
  std::cout << "  of given order.\n";
  std::cout << "\n";
  std::cout << "  Index             X                   W\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    clenshaw_curtis_set ( n, x, w );

    std::cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void fejer1_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    FEJER1_COMPUTE_TEST tests FEJER1_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int n_max = 10;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "FEJER1_COMPUTE_TEST\n";
  std::cout << "  FEJER1_COMPUTE computes a Fejer type 1 quadrature rule;\n";
  std::cout << "\n";
  std::cout << "     Order        W               X\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    fejer1_compute ( n, x, w );

    std::cout << "\n";
    std::cout << "  " << std::setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "          "
           << "  " << std::setw(14) << w[i]
           << "  " << std::setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void fejer1_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    FEJER1_SET_TEST tests FEJER1_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "FEJER1_SET_TEST\n";
  std::cout << "  FEJER1_SET sets a Fejer type 1 quadrature rule;\n";
  std::cout << "\n";
  std::cout << "     Order        W               X\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    fejer1_set ( n, x, w );

    std::cout << "\n";
    std::cout << "  " << std::setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "          "
           << "  " << std::setw(14) << w[i]
           << "  " << std::setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void fejer2_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    FEJER2_COMPUTE_TEST tests FEJER2_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "FEJER2_COMPUTE_TEST\n";
  std::cout << "  FEJER2_COMPUTE computes a Fejer type 2 quadrature rule;\n";
  std::cout << "\n";
  std::cout << "     Order        W               X\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    fejer2_compute ( n, x, w );

    std::cout << "\n";
    std::cout << "  " << std::setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "          "
           << "  " << std::setw(14) << w[i]
           << "  " << std::setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void fejer2_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    FEJER2_SET_TEST tests FEJER2_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "FEJER2_SET_TEST\n";
  std::cout << "  FEJER2_SET sets a Fejer type 2 quadrature rule;\n";
  std::cout << "\n";
  std::cout << "     Order        W               X\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    fejer2_set ( n, x, w );

    std::cout << "\n";
    std::cout << "  " << std::setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "          "
           << "  " << std::setw(14) << w[i]
           << "  " << std::setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void gegenbauer_integral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GEGENBAUER_INTEGRAL_TEST tests GEGENBAUER_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  int n;
  int prec;
  double value;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  alpha = 0.25;

  std::cout << "\n";
  std::cout << "GEGENBAUER_INTEGRAL_TEST\n";
  std::cout << "  GEGENBAUER_INTEGRAL evaluates\n";
  std::cout << "  Integral ( -1 < x < +1 ) x^n * (1-x*x)^alpha dx\n";
  std::cout << "\n";
  std::cout << "         N         Value\n";
  std::cout << "\n";

  for ( n = 0; n <= 10; n++ )
  {
    value = gegenbauer_integral ( n, alpha );
    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(24) << value << "\n";
  }

  std::cout.precision ( prec );
 
  return;
}
//****************************************************************************80

void gegenbauer_ek_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GEGENBAUER_EK_COMPUTE_TEST tests GEGENBAUER_EK_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double alpha;
  double b;
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  alpha = 0.5;
  a = -1.0;
  b = +1.0;

  std::cout << "\n";
  std::cout << "GEGENBAUER_EK_COMPUTE_TEST\n";
  std::cout << "  GEGENBAUER_EK_COMPUTE computes a Gauss-Gegenbauer rule;\n";
  std::cout << "\n";
  std::cout << "  with ALPHA = " << alpha << "\n";
  std::cout << "  and integration interval [" << a << ", " << b << "]\n";
  std::cout << "\n";
  std::cout << "                  W               X\n";

  for ( n = 1; n <= 10; n++ )
  {
    std::cout << "\n";

    w = new double[n];
    x = new double[n];

    gegenbauer_ek_compute ( n, alpha, a, b, x, w );

    for ( i = 0; i < n; i++ )
    {
      std::cout << "          "
           << "  " << std::setw(14) << w[i]
           << "  " << std::setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void gegenbauer_ss_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GEGENBAUER_SS_COMPUTE_TEST tests GEGENBAUER_SS_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  alpha = 0.5;

  std::cout << "\n";
  std::cout << "GEGENBAUER_SS_COMPUTE_TEST\n";
  std::cout << "  GEGENBAUER_SS_COMPUTE computes a Gauss-Gegenbauer rule;\n";
  std::cout << "\n";
  std::cout << "  with ALPHA = " << alpha << "\n";
  std::cout << "\n";
  std::cout << "                  W               X\n";

  for ( n = 1; n <= 10; n++ )
  {
    std::cout << "\n";

    w = new double[n];
    x = new double[n];

    gegenbauer_ss_compute ( n, alpha, x, w );

    for ( i = 0; i < n; i++ )
    {
      std::cout << "          "
           << "  " << std::setw(14) << w[i]
           << "  " << std::setw(14) << x[i] << "\n";
    }

    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void gen_hermite_ek_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GEN_HERMITE_EK_COMPUTE_TEST tests GEN_HERMITE_EK_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  alpha = 0.5;

  std::cout << "\n";
  std::cout << "GEN_HERMITE_EK_COMPUTE_TEST\n";
  std::cout << "  GEN_HERMITE_EK_COMPUTE computes \n";
  std::cout << "  a generalized Hermite quadrature rule\n";
  std::cout << "  using the Elhay-Kautsky algorithm.\n";
  std::cout << "\n";
  std::cout << "  Using ALPHA = " << alpha << "\n";
  std::cout << "\n";
  std::cout << "     Order        W               X\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    gen_hermite_ek_compute ( n, alpha, x, w );

    std::cout << "\n";
    std::cout << "  " << std::setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "          "
           << "  " << std::setw(14) << w[i]
           << "  " << std::setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void gen_hermite_integral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GEN_HERMITE_INTEGRAL_TEST tests GEN_HERMITE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  int n;
  int prec;
  double value;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  alpha = 0.5;

  std::cout << "\n";
  std::cout << "GEN_HERMITE_INTEGRAL_TEST\n";
  std::cout << "  GEN_HERMITE_INTEGRAL evaluates\n";
  std::cout << "  Integral ( -oo < x < +oo ) exp(-x^2) x^n |x|^alpha dx\n";
  std::cout << "\n";
  std::cout << "  Use ALPHA = " << alpha << "\n";
  std::cout << "\n";
  std::cout << "         N         Value\n";
  std::cout << "\n";

  for ( n = 0; n <= 10; n++ )
  {
    value = gen_hermite_integral ( n, alpha );
    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(24) << value << "\n";
  }
 
  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void gen_laguerre_ek_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GEN_LAGUERRE_EK_COMPUTE_TEST tests GEN_LAGUERRE_EK_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  alpha = 0.5;

  std::cout << "\n";
  std::cout << "GEN_LAGUERRE_EK_COMPUTE_TEST\n";
  std::cout << "  GEN_LAGUERRE_EK_COMPUTE computes \n";
  std::cout << "  a generalized Laguerre quadrature rule\n";
  std::cout << "  using the Elhay-Kautsky algorithm.\n";
  std::cout << "\n";
  std::cout << "  Using ALPHA = " << alpha << "\n";
  std::cout << "\n";
  std::cout << "     Order        W               X\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    gen_laguerre_ek_compute ( n, alpha, x, w );

    std::cout << "\n";
    std::cout << "  " << std::setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "          "
           << "  " << std::setw(14) << w[i]
           << "  " << std::setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void gen_laguerre_integral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GEN_LAGUERRE_INTEGRAL_TEST tests GEN_LAGUERRE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  int n;
  int prec;
  double value;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  alpha = 0.5;

  std::cout << "\n";
  std::cout << "GEN_LAGUERRE_INTEGRAL_TEST\n";
  std::cout << "  GEN_LAGUERRE_INTEGRAL evaluates\n";
  std::cout << "  Integral ( 0 < x < +oo ) exp(-x) x^n x^alpha dx\n";
  std::cout << "\n";
  std::cout << "  Use ALPHA = " << alpha << "\n";
  std::cout << "\n";
  std::cout << "         N         Value\n";
  std::cout << "\n";

  for ( n = 0; n <= 10; n++ )
  {
    value = gen_laguerre_integral ( n, alpha );
    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(24) << value << "\n";
  }
 
  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void gen_laguerre_ss_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GEN_LAGUERRE_SS_COMPUTE_TEST tests GEN_LAGUERRE_SS_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  alpha = 0.5;

  std::cout << "\n";
  std::cout << "GEN_LAGUERRE_SS_COMPUTE_TEST\n";
  std::cout << "  GEN_LAGUERRE_SS_COMPUTE computes \n";
  std::cout << "  a generalized Laguerre quadrature rule\n";
  std::cout << "  using the Stroud-Secrest algorithm.\n";
  std::cout << "\n";
  std::cout << "  Using ALPHA = " << alpha << "\n";
  std::cout << "\n";
  std::cout << "     Order        W               X\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    gen_laguerre_ss_compute ( n, alpha, x, w );

    std::cout << "\n";
    std::cout << "  " << std::setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "          "
           << "  " << std::setw(14) << w[i]
           << "  " << std::setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void hermite_ek_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_EK_COMPUTE_TEST tests HERMITE_EK_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "HERMITE_EK_COMPUTE_TEST\n";
  std::cout << "  HERMITE_EK_COMPUTE computes a Hermite quadrature rule\n";
  std::cout << "  using the Elhay-Kautsky algorithm.\n";
  std::cout << "\n";
  std::cout << "     Order        W               X\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    hermite_ek_compute ( n, x, w );

    std::cout << "\n";
    std::cout << "  " << std::setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "          "
           << "  " << std::setw(14) << w[i]
           << "  " << std::setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void hermite_integral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_INTEGRAL_TEST tests HERMITE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int prec;
  double value;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "HERMITE_INTEGRAL_TEST\n";
  std::cout << "  HERMITE_INTEGRAL evaluates\n";
  std::cout << "  Integral ( -oo < x < +oo ) exp(-x^2) x^n dx\n";
  std::cout << "\n";
  std::cout << "         N         Value\n";
  std::cout << "\n";

  for ( n = 0; n <= 10; n++ )
  {
    value = hermite_integral ( n );
    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(24) << value << "\n";
  }
 
  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void hermite_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_SET_TEST tests HERMITE_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "HERMITE_SET_TEST\n";
  std::cout << "  HERMITE_SET sets a Hermite quadrature rule over (-oo,+oo).\n";
  std::cout << "\n";
  std::cout << "  Index             X                   W\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    hermite_set ( n, x, w );

    std::cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void hermite_ss_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_SS_COMPUTE_TEST tests HERMITE_SS_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "HERMITE_SS_COMPUTE_TEST\n";
  std::cout << "  HERMITE_SS_COMPUTE computes a Hermite quadrature rule\n";
  std::cout << "  using the Stroud-Secrest algorithm.\n";
  std::cout << "\n";
  std::cout << "     Order        W               X\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    hermite_ss_compute ( n, x, w );

    std::cout << "\n";
    std::cout << "  " << std::setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "          "
           << "  " << std::setw(14) << w[i]
           << "  " << std::setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void hermite_gk16_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_GK16_SET_TEST tests HERMITE_GK16_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int l_max = 8;
  int l;
  int n;
  int n_list[9] = { 1, 3, 7, 9, 17, 19, 31, 33, 35 };
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "HERMITE_GK16_SET_TEST\n";
  std::cout << "  HERMITE_GK16_SET sets a nested rule\n";
  std::cout << "  for the Hermite integration problem.\n";
  std::cout << "\n";
  std::cout << "     Order        W               X\n";

  for ( l = 0; l <= l_max; l++ )
  {
    n = n_list[l];

    w = new double[n];
    x = new double[n];

    hermite_gk16_set ( n, x, w );

    std::cout << "\n";
    std::cout << "  " << std::setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "          "
           << "  " << std::setw(14) << w[i]
           << "  " << std::setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }
  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void hermite_gk18_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_GK18_SET_TEST tests HERMITE_GK18_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int l_max = 4;
  int l;
  int n;
  int n_list[9] = { 1, 3, 9, 19, 37 };
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "HERMITE_GK18_SET_TEST\n";
  std::cout << "  HERMITE_GK18_SET sets a nested rule\n";
  std::cout << "  for the Hermite integration problem.\n";
  std::cout << "\n";
  std::cout << "     Order        W               X\n";

  for ( l = 0; l <= l_max; l++ )
  {
    n = n_list[l];

    w = new double[n];
    x = new double[n];

    hermite_gk18_set ( n, x, w );

    std::cout << "\n";
    std::cout << "  " << std::setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "          "
           << "  " << std::setw(14) << w[i]
           << "  " << std::setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }
  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void hermite_gk22_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_GK22_SET_TEST tests HERMITE_GK22_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int l_max = 4;
  int l;
  int n;
  int n_list[9] = { 1, 3, 9, 19, 41 };
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "HERMITE_GK22_SET_TEST\n";
  std::cout << "  HERMITE_GK22_SET sets a nested rule\n";
  std::cout << "  for the Hermite integration problem.\n";
  std::cout << "\n";
  std::cout << "     Order        W               X\n";

  for ( l = 0; l <= l_max; l++ )
  {
    n = n_list[l];

    w = new double[n];
    x = new double[n];

    hermite_gk22_set ( n, x, w );

    std::cout << "\n";
    std::cout << "  " << std::setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "          "
           << "  " << std::setw(14) << w[i]
           << "  " << std::setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }
  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void hermite_gk24_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_GK24_SET_TEST tests HERMITE_GK24_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int l_max = 4;
  int l;
  int n;
  int n_list[9] = { 1, 3, 9, 19, 43 };
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "HERMITE_GK24_SET_TEST\n";
  std::cout << "  HERMITE_GK24_SET sets a nested rule\n";
  std::cout << "  for the Hermite integration problem.\n";
  std::cout << "\n";
  std::cout << "     Order        W               X\n";

  for ( l = 0; l <= l_max; l++ )
  {
    n = n_list[l];

    w = new double[n];
    x = new double[n];

    hermite_gk24_set ( n, x, w );

    std::cout << "\n";
    std::cout << "  " << std::setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "          "
           << "  " << std::setw(14) << w[i]
           << "  " << std::setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }
  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void hermite_1_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_1_SET_TEST tests HERMITE_1_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2013
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "HERMITE_1_SET_TEST\n";
  std::cout << "  HERMITE_1_SET sets a unit density Hermite quadrature rule.\n";
  std::cout << "  The integration interval is ( -oo, +oo ).\n";
  std::cout << "  The weight function is 1.\n";
  std::cout << "\n";
  std::cout << "  Index             X                   W\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    hermite_1_set ( n, x, w );

    std::cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void hermite_probabilist_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_PROBABILIST_SET_TEST tests HERMITE_PROBABILIST_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2013
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "HERMITE_PROBABILIST_SET_TEST\n";
  std::cout << "  HERMITE_PROBABILIST_SET sets a Hermite quadrature rule.\n";
  std::cout << "  The integration interval is ( -oo, +oo ).\n";
  std::cout << "  The weight function is exp ( - x * x / 2 ) / sqrt ( 2 * pi ).\n";
  std::cout << "\n";
  std::cout << "  Index             X                   W\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    hermite_probabilist_set ( n, x, w );

    std::cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void imtqlx_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    IMTQLX_TEST tests IMTQLX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 June 2015
//
//  Author:
//
//    John Burkardt.
//
{
  double angle;
  double d[5];
  double e[5];
  int i;
  double lam[5];
  double lam2[5];
  int n = 5;
  double qtz[5];
  double r8_pi = 3.141592653589793;
  double z[5];

  std::cout << "\n";
  std::cout << "IMTQLX_TEST\n";
  std::cout << "  IMTQLX takes a symmetric tridiagonal matrix A\n";
  std::cout << "  and computes its eigenvalues LAM.\n";
  std::cout << "  It also accepts a vector Z and computes Q'*Z,\n";
  std::cout << "  where Q is the matrix that diagonalizes A.\n";

  for ( i = 0; i < n; i++ )
  {
    d[i] = 2.0;
  }
  for ( i = 0; i < n - 1; i++ )
  {
    e[i] = -1.0;
  }
  e[n-1] = 0.0;
  for ( i = 0; i < n; i++ )
  {
    z[i] = 1.0;
  }
//
//  On input, LAM is D, and QTZ is Z.
//
  for ( i = 0; i < n; i++ )
  {
    lam[i] = d[i];
  }
  for ( i = 0; i < n; i++ )
  {
    qtz[i] = z[i];
  }

  imtqlx ( n, lam, e, qtz );

  r8vec_print ( n, lam, "  Computed eigenvalues:" );

  for ( i = 0; i < n; i++ )
  {
    angle = double( i + 1 ) * r8_pi / double( 2 * ( n + 1 ) );
    lam2[i] = 4.0 * pow ( sin ( angle ), 2 );
  }

  r8vec_print ( n, lam2, "  Exact eigenvalues:" );

  r8vec_print ( n, z, "  Vector Z:" );
  r8vec_print ( n, qtz, "  Vector Q'*Z:" );

  return;
}
//****************************************************************************80

void jacobi_ek_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_EK_COMPUTE_TEST tests JACOBI_EK_COMPUTE.
//
//  Discussion:
//
//    Compare with tabular values on page 178 of Stroud and Secrest.
//
//     In particular,
//
//             X              W
//
//     1  -0.9833999115   0.4615276287E-03
//     2  -0.9447138932   0.2732249104E-02
//     3  -0.8849310847   0.8045830455E-02
//    ..  .............   ................
//    19   0.9656375637   0.7613987785E-01
//    20   0.9934477866   0.3348337670E-01
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  double beta;
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  alpha = 1.5;
  beta = 0.5;

  std::cout << "\n";
  std::cout << "JACOBI_EK_COMPUTE_TEST\n";
  std::cout << "  JACOBI_EK_COMPUTE computes a Gauss-Jacobi rule;\n";
  std::cout << "\n";
  std::cout << "  ALPHA = " << alpha << "\n";
  std::cout << "  BETA =  " << beta << "\n";

  std::cout << "\n";
  std::cout << "     Order        W               X\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    jacobi_ek_compute ( n, alpha, beta, x, w );

    std::cout << "\n";
    std::cout << "  " << std::setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "          "
           << "  " << std::setw(14) << w[i]
           << "  " << std::setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }  

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void jacobi_integral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_INTEGRAL_TEST tests JACOBI_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  double beta;
  int n;
  int prec;
  double value;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  alpha = 1.5;
  beta = 0.5;

  std::cout << "\n";
  std::cout << "JACOBI_INTEGRAL_TEST\n";
  std::cout << "  JACOBI_INTEGRAL evaluates\n";
  std::cout << "  Integral ( -1 < x < +1 ) x^n (1-x)^alpha (1+x)^beta dx\n";
  std::cout << "\n";
  std::cout << "  Use ALPHA = " << alpha << "\n";
  std::cout << "      BETA  = " << beta << "\n";
  std::cout << "\n";
  std::cout << "         N         Value\n";
  std::cout << "\n";

  for ( n = 0; n <= 10; n++ )
  {
    value = jacobi_integral ( n, alpha, beta );
    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(24) << value << "\n";
  }
 
  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void jacobi_ss_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_SS_COMPUTE_TEST tests JACOBI_SS_COMPUTE.
//
//  Discussion:
//
//    Compare with tabular values on page 178 of Stroud and Secrest.
//
//     In particular,
//
//             X              W
//
//     1  -0.9833999115   0.4615276287E-03
//     2  -0.9447138932   0.2732249104E-02
//     3  -0.8849310847   0.8045830455E-02
//    ..  .............   ................
//    19   0.9656375637   0.7613987785E-01
//    20   0.9934477866   0.3348337670E-01
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  double beta;
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  alpha = 1.5;
  beta = 0.5;

  std::cout << "\n";
  std::cout << "JACOBI_SS_COMPUTE_TEST\n";
  std::cout << "  JACOBI_SS_COMPUTE computes a Gauss-Jacobi rule;\n";
  std::cout << "\n";
  std::cout << "  ALPHA = " << alpha << "\n";
  std::cout << "  BETA =  " << beta << "\n";

  std::cout << "\n";
  std::cout << "     Order        W               X\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    jacobi_ss_compute ( n, alpha, beta, x, w );

    std::cout << "\n";
    std::cout << "  " << std::setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "          "
           << "  " << std::setw(14) << w[i]
           << "  " << std::setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }  

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void kronrod_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    KRONROD_SET_TEST tests KRONROD_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int nk;
  int nl;
  int nl_test[4] = { 7, 10, 15, 20 };
  int test;
  double *wk;
  double *wl;
  double *xk;
  double *xl;

  std::cout << "\n";
  std::cout << "KRONROD_SET_TEST\n";
  std::cout << "  KRONROD_SET sets up a Kronrod quadrature rule;\n";
  std::cout << "  This is used following a lower order Legendre rule.\n";

  for ( test = 0; test < 4; test++ )
  {
    std::cout << "\n";
    std::cout << "  Legendre/Kronrod quadrature pair #" << test << "\n";
    std::cout << "                X                         W\n";
    std::cout << "\n";

    nl = nl_test[test];
    wl = new double[nl];
    xl = new double[nl];

    legendre_set ( nl, xl, wl );

    std::cout << "\n";

    for ( i = 0; i < nl; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << xl[i]
           << "  " << std::setw(24) << wl[i] << "\n";
    }
    delete [] wl;
    delete [] xl;

    std::cout << "\n";

    nk = 2 * nl + 1;
    wk = new double[nk];
    xk = new double[nk];

    kronrod_set ( nk, xk, wk );

    std::cout << "\n";

    for ( i = 0; i < nk; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << xk[i]
           << "  " << std::setw(24) << wk[i] << "\n";
    }
    delete [] wk;
    delete [] xk;
  }
  return;
}
//****************************************************************************80

void laguerre_ek_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_EK_COMPUTE_TEST tests LAGUERRE_EK_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "LAGUERRE_EK_COMPUTE_TEST\n";
  std::cout << "  LAGUERRE_EK_COMPUTE computes a Laguerre quadrature rule\n";
  std::cout << "  using the Elhay-Kautsky algorithm.\n";
  std::cout << "\n";
  std::cout << "     Order        W               X\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    laguerre_ek_compute ( n, x, w );

    std::cout << "\n";
    std::cout << "  " << std::setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "          "
           << "  " << std::setw(14) << w[i]
           << "  " << std::setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void laguerre_integral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_INTEGRAL_TEST tests LAGUERRE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int prec;
  double value;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "LAGUERRE_INTEGRAL_TEST\n";
  std::cout << "  LAGUERRE_INTEGRAL evaluates\n";
  std::cout << "  Integral ( 0 < x < oo ) x^n * exp(-x) dx\n";
  std::cout << "\n";
  std::cout << "         N         Value\n";
  std::cout << "\n";

  for ( n = 0; n <= 10; n++ )
  {
    value = laguerre_integral ( n );
    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(24) << value << "\n";
  }
 
  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void laguerre_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_SET_TEST tests LAGUERRE_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "LAGUERRE_SET_TEST\n";
  std::cout << "  LAGUERRE_SET sets a Laguerre rule.\n";
  std::cout << "\n";
  std::cout << "         I      X            W\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    laguerre_set ( n, x, w );

    std::cout << "\n";
    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(8)  << i
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void laguerre_ss_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_SS_COMPUTE_TEST tests LAGUERRE_SS_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "LAGUERRE_SS_COMPUTE_TEST\n";
  std::cout << "  LAGUERRE_SS_COMPUTE computes a Laguerre quadrature rule\n";
  std::cout << "  using the Stroud-Secrest algorithm.\n";
  std::cout << "\n";
  std::cout << "     Order        W               X\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    laguerre_ss_compute ( n, x, w );

    std::cout << "\n";
    std::cout << "  " << std::setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "          "
           << "  " << std::setw(14) << w[i]
           << "  " << std::setw(14) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void laguerre_1_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_1_SET_TEST tests LAGUERRE_1_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "LAGUERRE_1_SET_TEST\n";
  std::cout << "  LAGUERRE_1_SET sets a Laguerre rule.\n";
  std::cout << "  The density function is rho(x)=1.\n";
  std::cout << "\n";
  std::cout << "         I      X            W\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    laguerre_1_set ( n, x, w );

    std::cout << "\n";
    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(8)  << i
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void legendre_dr_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_DR_COMPUTE_TEST tests LEGENDRE_DR_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "LEGENDRE_DR_COMPUTE_TEST\n";
  std::cout << "  LEGENDRE_DR_COMPUTE computes a Legendre quadrature rule\n";
  std::cout << "  using the Davis-Rabinowitz algorithm.\n";
  std::cout << "\n";
  std::cout << "     Order        W               X\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

//  legendre_dr_compute ( n, x, w );

    std::cout << "\n";
    std::cout << "  " << std::setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "          "
           << "  " << std::setw(24) << w[i]
           << "  " << std::setw(24) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void legendre_ek_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_EK_COMPUTE_TEST tests LEGENDRE_EK_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "LEGENDRE_EK_COMPUTE_TEST\n";
  std::cout << "  LEGENDRE_EK_COMPUTE computes a Legendre quadrature rule\n";
  std::cout << "  using the Elhay-Kautsky algorithm.\n";
  std::cout << "\n";
  std::cout << "     Order        W               X\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    legendre_ek_compute ( n, x, w );

    std::cout << "\n";
    std::cout << "  " << std::setw(8) << n << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "          "
           << "  " << std::setw(24) << w[i]
           << "  " << std::setw(24) << x[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void legendre_integral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_INTEGRAL_TEST tests LEGENDRE_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int prec;
  double value;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "LEGENDRE_INTEGRAL_TEST\n";
  std::cout << "  LEGENDRE_INTEGRAL evaluates\n";
  std::cout << "  Integral ( -1 < x < +1 ) x^n dx\n";
  std::cout << "\n";
  std::cout << "         N         Value\n";
  std::cout << "\n";

  for ( n = 0; n <= 10; n++ )
  {
    value = legendre_integral ( n );
    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(24) << value << "\n";
  }
 
  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void legendre_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_SET_TEST tests LEGENDRE_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "LEGENDRE_SET_TEST\n";
  std::cout << "  LEGENDRE_SET sets a Legendre quadrature rule.\n";
  std::cout << "\n";
  std::cout << "         I      X            W\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    legendre_set ( n, x, w );

    std::cout << "\n";
    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(8)  << i
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void lobatto_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LOBATTO_COMPUTE_TEST tests LOBATTO_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "LOBATTO_COMPUTE_TEST\n";
  std::cout << "  LOBATTO_COMPUTE computes a Lobatto rule;\n";
  std::cout << "\n";
  std::cout << "         I      X             W\n";

  for ( n = 4; n <= 12; n = n + 3 )
  {
    w = new double[n];
    x = new double[n];

    lobatto_compute ( n, x, w );

    std::cout << "\n";
    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(8)  << i
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void lobatto_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LOBATTO_SET_TEST tests LOBATTO_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "LOBATTO_SET_TEST\n";
  std::cout << "  LOBATTO_SET sets a Lobatto rule;\n";
  std::cout << "\n";
  std::cout << "         I      X             W\n";

  for ( n = 4; n <= 12; n = n + 3 )
  {
    w = new double[n];
    x = new double[n];

    lobatto_set ( n, x, w );

    std::cout << "\n";
    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(8)  << i
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void nc_compute_weights_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NC_COMPUTE_WEIGHTS_TEST tests NCC_COMPUTE_WEIGHTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;
  double x_min;
  double x_max;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "NC_COMPUTE_WEIGHTS_TEST\n";
  std::cout << "  NC_COMPUTE_WEIGHTS computes weights for a Newton-Cotes rule;\n";

  x_min = 0.0;
  x_max = 1.0;
  
  std::cout << "\n";
  std::cout << "  Index             X                   W\n";

  for ( n = 1; n <= 10; n++ )
  {
    x = r8vec_linspace_new ( n, x_min, x_max );

    w = new double[n];
    nc_compute_weights ( n, x_min, x_max, x, w );

    std::cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }
 
  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void ncc_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NCC_COMPUTE_TEST tests NCC_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "NCC_COMPUTE_TEST\n";
  std::cout << "  NCC_COMPUTE computes a Newton-Cotes Closed rule;\n";
  std::cout << "\n";
  std::cout << "  Index             X                   W\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    ncc_compute ( n, x, w );

    std::cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void ncc_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NCC_SET_TEST tests NCC_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "NCC_SET_TEST\n";
  std::cout << "  NCC_SET sets up a Newton-Cotes Closed rule;\n";
  std::cout << "\n";
  std::cout << "  Index             X                   W\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    ncc_set ( n, x, w );

    std::cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void nco_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_COMPUTE_TEST tests NCO_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "NCO_COMPUTE_TEST\n";
  std::cout << "  NCO_COMPUTE computes a Newton-Cotes Open rule;\n";
  std::cout << "\n";
  std::cout << "  Index             X                   W\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    nco_compute ( n, x, w );

    std::cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void nco_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_SET_TEST tests NCO_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "NCO_SET_TEST\n";
  std::cout << "  NCO_SET sets up a Newton-Cotes Open rule;\n";
  std::cout << "\n";
  std::cout << "  Index             X                   W\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    nco_set ( n, x, w );

    std::cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void ncoh_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NCOH_COMPUTE_TEST tests NCOH_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "NCOH_COMPUTE_TEST\n";
  std::cout << "  NCOH_COMPUTE computes a Newton-Cotes Open Half rule;\n";
  std::cout << "\n";
  std::cout << "  Index             X                   W\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    ncoh_compute ( n, x, w );

    std::cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void ncoh_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NCOH_SET_TEST tests NCOH_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "NCOH_SET_TEST\n";
  std::cout << "  NCOH_SET sets up a Newton-Cotes Open Half rule;\n";
  std::cout << "\n";
  std::cout << "  Index             X                   W\n";

  for ( n = 1; n <= 10; n++ )
  {
    w = new double[n];
    x = new double[n];

    ncoh_set ( n, x, w );

    std::cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void patterson_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    PATTERSON_SET_TEST tests PATTERSON_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int n;
  int n_test[4] = { 1, 3, 7, 15 };
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "PATTERSON_SET_TEST\n";
  std::cout << "  PATTERSON_SET sets a Gauss-Patterson quadrature rule;\n";
  std::cout << "\n";
  std::cout << "  Index             X                   W\n";

  for ( j = 0; j < 4; j++ )
  {
    n = n_test[j];
    w = new double[n];
    x = new double[n];

    patterson_set ( n, x, w );

    std::cout << "\n";

    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(2) << i 
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }
 
  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void r8_psi_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_PSI_TEST tests R8_PSI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int n_data;
  int prec;
  double x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "R8_PSI_TEST:\n";
  std::cout << "  R8_PSI evaluates the Psi function.\n";
  std::cout << "\n";
  std::cout << "         X                  Psi(X)           " 
       << "         Psi(X)          DIFF\n";
  std::cout << "                         (Tabulated)         " 
       << "       (R8_PSI)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    psi_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r8_psi ( x );

    std::cout << "  " << std::setw(8)  << x   
         << "  " << std::setw(24) << fx  
         << "  " << std::setw(24) << fx2 
         << "  " << std::setw(10) << fabs ( fx - fx2 ) << "\n";

  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void radau_compute_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    RADAU_COMPUTE_TEST tests RADAU_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "RADAU_COMPUTE_TEST\n";
  std::cout << "  RADAU_COMPUTE computes a Radau rule;\n";
  std::cout << "\n";
  std::cout << "         I      X            W\n";

  for ( n = 4; n <= 12; n = n + 3 )
  {
    w = new double[n];
    x = new double[n];

    radau_compute ( n, x, w );

    std::cout << "\n";
    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(8)  << i
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
//****************************************************************************80

void radau_set_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    RADAU_SET_TEST tests RADAU_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 June 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prec;
  double *w;
  double *x;

  prec = std::cout.precision ( );
  std::cout.precision ( 16 );

  std::cout << "\n";
  std::cout << "RADAU_SET_TEST\n";
  std::cout << "  RADAU_SET sets a Radau rule from a table.\n";
  std::cout << "\n";
  std::cout << "         I      X            W\n";

  for ( n = 4; n <= 12; n = n + 3 )
  {
    w = new double[n];
    x = new double[n];

    radau_set ( n, x, w );

    std::cout << "\n";
    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(8)  << i
           << "  " << std::setw(24) << x[i]
           << "  " << std::setw(24) << w[i] << "\n";
    }
    delete [] w;
    delete [] x;
  }

  std::cout.precision ( prec );

  return;
}
