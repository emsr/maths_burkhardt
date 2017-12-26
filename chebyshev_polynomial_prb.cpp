# include <cmath>
# include <cstdlib>
# include <cstring>
# include <iomanip>
# include <iostream>

# include "chebyshev_polynomial.hpp"

int main ( );
void test01 ( );
void t_mass_matrix_test ( );
void t_moment_test ( );
void t_polynomial_test ( );
void t_polynomial_ab_test ( );
void t_polynomial_ab_value_test ( );
void t_polynomial_coefficients_test ( );
void t_polynomial_plot_test ( );
void t_polynomial_value_test ( );
void t_polynomial_zeros_test ( );
void t_quadrature_rule_test ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void tt_product_test ( );
void tt_product_integral_test ( );
void ttt_product_integral_test ( );
void tu_product_test ( );
void u_mass_matrix_test ( );
void u_moment_test ( );
void u_polynomial_test ( );
void u_polynomial_ab_test ( );
void u_polynomial_ab_value_test ( );
void u_polynomial_coefficients_test ( );
void u_polynomial_plot_test ( );
void u_polynomial_value_test ( );
void u_polynomial_zeros_test ( );
void u_quadrature_rule_test ( );
void uu_product_test ( );
void uu_product_integral_test ( );
void v_mass_matrix_test ( );
void v_moment_test ( );
void v_polynomial_test ( );
void v_polynomial_ab_test ( );
void v_polynomial_ab_value_test ( );
void v_polynomial_coefficients_test ( );
void v_polynomial_plot_test ( );
void v_polynomial_value_test ( );
void v_polynomial_zeros_test ( );
void v_quadrature_rule_test ( );
void vv_product_integral_test ( );
void w_mass_matrix_test ( );
void w_moment_test ( );
void w_polynomial_test ( );
void w_polynomial_ab_test ( );
void w_polynomial_ab_value_test ( );
void w_polynomial_coefficients_test ( );
void w_polynomial_plot_test ( );
void w_polynomial_value_test ( );
void w_polynomial_zeros_test ( );
void w_quadrature_rule_test ( );
void ww_product_integral_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CHEBYSHEV_POLYNOMIAL_PRB.
//
//  Discussion:
//
//    CHEBYSHEV_POLYNOMIAL_PRB tests the CHEBYSHEV_POLYNOMIAL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  std::cout << "\n";
  std::cout << "CHEBYSHEV_POLYNOMIAL_PRB\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the CHEBYSHEV_POLYNOMIAL library.\n";
 
  test01 ( );
  t_mass_matrix_test ( );
  t_moment_test ( );
  t_polynomial_test ( );
  t_polynomial_ab_test ( );
  t_polynomial_ab_value_test ( );
  t_polynomial_coefficients_test ( );
  t_polynomial_plot_test ( );
  t_polynomial_value_test ( );
  t_polynomial_zeros_test ( );
  t_quadrature_rule_test ( );
  test07 ( );
  test08 ( );
  test09 ( );
  test10 ( );
  tt_product_test ( );
  tt_product_integral_test ( );
  ttt_product_integral_test ( );
  tu_product_test ( );

  u_mass_matrix_test ( );
  u_moment_test ( );
  u_polynomial_test ( );
  u_polynomial_ab_test ( );
  u_polynomial_ab_value_test ( );
  u_polynomial_coefficients_test ( );
  u_polynomial_plot_test ( );
  u_polynomial_value_test ( );
  u_polynomial_zeros_test ( );
  u_quadrature_rule_test ( );
  uu_product_test ( );
  uu_product_integral_test ( );

  v_mass_matrix_test ( );
  v_moment_test ( );
  v_polynomial_test ( );
  v_polynomial_ab_test ( );
  v_polynomial_ab_value_test ( );
  v_polynomial_coefficients_test ( );
  v_polynomial_plot_test ( );
  v_polynomial_value_test ( );
  v_polynomial_zeros_test ( );
  v_quadrature_rule_test ( );
  vv_product_integral_test ( );

  w_mass_matrix_test ( );
  w_moment_test ( );
  w_polynomial_test ( );
  w_polynomial_ab_test ( );
  w_polynomial_ab_value_test ( );
  w_polynomial_coefficients_test ( );
  w_polynomial_plot_test ( );
  w_polynomial_value_test ( );
  w_polynomial_zeros_test ( );
  w_quadrature_rule_test ( );
  ww_product_integral_test ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "CHEBYSHEV_POLYNOMIAL_PRB\n";
  std::cout << "  Normal end of execution.\n";
  std::cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests T_PROJECT_COEFFICIENTS_DATA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2016
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double *c;
  double *d;
  double *d2;
  int i;
  int m;
  int n;
  int seed;
  double *x;

  std::cout << "\n";
  std::cout << "CHEBYSHEV_POLYNOMIAL_TEST01:\n";
  std::cout << "  T_PROJECT_COEFFICIENTS_DATA estimates the Chebyshev polynomial\n";
  std::cout << "  coefficients for a function given as data (x,fx).\n";
  std::cout << "\n";
  std::cout << "  Here, we use fx = f(x) = x^2 for the data.\n";
  std::cout << "\n";
  std::cout << "  Since T(0,x) = 1 and T(2,x) = 2*x^2 - 1, the correct expansion is\n";
  std::cout << "  f(x) = 1/2 T(0,x) + 0 T(1,x) + 1/2 T(2,x) + 0 * all other polys,\n";
  std::cout << "  if Chebyshev polynomials are based in [-1,+1].\n";
//
//  Data in [0,1];
//
  a = 0.0;
  b = 1.0;

  std::cout << "\n";
  std::cout << "  Chebyshev polynomial will be based in [" << a << "," << b << "]\n";
//
//  Compute sample data.
//
  m = 20;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( m, a, b, seed );
  d = new double[m];
  for ( i = 0; i < m; i++ )
  {
    d[i] = x[i] * x[i];
  }

  r8vec2_print ( m, x, d, "  Data ( X, D ):" );

  n = 4;
  c = t_project_coefficients_data ( a, b, m, n, x, d );
  
  r8vec_print ( n, c, "  Coefficients of Chebyshev expansion of degree 4." );
//
//  Compare Chebyshev expansion and original function.
//
  d2 = t_project_value_ab ( m, n, x, c, a, b );

  std::cout << "\n";
  std::cout << "   I      X(I)     Data(I)      Chebyshev(X(I))\n";
  std::cout << "\n";
  for ( i = 0; i < m; i++ )
  {
    std::cout << "  " << std::setw(2) << i
         << "  " << std::setw(12) << x[i]
         << "  " << std::setw(12) << d[i]
         << "  " << std::setw(12) << d2[i] << "\n";
  }
//
//  Free memory.
//
  delete [] c;
  delete [] d;
  delete [] d2;
  delete [] x;

  return;
}
//****************************************************************************80

void t_mass_matrix_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    T_MASS_MATRIX_TEST tests T_MASS_MATRIX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int n;

  std::cout << "\n";
  std::cout <<  "T_MASS_MATRIX_TEST:\n";
  std::cout <<  "  T_MASS_MATRIX computes the mass matrix for the\n";
  std::cout <<  "  Chebyshev polynomials T(i,x).\n";
  std::cout <<  "  A(I,J) = integral ( -1 <=x <= +1 ) T(i,x) T(j,x) / sqrt ( 1 - x^2 ) dx\n";
  std::cout <<  "  0    if i is not equal to j;\n";
  std::cout <<  "  pi   if i = j = 0;\n";
  std::cout <<  "  pi/2 if i = j =/= 0.\n";

  n = 3;
  a = t_mass_matrix ( n );

  r8mat_print ( n + 1, n + 1, a, "  T mass matrix:" );

  delete [] a;

  return;
}
//****************************************************************************80

void t_moment_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    T_MOMENT_TEST tests T_MOMENT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  int e;
  double value;

  std::cout << "\n";
  std::cout << "T_MOMENT_TEST:\n";
  std::cout << "  T_MOMENT returns the value of\n";
  std::cout << "  integral ( -1 <=x <= +1 ) x^e / sqrt ( 1 - x^2 ) dx\n";
  std::cout << "\n";
  std::cout << "   E       Integral\n";
  std::cout << "\n";
  for ( e = 0; e <= 10; e++ )
  {
    value = t_moment ( e );
    std::cout << "  " << std::setw(2) << e
         << "  " << std::setw(14) << value << "\n";
  }

  return;
}
//****************************************************************************80

void t_polynomial_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    T_POLYNOMIAL_TEST tests T_POLYNOMIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double *fx2;
  int n;
  int n_data;
  double x;
  double x_vec[1];

  std::cout << "\n";
  std::cout << "T_POLYNOMIAL_TEST:\n";
  std::cout << "  T_POLYNOMIAL evaluates the Chebyshev polynomial T(n,x).\n";
  std::cout << "\n";
  std::cout << "                   Tabulated      Computed\n";
  std::cout << "     N      X        T(n,x)        T(n,x)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    t_polynomial_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    if ( n < 0 )
    {
      continue;
    }

    x_vec[0] = x;
    fx2 = t_polynomial ( 1, n, x_vec );

    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(8)  << x
         << "  " << std::setw(14) << fx
         << "  " << std::setw(14) << fx2[n] << "\n";

    delete [] fx2;

  }

  return;
}
//****************************************************************************80

void t_polynomial_ab_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    T_POLYNOMIAL_AB_TEST tests T_POLYNOMIAL_AB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int m;
  int n;
  double *v;
  double *x;

  m = 11;
  n = 5;

  std::cout << "\n";
  std::cout << "T_POLYNOMIAL_AB_TEST:\n";
  std::cout << "  T_POLYNOMIAL_AB evaluates Chebyshev polynomials TAB(n,x)\n";
  std::cout << "  shifted from [-1,+1] to the domain [A,B].\n";
  std::cout << "\n";
  std::cout << "  Here, we will use the new domain [0,1]\n";
  std::cout << "  and the desired maximum polynomial degree will be N = 5.\n";

  a = 0.0;
  b = 1.0;
  x = r8vec_linspace_new ( m, a, b );
  
  v = t_polynomial_ab ( a, b, m, n, x );

  r8mat_print ( m, n + 1, v, "  Tables of T values:" );

  delete [] v;
  delete [] x;

  return;
}
//****************************************************************************80

void t_polynomial_ab_value_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    T_POLYNOMIAL_AB_VALUE_TEST tests T_POLYNOMIAL_AB_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  double fx2;
  int n;
  int n_data;
  double x01;

  std::cout << "\n";
  std::cout << "T_POLYNOMIAL_AB_VALUE_TEST:\n";
  std::cout << "  T_POLYNOMIAL_AB_VALUE evaluates Chebyshev polynomials TAB(n,x)\n";
  std::cout << "  shifted from [-1,+1] to the domain [A,B].\n";
  std::cout << "\n";
  std::cout << "  Here, we will use the new domain [0,1].\n";
  std::cout << "\n";
  std::cout << "                   Tabulated      Computed\n";
  std::cout << "     N      X01    T01(n,x)       T01(n,x)\n";
  std::cout << "\n";

  a = 0.0;
  b = 1.0;

  n_data = 0;

  for ( ; ; )
  {
    t_polynomial_01_values ( n_data, n, x01, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = t_polynomial_ab_value ( a, b, n, x01 );

    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(8)  << x01
         << "  " << std::setw(14) << fx
         << "  " << std::setw(14) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void t_polynomial_coefficients_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    T_POLYNOMIAL_COEFFICIENTS_TEST tests T_POLYNOMIAL_COEFFICIENTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  double *c2;
  int i;
  int j;
  int n;

  std::cout << "\n";
  std::cout << "T_POLYNOMIAL_COEFFICIENTS_TEST\n";
  std::cout << "  T_POLYNOMIAL_COEFFICIENTS determines the polynomial coefficients \n";
  std::cout << "  of T(n,x).\n";

  n = 5;

  c = t_polynomial_coefficients ( n );
 
  for ( i = 0; i <= n; i++ )
  {
    c2 = new double[i+1];
    for ( j = 0; j <= i; j++ )
    {
      c2[j] = c[i+j*(n+1)];
    }
    r8poly_print ( i, c2, "" );
    delete [] c2;
  }
  delete [] c;

  return;
}
//****************************************************************************80

void t_polynomial_plot_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    T_POLYNOMIAL_PLOT_TEST tests T_POLYNOMIAL_PLOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n_num = 6;
  int n_val[6];
  std::string output_filename;

  std::cout << "\n";
  std::cout << "T_POLYNOMIAL_PLOT_TEST\n";
  std::cout << "  T_POLYNOMIAL_PLOT plots selected\n";
  std::cout << "  Chebyshev polynomials T(n,x).\n";

  for ( i = 0; i <= 5; i++ )
  {
    n_val[i] = i;
  }

  output_filename = "t_polynomial_plot.png";

  t_polynomial_plot ( n_num, n_val, output_filename );

  return;
}
//****************************************************************************80

void t_polynomial_value_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    T_POLYNOMIAL_VALUE_TEST tests T_POLYNOMIAL_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "T_POLYNOMIAL_VALUE_TEST:\n";
  std::cout << "  T_POLYNOMIAL_VALUE evaluates the Chebyshev polynomial T(n,x).\n";
  std::cout << "\n";
  std::cout << "                   Tabulated      Computed\n";
  std::cout << "     N      X        T(n,x)        T(n,x)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    t_polynomial_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = t_polynomial_value ( n, x );

    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(8)  << x
         << "  " << std::setw(14) << fx
         << "  " << std::setw(14) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void t_polynomial_zeros_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    T_POLYNOMIAL_ZEROS_TEST tests T_POLYNOMIAL_ZEROS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *fx;
  int i;
  int n;
  int n_max = 5;
  double *z;

  std::cout << "\n";
  std::cout << "T_POLYNOMIAL_ZEROS_TEST:\n";
  std::cout << "  T_POLYNOMIAL_ZEROS returns zeroes of T(n,x).\n";
  std::cout << "\n";
  std::cout << "       N      X        T(n,x)\n";
  std::cout << "\n";

  for ( n = 1; n <= n_max; n++ )
  {
    z = t_polynomial_zeros ( n );
    fx = t_polynomial ( n, n, z );
    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(8) << n
           << "  " << std::setw(8) << z[i]
           << "  " << std::setw(14) << fx[i+n*n] << "\n";
    }
    std::cout << "\n";
    delete [] fx;
    delete [] z;
  }

  return;
}
//****************************************************************************80

void t_quadrature_rule_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    T_QUADRATURE_RULE_TEST tests T_QUADRATURE_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  int e;
  double *f;
  int i;
  int n;
  double q;
  double q_exact;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "T_QUADRATURE_RULE_TEST:\n";
  std::cout << "  T_QUADRATURE_RULE computes the quadrature rule\n";
  std::cout << "  associated with T(n,x);\n";

  n = 7;
  x = new double[n];
  w = new double[n];

  t_quadrature_rule ( n, x, w );

  r8vec2_print ( n, x, w, "    N      X            W" );

  std::cout << "\n";
  std::cout << "  Use the quadrature rule to estimate:\n";
  std::cout << "\n";
  std::cout << "    Q = Integral ( -1 <= X <= +1 ) X^E / sqrt ( 1-x^2) dx\n";
  std::cout << "\n";
  std::cout << "   E       Q_Estimate      Q_Exact\n";
  std::cout << "\n";

  f = new double[n];

  for ( e = 0; e <= 2 * n - 1; e++ )
  {
    if ( e == 0 )
    {
      for ( i = 0; i < n; i++ )
      {
        f[i] = 1.0;
      }
    }
    else
    {
      for ( i = 0; i < n; i++ )
      {
        f[i] = pow ( x[i], e );
      }
    }
    q = r8vec_dot_product ( n, w, f );
    q_exact = t_moment ( e );
    std::cout << "  " << std::setw(2) << e
         << "  " << std::setw(14) << q
         << "  " << std::setw(14) << q_exact << "\n";
  }

  delete [] f;
  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests T_PROJECT_COEFFICIENTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double *c;
  int n;

  std::cout << "\n";
  std::cout << "TEST07:\n";
  std::cout << "  T_PROJECT_COEFFICIENTS computes the Chebyshev coefficients\n";
  std::cout << "  of a function defined over [-1,+1].\n";
  std::cout << "  T_PROJECT_COEFFICIENTS_AB works in [A,B].\n";

  n = 3;
  c = new double[n+1];
  c = t_project_coefficients ( n, exp );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for exp(x) in [-1,+1]" );
  delete [] c;

  n = 5;
  c = new double[n+1];
  c = t_project_coefficients ( n, exp );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for exp(x) in [-1,+1]" );
  delete [] c;

  n = 5;
  c = new double[n+1];
  c = t_project_coefficients ( n, sin );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for sin(x) in [-1,+1]" );
  delete [] c;
//
//  Repeat calculation with T_PROJECT_COEFFICIENTS_AB.
//
  n = 5;
  c = new double[n+1];
  a = -1.0;
  b = +1.0;
  c = t_project_coefficients_ab ( n, sin, a, b );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for sin(x) in [-1,+1]" );
  delete [] c;
//
//  Now try a different interval.
//
  n = 5;
  c = new double[n+1];
  a = 0.0;
  b = 1.0;
  c = t_project_coefficients_ab ( n, sqrt, a, b );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for sqrt(x) in [0,+1]" );
  delete [] c;

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests T_PROJECT_COEFFICIENTS_DATA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double *c;
  double *d;
  int i;
  int m;
  int n;
  int seed;
  double *x;

  std::cout << "\n";
  std::cout << "TEST08:\n";
  std::cout << "  T_PROJECT_COEFFICIENTS_DATA computes the Chebyshev\n";
  std::cout << "  coefficients of a function defined by data.\n";
  std::cout << "\n";
  std::cout << "  We are looking for an approximation that is good in [-1,+1].\n";
  std::cout << "\n";
  std::cout << "  Begin by using equally spaced points in [-1,+1].\n";

  a = -1.0;
  b = +1.0;
  m = 10;
  x = r8vec_linspace_new ( m, a, b );
  d = new double[m];
  for ( i = 0; i < m; i++ )
  {
    d[i] = exp ( x[i] );
  }
  n = 3;
  c = t_project_coefficients_data ( a, b, m, n, x, d );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for exp(x) on [-1,+1]" );
  delete [] c;
  delete [] d;
  delete [] x;

  a = -1.0;
  b = +1.0;
  m = 10;
  x = r8vec_linspace_new ( m, a, b );
  d = new double[m];
  for ( i = 0; i < m; i++ )
  {
    d[i] = exp ( x[i] );
  }
  n = 5;
  c = t_project_coefficients_data ( a, b, m, n, x, d );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for exp(x) on [-1,+1]" );
  delete [] c;
  delete [] d;
  delete [] x;

  a = -1.0;
  b = +1.0;
  m = 10;
  x = r8vec_linspace_new ( m, a, b );
  d = new double[m];
  for ( i = 0; i < m; i++ )
  {
    d[i] = sin ( x[i] );
  }
  n = 5;
  c = t_project_coefficients_data ( a, b, m, n, x, d );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for sin(x) on [-1,+1]" );
  delete [] c;
  delete [] d;
  delete [] x;

  std::cout << "\n";
  std::cout << "  Now sample equally spaced points in [0,+1].\n";
  std::cout << "  The approximation still applies to the interval [-1,+1].\n";

  a = 0.0;
  b = +1.0;
  m = 10;
  x = r8vec_linspace_new ( m, a, b );
  d = new double[m];
  for ( i = 0; i < m; i++ )
  {
    d[i] = sin ( x[i] );
  }
  n = 5;
  c = t_project_coefficients_data ( a, b, m, n, x, d );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for sin(x) on [0,+1]" );
  delete [] c;
  delete [] d;
  delete [] x;

  a = 0.0;
  b = +1.0;
  m = 10;
  x = r8vec_linspace_new ( m, a, b );
  d = new double[m];
  for ( i = 0; i < m; i++ )
  {
    d[i] = sqrt ( x[i] );
  }
  n = 5;
  c = t_project_coefficients_data ( a, b, m, n, x, d );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for sqrt(x) on [0,+1]" );
  delete [] c;
  delete [] d;
  delete [] x;

  std::cout << "\n";
  std::cout << "  Now random points in [-1,+1].\n";

  a = -1.0;
  b = +1.0;
  m = 10;
  seed = 123456789;
  x = r8vec_uniform_ab_new ( m, a, b, seed );
  d = new double[m];
  for ( i = 0; i < m; i++ )
  {
    d[i] = sin ( x[i] );
  }
  n = 5;
  c = t_project_coefficients_data ( a, b, m, n, x, d );
  r8vec_print ( n + 1, c, "  Chebyshev coefficients for sin(x) on [-1,+1]" );
  delete [] c;
  delete [] d;
  delete [] x;

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 compares a function and projection over [-1,+1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double *c;
  int i;
  int m;
  int n;
  double r;
  double *v;
  double *x;

  std::cout << "\n";
  std::cout << "TEST09:\n";
  std::cout << "  T_PROJECT_COEFFICIENTS computes the Chebyshev interpolant C(F)(n,x)\n";
  std::cout << "  of a function F(x) defined over [-1,+1].\n";
  std::cout << "  T_PROJECT_VALUE evaluates that projection.\n";

  std::cout << "\n";
  std::cout << "  Compute projections of order N to exp(x) over [-1,+1],\n";
  std::cout << "\n";
  std::cout << "   N   Max||F(x)-C(F)(n,x)||\n";
  std::cout << "\n";

  a = -1.0;
  b = +1.0;

  for ( n = 0; n <= 10; n++ )
  {
    c = t_project_coefficients ( n, exp );
    m = 101;
    x = r8vec_linspace_new ( m, a, b );
    v = t_project_value ( m, n, x, c );
    r = 0.0;
    for ( i = 0; i < m; i++ )
    {
      r = r8_max ( r, fabs ( v[i] - exp ( x[i] ) ) );
    }
    std::cout << "  " << std::setw(2) << n
         << "  " << std::setw(14) << r << "\n";
    delete [] c;
    delete [] v;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 compares a function and projection over [A,B].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double *c;
  int i;
  int m;
  int n;
  double r;
  double *v;
  double *x;

  std::cout << "\n";
  std::cout << "TEST10:\n";
  std::cout << "  T_PROJECT_COEFFICIENTS_AB computes the Chebyshev interpolant C(F)(n,x)\n";
  std::cout << "  of a function F(x) defined over [A,B].\n";
  std::cout << "  T_PROJECT_VALUE_AB evaluates that projection.\n";

  a = 0.0;
  b = 1.5;

  std::cout << "\n";
  std::cout << "  Compute projections of order N to exp(x) over [" << a << "," << b << "]\n";
  std::cout << "\n";
  std::cout << "   N   Max||F(x)-C(F)(n,x)||\n";
  std::cout << "\n";

  for ( n = 0; n <= 10; n++ )
  {
    c = t_project_coefficients_ab ( n, exp, a, b );
    m = 101;
    x = r8vec_linspace_new ( m, a, b );
    v = t_project_value_ab ( m, n, x, c, a, b );
    r = 0.0;
    for ( i = 0; i < m; i++ )
    {
      r = r8_max ( r, fabs ( v[i] - exp ( x[i] ) ) );
    }
    std::cout << "  " << std::setw(2) << n
         << "  " << std::setw(14) << r << "\n";
    delete [] c;
    delete [] v;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void tt_product_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TT_PRODUCT_TEST tests TT_PRODUCT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  double r8_hi;
  double r8_lo;
  int seed;
  int test;
  double ti;
  double titj;
  double tj;
  double x;

  std::cout << "\n";
  std::cout << "TT_PRODUCT_TEST:\n";
  std::cout << "  TT_PRODUCT(I,J;X) = T(I,X) * T(J,X)\n";

  r8_lo = -1.0;
  r8_hi = +1.0;
  seed = 123456789;

  std::cout << "\n";
  std::cout << "   I   J      X               TI              TJ              TI*TJ       TT_PRODUCT\n";
  std::cout << "\n";
  for ( test = 1; test <= 10; test++ )
  {
    x = r8_uniform_ab ( r8_lo, r8_hi, seed );
    i = i4_uniform_ab ( 0, 6, seed );
    ti = t_polynomial_value ( i, x );
    j = i4_uniform_ab ( -1, 4, seed );
    tj = t_polynomial_value ( j, x );
    titj = tt_product ( i, j, x );
    std::cout << "  " << std::setw(2) << i
         << "  " << std::setw(2) << j
         << "  " << std::setw(14) << x
         << "  " << std::setw(14) << ti
         << "  " << std::setw(14) << tj
         << "  " << std::setw(14) << ti * tj
         << "  " << std::setw(14) << titj << "\n";
  }

  return;
}
//****************************************************************************80

void tt_product_integral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TT_PRODUCT_INTEGRAL_TEST tests TT_PRODUCT_INTEGRAL.
//
//  Discussion:
//
//    This process should match the T_MASS_MATRIX computation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int i;
  int j;
  int n;

  std::cout << "\n";
  std::cout << "TT_PRODUCT_INTEGRAL_TEST:\n";
  std::cout << "  TT_PRODUCT_INTEGRAL computes the product integral\n";
  std::cout << "  of a pair of Chebyshev T polynomials T(i,x) and T(j,x).\n";
  std::cout << "  A(I,J) = integral ( -1 <=x <= +1 ) T(i,x) T(j,x) / sqrt ( 1 - x^2 ) dx\n";
  std::cout << "  0    if i is not equal to j;\n";
  std::cout << "  pi   if i = j = 0;\n";
  std::cout << "  pi/2 if i = j =/= 0.\n";

  n = 4;
  a = new double[(n+1)*(n+1)];
  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      a[i+j*(n+1)] = tt_product_integral ( i, j );
    }
  }

  r8mat_print ( n + 1, n + 1, a, "  T(i,x)*T(j,x) integral matrix:" );

  delete [] a;

  return;
}
//****************************************************************************80

void ttt_product_integral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TTT_PRODUCT_INTEGRAL_TEST tests TTT_PRODUCT_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double fx1;
  double fx2;
  int i;
  int j;
  int k;
  int l;
  int n;
  int seed;
  int test;
  int test_num = 20;
  double ti;
  double tj;
  double tk;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TTT_PRODUCT_INTEGRAL_TEST:\n";
  std::cout << "  TTT_PRODUCT_INTEGRAL computes the triple integral\n";
  std::cout << "  Tijk = integral ( -1 <= x <= 1 ) T(i,x) T(j,x) T(k,x) / sqrt ( 1-x^2) dx\n";
  std::cout << "\n";
  std::cout << "   I   J   K     Tijk           Tijk\n";
  std::cout << "                 computed       exact\n";
  std::cout << "\n";

  n = 15;
  x = new double[n];
  w = new double[n];

  t_quadrature_rule ( n, x, w );

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    i = i4_uniform_ab ( 2, 6, seed );
    j = i4_uniform_ab ( 1, 3, seed );
    k = i4_uniform_ab ( 0, 4, seed );
    fx1 = ttt_product_integral ( i, j, k );
    fx2 = 0.0;
    for ( l = 0; l < n; l++ )
    {
      ti = t_polynomial_value ( i, x[l] );
      tj = t_polynomial_value ( j, x[l] );
      tk = t_polynomial_value ( k, x[l] );
      fx2 = fx2 + w[l] * ti * tj * tk;
    }
    std::cout << "  " << std::setw(2) << i
         << "  " << std::setw(2) << j
         << "  " << std::setw(2) << k
         << "  " << std::setw(14) << fx1
         << "  " << std::setw(14) << fx2 << "\n";
  }

  delete [] x;
  delete [] w;

  return;
}
//****************************************************************************80

void tu_product_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TU_PRODUCT_TEST tests TU_PRODUCT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  double r8_hi;
  double r8_lo;
  int seed;
  int test;
  double ti;
  double tiuj;
  double uj;
  double x;

  std::cout << "\n";
  std::cout << "TU_PRODUCT_TEST:\n";
  std::cout << "  TU_PRODUCT(I,J;X) = T(I,X) * U(J,X)\n";

  r8_lo = -1.0;
  r8_hi = +1.0;
  seed = 123456789;

  std::cout << "\n";
  std::cout << "   I   J      X               TI              UJ              TI*UJ       TU_PRODUCT\n";
  std::cout << "\n";
  for ( test = 1; test <= 10; test++ )
  {
    x = r8_uniform_ab ( r8_lo, r8_hi, seed );
    i = i4_uniform_ab ( 0, 6, seed );
    ti = t_polynomial_value ( i, x );
    j = i4_uniform_ab ( -1, 4, seed );
    uj = u_polynomial_value ( j, x );
    tiuj = tu_product ( i, j, x );
    std::cout << "  " << std::setw(2) << i
         << "  " << std::setw(2) << j
         << "  " << std::setw(14) << x
         << "  " << std::setw(14) << ti
         << "  " << std::setw(14) << uj
         << "  " << std::setw(14) << ti * uj
         << "  " << std::setw(14) << tiuj << "\n";
  }

  return;
}
//****************************************************************************80

void u_mass_matrix_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    U_MASS_MATRIX_TEST tests U_MASS_MATRIX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int n;

  std::cout << "\n";
  std::cout <<  "U_MASS_MATRIX_TEST:\n";
  std::cout <<  "  U_MASS_MATRIX computes the mass matrix for the\n";
  std::cout <<  "  Chebyshev U polynomials U(i,x).\n";
  std::cout <<  "  A(I,J) = integral ( -1 <=x <= +1 ) U(i,x) U(j,x) * sqrt ( 1 - x^2 ) dx\n";
  std::cout <<  "  0    if i is not equal to j;\n";
  std::cout <<  "  pi/2 if i = j.\n";

  n = 3;
  a = u_mass_matrix ( n );

  r8mat_print ( n + 1, n + 1, a, "  U mass matrix:" );

  delete [] a;

  return;
}
//****************************************************************************80

void u_moment_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    U_MOMENT_TEST tests U_MOMENT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  int e;
  double value;

  std::cout << "\n";
  std::cout << "U_MOMENT_TEST:\n";
  std::cout << "  U_MOMENT returns the value of\n";
  std::cout << "  integral ( -1 <=x <= +1 ) x^e * sqrt ( 1 - x^2 ) dx\n";
  std::cout << "\n";
  std::cout << "   E       Integral\n";
  std::cout << "\n";
  for ( e = 0; e <= 10; e++ )
  {
    value = u_moment ( e );
    std::cout << "  " << std::setw(2) << e
         << "  " << std::setw(14) << value << "\n";
  }

  return;
}
//****************************************************************************80

void u_polynomial_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    U_POLYNOMIAL_TEST tests U_POLYNOMIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double *fx2;
  int n;
  int n_data;
  double x;
  double x_vec[1];

  std::cout << "\n";
  std::cout << "U_POLYNOMIAL_TEST:\n";
  std::cout << "  U_POLYNOMIAL evaluates the Chebyshev polynomial U(n,x).\n";
  std::cout << "\n";
  std::cout << "                   Tabulated      Computed\n";
  std::cout << "     N      X        U(n,x)        U(n,x)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    u_polynomial_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    if ( n < 0 )
    {
      continue;
    }

    x_vec[0] = x;
    fx2 = u_polynomial ( 1, n, x_vec );

    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(8)  << x
         << "  " << std::setw(14) << fx
         << "  " << std::setw(14) << fx2[n] << "\n";

    delete [] fx2;

  }

  return;
}
//****************************************************************************80

void u_polynomial_ab_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    U_POLYNOMIAL_AB_TEST tests U_POLYNOMIAL_AB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int m;
  int n;
  double *v;
  double *x;

  m = 11;
  n = 5;

  std::cout << "\n";
  std::cout << "U_POLYNOMIAL_AB_TEST:\n";
  std::cout << "  U_POLYNOMIAL_AB evaluates Chebyshev polynomials UAB(n,x)\n";
  std::cout << "  shifted from [-1,+1] to the domain [A,B].\n";
  std::cout << "\n";
  std::cout << "  Here, we will use the new domain [0,1]\n";
  std::cout << "  and the desired maximum polynomial degree will be N = 5.\n";

  a = 0.0;
  b = 1.0;
  x = r8vec_linspace_new ( m, a, b );
  
  v = u_polynomial_ab ( a, b, m, n, x );

  r8mat_print ( m, n + 1, v, "  Tables of U values:" );

  delete [] v;
  delete [] x;

  return;
}
//****************************************************************************80

void u_polynomial_ab_value_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    U_POLYNOMIAL_AB_VALUE_TEST tests U_POLYNOMIAL_AB_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  double fx2;
  int n;
  int n_data;
  double x01;

  std::cout << "\n";
  std::cout << "U_POLYNOMIAL_AB_VALUE_TEST:\n";
  std::cout << "  U_POLYNOMIAL_AB_VALUE evaluates Chebyshev polynomials UAB(n,x)\n";
  std::cout << "  shifted from [-1,+1] to the domain [A,B].\n";
  std::cout << "\n";
  std::cout << "  Here, we will use the new domain [0,1].\n";
  std::cout << "\n";
  std::cout << "                   Tabulated      Computed\n";
  std::cout << "     N      X01    U01(n,x)       U01(n,x)\n";
  std::cout << "\n";

  a = 0.0;
  b = 1.0;

  n_data = 0;

  for ( ; ; )
  {
    u_polynomial_01_values ( n_data, n, x01, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = u_polynomial_ab_value ( a, b, n, x01 );

    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(8)  << x01
         << "  " << std::setw(14) << fx
         << "  " << std::setw(14) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void u_polynomial_coefficients_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    U_POLYNOMIAL_COEFFICIENTS_TEST tests U_POLYNOMIAL_COEFFICIENTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  double *c2;
  int i;
  int j;
  int n;

  std::cout << "\n";
  std::cout << "U_POLYNOMIAL_COEFFICIENTS_TEST\n";
  std::cout << "  U_POLYNOMIAL_COEFFICIENTS determines the polynomial coefficients \n";
  std::cout << "  of U(n,x).\n";

  n = 5;

  c = u_polynomial_coefficients ( n );
 
  for ( i = 0; i <= n; i++ )
  {
    c2 = new double[i+1];
    for ( j = 0; j <= i; j++ )
    {
      c2[j] = c[i+j*(n+1)];
    }
    r8poly_print ( i, c2, "" );
    delete [] c2;
  }
  delete [] c;

  return;
}
//****************************************************************************80

void u_polynomial_plot_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    U_POLYNOMIAL_PLOT_TEST tests U_POLYNOMIAL_PLOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n_num = 6;
  int n_val[6];
  std::string output_filename;

  std::cout << "\n";
  std::cout << "U_POLYNOMIAL_PLOT_TEST\n";
  std::cout << "  U_POLYNOMIAL_PLOT plots selected\n";
  std::cout << "  Chebyshev polynomials U(n,x).\n";

  for ( i = 0; i <= 5; i++ )
  {
    n_val[i] = i;
  }

  output_filename = "u_polynomial_plot.png";

  u_polynomial_plot ( n_num, n_val, output_filename );

  return;
}
//****************************************************************************80

void u_polynomial_value_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    U_POLYNOMIAL_VALUE_TEST tests U_POLYNOMIAL_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "U_POLYNOMIAL_VALUE_TEST:\n";
  std::cout << "  U_POLYNOMIAL_VALUE evaluates the Chebyshev polynomial U(n,x).\n";
  std::cout << "\n";
  std::cout << "                   Tabulated      Computed\n";
  std::cout << "     N      X        U(n,x)        U(n,x)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    u_polynomial_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = u_polynomial_value ( n, x );

    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(8)  << x
         << "  " << std::setw(14) << fx
         << "  " << std::setw(14) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void u_polynomial_zeros_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    U_POLYNOMIAL_ZEROS_TEST tests U_POLYNOMIAL_ZEROS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *fx;
  int i;
  int n;
  int n_max = 5;
  double *z;

  std::cout << "\n";
  std::cout << "U_POLYNOMIAL_ZEROS_TEST:\n";
  std::cout << "  U_POLYNOMIAL_ZEROS returns zeroes of U(n,x).\n";
  std::cout << "\n";
  std::cout << "       N      X        U(n,x)\n";
  std::cout << "\n";

  for ( n = 1; n <= n_max; n++ )
  {
    z = u_polynomial_zeros ( n );
    fx = u_polynomial ( n, n, z );
    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(8) << n
           << "  " << std::setw(8) << z[i]
           << "  " << std::setw(14) << fx[i+n*n] << "\n";
    }
    std::cout << "\n";
    delete [] fx;
    delete [] z;
  }

  return;
}
//****************************************************************************80

void u_quadrature_rule_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    U_QUADRATURE_RULE_TEST tests U_QUADRATURE_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  int e;
  double *f;
  int i;
  int n;
  double q;
  double q_exact;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "U_QUADRATURE_RULE_TEST:\n";
  std::cout << "  U_QUADRATURE_RULE computes the quadrature rule\n";
  std::cout << "  associated with U(n,x);\n";

  n = 7;
  x = new double[n];
  w = new double[n];

  u_quadrature_rule ( n, x, w );

  r8vec2_print ( n, x, w, "    N      X            W" );

  std::cout << "\n";
  std::cout << "  Use the quadrature rule to estimate:\n";
  std::cout << "\n";
  std::cout << "    Q = Integral ( -1 <= X <= +1 ) X^E * sqrt ( 1-x^2) dx\n";
  std::cout << "\n";
  std::cout << "   E       Q_Estimate      Q_Exact\n";
  std::cout << "\n";

  f = new double[n];

  for ( e = 0; e <= 2 * n - 1; e++ )
  {
    if ( e == 0 )
    {
      for ( i = 0; i < n; i++ )
      {
        f[i] = 1.0;
      }
    }
    else
    {
      for ( i = 0; i < n; i++ )
      {
        f[i] = pow ( x[i], e );
      }
    }
    q = r8vec_dot_product ( n, w, f );
    q_exact = u_moment ( e );
    std::cout << "  " << std::setw(2) << e
         << "  " << std::setw(14) << q
         << "  " << std::setw(14) << q_exact << "\n";
  }

  delete [] f;
  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void uu_product_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    UU_PRODUCT_TEST tests UU_PRODUCT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  double r8_hi;
  double r8_lo;
  int seed;
  int test;
  double ui;
  double uiuj;
  double uj;
  double x;

  std::cout << "\n";
  std::cout << "UU_PRODUCT_TEST:\n";
  std::cout << "  UU_PRODUCT(I,J;X) = U(I,X) * U(J,X)\n";

  r8_lo = -1.0;
  r8_hi = +1.0;
  seed = 123456789;

  std::cout << "\n";
  std::cout << "   I   J      X               UI              UJ              UI*UJ       UU_PRODUCT\n";
  std::cout << "\n";
  for ( test = 1; test <= 10; test++ )
  {
    x = r8_uniform_ab ( r8_lo, r8_hi, seed );
    i = i4_uniform_ab ( 0, 6, seed );
    ui = u_polynomial_value ( i, x );
    j = i4_uniform_ab ( -1, 4, seed );
    uj = u_polynomial_value ( j, x );
    uiuj = uu_product ( i, j, x );
    std::cout << "  " << std::setw(2) << i
         << "  " << std::setw(2) << j
         << "  " << std::setw(14) << x
         << "  " << std::setw(14) << ui
         << "  " << std::setw(14) << uj
         << "  " << std::setw(14) << ui * uj
         << "  " << std::setw(14) << uiuj << "\n";
  }

  return;
}
//****************************************************************************80

void uu_product_integral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    UU_PRODUCT_INTEGRAL_TEST tests UU_PRODUCT_INTEGRAL.
//
//  Discussion:
//
//    This process should match the U_MASS_MATRIX computation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int i;
  int j;
  int n;

  std::cout << "\n";
  std::cout << "UU_PRODUCT_INTEGRAL_TEST:\n";
  std::cout << "  UU_PRODUCT_INTEGRAL computes the product integral\n";
  std::cout << "  of a pair of Chebyshev U polynomials U(i,x) and U(j,x).\n";
  std::cout << "  A(I,J) = integral ( -1 <=x <= +1 ) U(i,x) U(j,x) * sqrt ( 1 - x^2 ) dx\n";
  std::cout << "  0    if i is not equal to j;\n";
  std::cout << "  pi/2 if i = j.\n";

  n = 4;
  a = new double[(n+1)*(n+1)];
  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      a[i+j*(n+1)] = uu_product_integral ( i, j );
    }
  }

  r8mat_print ( n + 1, n + 1, a, "  U(i,x)*U(j,x) integral matrix:" );

  delete [] a;

  return;
}
//****************************************************************************80

void v_mass_matrix_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    V_MASS_MATRIX_TEST tests V_MASS_MATRIX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int n;

  std::cout << "\n";
  std::cout <<  "V_MASS_MATRIX_TEST:\n";
  std::cout <<  "  V_MASS_MATRIX computes the mass matrix for the\n";
  std::cout <<  "  Chebyshev polynomials V(i,x).\n";
  std::cout <<  "  A(I,J) = integral ( -1 <=x <= +1 ) V(i,x) V(j,x) sqrt(1+x)/sqrt(1-x) dx\n";
  std::cout <<  "  0  if i is not equal to j;\n";
  std::cout <<  "  pi if i = j.\n";

  n = 3;
  a = v_mass_matrix ( n );

  r8mat_print ( n + 1, n + 1, a, "  V mass matrix:" );

  delete [] a;

  return;
}
//****************************************************************************80

void v_moment_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    V_MOMENT_TEST tests V_MOMENT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  int e;
  double value;

  std::cout << "\n";
  std::cout << "V_MOMENT_TEST:\n";
  std::cout << "  V_MOMENT returns the value of\n";
  std::cout << "  integral ( -1 <=x <= +1 ) x^e * sqrt(1+x)/sqrt(1-x) dx\n";
  std::cout << "\n";
  std::cout << "   E       Integral\n";
  std::cout << "\n";
  for ( e = 0; e <= 10; e++ )
  {
    value = v_moment ( e );
    std::cout << "  " << std::setw(2) << e
         << "  " << std::setw(14) << value << "\n";
  }

  return;
}
//****************************************************************************80

void v_polynomial_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    V_POLYNOMIAL_TEST tests V_POLYNOMIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double *fx2;
  int n;
  int n_data;
  double x;
  double x_vec[1];

  std::cout << "\n";
  std::cout << "V_POLYNOMIAL_TEST:\n";
  std::cout << "  V_POLYNOMIAL evaluates the Chebyshev polynomial V(n,x).\n";
  std::cout << "\n";
  std::cout << "                   Tabulated      Computed\n";
  std::cout << "     N      X        V(n,x)        V(n,x)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    v_polynomial_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    if ( n < 0 )
    {
      continue;
    }

    x_vec[0] = x;
    fx2 = v_polynomial ( 1, n, x_vec );

    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(8)  << x
         << "  " << std::setw(14) << fx
         << "  " << std::setw(14) << fx2[n] << "\n";

    delete [] fx2;

  }

  return;
}
//****************************************************************************80

void v_polynomial_ab_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    V_POLYNOMIAL_AB_TEST tests V_POLYNOMIAL_AB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int m;
  int n;
  double *v;
  double *x;

  m = 11;
  n = 5;

  std::cout << "\n";
  std::cout << "V_POLYNOMIAL_AB_TEST:\n";
  std::cout << "  V_POLYNOMIAL_AB evaluates Chebyshev polynomials VAB(n,x)\n";
  std::cout << "  shifted from [-1,+1] to the domain [A,B].\n";
  std::cout << "\n";
  std::cout << "  Here, we will use the new domain [0,1]\n";
  std::cout << "  and the desired maximum polynomial degree will be N = 5.\n";

  a = 0.0;
  b = 1.0;
  x = r8vec_linspace_new ( m, a, b );
  
  v = v_polynomial_ab ( a, b, m, n, x );

  r8mat_print ( m, n + 1, v, "  Tables of V values:" );

  delete [] v;
  delete [] x;

  return;
}
//****************************************************************************80

void v_polynomial_ab_value_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    V_POLYNOMIAL_AB_VALUE_TEST tests V_POLYNOMIAL_AB_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  double fx2;
  int n;
  int n_data;
  double x01;

  std::cout << "\n";
  std::cout << "V_POLYNOMIAL_AB_VALUE_TEST:\n";
  std::cout << "  V_POLYNOMIAL_AB_VALUE evaluates Chebyshev polynomials VAB(n,x)\n";
  std::cout << "  shifted from [-1,+1] to the domain [A,B].\n";
  std::cout << "\n";
  std::cout << "  Here, we will use the new domain [0,1].\n";
  std::cout << "\n";
  std::cout << "                   Tabulated      Computed\n";
  std::cout << "     N      X01    V01(n,x)       V01(n,x)\n";
  std::cout << "\n";

  a = 0.0;
  b = 1.0;

  n_data = 0;

  for ( ; ; )
  {
    v_polynomial_01_values ( n_data, n, x01, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = v_polynomial_ab_value ( a, b, n, x01 );

    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(8)  << x01
         << "  " << std::setw(14) << fx
         << "  " << std::setw(14) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void v_polynomial_coefficients_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    V_POLYNOMIAL_COEFFICIENTS_TEST tests V_POLYNOMIAL_COEFFICIENTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  double *c2;
  int i;
  int j;
  int n;

  std::cout << "\n";
  std::cout << "V_POLYNOMIAL_COEFFICIENTS_TEST\n";
  std::cout << "  V_POLYNOMIAL_COEFFICIENTS determines the polynomial coefficients \n";
  std::cout << "  of V(n,x).\n";

  n = 5;

  c = v_polynomial_coefficients ( n );
 
  for ( i = 0; i <= n; i++ )
  {
    c2 = new double[i+1];
    for ( j = 0; j <= i; j++ )
    {
      c2[j] = c[i+j*(n+1)];
    }
    r8poly_print ( i, c2, "" );
    delete [] c2;
  }
  delete [] c;

  return;
}
//****************************************************************************80

void v_polynomial_plot_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    V_POLYNOMIAL_PLOT_TEST tests V_POLYNOMIAL_PLOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n_num = 6;
  int n_val[6];
  std::string output_filename;

  std::cout << "\n";
  std::cout << "V_POLYNOMIAL_PLOT_TEST\n";
  std::cout << "  V_POLYNOMIAL_PLOT plots selected\n";
  std::cout << "  Chebyshev polynomials V(n,x).\n";

  for ( i = 0; i <= 5; i++ )
  {
    n_val[i] = i;
  }

  output_filename = "v_polynomial_plot.png";

  v_polynomial_plot ( n_num, n_val, output_filename );

  return;
}
//****************************************************************************80

void v_polynomial_value_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    V_POLYNOMIAL_VALUE_TEST tests V_POLYNOMIAL_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "V_POLYNOMIAL_VALUE_TEST:\n";
  std::cout << "  V_POLYNOMIAL_VALUE evaluates the Chebyshev polynomial V(n,x).\n";
  std::cout << "\n";
  std::cout << "                   Tabulated      Computed\n";
  std::cout << "     N      X        V(n,x)        V(n,x)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    v_polynomial_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = v_polynomial_value ( n, x );

    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(8)  << x
         << "  " << std::setw(14) << fx
         << "  " << std::setw(14) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void v_polynomial_zeros_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    V_POLYNOMIAL_ZEROS_TEST tests V_POLYNOMIAL_ZEROS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *fx;
  int i;
  int n;
  int n_max = 5;
  double *z;

  std::cout << "\n";
  std::cout << "V_POLYNOMIAL_ZEROS_TEST:\n";
  std::cout << "  V_POLYNOMIAL_ZEROS returns zeroes of V(n,x).\n";
  std::cout << "\n";
  std::cout << "       N      X        V(n,x)\n";
  std::cout << "\n";

  for ( n = 1; n <= n_max; n++ )
  {
    z = v_polynomial_zeros ( n );
    fx = v_polynomial ( n, n, z );
    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(8) << n
           << "  " << std::setw(8) << z[i]
           << "  " << std::setw(14) << fx[i+n*n] << "\n";
    }
    std::cout << "\n";
    delete [] fx;
    delete [] z;
  }

  return;
}
//****************************************************************************80

void v_quadrature_rule_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    V_QUADRATURE_RULE_TEST tests V_QUADRATURE_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  int e;
  double *f;
  int i;
  int n;
  double q;
  double q_exact;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "V_QUADRATURE_RULE_TEST:\n";
  std::cout << "  V_QUADRATURE_RULE computes the quadrature rule\n";
  std::cout << "  associated with V(n,x);\n";

  n = 7;
  x = new double[n];
  w = new double[n];

  v_quadrature_rule ( n, x, w );

  r8vec2_print ( n, x, w, "    N      X            W" );

  std::cout << "\n";
  std::cout << "  Use the quadrature rule to estimate:\n";
  std::cout << "\n";
  std::cout << "    Q = Integral ( -1 <= X <= +1 ) X^E * sqrt(1+x)/sqrt(1-x) dx\n";
  std::cout << "\n";
  std::cout << "   E       Q_Estimate      Q_Exact\n";
  std::cout << "\n";

  f = new double[n];

  for ( e = 0; e <= 2 * n - 1; e++ )
  {
    if ( e == 0 )
    {
      for ( i = 0; i < n; i++ )
      {
        f[i] = 1.0;
      }
    }
    else
    {
      for ( i = 0; i < n; i++ )
      {
        f[i] = pow ( x[i], e );
      }
    }
    q = r8vec_dot_product ( n, w, f );
    q_exact = v_moment ( e );
    std::cout << "  " << std::setw(2) << e
         << "  " << std::setw(14) << q
         << "  " << std::setw(14) << q_exact << "\n";
  }

  delete [] f;
  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void vv_product_integral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    VV_PRODUCT_INTEGRAL_TEST tests VV_PRODUCT_INTEGRAL.
//
//  Discussion:
//
//    This process should match the V_MASS_MATRIX computation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int i;
  int j;
  int n;

  std::cout << "\n";
  std::cout << "VV_PRODUCT_INTEGRAL_TEST:\n";
  std::cout << "  VV_PRODUCT_INTEGRAL computes the product integral\n";
  std::cout << "  of a pair of Chebyshev V polynomials V(i,x) and V(j,x).\n";
  std::cout << "  A(I,J) = integral ( -1 <=x <= +1 ) V(i,x) V(j,x) * sqrt ( 1 + x ) / sqrt ( 1 - x ) dx\n";
  std::cout << "  0  if i is not equal to j;\n";
  std::cout << "  pi if i = j.\n";

  n = 4;
  a = new double[(n+1)*(n+1)];
  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      a[i+j*(n+1)] = vv_product_integral ( i, j );
    }
  }

  r8mat_print ( n + 1, n + 1, a, "  V(i,x)*V(j,x) integral matrix:" );

  delete [] a;

  return;
}
//****************************************************************************80

void w_mass_matrix_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    W_MASS_MATRIX_TEST tests W_MASS_MATRIX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int n;

  std::cout << "\n";
  std::cout <<  "W_MASS_MATRIX_TEST:\n";
  std::cout <<  "  W_MASS_MATRIX computes the mass matrix for the\n";
  std::cout <<  "  Chebyshev polynomials W(i,x).\n";
  std::cout <<  "  A(I,J) = integral ( -1 <=x <= +1 ) W(i,x) W(j,x) sqrt(1-x)/sqrt(1+x) dx\n";
  std::cout <<  "  0  if i is not equal to j;\n";
  std::cout <<  "  pi if i = j.\n";

  n = 3;
  a = w_mass_matrix ( n );

  r8mat_print ( n + 1, n + 1, a, "  W mass matrix:" );

  delete [] a;

  return;
}
//****************************************************************************80

void w_moment_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    W_MOMENT_TEST tests W_MOMENT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  int e;
  double value;

  std::cout << "\n";
  std::cout << "W_MOMENT_TEST:\n";
  std::cout << "  W_MOMENT returns the value of\n";
  std::cout << "  integral ( -1 <=x <= +1 ) x^e * sqrt(1-x)/sqrt(1+x) dx\n";
  std::cout << "\n";
  std::cout << "   E       Integral\n";
  std::cout << "\n";
  for ( e = 0; e <= 10; e++ )
  {
    value = w_moment ( e );
    std::cout << "  " << std::setw(2) << e
         << "  " << std::setw(14) << value << "\n";
  }

  return;
}
//****************************************************************************80

void w_polynomial_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    W_POLYNOMIAL_TEST tests W_POLYNOMIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double *fx2;
  int n;
  int n_data;
  double x;
  double x_vec[1];

  std::cout << "\n";
  std::cout << "W_POLYNOMIAL_TEST:\n";
  std::cout << "  W_POLYNOMIAL evaluates the Chebyshev polynomial W(n,x).\n";
  std::cout << "\n";
  std::cout << "                   Tabulated      Computed\n";
  std::cout << "     N      X        W(n,x)        W(n,x)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    w_polynomial_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    if ( n < 0 )
    {
      continue;
    }

    x_vec[0] = x;
    fx2 = w_polynomial ( 1, n, x_vec );

    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(8)  << x
         << "  " << std::setw(14) << fx
         << "  " << std::setw(14) << fx2[n] << "\n";

    delete [] fx2;

  }

  return;
}
//****************************************************************************80

void w_polynomial_ab_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    W_POLYNOMIAL_AB_TEST tests W_POLYNOMIAL_AB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int m;
  int n;
  double *v;
  double *x;

  m = 11;
  n = 5;

  std::cout << "\n";
  std::cout << "W_POLYNOMIAL_AB_TEST:\n";
  std::cout << "  W_POLYNOMIAL_AB evaluates Chebyshev polynomials WAB(n,x)\n";
  std::cout << "  shifted from [-1,+1] to the domain [A,B].\n";
  std::cout << "\n";
  std::cout << "  Here, we will use the new domain [0,1]\n";
  std::cout << "  and the desired maximum polynomial degree will be N = 5.\n";

  a = 0.0;
  b = 1.0;
  x = r8vec_linspace_new ( m, a, b );
  
  v = w_polynomial_ab ( a, b, m, n, x );

  r8mat_print ( m, n + 1, v, "  Tables of W values:" );

  delete [] v;
  delete [] x;

  return;
}
//****************************************************************************80

void w_polynomial_ab_value_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    W_POLYNOMIAL_AB_VALUE_TEST tests W_POLYNOMIAL_AB_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  double fx2;
  int n;
  int n_data;
  double x01;

  std::cout << "\n";
  std::cout << "W_POLYNOMIAL_AB_VALUE_TEST:\n";
  std::cout << "  W_POLYNOMIAL_AB_VALUE evaluates Chebyshev polynomials WAB(n,x)\n";
  std::cout << "  shifted from [-1,+1] to the domain [A,B].\n";
  std::cout << "\n";
  std::cout << "  Here, we will use the new domain [0,1].\n";
  std::cout << "\n";
  std::cout << "                   Tabulated      Computed\n";
  std::cout << "     N      X01    W01(n,x)       W01(n,x)\n";
  std::cout << "\n";

  a = 0.0;
  b = 1.0;

  n_data = 0;

  for ( ; ; )
  {
    w_polynomial_01_values ( n_data, n, x01, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = w_polynomial_ab_value ( a, b, n, x01 );

    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(8)  << x01
         << "  " << std::setw(14) << fx
         << "  " << std::setw(14) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void w_polynomial_coefficients_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    W_POLYNOMIAL_COEFFICIENTS_TEST tests W_POLYNOMIAL_COEFFICIENTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  double *c2;
  int i;
  int j;
  int n;

  std::cout << "\n";
  std::cout << "W_POLYNOMIAL_COEFFICIENTS_TEST\n";
  std::cout << "  W_POLYNOMIAL_COEFFICIENTS determines the polynomial coefficients \n";
  std::cout << "  of W(n,x).\n";

  n = 5;

  c = w_polynomial_coefficients ( n );
 
  for ( i = 0; i <= n; i++ )
  {
    c2 = new double[i+1];
    for ( j = 0; j <= i; j++ )
    {
      c2[j] = c[i+j*(n+1)];
    }
    r8poly_print ( i, c2, "" );
    delete [] c2;
  }
  delete [] c;

  return;
}
//****************************************************************************80

void w_polynomial_plot_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    W_POLYNOMIAL_PLOT_TEST tests W_POLYNOMIAL_PLOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n_num = 6;
  int n_val[6];
  std::string output_filename;

  std::cout << "\n";
  std::cout << "W_POLYNOMIAL_PLOT_TEST\n";
  std::cout << "  W_POLYNOMIAL_PLOT plots selected\n";
  std::cout << "  Chebyshev polynomials W(n,x).\n";

  for ( i = 0; i <= 5; i++ )
  {
    n_val[i] = i;
  }

  output_filename = "w_polynomial_plot.png";

  w_polynomial_plot ( n_num, n_val, output_filename );

  return;
}
//****************************************************************************80

void w_polynomial_value_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    W_POLYNOMIAL_VALUE_TEST tests W_POLYNOMIAL_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "W_POLYNOMIAL_VALUE_TEST:\n";
  std::cout << "  W_POLYNOMIAL_VALUE evaluates the Chebyshev polynomial W(n,x).\n";
  std::cout << "\n";
  std::cout << "                   Tabulated      Computed\n";
  std::cout << "     N      X        W(n,x)        W(n,x)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    w_polynomial_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = w_polynomial_value ( n, x );

    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(8)  << x
         << "  " << std::setw(14) << fx
         << "  " << std::setw(14) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void w_polynomial_zeros_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    W_POLYNOMIAL_ZEROS_TEST tests W_POLYNOMIAL_ZEROS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *fx;
  int i;
  int n;
  int n_max = 5;
  double *z;

  std::cout << "\n";
  std::cout << "W_POLYNOMIAL_ZEROS_TEST:\n";
  std::cout << "  W_POLYNOMIAL_ZEROS returns zeroes of W(n,x).\n";
  std::cout << "\n";
  std::cout << "       N      X        W(n,x)\n";
  std::cout << "\n";

  for ( n = 1; n <= n_max; n++ )
  {
    z = w_polynomial_zeros ( n );
    fx = w_polynomial ( n, n, z );
    for ( i = 0; i < n; i++ )
    {
      std::cout << "  " << std::setw(8) << n
           << "  " << std::setw(8) << z[i]
           << "  " << std::setw(14) << fx[i+n*n] << "\n";
    }
    std::cout << "\n";
    delete [] fx;
    delete [] z;
  }

  return;
}
//****************************************************************************80

void w_quadrature_rule_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    W_QUADRATURE_RULE_TEST tests W_QUADRATURE_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  int e;
  double *f;
  int i;
  int n;
  double q;
  double q_exact;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "W_QUADRATURE_RULE_TEST:\n";
  std::cout << "  W_QUADRATURE_RULE computes the quadrature rule\n";
  std::cout << "  associated with W(n,x);\n";

  n = 7;
  x = new double[n];
  w = new double[n];

  w_quadrature_rule ( n, x, w );

  r8vec2_print ( n, x, w, "    N      X            W" );

  std::cout << "\n";
  std::cout << "  Use the quadrature rule to estimate:\n";
  std::cout << "\n";
  std::cout << "    Q = Integral ( -1 <= X <= +1 ) X^E * sqrt(1-x)/sqrt(1+x) dx\n";
  std::cout << "\n";
  std::cout << "   E       Q_Estimate      Q_Exact\n";
  std::cout << "\n";

  f = new double[n];

  for ( e = 0; e <= 2 * n - 1; e++ )
  {
    if ( e == 0 )
    {
      for ( i = 0; i < n; i++ )
      {
        f[i] = 1.0;
      }
    }
    else
    {
      for ( i = 0; i < n; i++ )
      {
        f[i] = pow ( x[i], e );
      }
    }
    q = r8vec_dot_product ( n, w, f );
    q_exact = w_moment ( e );
    std::cout << "  " << std::setw(2) << e
         << "  " << std::setw(14) << q
         << "  " << std::setw(14) << q_exact << "\n";
  }

  delete [] f;
  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void ww_product_integral_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    WW_PRODUCT_INTEGRAL_TEST tests WW_PRODUCT_INTEGRAL.
//
//  Discussion:
//
//    This process should match the W_MASS_MATRIX computation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int i;
  int j;
  int n;

  std::cout << "\n";
  std::cout << "WW_PRODUCT_INTEGRAL_TEST:\n";
  std::cout << "  WW_PRODUCT_INTEGRAL computes the product integral\n";
  std::cout << "  of a pair of Chebyshev W polynomials W(i,x) and W(j,x).\n";
  std::cout << "  A(I,J) = integral ( -1 <=x <= +1 ) W(i,x) W(j,x) * sqrt ( 1 - x ) / sqrt ( 1 + x ) dx\n";
  std::cout << "  0  if i is not equal to j;\n";
  std::cout << "  pi if i = j.\n";

  n = 4;
  a = new double[(n+1)*(n+1)];
  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      a[i+j*(n+1)] = ww_product_integral ( i, j );
    }
  }

  r8mat_print ( n + 1, n + 1, a, "  W(i,x)*W(j,x) integral matrix:" );

  delete [] a;

  return;
}
