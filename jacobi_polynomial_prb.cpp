# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "jacobi_polynomial.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for JACOBI_POLYNOMIAL_PRB.
//
//  Discussion:
//
//    JACOBI_POLYNOMIAL_PRB tests the JACOBI_POLYNOMIAL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  std::cout << "\n";
  std::cout << "JACOBI_POLYNOMIAL_PRB\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the JACOBI_POLYNOMIAL library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "JACOBI_POLYNOMIAL_PRB\n";
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
//    TEST01 tests J_POLYNOMIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double e;
  double fx1;
  double fx2;
  double *fx2_vec;
  int m;
  int n;
  int n_data;
  double x;
  double x_vec[1];

  std::cout << "\n";
  std::cout << "TEST01:\n";
  std::cout << "  J_POLYNOMIAL_VALUES stores values of\n";
  std::cout << "  the Jacobi polynomials.\n";
  std::cout << "  J_POLYNOMIAL evaluates the polynomial.\n";
  std::cout << "\n";
  std::cout << "                                    Tabulated                 Computed\n";
  std::cout << "     N     A     B        X           J(N,A,B,X)                    J(N,A,B,X)                     Error\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    j_polynomial_values ( n_data, n, a, b, x, fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    m = 1;
    x_vec[0] = x;
    fx2_vec = j_polynomial ( m, n, a, b, x_vec );
    fx2 = fx2_vec[0+n*1];
    e = fx1 - fx2;

    std::cout << "  " << std::setw(4) << n
         << "  " << std::setw(6) << a
         << "  " << std::setw(6) << b
         << "  " << std::setw(6) << x
         << "  " << std::setw(24) << fx1
         << "  " << std::setw(24) << fx2
         << "  " << std::setw(8) << e << "\n";
    
    delete [] fx2_vec;
  }
  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests J_POLYNOMIAL_ZEROS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double a_test[3] = { 0.5, 1.0, 2.0 };
  double b;
  double b_test[3] = { 0.5, 1.5, 0.5 };
  int degree;
  double *hz;
  int test;
  int test_num = 3;
  std::string title;
  double *z;

  std::cout << "\n";
  std::cout << "TEST02:\n";
  std::cout << "  J_POLYNOMIAL_ZEROS computes the zeros of J(n,a,b,x);\n";
  std::cout << "  Check by calling J_POLYNOMIAL there.\n";

  for ( test = 0; test < test_num; test++ )
  {
    a = a_test[test];
    b = b_test[test];

    for ( degree = 1; degree <= 5; degree++ )
    {
      z = j_polynomial_zeros ( degree, a, b );
      title = "Zeros for J(" + i4_to_string ( degree ) + "," 
        + r8_to_string ( a, "%f" ) + "," + r8_to_string ( b, "%f" ) + ")";
      r8vec_print ( degree, z, title );

      hz = j_polynomial ( degree, degree, a, b, z );
      title = "Evaluate J(" + i4_to_string ( degree ) + "," 
        + r8_to_string ( a, "%f" ) + "," + r8_to_string ( b, "%f" ) + ")";
      r8vec_print ( degree, hz + degree * degree, title );

      delete [] hz;
      delete [] z;
    }
  }
  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests J_QUADRATURE_RULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int i;
  int j;
  double *ji;
  double *jj;
  int k;
  int n;
  double q;
  double q_exact;
  double *w;
  double *x;

  std::cout << "\n";
  std::cout << "TEST03:\n";
  std::cout << "  J_QUADRATURE_RULE computes the quadrature rule\n";
  std::cout << "  associated with J(n,a,b,x);\n";

  n = 7;
  a = 1.0;
  b = 2.5;

  x = new double[n];
  w = new double[n];

  j_quadrature_rule ( n, a, b, x, w );

  r8vec2_print ( n, x, w, "      X            W" );

  std::cout << "\n";
  std::cout << "  Use the quadrature rule to estimate:\n";
  std::cout << "\n";
  std::cout << "    Q = Integral (-1<x<+1) J(i,a,b,x) J(j,a,b,x) (1-x)^a (1+x)^b dx\n";
  std::cout << "\n";
  std::cout << "   I   J      Q_Estimate         Q_Exact\n";
  std::cout << "\n";

  for ( i = 0; i <= 5; i++ )
  {
    ji = j_polynomial ( n, i, a, b, x );
    for ( j = i; j <= 5; j++ )
    {
      jj = j_polynomial ( n, j, a, b, x );
      q = 0.0;
      for ( k = 0; k < n; k++ )
      {
        q = q + w[k] * ji[k+i*n] * jj[k+j*n];
      }
      q_exact = j_double_product_integral ( i, j, a, b );
      std::cout << "  " << std::setw(2) << i
           << "  " << std::setw(2) << j
           << "  " << std::setw(14) << q
           << "  " << std::setw(14) << q_exact << "\n";
      delete [] jj;
    }
    delete [] ji;
  }
  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests J_DOUBLE_PRODUCT_INTEGRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int i;
  int j;
  double q;

  std::cout << "\n";
  std::cout << "TEST04:\n";
  std::cout << "  J_DOUBLE_PRODUCT_INTEGRAL returns the weighted integral of\n";
  std::cout << "  J(i,a,b,x) * J(j,a,b,x);\n";

  a = 1.0;
  b = 2.5;

  std::cout << "\n";
  std::cout << "    Q = Integral (-1<x<+1) J(i,a,b,x) J(j,a,b,x) (1-x)^a (1+x)^b dx\n";
  std::cout << "\n";
  std::cout << "   I   J      Q\n";
  std::cout << "\n";

  for ( i = 0; i <= 5; i++ )
  {
    for ( j = i; j <= 5; j++ )
    {
      q = j_double_product_integral ( i, j, a, b );
      std::cout << "  " << std::setw(2) << i
           << "  " << std::setw(2) << j
           << "  " << std::setw(14) << q << "\n";
    }
  }
  return;
}

