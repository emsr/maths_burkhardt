# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "asa109.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA109_PRB.
//
//  Discussion:
//
//    ASA109_PRB tests the ASA109 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 September 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  std::cout << "\n";
  std::cout << "ASA109_PRB:\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the ASA109 library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "ASA109_PRB:\n";
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
//    TEST01 demonstrates the use of XINBTA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2013
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double beta_log;
  double fx;
  int ifault;
  int n_data;
  double x;
  double x2;

  std::cout << "\n";
  std::cout << "TEST01:\n";
  std::cout << "  XINBTA inverts the incomplete Beta function.\n";
  std::cout << "  Given CDF, it computes an X.\n";
  std::cout << "\n";
  std::cout << "           A           B           CDF    "
       << "    X                         X\n";
  std::cout << "                                          "
       << "    (Tabulated)               (XINBTA)            DIFF\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( n_data, a, b, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    beta_log = lgamma ( a )
             + lgamma ( b )
             - lgamma ( a + b );

    x2 = xinbta ( a, b, beta_log, fx, ifault );

    std::cout << "  " << std::setprecision(4) << std::setw(10) << a
         << "  " << std::setprecision(4) << std::setw(10) << b
         << "  " << std::setprecision(4) << std::setw(10) << fx
         << "  " << std::setprecision(16) << std::setw(24) << x
         << "  " << std::setprecision(16) << std::setw(24) << x2
         << "  " << std::setprecision(4) << std::setw(10) << fabs ( x - x2 ) << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests BETA_INC_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 September 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double beta_log;
  double fx;
  double fx2;
  int ifault;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TEST02:\n";
  std::cout << "  BETA_INC_VALUES stores values of\n";
  std::cout << "  the incomplete Beta function.\n";
  std::cout << "\n";
  std::cout << "      A            B            X            BETA_INC(A,B)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( n_data, a, b, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    beta_log = lgamma ( a )
             + lgamma ( b )
             - lgamma ( a + b );

    ifault = 0;
    fx2 = betain ( x, a, b, beta_log, ifault );

    std::cout << "  " << std::setprecision(4) << std::setw(10) << a
         << "  " << std::setprecision(4) << std::setw(10) << b
         << "  " << std::setprecision(4) << std::setw(10) << x
         << "  " << std::setprecision(16) << std::setw(24) << fx
         << "  " << std::setprecision(16) << std::setw(24) << fx2
         << "  " << std::setprecision(4) << std::setw(10) << fabs ( fx - fx2 ) << "\n";
  }
  return;
}
