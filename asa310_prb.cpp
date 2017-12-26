# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "asa310.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA310_PRB.
//
//  Discussion:
//
//    ASA310_PRB tests the ASA310 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  std::cout << "\n";
  std::cout << "ASA310_PRB:\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the ASA310 library.\n";

  test01 ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "ASA310_PRB:\n";
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
//    TEST01 demonstrates the use of NCBETA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double errmax = 1.0E-10;
  double fx;
  double fx2;
  int ifault;
  double lambda;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TEST01:\n";
  std::cout << "  NCBETA computes the noncentral incomplete Beta function.\n";
  std::cout << "  Compare to tabulated values.\n";
  std::cout << "\n";
  std::cout << "      A        B     LAMBDA        X      "
       << "    FX                        FX2\n";
  std::cout << "                                          "
       << "    (Tabulated)               (NCBETA)            DIFF\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_noncentral_cdf_values ( &n_data, &a, &b, &lambda, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = ncbeta ( a, b, lambda, x, errmax, &ifault );

    std::cout << "  " << std::setprecision(2) << std::setw(7) << a
         << "  " << std::setprecision(2) << std::setw(7) << b
         << "  " << std::setprecision(3) << std::setw(7) << lambda
         << "  " << std::setprecision(4) << std::setw(10) << x
         << "  " << std::setprecision(16) << std::setw(24) << fx
         << "  " << std::setprecision(16) << std::setw(24) << fx2
         << "  " << std::setprecision(4) << std::setw(10) << fabs ( fx - fx2 ) << "\n";
  }

  return;
}
