# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "beta_nc.H"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BETA_NC_PRB.
//
//  Discussion:
//
//    BETA_NC_PRB tests the BETA_NC library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  std::cout << "\n";
  std::cout << "BETA_NC_PRB:\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the BETA_NC library.\n";

  test01 ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "BETA_NC_PRB:\n";
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
//    TEST01 demonstrates the use of BETA_NONCENTRAL_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double error_max;
  double fx;
  double fx2;
  double lambda;
  int n_data;
  double x;

  error_max = 1.0E-10;

  std::cout << "\n";
  std::cout << "TEST01:\n";
  std::cout << "  BETA_NONCENTRAL_CDF computes the noncentral incomplete \n";
  std::cout << "  Beta function.\n";
  std::cout << "  Compare to tabulated values.\n";
  std::cout << "\n";
  std::cout << "      A      B     LAMBDA        X      "
       << "    FX                        FX2\n";
  std::cout << "                                        "
       << "    (Tabulated)               (computed)          DIFF\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_noncentral_cdf_values ( &n_data, &a, &b, &lambda, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = beta_noncentral_cdf ( a, b, lambda, x, error_max );

    std::cout << "  " << std::setprecision(2) << std::setw(5) << a
         << "  " << std::setprecision(2) << std::setw(5) << b
         << "  " << std::setprecision(3) << std::setw(7) << lambda
         << "  " << std::setprecision(4) << std::setw(10) << x
         << "  " << std::setprecision(16) << std::setw(24) << fx
         << "  " << std::setprecision(16) << std::setw(24) << fx2
         << "  " << std::setprecision(4) << std::setw(10) << r8_abs ( fx - fx2 ) << "\n";
  }

  return;
}
