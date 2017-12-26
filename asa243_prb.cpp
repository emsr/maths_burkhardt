# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

# include "asa243.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA243_PRB.
//
//  Discussion:
//
//    ASA243_PRB tests the ASA243 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  std::cout << "\n";
  std::cout << "ASA243_PRB:\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the ASA243 library.\n";

  test01 ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "ASA243_PRB:\n";
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
//    TEST01 demonstrates the use of TNC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  double delta;
  int df;
  double df_real;
  double fx;
  double fx2;
  int ifault;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TEST01:\n";
  std::cout << "  TNC computes the noncentral Student T\n";
  std::cout << "  Cumulative Density Function.\n";
  std::cout << "  Compare with tabulated values.\n";
  std::cout << "\n";
  std::cout << "        X         LAMBDA        DF     "
       << " CDF             CDF           DIFF\n";
  std::cout << "                                       "
       << " Tabulated       PRNCST\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    student_noncentral_cdf_values ( &n_data, &df, &delta, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    df_real = double( df );

    fx2 = tnc ( x, df_real, delta, &ifault );

    std::cout << "  " << std::setw(10) << std::setprecision(4) << x
         << "  " << std::setw(10) << std::setprecision(4) << delta
         << "  " << std::setw(8)                     << df
         << "  " << std::setw(24) << std::setprecision(16) << fx
         << "  " << std::setw(24) << std::setprecision(16) << fx2
         << "  " << std::setw(10) << std::setprecision(4) << fabs ( fx - fx2 ) << "\n";
  }

  return;
}
