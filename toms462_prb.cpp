# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "toms462.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TOMS462_PRB.
//
//  Discussion:
//
//    TOMS462_PRB tests the TOMS462 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  std::cout << "\n";
  std::cout << "TOMS462_PRB\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the TOMS462 library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "TOMS462_PRB\n";
  std::cout << "  Normal end of execution.\n";
  std::cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************807
//
//  Purpose:
//
//    TEST01 tests BIVNOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double cdf;
  double expect;
  double r;
  double x;
  double y;

  std::cout << "\n";
  std::cout << "TEST01\n";
  std::cout << "  Compare BIVNOR with some simple data\n";
  std::cout << "  with 3 digit accuracy.\n";
  std::cout << "\n";
  std::cout << "       X         Y          R          P               P\n";
  std::cout << "                                      (Tabulated)     (BIVNOR)\n";
  std::cout << "\n";

  x =  0.8;
  y = -1.5;
  r =  -0.9;
  expect = 0.148;
  cdf = bivnor ( x, y, r );
  std::cout << "  " << std::setw(9) << x
       << "  " << std::setw(8) << y
       << "  " << std::setw(8) << r
       << "  " << std::setw(14) << expect
       << "  " << std::setw(14) << cdf << "\n";

  x =  0.6;
  y = -1.4;
  r =  -0.7;
  expect = 0.208;
  cdf = bivnor ( x, y, r );
  std::cout << "  " << std::setw(9) << x
       << "  " << std::setw(8) << y
       << "  " << std::setw(8) << r
       << "  " << std::setw(14) << expect
       << "  " << std::setw(14) << cdf << "\n";

  x =  0.2;
  y = -1.0;
  r =  -0.5;
  expect = 0.304;
  cdf = bivnor ( x, y, r );
  std::cout << "  " << std::setw(9) << x
       << "  " << std::setw(8) << y
       << "  " << std::setw(8) << r
       << "  " << std::setw(14) << expect
       << "  " << std::setw(14) << cdf << "\n";

  x = -1.2;
  y =  0.1;
  r =   0.0;
  expect = 0.407;
  cdf = bivnor ( x, y, r );
  std::cout << "  " << std::setw(9) << x
       << "  " << std::setw(8) << y
       << "  " << std::setw(8) << r
       << "  " << std::setw(14) << expect
       << "  " << std::setw(14) << cdf << "\n";

  x = -1.2;
  y = -0.1;
  r =   0.3;
  expect = 0.501;
  cdf = bivnor ( x, y, r );
  std::cout << "  " << std::setw(9) << x
       << "  " << std::setw(8) << y
       << "  " << std::setw(8) << r
       << "  " << std::setw(14) << expect
       << "  " << std::setw(14) << cdf << "\n";

  x = -0.4;
  y = -0.9;
  r =   0.6;
  expect = 0.601;
  cdf = bivnor ( x, y, r );
  std::cout << "  " << std::setw(9) << x
       << "  " << std::setw(8) << y
       << "  " << std::setw(8) << r
       << "  " << std::setw(14) << expect
       << "  " << std::setw(14) << cdf << "\n";

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests BIVNOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double fxy1;
  double fxy2;
  int n_data;
  double r;
  double x;
  double y;

  std::cout << "\n";
  std::cout << "TEST02\n";
  std::cout << "  Compare BIVNOR with some tabulated data.\n";
  std::cout << "\n";
  std::cout << "      X          Y          ";
  std::cout << "R           P                         P";
  std::cout << "                      DIFF\n";
  std::cout << "                                ";
  std::cout << "       (Tabulated)               (BIVNOR)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bivariate_normal_cdf_values ( n_data, x, y, r, fxy1 );

    if ( n_data == 0 )
    {
      break;
    }
//
//  BIVNOR computes the "tail" of the probability, and we want the
//  initial part//
//
    fxy2 = bivnor ( - x, - y, r );

    std::cout << "  " << std::setw(8) << x
         << "  " << std::setw(8) << y
         << "  " << std::setw(8) << r
         << "  " << std::setw(24) << fxy1
         << "  " << std::setw(24) << fxy2
         << "  " << std::setw(10) << r8_abs ( fxy1 - fxy2 ) << "\n";
  }
  return;
}
