# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "toms443.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TOMS443_PRB.
//
//  Discussion:
//
//    TOMS443_PRB tests the TOMS443 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  std::cout << "\n";
  std::cout << "TOMS443_PRB\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the TOMS443 library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "TOMS433_PRB\n";
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
//    TEST01 tests WEW_A
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  double en;
  int n_data;
  double w1;
  double w2;
  double x;

  std::cout << "\n";
  std::cout << "TEST01\n";
  std::cout << "  Test WEW_A to evaluate\n";
  std::cout << "  Lambert's W function.\n";
  std::cout << "\n";
  std::cout << "          X             Exact             Computed      Error\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    lambert_w_values ( n_data, x, w1 );

    if ( n_data <= 0 )
    {
      break;
    }

    if ( x == 0.0 )
    {
      w2 = 0.0;
    }
    else
    {
      w2 = wew_a ( x, en );
    }

    std::cout << std::setw(14) << x << "  "
         << std::setw(16) << w1 << "  "
         << std::setw(16) << w2 << "  "
         << std::setw(10) << fabs ( w1 - w2 ) << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests WEW_B
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  double en;
  int n_data;
  double w1;
  double w2;
  double x;

  std::cout << "\n";
  std::cout << "TEST02\n";
  std::cout << "  Test WEW_B to evaluate\n";
  std::cout << "  Lambert's W function.\n";
  std::cout << "\n";
  std::cout << "          X             Exact             Computed      Error\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    lambert_w_values ( n_data, x, w1 );

    if ( n_data <= 0 )
    {
      break;
    }

    if ( x == 0.0 )
    {
      w2 = 0.0;
    }
    else
    {
      w2 = wew_b ( x, en );
    }

    std::cout << std::setw(14) << x << "  "
         << std::setw(16) << w1 << "  "
         << std::setw(16) << w2 << "  "
         << std::setw(10) << fabs ( w1 - w2 ) << "\n";
  }

  return;
}

