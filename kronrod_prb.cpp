# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "kronrod.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
double f ( double x );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for KRONROD_PRB.
//
//  Discussion:
//
//    KRONROD_PRB tests the KRONROD library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  std::cout << "\n";
  std::cout << "KRONROD_PRB:\n";
  std::cout << "  C++ version.\n";
  std::cout << "  Test the KRONROD library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "KRONROD_PRB:\n";
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
//    TEST01 tests the code for the odd case N = 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  double eps;
  int i;
  int i2;
  int n = 3;
  double s;
  double *w1;
  double *w2;
  double wg[3] = {
    0.555555555555555555556, 
    0.888888888888888888889, 
    0.555555555555555555556 };
  double wk[7] = {
    0.104656226026467265194, 
    0.268488089868333440729, 
    0.401397414775962222905, 
    0.450916538658474142345, 
    0.401397414775962222905, 
    0.268488089868333440729, 
    0.104656226026467265194 };
  double *x;
  double xg[3] = {
   -0.77459666924148337704, 
    0.0, 
    0.77459666924148337704 };
  double xk[7]= {
   -0.96049126870802028342, 
   -0.77459666924148337704, 
   -0.43424374934680255800, 
    0.0, 
    0.43424374934680255800, 
    0.77459666924148337704, 
    0.96049126870802028342 };

  std::cout << "\n";
  std::cout << "TEST01\n";
  std::cout << "  Request KRONROD to compute the Gauss rule\n";
  std::cout << "  of order 3, and the Kronrod extension of\n";
  std::cout << "  order 3+4=7.\n";
  std::cout << "\n";
  std::cout << "  Compare to exact data.\n";

  eps = 0.000001;
  w1 = new double[n+1];
  w2 = new double[n+1];
  x = new double[n+1];

  kronrod ( n, eps, x, w1, w2 );

  std::cout << "\n";
  std::cout << "  KRONROD returns 3 vectors of length " << n + 1 << "\n";
  std::cout << "\n";
  std::cout << "     I      X               WK              WG\n";
  std::cout << "\n";
  for ( i = 1; i <= n + 1; i++ )
  {
    std::cout << "  " << std::setw(4) << i
         << "  " << std::setw(14) << x[i-1]
         << "  " << std::setw(14) << w1[i-1]
         << "  " << std::setw(14) << w2[i-1] << "\n";
  }

  std::cout << "\n";
  std::cout << "               Gauss Abscissas\n";
  std::cout << "            Exact           Computed\n";
  std::cout << "\n";
  for ( i = 1; i <= n; i++ )
  {
    if ( 2 * i <= n + 1 )
    {
      i2 = 2 * i;
      s = -1.0;
    }
    else
    {
      i2 = 2 * ( n + 1 ) - 2 * i;
      s = +1.0;
    }
    std::cout << "  " << std::setw(4) << i
         << "  " << std::setw(14) << xg[i-1]
         << "  " << std::setw(14) << s * x[i2-1] << "\n";
  }
  std::cout << "\n";
  std::cout << "               Gauss Weights\n";
  std::cout << "            Exact           Computed\n";
  std::cout << "\n";
  for ( i = 1; i <= n; i++ )
  {
    if ( 2 * i <= n + 1 )
    {
      i2 = 2 * i;
    }
    else
    {
      i2 = 2 * ( n + 1 ) - 2 * i;
    }
    std::cout << "  " << std::setw(4) << i
         << "  " << std::setw(14) << wg[i-1]
         << "  " << std::setw(14) << w2[i2-1] << "\n";
  }

  std::cout << "\n";
  std::cout << "             Gauss Kronrod Abscissas\n";
  std::cout << "            Exact           Computed\n";
  std::cout << "\n";
  for ( i = 1; i <= 2 * n + 1; i++ )
  {
    if ( i <= n + 1 )
    {
      i2 = i;
      s = -1.0;
    }
    else
    {
      i2 = 2 * ( n + 1 ) - i;
      s = +1.0;
    }
    std::cout << "  " << std::setw(4) << i
         << "  " << std::setw(14) << xk[i-1]
         << "  " << std::setw(14) << s * x[i2-1] << "\n";
  }
  std::cout << "\n";
  std::cout << "             Gauss Kronrod Weights\n";
  std::cout << "            Exact           Computed\n";
  std::cout << "\n";
  for ( i = 1; i <= 2 * n + 1; i++ )
  {
    if ( i <= n + 1 )
    {
      i2 = i;
    }
    else
    {
      i2 = 2 * ( n + 1 ) - i;
    }
    std::cout << "  " << std::setw(4) << i
         << "  " << std::setw(14) << wk[i-1]
         << "  " << std::setw(14) << w1[i2-1] << "\n";
  }

  delete [] w1;
  delete [] w2;
  delete [] x;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests the code for the even case N = 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  double eps;
  int i;
  int i2;
  int n = 4;
  double s;
  double *w1;
  double *w2;
  double *x;

  std::cout << "\n";
  std::cout << "TEST02\n";
  std::cout << "  Request KRONROD to compute the Gauss rule\n";
  std::cout << "  of order 4, and the Kronrod extension of\n";
  std::cout << "  order 4+5=9.\n";

  eps = 0.000001;
  w1 = new double[n+1];
  w2 = new double[n+1];
  x = new double[n+1];

  kronrod ( n, eps, x, w1, w2 );

  std::cout << "\n";
  std::cout << "  KRONROD returns 3 vectors of length " << n + 1 <<"\n";
  std::cout << "\n";
  std::cout << "     I      X               WK              WG\n";
  std::cout << "\n";
  for ( i = 1; i <= n + 1; i++ )
  {
    std::cout << "  " << std::setw(4) << i
         << "  " << std::setw(14) << x[i-1]
         << "  " << std::setw(14) << w1[i-1]
         << "  " << std::setw(14) << w2[i-1] << "\n";
  }

  delete [] w1;
  delete [] w2;
  delete [] x;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 uses the program to estimate an integral.
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
  double eps;
  double exact = 1.5643964440690497731;
  int i;
  double i1;
  double i2;
  int n;
  double *w1;
  double *w2;
  double *x;

  std::cout << "\n";
  std::cout << "TEST03\n";
  std::cout << "  Call Kronrod to estimate the integral of a function.\n";
  std::cout << "  Keep trying until the error is small.\n";
//
//  EPS just tells KRONROD how carefully it must compute X, W1 and W2.
//  It is NOT a statement about the accuracy of your integral estimate!
//
  eps = 0.000001;
//
//  Start the process with a 1 point rule.
//
  n = 1;

  for ( ; ; )
  {
//
//  Make space.
//
    w1 = new double[n+1];
    w2 = new double[n+1];
    x = new double[n+1];

    kronrod ( n, eps, x, w1, w2 );
//
//  Compute the estimates.
//  There are two complications here:
//
//  1) Both rules use all the points.  However, the lower order rule uses
//     a zero weight for the points it doesn't need.
//
//  2) The points X are all positive, and are listed in descending order.
//     this means that 0 is always in the list, and always occurs as the
//     last member.  Therefore, the integral estimates should use the
//     function value at 0 once, and the function values at the other
//     X values "twice", that is, once at X and once at -X.
//
    i1 = w1[n] * f ( x[n] );
    i2 = w2[n] * f ( x[n] );

    for ( i = 0; i < n; i++ )
    {
      i1 = i1 + w1[i] * ( f ( - x[i] ) + f ( x[i] ) );
      i2 = i2 + w2[i] * ( f ( - x[i] ) + f ( x[i] ) );
    }

    delete [] w1;
    delete [] w2;
    delete [] x;

    if ( fabs ( i1 - i2 ) < 0.0001 )
    {
      std::cout << "\n";
      std::cout << "  Error tolerance satisfied with N = " << n << "\n";
      std::cout << "  Coarse integral estimate = " << std::setprecision ( 8 ) << i1 << "\n";
      std::cout << "  Fine   integral estimate = " << i2 << "\n";
      std::cout << "  Error estimate = " << fabs ( i2 - i1 ) << "\n";
      std::cout << "  Actual error = " << fabs ( exact - i2 ) << "\n";
      break;
    }

    if ( 25 < n )
    {
      std::cout << "\n";
      std::cout << "  Error tolerance failed even for n = " << n << "\n";
      std::cout << "  Canceling iteration, and accepting bad estimates!\n";
      std::cout << "  Coarse integral estimate = " << i1 << "\n";
      std::cout << "  Fine   integral estimate = " << i2 << "\n";
      std::cout << "  Error estimate = " << fabs ( i2 - i1 ) << "\n";
      std::cout << "  Actual error = " << fabs ( exact - i2 ) << "\n";
      break;
    }
    n = 2 * n + 1;
  }

  return;
}
//****************************************************************************80

double f ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    F is a function whose integral from -1 to +1 is to be estimated.
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
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double F, the value.
//
{
  double value;

  value = 1.0 / ( x * x + 1.005 );

  return value;
}
