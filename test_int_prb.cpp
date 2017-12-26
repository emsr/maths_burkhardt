# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "test_int.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TEST_INT_PRB.
//
//  Discussion:
//
//    TEST_INT_PRB tests the TEST_INT library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  std::cout << "\n";
  std::cout << "TEST_INT_PRB\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the TEST_INT library.\n";
 
  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "TEST_INT_PRB\n";
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
//    TEST01 applies a composite midpoint rule to finite interval 1D problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  double error;
  double exact;
  int int_log;
  int int_num;
  int prob;
  int prob_num;
  double result;

  std::cout << "\n";
  std::cout << "TEST01\n";
  std::cout << "  Composite midpoint rule,\n";
  std::cout << "  for 1D finite interval problems.\n";

  prob_num = p00_prob_num ( );

  std::cout << "\n";
  std::cout << "  Prob   Ints         Exact       Error\n";
  std::cout << "                      Approx\n";
//
//  Pick a problem.
//
  for ( prob = 1; prob <= prob_num; prob++ )
  {
    exact = p00_exact ( prob );

    std::cout << "\n";
    std::cout << "  " << std::setw(4) << prob
         << "  " << "    "
         << "  " << std::setw(14) << exact << "\n";
//
//  Pick a number of subintervals.
//
    for ( int_log = 0; int_log <= 7; int_log++ )
    {
      int_num = i4_power ( 2, int_log );

      result = p00_midpoint ( prob, int_num );

      error = r8_abs ( exact - result );

      std::cout << "  " << "    "
           << "  " << std::setw(4) << int_num
           << "  " << std::setw(14) << result
           << "  " << std::setw(14) << error << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 applies a composite Simpson rule to finite interval 1D problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  double error;
  double exact;
  int int_log;
  int int_num;
  int prob;
  int prob_num;
  double result;

  std::cout << "\n";
  std::cout << "TEST02\n";
  std::cout << "  Composite Simpson rule,\n";
  std::cout << "  for 1D finite interval problems.\n";

  prob_num = p00_prob_num ( );

  std::cout << "\n";
  std::cout << "  Prob   Ints         Exact       Error\n";
  std::cout << "                      Approx\n";
//
//  Pick a problem.
//
  for ( prob = 1; prob <= prob_num; prob++ )
  {
//
//  Some problems have singularities that kill the calculation.
//
    exact = p00_exact ( prob );

    std::cout << "\n";
    std::cout << "  " << std::setw(4) << prob
         << "  " << "    "
         << "  " << std::setw(14) << exact << "\n";
//
//  Pick a number of subintervals.
//
    for ( int_log = 0; int_log <= 10; int_log++ )
    {
      int_num = i4_power ( 2, int_log );

      result = p00_simpson ( prob, int_num );

      error = r8_abs ( exact - result );

      std::cout << "  " << "    "
           << "  " << std::setw(4) << int_num
           << "  " << std::setw(14) << result
           << "  " << std::setw(14) << error << "\n";
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
//    TEST03 applies a Monte Carlo rule to finite interval 1D problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2006
//
//  Author:
//
//    John Burkardt
//
{
  double error;
  double exact;
  int int_log;
  int int_num;
  int prob;
  int prob_num;
  double result;

  std::cout << "\n";
  std::cout << "TEST03\n";
  std::cout << "  Monte Carlo rule,\n";
  std::cout << "  for 1D finite interval problems.\n";

  prob_num = p00_prob_num ( );

  std::cout << "\n";
  std::cout << "  Prob   Ints         Exact       Error\n";
  std::cout << "                      Approx\n";
//
//  Pick a problem.
//
  for ( prob = 1; prob <= prob_num; prob++ )
  {
    exact = p00_exact ( prob );

    std::cout << "\n";
    std::cout << "  " << std::setw(4) << prob
         << "  " << "    "
         << "  " << std::setw(14) << exact << "\n";
//
//  Pick a number of points.
//
    for ( int_log = 0; int_log <= 10; int_log++ )
    {
      int_num = i4_power ( 2, int_log );

      result = p00_montecarlo ( prob, int_num );

      error = r8_abs ( exact - result );

      std::cout << "  " << "    "
           << "  " << std::setw(4) << int_num
           << "  " << std::setw(14) << result
           << "  " << std::setw(14) << error << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 applies a composite Gauss-Legendre rule to finite interval 1D problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  double error;
  double exact;
  int int_log;
  int int_num;
  int prob;
  int prob_num;
  double result;

  std::cout << "\n";
  std::cout << "TEST04\n";
  std::cout << "  Use a composite 4 point Gauss-Legendre rule,\n";
  std::cout << "  for 1D finite interval problems.\n";

  prob_num = p00_prob_num ( );

  std::cout << "\n";
  std::cout << "  Prob   Ints         Exact       Error\n";
  std::cout << "                      Approx\n";
//
//  Pick a problem.
//
  for ( prob = 1; prob <= prob_num; prob++ )
  {
    exact = p00_exact ( prob );

    std::cout << "\n";
    std::cout << "  " << std::setw(4) << prob
         << "  " << "    "
         << "  " << std::setw(14) << exact << "\n";
//
//  Pick a number of subintervals.
//
    for ( int_log = 0; int_log <= 10; int_log++ )
    {
      int_num = i4_power ( 2, int_log );

      result = p00_gauss_legendre ( prob, int_num );

      error = r8_abs ( exact - result );

      std::cout << "  " << "    "
           << "  " << std::setw(4) << int_num
           << "  " << std::setw(14) << result
           << "  " << std::setw(14) << error << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 applies a composite trapezoid rule to finite interval 1D problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  double error;
  double exact;
  int int_log;
  int int_num;
  int prob;
  int prob_num;
  double result;

  std::cout << "\n";
  std::cout << "TEST05\n";
  std::cout << "  Composite trapezoid rule,\n";
  std::cout << "  for 1D finite interval problems.\n";

  prob_num = p00_prob_num ( );

  std::cout << "\n";
  std::cout << "  Prob   Ints         Exact       Error\n";
  std::cout << "                      Approx\n";
//
//  Pick a problem.
//
  for ( prob = 1; prob <= prob_num; prob++ )
  {
    exact = p00_exact ( prob );

    std::cout << "\n";
    std::cout << "  " << std::setw(4) << prob
         << "  " << "    "
         << "  " << std::setw(14) << exact << "\n";
//
//  Pick a number of subintervals.
//
    for ( int_log = 0; int_log <= 10; int_log++ )
    {
      int_num = i4_power ( 2, int_log );

      result = p00_trapezoid ( prob, int_num );

      error = r8_abs ( exact - result );

      std::cout << "  " << "    "
           << "  " << std::setw(4) << int_num
           << "  " << std::setw(14) << result
           << "  " << std::setw(14) << error << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 applies a Halton sequence rule to finite interval 1D problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  double error;
  double exact;
  int int_log;
  int int_num;
  int prob;
  int prob_num;
  double result;

  std::cout << "\n";
  std::cout << "TEST06\n";
  std::cout << "  Halton sequence rule,\n";
  std::cout << "  for 1D finite interval problems.\n";

  prob_num = p00_prob_num ( );

  std::cout << "\n";
  std::cout << "  Prob   Ints         Exact       Error\n";
  std::cout << "                      Approx\n";
//
//  Pick a problem.
//
  for ( prob = 1; prob <= prob_num; prob++ )
  {
    exact = p00_exact ( prob );

    std::cout << "\n";
    std::cout << "  " << std::setw(4) << prob
         << "  " << "    "
         << "  " << std::setw(14) << exact << "\n";
//
//  Pick a number of points.
//
    for ( int_log = 0; int_log <= 10; int_log++ )
    {
      int_num = i4_power ( 2, int_log );

      result = p00_halton ( prob, int_num );

      error = r8_abs ( exact - result );

      std::cout << "  " << "    "
           << "  " << std::setw(4) << int_num
           << "  " << std::setw(14) << result
           << "  " << std::setw(14) << error << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 applies an evenly spaced point rule to finite interval 1D problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  double error;
  double exact;
  int int_log;
  int int_num;
  int prob;
  int prob_num;
  double result;

  std::cout << "\n";
  std::cout << "TEST07\n";
  std::cout << "  Evenly spaced point sequence rule,\n";
  std::cout << "  for 1D finite interval problems.\n";

  prob_num = p00_prob_num ( );

  std::cout << "\n";
  std::cout << "  Prob   Ints         Exact       Error\n";
  std::cout << "                      Approx\n";
//
//  Pick a problem.
//
  for ( prob = 1; prob <= prob_num; prob++ )
  {
    exact = p00_exact ( prob );

    std::cout << "\n";
    std::cout << "  " << std::setw(4) << prob
         << "  " << "    "
         << "  " << std::setw(14) << exact << "\n";
//
//  Pick a number of points.
//
    for ( int_log = 0; int_log <= 10; int_log++ )
    {
      int_num = i4_power ( 2, int_log );

      result = p00_even ( prob, int_num );

      error = r8_abs ( exact - result );

      std::cout << "  " << "    "
           << "  " << std::setw(4) << int_num
           << "  " << std::setw(14) << result
           << "  " << std::setw(14) << error << "\n";
    }
  }
  return;
}
