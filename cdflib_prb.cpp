# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

# include "cdflib.hpp"

int main ( );
void test005 ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );

void test10 ( );
void test11 ( );
void test12 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test17 ( );
void test18 ( );
void test19 ( );

void test20 ( );
void test21 ( );
void test22 ( );
void test23 ( );
void test24 ( );
void test25 ( );
void test26 ( );
void test27 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CDFLIB_PRB.
//
//  Discussion:
//
//    CDFLIB_PRB tests the CDFLIB library.
//
//  Modified:
//
//    17 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  std::cout << "\n";
  std::cout << "CDFLIB_PRB\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the CDFLIB library.\n";

  test005 ( );
  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test17 ( );
  test18 ( );
  test19 ( );

  test20 ( );
  test21 ( );
  test22 ( );
  test23 ( );
  test24 ( );
  test25 ( );
  test26 ( );
  test27 ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "CDFLIB_PRB\n";
  std::cout << "  Normal end of execution.\n";
  std::cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test005 ( )

//****************************************************************************80
//
//  Purpose:
//
//   TEST005 tests BETA_INC and BETA_INC_VALUES.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int ierror;
  int n_data;
  double x;
  double y;

  std::cout << "\n";
  std::cout << "TEST005\n";
  std::cout << "  BETA_INC computes the incomplete Beta ratio.\n";
  std::cout << "  BETA_INC_CDF_VALUES looks up some values.\n";
  std::cout << "\n";
  std::cout << "    X         Y         A         B         CDF           CDF\n";
  std::cout <<
    "                                           (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( &n_data, &a, &b, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    y = 1.0 - x;

    beta_inc ( &a, &b, &x, &y, &cdf_compute, &ccdf_compute, &ierror );

    std::cout << std::setw(10) << x
         << std::setw(10) << y
         << std::setw(10) << a
         << std::setw(10) << b
         << std::setw(14) << cdf_lookup
         << std::setw(14) << cdf_compute << "\n";
  }

  std::cout << "\n";
  std::cout << "    X         Y         A         B         1-CDF         CCDF\n";
  std::cout <<
    "                                           (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( &n_data, &a, &b, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    y = 1.0 - x;

    beta_inc ( &a, &b, &x, &y, &cdf_compute, &ccdf_compute, &ierror );

    std::cout << std::setw(10) << x
         << std::setw(10) << y
         << std::setw(10) << a
         << std::setw(10) << b
         << std::setw(14) << ccdf_lookup
         << std::setw(14) << ccdf_compute << "\n";

  }

  return;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests CDFBET.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double bound;
  double p;
  double q;
  int status;
  int which;
  double x;
  double y;

  std::cout << "\n";
  std::cout << "TEST01\n";
  std::cout << "  CDFBET computes one missing parameter from the\n";
  std::cout << "    BETA CDF:\n";
  std::cout << "\n";
  std::cout << "   BETA_CDF ( (P,Q), (X,Y), A, B )\n";
  std::cout << "\n";
  std::cout << "      P           Q               X           Y"
       << "            A           B\n";
  std::cout << "\n";

  for ( which = 1; which <= 4; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      x = 0.25;
      y = 1.0 - x;
      a = 2.0;
      b = 3.0;
    }
    else if ( which == 2 )
    {
      p = 0.261719;
      q = 1.0 - p;
      x = -1.0;
      y = -1.0;
      a = 2.0;
      b = 3.0;
    }
    else if ( which == 3 )
    {
      p = 0.261719;
      q = 1.0 - p;
      x = 0.25;
      y = 1.0 - x;
      a = -1.0;
      b = 3.0;
    }
    else if ( which == 4 )
    {
      p = 0.261719;
      q = 1.0 - p;
      x = 0.25;
      y = 1.0 - x;
      a = 2.0;
      b = -1.0;
    }

    cdfbet ( &which, &p, &q, &x, &y, &a, &b, &status, &bound );

    if ( status != 0 )
    {
      std::cout << "\n";
      std::cout << "  CDFBET returned STATUS = " << status << "\n";
      continue;
    }
    std::cout                  << "  "
         << std::setw(10) << p << "  "
         << std::setw(10) << q << "  "
         << std::setw(10) << x << "  "
         << std::setw(10) << y << "  "
         << std::setw(10) << a << "  "
         << std::setw(10) << b << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests CDFBIN.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bound;
  double ompr;
  double p;
  double pr;
  double q;
  double s;
  int status;
  int which;
  double xn;

  std::cout << "\n";
  std::cout << "TEST02\n";
  std::cout << "  CDFBIN computes one missing parameter from the\n";
  std::cout << "    Binomial CDF:\n";
  std::cout << "\n";
  std::cout << "   BINOMIAL_CDF ( (P,Q), S, XN, (PR,OMPR) )\n";
  std::cout << "\n";
  std::cout << "      P           Q                S          "
       << "XN         PR         OMPR\n";
  std::cout << "\n";

  for ( which = 1; which <= 4; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      s = 5.0;
      xn = 8.0;
      pr = 0.875;
      ompr = 1.0 - pr;
    }
    else if ( which == 2 )
    {
      p = 0.067347;
      q = 1.0 - p;
      s = -1.0;
      xn = 8.0;
      pr = 0.875;
      ompr = 1.0 - pr;
    }
    else if ( which == 3 )
    {
      p = 0.067347;
      q = 1.0 - p;
      s = 5.0;
      xn = -1.0;
      pr = 0.875;
      ompr = 1.0 - pr;
    }
    else if ( which == 4 )
    {
      p = 0.067347;
      q = 1.0 - p;
      s = 5.0;
      xn = 8.0;
      pr = -1.0;
      ompr = -1.0;
    }

    cdfbin ( &which, &p, &q, &s, &xn, &pr, &ompr, &status, &bound );

    if ( status != 0 )
    {
      std::cout << "\n";
      std::cout << "  CDFBIN returned STATUS = " << status << "\n";
      continue;
    }
    std::cout                     << "  "
         << std::setw(10) << p    << "  "
         << std::setw(10) << q    << "  "
         << std::setw(10) << s    << "  "
         << std::setw(10) << xn   << "  "
         << std::setw(10) << pr   << "  "
         << std::setw(10) << ompr << "\n";
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests CDFCHI.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bound;
  double df;
  double p;
  double q;
  int status;
  int which;
  double x;

  std::cout << "\n";
  std::cout << "TEST03\n";
  std::cout << "  CDFCHI computes one missing parameter from the\n";
  std::cout << "    Chi Square CDF:\n";
  std::cout << "\n";
  std::cout << "   CHI_CDF ( (P,Q), X, DF )\n";
  std::cout << "\n";
  std::cout << "      P           Q                X          DF\n";
  std::cout << "\n";

  for ( which = 1; which <= 3; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      x = 5.0;
      df = 8.0;
    }
    else if ( which == 2 )
    {
      p = 0.242424;
      q = 1.0 - p;
      x = -1.0;
      df = 8.0;
    }
    else if ( which == 3 )
    {
      p = 0.242424;
      q = 1.0 - p;
      x = 5.0;
      df = -1.0;
    }

    cdfchi ( &which, &p, &q, &x, &df, &status, &bound );

    if ( status != 0 )
    {
      std::cout << "\n";
      std::cout << "  CDFCHI returned STATUS = " << status << "\n";
      continue;
    }
    std::cout                     << "  "
         << std::setw(10) << p    << "  "
         << std::setw(10) << q    << "  "
         << std::setw(10) << x    << "  "
         << std::setw(10) << df   << "\n";
  }
  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests CDFCHN.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bound;
  double df;
  double p;
  double pnonc;
  double q;
  int status;
  int which;
  double x;

  std::cout << "\n";
  std::cout << "TEST04\n";
  std::cout << "  CDFCHN computes one missing parameter from the\n";
  std::cout << "    Chi Square CDF:\n";
  std::cout << "\n";
  std::cout << "   CHI_Noncentral_CDF ( (P,Q), X, DF, PNONC )\n";
  std::cout << "\n";
  std::cout << "     P         Q             X        DF     PNONC\n";
  std::cout << "\n";

  for ( which = 1; which <= 4; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      x = 5.0;
      df = 8.0;
      pnonc = 0.5;
    }
    else if ( which == 2 )
    {
      p = 0.211040;
      q = 1.0 - p;
      x = -1.0;
      df = 8.0;
      pnonc = 0.5;
    }
    else if ( which == 3 )
    {
      p = 0.211040;
      q = 1.0 - p;
      x = 5.0;
      df = -1.0;
      pnonc = 0.5;
    }
    else if ( which == 4 )
    {
      p = 0.211040;
      q = 1.0 - p;
      x = 5.0;
      df = 8.0;
      pnonc = -1.0;
    }

    cdfchn ( &which, &p, &q, &x, &df, &pnonc, &status, &bound );

    if ( status != 0 )
    {
      std::cout << "\n";
      std::cout << "  CDFCHN returned STATUS = " << status << "\n";
      continue;
    }

    std::cout <<                     "  "
         << std::setw(8) << p     << "  "
         << std::setw(8) << q     << "  "
         << std::setw(8) << x     << "  "
         << std::setw(8) << df    << "  "
         << std::setw(8) << pnonc << "\n";
  }
  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests CDFF.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bound;
  double dfd;
  double dfn;
  double f;
  double p;
  double q;
  int status;
  int which;

  std::cout << "\n";
  std::cout << "TEST05\n";
  std::cout << "  CDFF computes one missing parameter from the\n";
  std::cout << "    F CDF:\n";
  std::cout << "\n";
  std::cout << "   F_CDF ( (P,Q), F, DFN, DFD )\n";
  std::cout << "\n";
  std::cout << "     P         Q             F       DFN       DFD\n";
  std::cout << "\n";

  for ( which = 1; which <= 4; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      f = 5.0;
      dfn = 8.0;
      dfd = 3.0;
    }
    else if ( which == 2 )
    {
      p = 0.893510;
      q = 1.0 - p;
      f = -1.0;
      dfn = 8.0;
      dfd = 3.0;
    }
    else if ( which == 3 )
    {
      p = 0.893510;
      q = 1.0 - p;
      f = 5.0;
      dfn = -1.0;
      dfd = 3.0;
    }
    else if ( which == 4 )
    {
      p = 0.893510;
      q = 1.0 - p;
      f = 5.0;
      dfn = 8.0;
      dfd = -1.0;
    }

    cdff ( &which, &p, &q, &f, &dfn, &dfd, &status, &bound );

    if ( status != 0 )
    {
      std::cout << "\n";
      std::cout << "  CDFF returned STATUS = " << status << "\n";
      continue;
    }

    std::cout <<                   "  "
         << std::setw(8) << p   << "  "
         << std::setw(8) << q   << "  "
         << std::setw(8) << f   << "  "
         << std::setw(8) << dfn << "  "
         << std::setw(8) << dfd << "\n";
  }
  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests CDFFNC.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bound;
  double dfd;
  double dfn;
  double f;
  double p;
  double pnonc;
  double q;
  int status;
  int which;

  std::cout << "\n";
  std::cout << "TEST06\n";
  std::cout << "  CDFFNC computes one missing parameter from the\n";
  std::cout << "    noncentral F CDF:\n";
  std::cout << "\n";
  std::cout << "   F_noncentral_CDF ( (P,Q), F, DFN, DFD, PNONC )\n";
  std::cout << "\n";
  std::cout << "         P         Q         F       DFN       DFD     PNONC\n";
  std::cout << "\n";

  for ( which = 1; which <= 5; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      f = 5.0;
      dfn = 8.0;
      dfd = 3.0;
      pnonc = 17.648016;
    }
    else if ( which == 2 )
    {
      p = 0.60;
      q = 1.0 - p;
      f = -1.0;
      dfn = 8.0;
      dfd = 3.0;
      pnonc = 17.648016;
    }
    else if ( which == 3 )
    {
      p = 0.60;
      q = 1.0 - p;
      f = 5.0;
      dfn = -1.0;
      dfd = 3.0;
      pnonc = 17.648016;
    }
    else if ( which == 4 )
    {
      p = 0.60;
      q = 1.0 - p;
      f = 5.0;
      dfn = 8.0;
      dfd = -1.0;
      pnonc = 17.648016;
    }
    else if ( which == 5 )
    {
      p = 0.60;
      q = 1.0 - p;
      f = 5.0;
      dfn = 8.0;
      dfd = 3.0;
      pnonc = -1.0;
    }

    cdffnc ( &which, &p, &q, &f, &dfn, &dfd, &pnonc, &status, &bound );

    if ( status != 0 )
    {
      std::cout << "\n";
      std::cout << "  CDFFNC returned STATUS = " << status << "\n";
      continue;
    }

    std::cout <<                   "  "
         << std::setw(8) << p     << "  "
         << std::setw(8) << q     << "  "
         << std::setw(8) << f     << "  "
         << std::setw(8) << dfn   << "  "
         << std::setw(8) << dfd   << "  "
         << std::setw(8) << pnonc << "\n";
  }

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests CDFGAM.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bound;
  double p;
  double q;
  double scale;
  double shape;
  int status;
  int which;
  double x;

  std::cout << "\n";
  std::cout << "TEST07\n";
  std::cout << "  CDFGAM computes one missing parameter from the\n";
  std::cout << "    Gamma CDF:\n";
  std::cout << "\n";
  std::cout << "   Gamma_CDF ( (P,Q), X, SHAPE, SCALE )\n";
  std::cout << "\n";
  std::cout << "    P         Q              X     SHAPE     SCALE\n";
  std::cout << "\n";

  for ( which = 1; which <= 4; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      x = 5.0;
      shape = 8.0;
      scale = 3.0;
    }
    else if ( which == 2 )
    {
      p = 0.981998;
      q = 1.0 - p;
      x = -1.0;
      shape = 8.0;
      scale = 3.0;
    }
    else if ( which == 3 )
    {
      p = 0.981998;
      q = 1.0 - p;
      x = 5.0;
      shape = -1.0;
      scale = 3.0;
    }
    else if ( which == 4 )
    {
      p = 0.981998;
      q = 1.0 - p;
      x = 5.0;
      shape = 8.0;
      scale = -1.0;
    }

    cdfgam ( &which, &p, &q, &x, &shape, &scale, &status, &bound );

    if ( status != 0 )
    {
      std::cout << "\n";
      std::cout << "  CDFGAM returned STATUS = " << status << "\n";
      continue;
    }

    std::cout <<                     "  "
         << std::setw(8) << p     << "  "
         << std::setw(9) << q     << "  "
         << std::setw(8) << x     << "  "
         << std::setw(8) << shape << "  "
         << std::setw(8) << scale << "\n";
  }

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests CDFNBN.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bound;
  double f;
  double ompr;
  double p;
  double pr;
  double q;
  double s;
  int status;
  int which;

  std::cout << "\n";
  std::cout << "TEST08\n";
  std::cout << "  CDFNBN computes one missing parameter from the\n";
  std::cout << "    Negative_Binomial CDF:\n";
  std::cout << "\n";
  std::cout << "   Negative_BINOMIAL_CDF ( (P,Q), F, S, (PR,OMPR) )\n";
  std::cout << "\n";
  std::cout << "    P         Q               F         S       PR        OMPR\n";
  std::cout << "\n";

  for ( which = 1; which <= 4; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      f = 3.0;
      s = 5.0;
      pr = 0.875;
      ompr = 1.0 - pr;
    }
    else if ( which == 2 )
    {
      p = 0.988752;
      q = 1.0 - p;
      f = -1.0;
      s = 5.0;
      pr = 0.875;
      ompr = 1.0 - pr;
    }
    else if ( which == 3 )
    {
      p = 0.988752;
      q = 1.0 - p;
      f = 3.0;
      s = -1.0;
      pr = 0.875;
      ompr = 1.0 - pr;
    }
    else if ( which == 4 )
    {
      p = 0.988752;
      q = 1.0 - p;
      f = 3.0;
      s = 5.0;
      pr = -1.0;
      ompr = -1.0;
    }

    cdfnbn ( &which, &p, &q, &f, &s, &pr, &ompr, &status, &bound );

    if ( status != 0 )
    {
      std::cout << "\n";
      std::cout << "  CDFNBN returned STATUS = " << status << "\n";
      continue;
    }
    std::cout <<                    "  "
         << std::setw(8) << p    << "  "
         << std::setw(9) << q    << "  "
         << std::setw(8) << f    << "  "
         << std::setw(8) << s    << "  "
         << std::setw(8) << pr   << "  "
         << std::setw(8) << ompr << "\n";
  }

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests CDFNOR.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bound;
  double mean;
  double p;
  double q;
  double sd;
  int status;
  int which;
  double x;

  std::cout << "\n";
  std::cout << "TEST09\n";
  std::cout << "  CDFNOR computes one missing parameter from the\n";
  std::cout << "    Normal CDF:\n";
  std::cout << "\n";
  std::cout << "   Normal_CDF ( (P,Q), X, MEAN, SD )\n";
  std::cout << "\n";
  std::cout << "    P         Q               X      MEAN       SD\n";
  std::cout << "\n";

  for ( which = 1; which <= 4; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      x = 3.0;
      mean = 5.0;
      sd = 0.875;
    }
    else if ( which == 2 )
    {
      p = 0.011135;
      q = 1.0 - p;
      x = -1.0;
      mean = 5.0;
      sd = 0.875;
    }
    else if ( which == 3 )
    {
      p = 0.011135;
      q = 1.0 - p;
      x = 3.0;
      mean = -1.0;
      sd = 0.875;
    }
    else if ( which == 4 )
    {
      p = 0.011135;
      q = 1.0 - p;
      x = 3.0;
      mean = 5.0;
      sd = -1.0;
    }

    cdfnor ( &which, &p, &q, &x, &mean, &sd, &status, &bound );

    if ( status != 0 )
    {
      std::cout << "\n";
      std::cout << "  CDFNOR returned STATUS = " << status << "\n";
      continue;
    }
    std::cout <<                    "  "
         << std::setw(9) << p    << "  "
         << std::setw(8) << q    << "  "
         << std::setw(8) << x    << "  "
         << std::setw(8) << mean << "  "
         << std::setw(8) << sd   << "\n";
  }

  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests CDFPOI.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bound;
  double p;
  double q;
  double s;
  int status;
  int which;
  double xlam;

  std::cout << "\n";
  std::cout << "TEST10\n";
  std::cout << "  CDFPOI computes one missing parameter from the\n";
  std::cout << "    Poisson CDF:\n";
  std::cout << "\n";
  std::cout << "   POISSON_CDF ( (P,Q), S, XLAM )\n";
  std::cout << "\n";
  std::cout << "     P         Q         S         XLAM\n";
  std::cout << "\n";

  for ( which = 1; which <= 3; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      s = 3.0;
      xlam = 5.0;
    }
    else if ( which == 2 )
    {
      p = 0.265026;
      q = 1.0 - p;
      s = -1.0;
      xlam = 5.0;
    }
    else if ( which == 3 )
    {
      p = 0.265026;
      q = 1.0 - p;
      s = 3.0;
      xlam = -1.0;
    }

    cdfpoi ( &which, &p, &q, &s, &xlam, &status, &bound );

    if ( status != 0 )
    {
      std::cout << "\n";
      std::cout << "  CDFPOI returned STATUS = " << status << "\n";
      continue;
    }
    std::cout <<                    "  "
         << std::setw(9) << p    << "  "
         << std::setw(9) << q    << "  "
         << std::setw(9) << s    << "  "
         << std::setw(9) << xlam << "\n";
  }

  return;
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests CDFT.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bound;
  double df;
  double p;
  double q;
  int status;
  double t;
  int which;

  std::cout << "\n";
  std::cout << "TEST11\n";
  std::cout << "  CDFT computes one missing parameter from the\n";
  std::cout << "    T CDF:\n";
  std::cout << "\n";
  std::cout << "   T_CDF ( (P,Q), T, DF )\n";
  std::cout << "\n";
  std::cout << "    P         Q         T         DF\n";
  std::cout << "\n";

  for ( which = 1; which <= 3; which++ )
  {
    if ( which == 1 )
    {
      p = -1.0;
      q = -1.0;
      t = 3.0;
      df = 5.0;
    }
    else if ( which == 2 )
    {
      p = 0.984950;
      q = 1.0 - p;
      t = -1.0;
      df = 5.0;
    }
    else if ( which == 3 )
    {
      p = 0.984950;
      q = 1.0 - p;
      t = 3.0;
      df = -1.0;
    }

    cdft ( &which, &p, &q, &t, &df, &status, &bound );

    if ( status != 0 )
    {
      std::cout << "\n";
      std::cout << "  CDFT returned STATUS = " << status << "\n";
      continue;
    }
    std::cout <<                  "  "
         << std::setw(9) << p  << "  "
         << std::setw(9) << q  << "  "
         << std::setw(9) << t  << "  "
         << std::setw(9) << df << "\n";
  }

  return;
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests CUMBET, BETA_INC_VALUES.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int n_data;
  double x;
  double y;

  std::cout << "\n";
  std::cout << "TEST12\n";
  std::cout << "  CUMBET computes the Beta CDF\n";
  std::cout << "    and the complementary CDF.\n";
  std::cout << "  BETA_INC_CDF_VALUES looks up some values.\n";
  std::cout << "\n";
  std::cout << "    X         Y         A         B         CDF           CDF\n";
  std::cout <<
    "                                           (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( &n_data, &a, &b, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    y = 1.0 - x;

    cumbet ( &x, &y, &a, &b, &cdf_compute, &ccdf_compute );

    std::cout << " "
         << std::setw(9) << x           << "  "
         << std::setw(9) << y           << "  "
         << std::setw(9) << a           << "  "
         << std::setw(9) << b           << "  "
         << std::setw(9) << cdf_lookup  << "  "
         << std::setw(9) << cdf_compute << "\n";
  }

  std::cout << "\n";
  std::cout << "    X         Y         A         B         1-CDF         CCDF\n";
  std::cout <<
    "                                           (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( &n_data, &a, &b, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    y = 1.0 - x;

    cumbet ( &x, &y, &a, &b, &cdf_compute, &ccdf_compute );

    std::cout << "  "
         << std::setw(9) << x            << "  "
         << std::setw(9) << y            << "  "
         << std::setw(9) << a            << "  "
         << std::setw(9) << b            << "  "
         << std::setw(9) << ccdf_lookup  << "  "
         << std::setw(9) << ccdf_compute << "\n";
  }

  return;
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests CUMBIN, BINOMIAL_CDF_VALUES.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int n_data;
  double ompr;
  int s;
  double s_double;
  double pr;
  int x;
  double x_double;

  std::cout << "\n";
  std::cout << "TEST13\n";
  std::cout << "  CUMBIN computes the Binomial CDF\n";
  std::cout << "    and the complementary CDF.\n";
  std::cout << "  BINOMIAL_CDF_VALUES looks up some values.\n";
  std::cout << "\n";
  std::cout << "   X   S    Pr       CDF           CDF\n";
  std::cout << "                    (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    binomial_cdf_values ( &n_data, &x, &pr, &s, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ompr = 1.0 - pr;

    s_double = ( double ) s;
    x_double = ( double ) x;

    cumbin ( &s_double, &x_double, &pr, &ompr, &cdf_compute, &ccdf_compute );

    std::cout << "  "
         << std::setw(2)  << s           << "  "
         << std::setw(2)  << x           << "  "
         << std::setw(8)  << pr          << "  "
         << std::setw(12) << cdf_lookup  << "  "
         << std::setw(12) << cdf_compute << "\n";
  }

  std::cout << "\n";
  std::cout << "   X   S    Pr       1-CDF         CCDF\n";
  std::cout << "                    (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    binomial_cdf_values ( &n_data, &x, &pr, &s, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    ompr = 1.0 - pr;

    s_double = ( double ) s;
    x_double = ( double ) x;

    cumbin ( &s_double, &x_double, &pr, &ompr, &cdf_compute, &ccdf_compute );

    std::cout << "  "
         << std::setw(2)  << s            << "  "
         << std::setw(2)  << x            << "  "
         << std::setw(8)  << pr           << "  "
         << std::setw(12) << ccdf_lookup  << "  "
         << std::setw(12) << ccdf_compute << "\n";
  }

  return;
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests CUMCHI, CHI_SQUARE_CDF_VALUES.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int df;
  double df_double;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TEST14\n";
  std::cout << "  CUMCHI computes the chi square CDF\n";
  std::cout << "    and the complementary CDF.\n";
  std::cout << "  CHI_SQUARE_CDF_VALUES looks up some values.\n";
  std::cout << "\n";
  std::cout << "    X       DF    CDF           CDF\n";
  std::cout << "                 (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    chi_square_cdf_values ( &n_data, &df, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    df_double = ( double ) df;

    cumchi ( &x, &df_double, &cdf_compute, &ccdf_compute );

    std::cout << "  "
         << std::setw(8)  << x           << "  "
         << std::setw(2)  << df          << "  "
         << std::setw(12) << cdf_lookup  << "  "
         << std::setw(12) << cdf_compute << "\n";
  }

  std::cout << "\n";
  std::cout << "    X       DF    1-CDF         CCDF\n";
  std::cout << "                 (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    chi_square_cdf_values ( &n_data, &df, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    df_double = ( double ) df;

    cumchi ( &x, &df_double, &cdf_compute, &ccdf_compute );

    std::cout << "  "
         << std::setw(8)  << x            << "  "
         << std::setw(2)  << df           << "  "
         << std::setw(12) << ccdf_lookup  << "  "
         << std::setw(12) << ccdf_compute << "\n";
  }

  return;
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests CUMCHN, CHI_NONCENTRAL_CDF_VALUES.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int df;
  double df_double;
  double lambda;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TEST15\n";
  std::cout << "  CUMCHN computes the cumulative density\n";
  std::cout << "    function for the noncentral chi-squared\n";
  std::cout << "    distribution.\n";
  std::cout << "  CHI_NONCENTRAL_CDF_VALUES looks up some values.\n";
  std::cout << "\n";
  std::cout << "    DF    Lambda    X         CDF           CDF\n";
  std::cout << "                             (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    chi_noncentral_cdf_values ( &n_data, &x, &lambda, &df, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    df_double = ( double ) df;

    cumchn ( &x, &df_double, &lambda, &cdf_compute, &ccdf_compute );

    std::cout << "  "
         << std::setw(6)  << df          << "  "
         << std::setw(8)  << lambda      << "  "
         << std::setw(8)  << x           << "  "
         << std::setw(12) << cdf_lookup  << "  "
         << std::setw(12) << cdf_compute << "\n";
  }

  std::cout << "\n";
  std::cout << "    DF    Lambda    X         1-CDF         CCDF\n";
  std::cout << "                             (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    chi_noncentral_cdf_values ( &n_data, &x, &lambda, &df, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    df_double = ( double ) df;

    cumchn ( &x, &df_double, &lambda, &cdf_compute, &ccdf_compute );

    std::cout << "  "
         << std::setw(6)  << df           << "  "
         << std::setw(8)  << lambda       << "  "
         << std::setw(8)  << x            << "  "
         << std::setw(12) << ccdf_lookup  << "  "
         << std::setw(12) << ccdf_compute << "\n";
  }

  return;
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 tests CUMF, F_CDF_VALUES.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int dfd;
  double dfd_double;
  int dfn;
  double dfn_double;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TEST16\n";
  std::cout << "  CUMF computes the F CDF\n";
  std::cout << "    and the complementary CDF.\n";
  std::cout << "  F_CDF_VALUES looks up some values.\n";
  std::cout << "\n";
  std::cout << "    X      DFN DFD    CDF           CDF\n";
  std::cout << "                     (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    f_cdf_values ( &n_data, &dfn, &dfd, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    dfn_double = ( double ) dfn;
    dfd_double = ( double ) dfd;

    cumf ( &x, &dfn_double, &dfd_double, &cdf_compute, &ccdf_compute );

    std::cout << "  "
         << std::setw(8)  << x           << "  "
         << std::setw(2)  << dfn         << "  "
         << std::setw(2)  << dfd         << "  "
         << std::setw(12) << cdf_lookup  << "  "
         << std::setw(12) << cdf_compute << "\n";

  }

  std::cout << "\n";
  std::cout << "    X      DFN DFD    1-CDF         CCDF\n";
  std::cout << "                     (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    f_cdf_values ( &n_data, &dfn, &dfd, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    dfn_double = ( double ) dfn;
    dfd_double = ( double ) dfd;

    cumf ( &x, &dfn_double, &dfd_double, &cdf_compute, &ccdf_compute );

    std::cout << "  "
         << std::setw(8)  << x            << "  "
         << std::setw(2)  << dfn          << "  "
         << std::setw(2)  << dfd          << "  "
         << std::setw(12) << ccdf_lookup  << "  "
         << std::setw(12) << ccdf_compute << "\n";
  }

  return;
}
//****************************************************************************80

void test17 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17 tests CUMFNC, F_NONCENTRAL_CDF_VALUES.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int dfd;
  double dfd_double;
  int dfn;
  double dfn_double;
  double lambda;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TEST17\n";
  std::cout << "  CUMFNC computes the noncentral F CDF\n";
  std::cout << "    and the complementary CDF.\n";
  std::cout << "  F_NONCENTRAL_CDF_VALUES looks up some values.\n";
  std::cout << "\n";
  std::cout << "    X      DFN DFD    LAMBDA    CDF           CDF\n";
  std::cout << "                               (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    f_noncentral_cdf_values ( &n_data, &dfn, &dfd, &lambda, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    dfn_double = ( double ) dfn;
    dfd_double = ( double ) dfd;

    cumfnc ( &x, &dfn_double, &dfd_double, &lambda, &cdf_compute,
      &ccdf_compute );

    std::cout << "  "
         << std::setw(8)  << x           << "  "
         << std::setw(2)  << dfn         << "  "
         << std::setw(2)  << dfd         << "  "
         << std::setw(8)  << lambda      << "  "
         << std::setw(12) << cdf_lookup  << "  "
         << std::setw(12) << cdf_compute << "\n";
  }

  std::cout << "\n";
  std::cout << "    X      DFN DFD    LAMBDA    1-CDF         CCDF\n";
  std::cout << "                               (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    f_noncentral_cdf_values ( &n_data, &dfn, &dfd, &lambda, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    dfn_double = ( double ) dfn;
    dfd_double = ( double ) dfd;

    cumfnc ( &x, &dfn_double, &dfd_double, &lambda, &cdf_compute,
      &ccdf_compute );

    std::cout << "  "
         << std::setw(8)  << x            << "  "
         << std::setw(2)  << dfn          << "  "
         << std::setw(2)  << dfd          << "  "
         << std::setw(8)  << lambda       << "  "
         << std::setw(12) << ccdf_lookup  << "  "
         << std::setw(12) << ccdf_compute << "\n";
  }

  return;
}
//****************************************************************************80

void test18 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST18 tests CUMGAM, GAMMA_INC_VALUES.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TEST18\n";
  std::cout << "  CUMGAM computes the Gamma CDF\n";
  std::cout << "    and the complementary CDF.\n";
  std::cout << "  GAMMA_INC_VALUES looks up some values.\n";
  std::cout << "\n";
  std::cout << "    A         X         CDF           CDF\n";
  std::cout << "                        (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_values ( &n_data, &a, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    cumgam ( &x, &a, &cdf_compute, &ccdf_compute );

    std::cout << "  "
         << std::setw(8)  << a           << "  "
         << std::setw(8)  << x           << "  "
         << std::setw(12) << cdf_lookup  << "  "
         << std::setw(12) << cdf_compute << "\n";
  }

  std::cout << "\n";
  std::cout << "    A         X         CDF           CDF\n";
  std::cout << "                        (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_values ( &n_data, &a, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    cumgam ( &x, &a, &cdf_compute, &ccdf_compute );

    std::cout << "  "
         << std::setw(8)  << a            << "  "
         << std::setw(8)  << x            << "  "
         << std::setw(12) << ccdf_lookup  << "  "
         << std::setw(12) << ccdf_compute << "\n";
  }

  return;
}
//****************************************************************************80

void test19 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST19 tests CUMNBN, NEGATIVE_BINOMIAL_CDF_VALUES.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int f;
  double f_double;
  int n_data;
  double ompr;
  int s;
  double s_double;
  double pr;

  std::cout << "\n";
  std::cout << "TEST19\n";
  std::cout << "  CUMNBN computes the Negative Binomial CDF\n";
  std::cout << "    and the complementary CDF.\n";
  std::cout << "  NEGATIVE_BINOMIAL_CDF_VALUES looks up some values.\n";
  std::cout << "\n";
  std::cout << "   F   S    Pr       CDF           CDF\n";
  std::cout << "                     (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    negative_binomial_cdf_values ( &n_data, &f, &s, &pr, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ompr = 1.0 - pr;

    f_double = ( double ) f;
    s_double = ( double ) s;

    cumnbn ( &f_double, &s_double, &pr, &ompr, &cdf_compute, &ccdf_compute );

    std::cout << "  "
         << std::setw(2)  << f           << "  "
         << std::setw(2)  << s           << "  "
         << std::setw(8)  << pr          << "  "
         << std::setw(12) << cdf_lookup  << "  "
         << std::setw(12) << cdf_compute << "\n";
  }

  std::cout << "\n";
  std::cout << "   F   S    Pr       1-CDF         CCDF\n";
  std::cout << "                     (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    negative_binomial_cdf_values ( &n_data, &f, &s, &pr, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    ompr = 1.0 - pr;

    f_double = ( double ) f;
    s_double = ( double ) s;

    cumnbn ( &f_double, &s_double, &pr, &ompr, &cdf_compute, &ccdf_compute );

    std::cout << "  "
         << std::setw(2)  << f            << "  "
         << std::setw(2)  << s            << "  "
         << std::setw(8)  << pr           << "  "
         << std::setw(12) << ccdf_lookup  << "  "
         << std::setw(12) << ccdf_compute << "\n";
  }

  return;
}
//****************************************************************************80

void test20 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20 tests CUMNOR, NORMAL_CDF_VALUES.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TEST20\n";
  std::cout << "  CUMNOR computes the Normal CDF\n";
  std::cout << "    and the complementary CDF.\n";
  std::cout << "  NORMAL_CDF_VALUES looks up some values.\n";
  std::cout << "\n";
  std::cout << "    X         CDF           CDF\n";
  std::cout << "              (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    normal_cdf_values ( &n_data, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    cumnor ( &x, &cdf_compute, &ccdf_compute );

    std::cout << "  "
         << std::setw(8)  << x           << "  "
         << std::setw(12) << cdf_lookup  << "  "
         << std::setw(12) << cdf_compute << "\n";
  }

  std::cout << "\n";
  std::cout << "    X         1-CDF         CCDF\n";
  std::cout << "              (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    normal_cdf_values ( &n_data, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    cumnor ( &x, &cdf_compute, &ccdf_compute );

    std::cout << "  "
         << std::setw(8)  << x            << "  "
         << std::setw(12) << ccdf_lookup  << "  "
         << std::setw(12) << ccdf_compute << "\n";
  }

  return;
}
//****************************************************************************80

void test21 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST21 tests CUMPOI, POISSON_CDF_VALUES.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double  cdf_lookup;
  double lambda;
  int n_data;
  int x;
  double x_double;

  std::cout << "\n";
  std::cout << "TEST21\n";
  std::cout << "  CUMPOI computes the Poisson CDF\n";
  std::cout << "    and the complementary CDF.\n";
  std::cout << "  POISSON_CDF_VALUES looks up some values.\n";
  std::cout << "\n";
  std::cout << "     X    LAMBDA    CDF           CDF\n";
  std::cout << "                   (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    poisson_cdf_values ( &n_data, &lambda, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    x_double = ( double ) x;

    cumpoi ( &x_double, &lambda, &cdf_compute, &ccdf_compute );

    std::cout << "  "
         << std::setw(4)  << x           << "  "
         << std::setw(8)  << lambda      << "  "
         << std::setw(12) << cdf_lookup  << "  "
         << std::setw(12) << cdf_compute << "\n";
  }

  std::cout << "\n";
  std::cout << "     X    LAMBDA    1-CDF         CCDF\n";
  std::cout << "                   (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    poisson_cdf_values ( &n_data, &lambda, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    x_double = ( double ) x;
    ccdf_lookup = 1.0 - cdf_lookup;

    cumpoi ( &x_double, &lambda, &cdf_compute, &ccdf_compute );

    std::cout << "  "
         << std::setw(4)  << x            << "  "
         << std::setw(8)  << lambda       << "  "
         << std::setw(12) << ccdf_lookup  << "  "
         << std::setw(12) << ccdf_compute << "\n";
  }

  return;
}
//****************************************************************************80

void test22 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST22 tests CUMT, STUDENT_CDF_VALUES.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double ccdf_compute;
  double ccdf_lookup;
  double cdf_compute;
  double cdf_lookup;
  int df;
  double df_double;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TEST22\n";
  std::cout << "  CUMT computes the Student T CDF\n";
  std::cout << "    and the complementary CDF.\n";
  std::cout << "  STUDENT_CDF_VALUES looks up some values.\n";
  std::cout << "\n";
  std::cout << "    X       DF    CDF           CDF\n";
  std::cout << "                 (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    student_cdf_values ( &n_data, &df, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }
    df_double = ( double ) df;

    cumt ( &x, &df_double, &cdf_compute, &ccdf_compute );

    std::cout << "  "
         << std::setw(8)  << x           << "  "
         << std::setw(2)  << df          << "  "
         << std::setw(12) << cdf_lookup  << "  "
         << std::setw(12) << cdf_compute << "\n";
  }

  std::cout << "\n";
  std::cout << "    X       DF    1-CDF         CCDF\n";
  std::cout << "                 (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    student_cdf_values ( &n_data, &df, &x, &cdf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    ccdf_lookup = 1.0 - cdf_lookup;

    df_double = ( double ) df;

    cumt ( &x, &df_double, &cdf_compute, &ccdf_compute );

    std::cout << "  "
         << std::setw(8)  << x            << "  "
         << std::setw(2)  << df           << "  "
         << std::setw(12) << ccdf_lookup  << "  "
         << std::setw(12) << ccdf_compute << "\n";
  }

  return;
}
//****************************************************************************80

void test23 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST23 tests BETA, GAMMA_X.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double apb;
  double beta1;
  double beta2;

  std::cout << "\n";
  std::cout << "TEST23\n";
  std::cout << "  BETA evaluates the Beta function;\n";
  std::cout << "  GAMMA_X evaluates the Gamma function.\n";

  a = 2.2;
  b = 3.7;
  apb = a + b;

  beta1 = beta ( a, b );
  beta2 = gamma_x ( &a ) * gamma_x ( &b ) / gamma_x ( &apb );

  std::cout << "\n";
  std::cout << "  Argument A =                   " << a << "\n";
  std::cout << "  Argument B =                   " << b << "\n";
  std::cout << "  Beta(A,B) =                    " << beta1 << "\n";
  std::cout << "  (Expected value = 0.0454 )\n";
  std::cout << "\n";
  std::cout << "  Gamma(A)*Gamma(B)/Gamma(A+B) = " << beta2 << "\n";

  return;
}
//****************************************************************************80

void test24 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST24 tests ERROR_F, ERROR_FC, ERF_VALUES..
//
//  Modified:
//
//    17 November 2006
//
//  Author:
//
//    John Burkardt
//
{
  double erf_compute;
  double erf_lookup;
  double erfc_compute;
  double erfc_lookup;
  int ind;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TEST24\n";
  std::cout << "  ERROR_F computes the error function ERF;\n";
  std::cout << "  ERROR_FC the complementary error function ERFC.\n";
  std::cout << "  ERF_VALUES looks up some values.\n";
  std::cout << "\n";
  std::cout << "    X         ERF           ERF\n";
  std::cout << "              (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    erf_values ( &n_data, &x, &erf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    erf_compute = error_f ( &x );

    std::cout << "  "
         << std::setw(8)  << x           << "  "
         << std::setw(12) << erf_lookup  << "  "
         << std::setw(12) << erf_compute << "\n";
  }

  std::cout << "\n";
  std::cout << "    X         ERFC          ERFC\n";
  std::cout << "              (Lookup)      (Computed)\n";
  std::cout << "\n";

  ind = 0;
  n_data = 0;

  for ( ; ; )
  {
    erf_values ( &n_data, &x, &erf_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    erfc_lookup = 1.0 - erf_lookup;
    erfc_compute = error_fc ( &ind, &x );

    std::cout << "  "
         << std::setw(8)  << x            << "  "
         << std::setw(12) << erfc_lookup  << "  "
         << std::setw(12) << erfc_compute << "\n";
  }

  return;
}
//****************************************************************************80

void test25 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST25 tests XGAMM, GAMMA_VALUES.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double gamma_compute;
  double gamma_lookup;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TEST25\n";
  std::cout << "  XGAMM computes the Gamma function;\n";
  std::cout << "  GAMMA_VALUES looks up some values.\n";
  std::cout << "\n";
  std::cout << "    X         GAMMA         GAMMA\n";
  std::cout << "              (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_values ( &n_data, &x, &gamma_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    gamma_compute = gamma_x ( &x );

    std::cout << "  "
         << std::setw(8)  << x             << "  "
         << std::setw(12) << gamma_lookup  << "  "
         << std::setw(12) << gamma_compute << "\n";
  }

  return;
}
//****************************************************************************80

void test26 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST26 tests GAMMA_INC, GAMMA_INC_INV.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  int i;
  int ierror;
  int ind;
  double p;
  double q;
  int test_num = 10;
  double x;
  double x0;
  double x2;

  a = 3.0;
  ind = 1;
  x0 = 0;

  std::cout << "\n";
  std::cout << "TEST26\n";
  std::cout << "  GAMMA_INC evaluates the incomplete Gamma ratio;\n";
  std::cout << "  GAMMA_INC_INV inverts it.\n";
  std::cout << "\n";
  std::cout << "  Parameters:\n";
  std::cout << "\n";
  std::cout << "    A = " << a << "\n";
  std::cout << "\n";
  std::cout << "    X             P             Q             Inverse\n";
  std::cout << "\n";

  for ( i = 0; i <= test_num; i++ )
  {
    x = ( double ) i / ( double ) test_num;

    gamma_inc ( &a, &x, &p, &q, &ind );

    gamma_inc_inv ( &a, &x2, &x0, &p, &q, &ierror );

    std::cout << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << p  << "  "
         << std::setw(12) << q  << "  "
         << std::setw(12) << x2 << "\n";
  }

  return;
}
//****************************************************************************80

void test27 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST27 tests PSI, PSI_VALUES.
//
//  Modified:
//
//    14 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double psi_compute;
  double psi_lookup;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TEST27\n";
  std::cout << "  PSI computes the Psi function;\n";
  std::cout << "  PSI_VALUES looks up some values.\n";
  std::cout << "\n";
  std::cout << "    X         PSI           PSI\n";
  std::cout << "              (Lookup)      (Computed)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    psi_values ( &n_data, &x, &psi_lookup );

    if ( n_data == 0 )
    {
      break;
    }

    psi_compute = psi ( &x );

    std::cout << "  "
         << std::setw(8)  << x           << "  "
         << std::setw(12) << psi_lookup  << "  "
         << std::setw(12) << psi_compute << "\n";
  }

  return;
}
