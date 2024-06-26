# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <iomanip>
# include <iostream>

using namespace std;

# include "blas0.hpp"

//****************************************************************************80

complex <float> c4_uniform_01 ( int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    C4_UNIFORM_01 returns a unit complex pseudorandom number.
//
//  Discussion:
//
//    The angle should be uniformly distributed between 0 and 2 * PI,
//    the square root of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <float> C4_UNIFORM_01, a pseudorandom complex value.
//
{
  float r;
  int k;
  float pi = 3.1415926;
  float theta;
  complex <float> value;

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + 2147483647;
  }

  r = sqrt ( static_cast<float>(static_cast<double>(seed) * 4.656612875E-10) );

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + 2147483647;
  }

  theta = 2.0 * pi * static_cast<float>(static_cast<double>(seed) * 4.656612875E-10);

  value = complex <float> ( r * cos ( theta ), r * sin ( theta ) );

  return value;
}
//****************************************************************************80

void c4mat_print ( int m, int n, complex <float> a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_PRINT prints a C4MAT.
//
//  Discussion:
//
//    A C4MAT is an array of complex <float> values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input, complex <float> A[M*N], the matrix.
//
//    Input, string TITLE, a title.
//
{
  c4mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void c4mat_print_some ( int m, int n, complex <float> a[], int ilo, int jlo, 
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_PRINT_SOME prints some of a C4MAT.
//
//  Discussion:
//
//    A C4MAT is an array of complex <float> values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input, complex <float> A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
  complex <float> c;
  int i;
  int i2hi;
  int i2lo;
  int inc;
  int incx = 4;
  int j;
  int j2;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }
//
//  Print the columns of the matrix, in strips of INCX.
//
  for ( j2lo = jlo; j2lo <= i4_min ( jhi, n ); j2lo = j2lo + incx )
  {
    j2hi = j2lo + incx - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    inc = j2hi + 1 - j2lo;

    cout << "\n";
    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      j2 = j + 1 - j2lo;
      cout << "     " << setw(10) << j << "     ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(5) << i << ":";
//
//  Print out (up to) INCX entries in row I, that lie in the current strip.
//
      for ( j2 = 1; j2 <= inc; j2++ )
      {
        j = j2lo - 1 + j2;
        c = a[i-1+(j-1)*m];
        cout << "  " << setw(8) << real ( c )
             << "  " << setw(8) << imag ( c );
      }
      cout << "\n";
    }
  }

  return;
}
//****************************************************************************80

complex <float> *c4mat_test ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_TEST returns a test matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Output, complex <float> C4MAT_TEST[N*N], the matrix.
//
{
  complex <float> *a;
  float angle;
  int i;
  complex <float> I = complex <float> ( 0.0, 1.0 );
  int j;
  float pi = 3.141592653589793;

  a = new complex <float>[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      angle = 2.0 * pi * static_cast<float>(i * j) / static_cast<float>(n);

      a[i+j*n] = exp ( I * angle ) / sqrt ( static_cast<float>(n) );
    }
  }
  return a;
}
//****************************************************************************80

complex <float> *c4mat_test_inverse ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_TEST_INVERSE returns the inverse of a test matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Output, complex <float> C4MAT_TEST_INVERSE[N*N], the matrix.
//
{
  complex <float> *a;
  int i;
  int j;
  complex <float> t;

  a = c4mat_test ( n );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      t        = conj ( a[i+j*n] );
      a[i+j*n] = conj ( a[j+i*n] );
      a[j+i*n] = t;
    }
  }
  return a;
}
//****************************************************************************80

complex <float> *c4mat_uniform_01_new ( int m, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    C4MAT_UNIFORM_01_NEW returns a new unit pseudorandom C4MAT.
//
//  Discussion:
//
//    The angles should be uniformly distributed between 0 and 2 * PI,
//    the square roots of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 June 2015
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <float> C4MAT_UNIFORM_01[M*N], the pseudorandom complex matrix.
//
{
  complex <float> *c;
  int i;
  const int i4_huge = 2147483647;
  int j;
  float r;
  int k;
  const float r4_pi = 3.1415926;
  float theta;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "C4MAT_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  c = new complex <float>[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + i4_huge;
      }

      r = sqrt ( static_cast<float>(seed) * 4.656612875E-10 );

      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + i4_huge;
      }

      theta = 2.0 * r4_pi * ( static_cast<float>(seed) * 4.656612875E-10 );

      c[i+j*m] = r * complex <float> ( cos ( theta ), sin ( theta ) );
    }
  }

  return c;
}
//****************************************************************************80

complex <float> *c4vec_uniform_01_new ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    C4VEC_UNIFORM_01_NEW returns a unit pseudorandom C4VEC.
//
//  Discussion:
//
//    The angles should be uniformly distributed between 0 and 2 * PI,
//    the square roots of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 June 2015
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int N, the number of values to compute.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <float> C4VEC_UNIFORM_01_NEW[N], the pseudorandom 
//    complex vector.
//
{
  complex <float> *c;
  int i;
  const int i4_huge = 2147483647;
  int k;
  const float r4_pi = 3.1415926;
  float r;
  float theta;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "C4VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  c = new complex <float> [n];

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r = sqrt ( static_cast<float>(seed) * 4.656612875E-10 );

    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    theta = 2.0 * r4_pi * ( static_cast<float>(seed) * 4.656612875E-10 );

    c[i] = r * complex <float> ( cos ( theta ), sin ( theta ) );
  }

  return c;
}
//****************************************************************************80

complex <double> c8_uniform_01 ( int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    C8_UNIFORM_01 returns a unit double complex pseudorandom number.
//
//  Discussion:
//
//    The angle should be uniformly distributed between 0 and 2 * PI,
//    the square root of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <double> C8_UNIFORM_01, a pseudorandom complex value.
//
{
  double r;
  int k;
  double pi = 3.141592653589793;
  double theta;
  complex <double> value;

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + 2147483647;
  }

  r = sqrt ( ( static_cast<double>(seed) * 4.656612875E-10 ) );

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + 2147483647;
  }

  theta = 2.0 * pi * ( static_cast<double>(seed) * 4.656612875E-10 );

  value = complex <double> ( r * cos ( theta ), r * sin ( theta ) );

  return value;
}
//****************************************************************************80

void c8mat_print ( int m, int n, complex <double> a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_PRINT prints a C8MAT.
//
//  Discussion:
//
//    A C8MAT is an array of complex <double> values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input, complex <double> A[M*N], the matrix.
//
//    Input, string TITLE, a title.
//
{
  c8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void c8mat_print_some ( int m, int n, complex <double> a[], int ilo, int jlo, 
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_PRINT_SOME prints some of a C8MAT.
//
//  Discussion:
//
//    A C8MAT is an array of complex <double> values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input, complex <double> A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
  complex <double> c;
  int i;
  int i2hi;
  int i2lo;
  int inc;
  int incx = 4;
  int j;
  int j2;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }
//
//  Print the columns of the matrix, in strips of INCX.
//
  for ( j2lo = jlo; j2lo <= i4_min ( jhi, n ); j2lo = j2lo + incx )
  {
    j2hi = j2lo + incx - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    inc = j2hi + 1 - j2lo;

    cout << "\n";
    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      j2 = j + 1 - j2lo;
      cout << "     " << setw(10) << j << "     ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(5) << i << ":";
//
//  Print out (up to) INCX entries in row I, that lie in the current strip.
//
      for ( j2 = 1; j2 <= inc; j2++ )
      {
        j = j2lo - 1 + j2;
        c = a[i-1+(j-1)*m];
        cout << "  " << setw(8) << real ( c )
             << "  " << setw(8) << imag ( c );
      }
      cout << "\n";
    }
  }

  return;
}
//****************************************************************************80

complex <double> *c8mat_test ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_TEST returns a test matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Output, complex <double> C8MAT_TEST[N*N], the matrix.
//
{
  complex <double> *a;
  double angle;
  complex <double> I = complex <double> ( 0.0, 1.0 );
  int i;
  int j;
  double pi = 3.141592653589793;

  a = new complex <double>[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      angle = 2.0 * pi * static_cast<double>(i * j) / static_cast<double>(n);

      a[i+j*n] = exp ( I * angle ) / sqrt ( static_cast<double>(n) );
    }
  }
  return a;
}
//****************************************************************************80

complex <double> *c8mat_test_inverse ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_TEST_INVERSE returns the inverse of a test matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Output, complex <double> C8MAT_TEST_INVERSE[N*N], the matrix.
//
{
  complex <double> *a;
  int i;
  int j;
  complex <double> t;

  a = c8mat_test ( n );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      t        = conj ( a[i+j*n] );
      a[i+j*n] = conj ( a[j+i*n] );
      a[j+i*n] = t;
    }
  }
  return a;
}
//****************************************************************************80

complex <double> *c8mat_uniform_01_new ( int m, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    C8MAT_UNIFORM_01_NEW returns a unit pseudorandom C8MAT.
//
//  Discussion:
//
//    A C8MAT is an array of complex <double> values.
//
//    The angles should be uniformly distributed between 0 and 2 * PI,
//    the square roots of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <double> C8MAT_UNIFORM_01_NEW[M*N], the pseudorandom 
//    complex matrix.
//
{
  complex <double> *c;
  int i;
  int j;
  double r;
  int k;
  const double r8_pi = 3.141592653589793;
  double theta;

  c = new complex <double> [m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + 2147483647;
      }

      r = sqrt ( static_cast<double>(seed) * 4.656612875E-10 );

      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + 2147483647;
      }

      theta = 2.0 * r8_pi * ( static_cast<double>(seed) * 4.656612875E-10 );

      c[i+j*m] = r * complex <double> ( cos ( theta ), sin ( theta ) );
    }
  }
  return c;
}
//****************************************************************************80

complex <double> *c8vec_uniform_01_new ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    C8VEC_UNIFORM_01_NEW returns a unit pseudorandom C8VEC.
//
//  Discussion:
//
//    A C8VEC is a vector of complex <double> values.
//
//    The angles should be uniformly distributed between 0 and 2 * PI,
//    the square roots of the radius uniformly distributed between 0 and 1.
//
//    This results in a uniform distribution of values in the unit circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values to compute.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, complex <double> C8VEC_UNIFORM_01_NEW[N], the pseudorandom 
//    complex vector.
//
{
  complex <double> *c;
  int i;
  double r;
  int k;
  const double r8_pi = 3.141592653589793;
  double theta;

  c = new complex <double> [n];

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + 2147483647;
    }

    r = sqrt ( static_cast<double>(seed) * 4.656612875E-10 );

    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + 2147483647;
    }

    theta = 2.0 * r8_pi * ( static_cast<double>(seed) * 4.656612875E-10 );

    c[i] = r * complex <double> ( cos ( theta ), sin ( theta ) );
  }

  return c;
}
//****************************************************************************80

float cabs1 ( complex <float> z )

//****************************************************************************80
//
//  Purpose:
//
//    CABS1 returns the L1 norm of a number.
//
//  Discussion:
//
//    This routine uses single precision complex arithmetic.
//
//    The L1 norm of a complex number is the sum of the absolute values
//    of the real and imaginary components.
//
//    CABS1 ( Z ) = abs ( real ( Z ) ) + abs ( imaginary ( Z ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2007
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for FORTRAN usage,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, pages 308-323, 1979.
//
//  Parameters:
//
//    Input, complex <float> Z, the number whose norm is desired.
//
//    Output, float CABS1, the L1 norm of Z.
//
{
  float value;

  value = r4_abs ( real ( z ) ) + r4_abs ( imag ( z ) );

  return value;
}
//****************************************************************************80

float cabs2 ( complex <float> z )

//****************************************************************************80
//
//  Purpose:
//
//    CABS2 returns the L2 norm of a number.
//
//  Discussion:
//
//    This routine uses single precision complex arithmetic.
//
//    The L2 norm of a complex number is the square root of the sum
//    of the squares of the real and imaginary components.
//
//    CABS2 ( Z ) = sqrt ( ( real ( Z ) )**2 + ( imaginary ( Z ) )**2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for FORTRAN usage,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, pages 308-323, 1979.
//
//  Parameters:
//
//    Input, complex <float> Z, the number whose norm is desired.
//
//    Output, float CABS2, the L2 norm of Z.
//
{
  float value;

  value = sqrt ( pow ( real ( z ), 2 )
               + pow ( imag ( z ), 2 ) );

  return value;
}
//****************************************************************************80

float cmach ( int job )

//****************************************************************************80
//
//  Purpose:
//
//    CMACH computes machine parameters for complex arithmetic.
//
//  Discussion:
//
//    Assume the computer has
//
//      B = base of arithmetic;
//      T = number of base B digits;
//      L = smallest possible exponent;
//      U = largest possible exponent;
//
//    then
//
//      EPS = B^(1-T)
//      TINY = 100.0 * B^(-L+T)
//      HUGE = 0.01 * B^(U-T)
//
//    If complex division is done by
//
//      1 / (X+i*Y) = (X-i*Y) / (X^2+Y^2)
//
//    then
//
//      TINY = sqrt ( TINY )
//      HUGE = sqrt ( HUGE )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for FORTRAN usage,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, pages 308-323, 1979.
//
//  Parameters:
//
//    Input, int JOB:
//    1, EPS is desired;
//    2, TINY is desired;
//    3, HUGE is desired.
//
//    Output, float CMACH, the requested value.
//
{
  float eps;
  float huge;
  float s;
  complex <float> temp1;
  complex <float> temp2;
  complex <float> temp3;
  float tiny;
  float value;

  eps = 1.0;

  for ( ; ; )
  {
    eps = eps / 2.0;
    s = 1.0 + eps;
    if ( s <= 1.0 )
    {
      break;
    }
  }

  eps = 2.0 * eps;

  s = 1.0;

  for ( ; ; )
  {
    tiny = s;
    s = s / 16.0;

    if ( s * 1.0 == 0.0 )
    {
      break;
    }
  }

  tiny = ( tiny / eps ) * 100.0;
//
//  Had to insert this manually!
//
  tiny = sqrt ( tiny );

  if ( false )
  {
    temp1 = complex <float> ( 1.0,  0.0 );
    temp2 = complex <float> ( tiny, 0.0 );
    temp3 = temp1 / temp2;

    s = real ( temp3 );

    if ( s != 1.0 / tiny )
    {
      tiny = sqrt ( tiny );
    }
  }

  huge = 1.0 / tiny;

  if ( job == 1 )
  {
    value = eps;
  }
  else if ( job == 2 )
  {
    value = tiny;
  }
  else if ( job == 3 )
  {
    value = huge;
  }
  else
  {
    value = 0.0;
  }

  return value;
}
//****************************************************************************80

complex <float> csign1 ( complex <float> z1, complex <float> z2 )

//****************************************************************************80
//
//  Purpose:
//
//    CSIGN1 is a transfer-of-sign function.
//
//  Discussion:
//
//    This routine uses single precision complex arithmetic.
//
//    The L1 norm is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> Z1, Z2, the arguments.
//
//    Output, complex <float> CSIGN1,  a complex value, with the magnitude of
//    Z1, and the argument of Z2.
//
{
  complex <float> value;

  if ( cabs1 ( z2 ) == 0.0 )
  {
    value = complex <float> ( 0.0, 0.0 );
  }
  else
  {
    value = cabs1 ( z1 ) * ( z2 / cabs1 ( z2 ) );
  }

  return value;
}
//****************************************************************************80

complex <float> csign2 ( complex <float> z1, complex <float> z2 )

//****************************************************************************80
//
//  Purpose:
//
//    CSIGN2 is a transfer-of-sign function.
//
//  Discussion:
//
//    This routine uses single precision complex arithmetic.
//
//    The L2 norm is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <float> Z1, Z2, the arguments.
//
//    Output, complex <float> CSIGN2,  a complex value, with the magnitude of
//    Z1, and the argument of Z2.
//
{
  complex <float> value;

  if ( cabs2 ( z2 ) == 0.0 )
  {
    value = complex <float> ( 0.0, 0.0 );
  }
  else
  {
    value = cabs2 ( z1 ) * ( z2 / cabs2 ( z2 ) );
  }

  return value;
}
//****************************************************************************80

double dmach ( int job )

//****************************************************************************80
//
//  Purpose:
//
//    DMACH computes machine parameters of double precision real arithmetic.
//
//  Discussion:
//
//    This routine is for testing only.  It is not required by LINPACK.
//
//    If there is trouble with the automatic computation of these quantities,
//    they can be set by direct assignment statements.
//
//    We assume the computer has
//
//      B = base of arithmetic;
//      T = number of base B digits;
//      L = smallest possible exponent;
//      U = largest possible exponent.
//
//    then
//
//      EPS = B^(1-T)
//      TINY = 100.0 * B^(-L+T)
//      HUGE = 0.01 * B^(U-T)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
//    David Kincaid, Fred Krogh.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int JOB:
//    1: requests EPS;
//    2: requests TINY;
//    3: requests HUGE.
//
//    Output, double DMACH, the requested value.
//
{
  double eps;
  double huge;
  double s;
  double tiny;
  double value;

  eps = 1.0;
  for ( ; ; )
  {
    value = 1.0 + ( eps / 2.0 );
    if ( value <= 1.0 )
    {
      break;
    }
    eps = eps / 2.0;
  }

  s = 1.0;

  for ( ; ; )
  {
    tiny = s;
    s = s / 16.0;

    if ( s * 1.0 == 0.0 )
    {
      break;
    }

  }

  tiny = ( tiny / eps ) * 100.0;
  huge = 1.0 / tiny;

  if ( job == 1 )
  {
    value = eps;
  }
  else if ( job == 2 )
  {
    value = tiny;
  }
  else if ( job == 3 )
  {
    value = huge;
  }
  else
  {
    cout << "\n";
    cout << "DMACH - Fatal error!\n";
    cout << "  Illegal input value of JOB = " << job << "\n";
    exit ( 1 );
  }

  return value;
}
//****************************************************************************80

bool lsame ( char ca, char cb )

//****************************************************************************80
//
//  Purpose:
//
//    LSAME returns TRUE if CA is the same letter as CB regardless of case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, char CA, CB, the characters to compare.
//
//    Output, bool LSAME, is TRUE if the characters are equal,
//    disregarding case.
//
{
  if ( ca == cb )
  {
    return true;
  }

  if ( 'A' <= ca && ca <= 'Z' )
  {
    if ( ca - 'A' == cb - 'a' )
    {
      return true;
    }
  }
  else if ( 'a' <= ca && ca <= 'z' )
  {
    if ( ca - 'a' == cb - 'A' )
    {
      return true;
    }
  }

  return false;
}
//****************************************************************************80

float r4_abs ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ABS returns the absolute value of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the quantity whose absolute value is desired.
//
//    Output, float R4_ABS, the absolute value of X.
//
{
  float value;

  if ( 0.0 <= x )
  {
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

float r4_max ( float x, float y )

//****************************************************************************80
//
//  Purpose:
//
//    R4_MAX returns the maximum of two R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, Y, the quantities to compare.
//
//    Output, float R4_MAX, the maximum of X and Y.
//
{
  float value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

float r4_sign ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_SIGN returns the sign of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the number whose sign is desired.
//
//    Output, float R4_SIGN, the sign of X.
//
{
  float value;

  if ( x < 0.0 )
  {
    value = -1.0;
  }
  else
  {
    value = 1.0;
  }
  return value;
}
//****************************************************************************80

float r4_uniform_01 ( int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4_UNIFORM_01 returns a unit pseudorandom R4.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2^31 - 1 )
//      r4_uniform_01 = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R4_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, float R4_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  float value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R4_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + 2147483647;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number.
//
  value = static_cast<float>(seed) * 4.656612875E-10;

  return value;
}
//****************************************************************************80

float r4_uniform_ab ( float a, float b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4_UNIFORM_AB returns a scaled pseudorandom R4.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, float R4_UNIFORM_AB, a number strictly between A and B.
//
{
  int i4_huge = 2147483647;
  int k;
  float value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R4_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  value = static_cast<float>(seed) * 4.656612875E-10;

  value = a + ( b - a ) * value;

  return value;
}
//****************************************************************************80

void r4mat_print ( int m, int n, float a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_PRINT prints an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
//    in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*M]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, float A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  r4mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r4mat_print_some ( int m, int n, float a[], int ilo, int jlo, int ihi,
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_PRINT_SOME prints some of an R4MAT.
//
//  Discussion:
//
//    An R4MAT is a doubly dimensioned array of R4 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, float A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j - 1 << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(5) << i - 1 << ": ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

float *r4mat_test ( char trans, int lda, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_TEST sets up a test matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
//
//    10 February 2014
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, char TRANS, indicates whether matrix is to be transposed.
//    'N', no transpose.
//    'T', transpose the matrix.
//
//    Input, int LDA, the leading dimension of the matrix.
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Output, float R4MAT_TEST[LDA*?], the matrix.
//    if TRANS is 'N', then the matrix is stored in LDA*N entries,
//    as an M x N matrix;
//    if TRANS is 'T', then the matrix is stored in LDA*M entries,
//    as an N x M matrix.
//
{
  float *a;
  int i;
  int j;

  if ( trans == 'N' )
  {
    a = new float[lda*n];

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        a[i+j*lda] = static_cast<float> ( 10 * ( i + 1 ) + ( j + 1 ) );
      }
    }
  }
  else
  {
    a = new float[lda*m];

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        a[j+i*lda] = static_cast<float> ( 10 * ( i + 1 ) + ( j + 1 ) );
      }
    }
  }
  return a;
}
//****************************************************************************80

void r4mat_uniform_01 ( int m, int n, int &seed, float r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_UNIFORM_01 returns a unit pseudorandom R4MAT.
//
//  Discussion:
//
//    An R4MAT is an array of R4's.
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has
//    been updated.
//
//    Output, float R[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int j;
  int k;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R4MAT_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + i4_huge;
      }

      r[i+j*m] = static_cast<float>(seed) * 4.656612875E-10;
    }
  }
  return;
}
//****************************************************************************80

void r4vec_print ( int n, float a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_PRINT prints an R4VEC.
//
//  Discussion:
//
//    An R4VEC is a vector of R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, float A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(8)  << i
         << ": " << setw(14) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

double *r8mat_test ( char trans, int lda, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TEST sets up a test matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//  
//  Modified:
//
//    10 February 2014
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, char TRANS, indicates whether matrix is to be transposed.
//    'N', no transpose.
//    'T', transpose the matrix.
//
//    Input, int LDA, the leading dimension of the matrix.
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Output, double R8MAT_TEST[?], the matrix.
//    if TRANS is 'N', then the matrix is stored in LDA*N entries,
//    as an M x N matrix;
//    if TRANS is 'T', then the matrix is stored in LDA*M entries,
//    as an N x M matrix.
//
{
  double *a;
  int i;
  int j;

  if ( trans == 'N' )
  {
    a = new double[lda*n];

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        a[i+j*lda] = static_cast<double> ( 10 * ( i + 1 ) + ( j + 1 ) );
      }
    }
  }
  else
  {
    a = new double[lda*m];

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        a[j+i*lda] = static_cast<double> ( 10 * ( i + 1 ) + ( j + 1 ) );
      }
    }
  }
  return a;
}
//****************************************************************************80

void r8mat_uniform_01 ( int m, int n, int &seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has
//    been updated.
//
//    Output, double R[M*N], a matrix of pseudorandom values.
//
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + i4_huge;
      }

      r[i+j*m] = static_cast<double>(seed) * 4.656612875E-10;
    }
  }
  return;
}
//****************************************************************************80

float smach ( int job )

//****************************************************************************80
//
//  Purpose:
//
//    SMACH computes machine parameters of float arithmetic.
//
//  Discussion:
//
//    This routine is for testing only.  It is not required by LINPACK.
//
//    If there is trouble with the automatic computation of these quantities,
//    they can be set by direct assignment statements.
//
//    We assume the computer has
//
//      B = base of arithmetic;
//      T = number of base B digits;
//      L = smallest possible exponent;
//      U = largest possible exponent.
//
//    then
//
//      EPS = B**(1-T)
//      TINY = 100.0 * B**(-L+T)
//      HUGE = 0.01 * B**(U-T)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int JOB:
//    1: requests EPS;
//    2: requests TINY;
//    3: requests HUGE.
//
//    Output, float SMACH, the requested value.
//
{
  float eps;
  float huge;
  float s;
  float tiny;
  float value;

  eps = 1.0;
  for ( ; ; )
  {
    value = 1.0 + ( eps / 2.0 );
    if ( value <= 1.0 )
    {
      break;
    }
    eps = eps / 2.0;
  }

  s = 1.0;

  for ( ; ; )
  {
    tiny = s;
    s = s / 16.0;

    if ( s * 1.0 == 0.0 )
    {
      break;
    }

  }

  tiny = ( tiny / eps ) * 100.0;
  huge = 1.0 / tiny;

  if ( job == 1 )
  {
    value = eps;
  }
  else if ( job == 2 )
  {
    value = tiny;
  }
  else if ( job == 3 )
  {
    value = huge;
  }
  else
  {
    xerbla ( "SMACH", 1 );
  }

  return value;
}
//****************************************************************************80

int triangle_upper_to_i4 ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UPPER_TO_I4 converts an upper triangular coordinate to an integer.
//
//  Discussion:
//
//    Triangular coordinates are handy when storing a naturally triangular
//    array (such as the upper half of a matrix) in a linear array.
//
//    Thus, for example, we might consider storing
//
//    (0,0) (0,1) (0,2) (0,3)
//          (1,1) (1,2) (1,3)
//                (2,2) (2,3)
//                      (3,3)
//
//    as the linear array
//
//    (0,0) (0,1) (1,1) (0,2) (1,2) (2,2) (0,3) (1,3) (2,3) (3,3)
//
//    Here, the quantities in parenthesis represent the natural row and
//    column indices of a single number when stored in a rectangular array.
//
//    Thus, our goal is, given the row I and column J of the data,
//    to produce the value K which indicates its position in the linear
//    array.
//
//    The triangular numbers are the indices associated with the
//    diagonal elements of the original array, T(0,0), T(1,1), T(2,2), 
//    T(3,3) and so on.
//
//    The formula is:
//
//      K = I + ( ( J * ( J + 1 ) ) / 2
//
//  Example:
//
//    I  J  K
//
//    0  0  0
//    0  1  1
//    1  1  2
//    0  2  3
//    1  2  4
//    2  2  5
//    0  3  6
//    1  3  7
//    2  3  8
//    3  3  9
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2017
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the row and column indices.  I and J must
//    be nonnegative, and I must not be greater than J.
//
//    Output, int TRIANGLE_UPPER_TO_I4, the linear index of the (I,J) element.
//
{
  int value;

  if ( i < 0 )
  {
    cerr << "\n";
    cerr << "TRIANGLE_UPPER_TO_I4 - Fatal error!\n";
    cerr << "  I < 0.\n";
    cerr << "  I = " << i << "\n";
    exit ( 1 );
  }
  else if ( j < 0 )
  {
    cerr << "\n";
    cerr << "TRIANGLE_UPPER_TO_I4 - Fatal error!\n";
    cerr << "  J < 0.\n";
    cerr << "  J = " << j << "\n";
    exit ( 1 );
  }
  else if ( j < i )
  {
    cerr << "\n";
    cerr << "TRIANGLE_UPPER_TO_I4 - Fatal error!\n";
    cerr << "  J < I.\n";
    cerr << "  I = " << i << "\n";
    cerr << "  J = " << j << "\n";
    exit ( 1 );
  }

  value = i + (  j * ( j + 1 ) ) / 2;

  return value;
}
//****************************************************************************80

void xerbla ( string srname, int info )

//****************************************************************************80
//
//  Purpose:
//
//    XERBLA is an error handler for the LAPACK routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, string SRNAME, the name of the routine
//    which called XERBLA.
//
//    Input, int INFO, the position of the invalid parameter in
//    the parameter list of the calling routine.
//
{
  cerr << "\n";
  cerr << "XERBLA - Fatal error!\n";
  cerr << "  On entry to routine '" << srname << "',\n";
  cerr << "  input parameter number " << info << " had an illegal value.\n";

  exit ( 1 );
}
//****************************************************************************80

double zabs1 ( complex <double> z )

//****************************************************************************80
//
//  Purpose:
//
//    ZABS1 returns the L1 norm of a complex <double>.
//
//  Discussion:
//
//    This routine uses double precision complex arithmetic.
//
//    The L1 norm of a complex number is the sum of the absolute values
//    of the real and imaginary components.
//
//    ZABS1 ( Z ) = abs ( real ( Z ) ) + abs ( imaginary ( Z ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for FORTRAN usage,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, pages 308-323, 1979.
//
//  Parameters:
//
//    Input, complex <double> Z, the number whose norm is desired.
//
//    Output, double ZABS1, the L1 norm of Z.
//
{
  double value;

  value = fabs ( real ( z ) ) + fabs ( imag ( z ) );

  return value;
}
//****************************************************************************80

double zabs2 ( complex <double> z )

//****************************************************************************80
//
//  Purpose:
//
//    ZABS2 returns the L2 norm of a complex <double>.
//
//  Discussion:
//
//    This routine uses double precision complex arithmetic.
//
//    The L2 norm of a complex number is the square root of the sum
//    of the squares of the real and imaginary components.
//
//    ZABS2 ( Z ) = sqrt ( ( real ( Z ) )**2 + ( imaginary ( Z ) )**2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for FORTRAN usage,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, pages 308-323, 1979.
//
//  Parameters:
//
//    Input, complex <double> Z, the number whose norm is desired.
//
//    Output, float ZABS2, the L2 norm of Z.
//
{
  double value;
  double zi;
  double zr;

  zr = real ( z );
  zi = imag ( z );
  value = sqrt ( zr * zr + zi * zi );

  return value;
}
//****************************************************************************80

double zmach ( int job )

//****************************************************************************80
//
//  Purpose:
//
//    ZMACH computes machine parameters for complex <double> arithmetic.
//
//  Discussion:
//
//    Assume the computer has
//
//      B = base of arithmetic;
//      T = number of base B digits;
//      L = smallest possible exponent;
//      U = largest possible exponent;
//
//    then
//
//      EPS = B^(1-T)
//      TINY = 100.0 * B^(-L+T)
//      HUGE = 0.01 * B^(U-T)
//
//    If complex division is done by
//
//      1 / (X+i*Y) = (X-i*Y) / (X^2+Y^2)
//
//    then
//
//      TINY = sqrt ( TINY )
//      HUGE = sqrt ( HUGE )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for FORTRAN usage,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, pages 308-323, 1979.
//
//  Parameters:
//
//    Input, int JOB:
//    1, EPS is desired;
//    2, TINY is desired;
//    3, HUGE is desired.
//
//    Output, double ZMACH, the requested value.
//
{
  double eps;
  double huge;
  double s;
  complex <double> temp1;
  complex <double> temp2;
  complex <double> temp3;
  double tiny;
  double value;

  eps = 1.0;

  for ( ; ; )
  {
    eps = eps / 2.0;
    s = 1.0 + eps;
    if ( s <= 1.0 )
    {
      break;
    }
  }

  eps = 2.0 * eps;

  s = 1.0;

  for ( ; ; )
  {
    tiny = s;
    s = s / 16.0;

    if ( s * 1.0 == 0.0 )
    {
      break;
    }
  }

  tiny = ( tiny / eps ) * 100.0;
//
//  Had to insert this manually!
//
  tiny = sqrt ( tiny );

  if ( false )
  {
    temp1 = complex <double> ( 1.0,  0.0 );
    temp2 = complex <double> ( tiny, 0.0 );
    temp3 = temp1 / temp2;

    s = real ( temp3 );

    if ( s != 1.0 / tiny )
    {
      tiny = sqrt ( tiny );
    }
  }

  huge = 1.0 / tiny;

  if ( job == 1 )
  {
    value = eps;
  }
  else if ( job == 2 )
  {
    value = tiny;
  }
  else if ( job == 3 )
  {
    value = huge;
  }
  else
  {
    value = 0.0;
  }

  return value;
}
//****************************************************************************80

complex <double> zsign1 ( complex <double> z1, complex <double> z2 )

//****************************************************************************80
//
//  Purpose:
//
//    ZSIGN1 is a complex <double> transfer-of-sign function.
//
//  Discussion:
//
//    This routine uses double precision complex arithmetic.
//
//    The L1 norm is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <double> Z1, Z2, the arguments.
//
//    Output, complex <double> ZSIGN1,  a complex value, with the magnitude of
//    Z1, and the argument of Z2.
//
{
  complex <double> value;

  if ( zabs1 ( z2 ) == 0.0 )
  {
    value = complex <double> ( 0.0, 0.0 );
  }
  else
  {
    value = zabs1 ( z1 ) * ( z2 / zabs1 ( z2 ) );
  }

  return value;
}
//****************************************************************************80

complex <double> zsign2 ( complex <double> z1, complex <double> z2 )

//****************************************************************************80
//
//  Purpose:
//
//    ZSIGN2 is a complex <double> transfer-of-sign function.
//
//  Discussion:
//
//    This routine uses double precision complex arithmetic.
//
//    The L2 norm is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, complex <double> Z1, Z2, the arguments.
//
//    Output, complex <double> ZSIGN2,  a complex value, with the magnitude of
//    Z1, and the argument of Z2.
//
{
  complex <double> value;

  if ( zabs2 ( z2 ) == 0.0 )
  {
    value = complex <double> ( 0.0, 0.0 );
  }
  else
  {
    value = zabs2 ( z1 ) * ( z2 / zabs2 ( z2 ) );
  }

  return value;
}
