# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "qwv_2d.hpp"

int main ( );
void test01 ( );
void test02 ( );
void padua_point_set ( int l, double x[], double y[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for QWV_2D_PRB.
//
//  Discussion:
//
//    QWV_2D_PRB tests the QWV_2D library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 May 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  std::cout << "\n";
  std::cout << "QWV_2D_PRB:\n";
  std::cout << "  C++ version\n";
  std::cout << "  Test the QWV_2D library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "QWV_2D_PRB:\n";
  std::cout << "  Normal end of execution.\n";
  std::cout << "\n";
  timestamp ( );

  return 0;
}
//*****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests QWV_2D for a trio of points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 May 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double d;
  int n;
  int t;
  double *w;
  double *x;
  double *y;

  a = -1.0;
  b = +1.0;
  c = -1.0;
  d = +1.0;

  t = 1;
  n = ( ( t + 1 ) * ( t + 2 ) ) / 2;

  std::cout << "\n";
  std::cout << "\n";
  std::cout << "TEST01:\n";
  std::cout << "  Compute the weights associated with an interpolatory\n";
  std::cout << "  quadrature rule defined by N=(T+1)*(T+2)/2 points,\n";
  std::cout << "  exact for polynomials of total degree T or less.\n";
  std::cout << "\n";
  std::cout << "  Degree T = " << t << "\n";
  std::cout << "  Number of points N = " << n << "\n";
  std::cout << "  X Interval = [" << a << "," << b << "]\n";
  std::cout << "  Y Interval = [" << c << "," << d << "]\n";
//
//  Set the points.
//
  x = new double[n];
  y = new double[n];

  x[0] =  1.0;
  y[0] = -1.0;
  x[1] =  1.0;
  y[1] =  1.0;
  x[2] = -1.0;
  y[2] =  1.0;
  r8vec2_print ( n, x, y, "  Abscissas:" );
//
//  Compute the weights.
//
  w = qwv_2d ( t, n, a, b, c, d, x, y );

  r8vec_print ( n, w, "  Weights:" );
//
//  Free memory.
//
  delete [] w;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests QWV_2D for Padua points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 May 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double d;
  int n;
  int t;
  double *w;
  double *x;
  double *y;

  a = -1.0;
  b = +1.0;
  c = -1.0;
  d = +1.0;

  std::cout << "\n";
  std::cout << "TEST02:\n";
  std::cout << "  Compute the weights associated with an interpolatory\n";
  std::cout << "  quadrature rule defined by N=(T+1)*(T+2)/2 points,\n";
  std::cout << "  exact for polynomials of total degree T or less.\n";
  std::cout << "  X Interval = [" << a << "," << b << "]\n";
  std::cout << "  Y Interval = [" << c << "," << d << "]\n";

  for ( t = 0; t <= 10; t++ )
  {
    n = ( ( t + 1 ) * ( t + 2 ) ) / 2;

    std::cout << "\n";
    std::cout << "  Degree T = " << t << "\n";
    std::cout << "  Number of points N = " << n << "\n";

    x = new double[n];
    y = new double[n];

    padua_point_set ( t, x, y );
//
//  Compute the weights.
//
    w = qwv_2d ( t, n, a, b, c, d, x, y );

    r8vec_print_16 ( n, w, "  Weights:" );
//
//  Free memory.
//
    delete [] w;
    delete [] x;
    delete [] y;
  }
  return;
}
//****************************************************************************80

void padua_point_set ( int l, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    PADUA_POINT_SET sets the Padua points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Marco Caliari, Stefano de Marchi, Marco Vianello,
//    Bivariate interpolation on the square at new nodal sets,
//    Applied Mathematics and Computation,
//    Volume 165, Number 2, 2005, pages 261-274.
//
//  Parameters:
//
//    Input, int L, the level.
//    0 <= L <= 10.
//
//    Output, double X[N], Y[N], the Padua points.
//
{
  if ( l == 0 )
  {
    x[ 0] =  0.000000000000000;
    y[ 0] =  0.000000000000000;    
  }
  else if ( l == 1 )
  {
    x[ 0] = -1.000000000000000;
    y[ 0] = -1.000000000000000;    
    x[ 1] = -1.000000000000000;
    y[ 1] =  1.000000000000000;    
    x[ 2] =  1.000000000000000;
    y[ 2] =  0.000000000000000;    
  }
  else if ( l == 2 )
  {
    x[ 0] = -1.000000000000000;
    y[ 0] = -1.000000000000000;    
    x[ 1] = -1.000000000000000;
    y[ 1] =  0.5000000000000001;
    x[ 2] =  0.000000000000000;
    y[ 2] = -0.4999999999999998;    
    x[ 3] =  0.000000000000000;
    y[ 3] =  1.000000000000000;    
    x[ 4] =  1.000000000000000;
    y[ 4] = -1.000000000000000;
    x[ 5] =  1.000000000000000;
    y[ 5] =  0.5000000000000001;    
  }
  else if ( l == 3 )
  {
    x[ 0] = -1.000000000000000;
    y[ 0] = -1.000000000000000;    
    x[ 1] = -1.000000000000000;
    y[ 1] =  0.000000000000000;
    x[ 2] = -1.000000000000000;
    y[ 2] =  1.000000000000000;    
    x[ 3] = -0.4999999999999998;
    y[ 3] = -0.7071067811865475;    
    x[ 4] = -0.4999999999999998;
    y[ 4] =  0.7071067811865476;    
    x[ 5] =  0.5000000000000001;
    y[ 5] = -1.000000000000000;    
    x[ 6] =  0.5000000000000001;
    y[ 6] =  0.000000000000000;    
    x[ 7] =  0.5000000000000001;
    y[ 7] =  1.000000000000000;    
    x[ 8] =  1.000000000000000;
    y[ 8] = -0.7071067811865475;    
    x[ 9] =  1.000000000000000;
    y[ 9] =  0.7071067811865476;    
  }
  else if ( l == 4 )
  {
    x[ 0] = -1.000000000000000;
    y[ 0] = -1.000000000000000;   
    x[ 1] = -1.000000000000000;
    y[ 1] = -0.3090169943749473;    
    x[ 2] = -1.000000000000000;
    y[ 2] =  0.8090169943749475;    
    x[ 3] = -0.7071067811865475;
    y[ 3] = -0.8090169943749473;    
    x[ 4] = -0.7071067811865475;
    y[ 4] =  0.3090169943749475;
    x[ 5] = -0.7071067811865475;
    y[ 5] =  1.000000000000000;    
    x[ 6] =  0.000000000000000;
    y[ 6] = -1.000000000000000;   
    x[ 7] =  0.000000000000000;
    y[ 7] = -0.3090169943749473;    
    x[ 8] =  0.000000000000000;
    y[ 8] =  0.8090169943749475;    
    x[ 9] =  0.7071067811865476;
    y[ 9] = -0.8090169943749473;    
    x[10] =  0.7071067811865476;
    y[10] =  0.3090169943749475; 
    x[11] =  0.7071067811865476;
    y[11] =  1.000000000000000;    
    x[12] =  1.000000000000000;
    y[12] = -1.000000000000000;   
    x[13] =  1.000000000000000;
    y[13] = -0.3090169943749473;    
    x[14] =  1.000000000000000;
    y[14] =  0.8090169943749475;    
  }
  else if ( l == 5 )
  {
    x[ 0] = -1.000000000000000;
    y[ 0] = -1.000000000000000;
    x[ 1] = -1.000000000000000;
    y[ 1] = -0.4999999999999998;    
    x[ 2] = -1.000000000000000;
    y[ 2] =  0.5000000000000001;    
    x[ 3] = -1.000000000000000;
    y[ 3] =  1.000000000000000;    
    x[ 4] = -0.8090169943749473;
    y[ 4] = -0.8660254037844387;
    x[ 5] = -0.8090169943749473;
    y[ 5] =  0.000000000000000;   
    x[ 6] = -0.8090169943749473;
    y[ 6] =  0.8660254037844387;
    x[ 7] = -0.3090169943749473;
    y[ 7] = -1.000000000000000;   
    x[ 8] = -0.3090169943749473;
    y[ 8] = -0.4999999999999998;    
    x[ 9] = -0.3090169943749473;
    y[ 9] =  0.5000000000000001;
    x[10] = -0.3090169943749473;
    y[10] =  1.000000000000000;   
    x[11] =  0.3090169943749475;
    y[11] = -0.8660254037844387;
    x[12] =  0.3090169943749475;
    y[12] =  0.000000000000000;    
    x[13] =  0.3090169943749475;
    y[13] =  0.8660254037844387;    
    x[14] =  0.8090169943749475;
    y[14] = -1.000000000000000;    
    x[15] =  0.8090169943749475;
    y[15] = -0.4999999999999998;    
    x[16] =  0.8090169943749475;
    y[16] =  0.5000000000000001;
    x[17] =  0.8090169943749475;
    y[17] =  1.000000000000000;
    x[18] =  1.000000000000000;
    y[18] = -0.8660254037844387;    
    x[19] =  1.000000000000000;
    y[19] =  0.000000000000000;
    x[20] =  1.000000000000000;
    y[20] =  0.8660254037844387; 
  }
  else if ( l == 6 )
  {   
    x[ 0] = -1.000000000000000;
    y[ 0] = -1.000000000000000;   
    x[ 1] = -1.000000000000000;
    y[ 1] = -0.6234898018587335;    
    x[ 2] = -1.000000000000000;
    y[ 2] =  0.2225209339563144;    
    x[ 3] = -1.000000000000000;
    y[ 3] =  0.9009688679024191;    
    x[ 4] = -0.8660254037844387;
    y[ 4] = -0.9009688679024190;    
    x[ 5] = -0.8660254037844387;
    y[ 5] = -0.2225209339563143;    
    x[ 6] = -0.8660254037844387;
    y[ 6] =  0.6234898018587336;
    x[ 7] = -0.8660254037844387;
    y[ 7] =  1.000000000000000;
    x[ 8] = -0.4999999999999998;
    y[ 8] = -1.000000000000000;   
    x[ 9] = -0.4999999999999998;
    y[ 9] = -0.6234898018587335;    
    x[10] = -0.4999999999999998;
    y[10] =  0.2225209339563144;    
    x[11] = -0.4999999999999998;
    y[11] =  0.9009688679024191;    
    x[12] =  0.000000000000000;
    y[12] = -0.9009688679024190;    
    x[13] =  0.000000000000000;
    y[13] = -0.2225209339563143;
    x[14] =  0.000000000000000;
    y[14] =  0.6234898018587336;    
    x[15] =  0.000000000000000;
    y[15] =  1.000000000000000;
    x[16] =  0.5000000000000001;
    y[16] = -1.000000000000000;    
    x[17] =  0.5000000000000001;
    y[17] = -0.6234898018587335;    
    x[18] =  0.5000000000000001;
    y[18] =  0.2225209339563144;    
    x[19] =  0.5000000000000001;
    y[19] =  0.9009688679024191;    
    x[20] =  0.8660254037844387;
    y[20] = -0.9009688679024190;    
    x[21] =  0.8660254037844387;
    y[21] = -0.2225209339563143;    
    x[22] =  0.8660254037844387;
    y[22] =  0.6234898018587336;
    x[23] =  0.8660254037844387;
    y[23] =  1.000000000000000;    
    x[24] =  1.000000000000000;
    y[24] = -1.000000000000000;   
    x[25] =  1.000000000000000;
    y[25] = -0.6234898018587335;    
    x[26] =  1.000000000000000;
    y[26] =  0.2225209339563144;
    x[27] =  1.000000000000000;
    y[27] =  0.9009688679024191;
  }
  else if ( l == 7 )
  {    
    x[ 0] = -1.000000000000000;
    y[ 0] = -1.000000000000000;
    x[ 1] = -1.000000000000000;
    y[ 1] = -0.7071067811865475;    
    x[ 2] = -1.000000000000000;
    y[ 2] =  0.000000000000000;   
    x[ 3] = -1.000000000000000;
    y[ 3] =  0.7071067811865476; 
    x[ 4] = -1.000000000000000;
    y[ 4] =  1.000000000000000;    
    x[ 5] = -0.9009688679024190;
    y[ 5] = -0.9238795325112867;    
    x[ 6] = -0.9009688679024190;
    y[ 6] = -0.3826834323650897;    
    x[ 7] = -0.9009688679024190;
    y[ 7] =  0.3826834323650898;    
    x[ 8] = -0.9009688679024190;
    y[ 8] =  0.9238795325112867;
    x[ 9] = -0.6234898018587335;
    y[ 9] = -1.000000000000000;   
    x[10] = -0.6234898018587335;
    y[10] = -0.7071067811865475;
    x[11] = -0.6234898018587335;
    y[11] =  0.000000000000000;   
    x[12] = -0.6234898018587335;
    y[12] =  0.7071067811865476;
    x[13] = -0.6234898018587335;
    y[13] =  1.000000000000000;   
    x[14] = -0.2225209339563143;
    y[14] = -0.9238795325112867;    
    x[15] = -0.2225209339563143;
    y[15] = -0.3826834323650897;    
    x[16] = -0.2225209339563143;
    y[16] =  0.3826834323650898;    
    x[17] = -0.2225209339563143;
    y[17] =  0.9238795325112867;
    x[18] =  0.2225209339563144;
    y[18] = -1.000000000000000;   
    x[19] =  0.2225209339563144;
    y[19] = -0.7071067811865475;
    x[20] =  0.2225209339563144;
    y[20] =  0.000000000000000;
    x[21] =  0.2225209339563144;
    y[21] =  0.7071067811865476;
    x[22] =  0.2225209339563144;
    y[22] =  1.000000000000000;   
    x[23] =  0.6234898018587336;
    y[23] = -0.9238795325112867;    
    x[24] =  0.6234898018587336;
    y[24] = -0.3826834323650897;    
    x[25] =  0.6234898018587336;
    y[25] =  0.3826834323650898;    
    x[26] =  0.6234898018587336;
    y[26] =  0.9238795325112867;
    x[27] =  0.9009688679024191;
    y[27] = -1.000000000000000;   
    x[28] =  0.9009688679024191;
    y[28] = -0.7071067811865475;
    x[29] =  0.9009688679024191;
    y[29] =  0.000000000000000;   
    x[30] =  0.9009688679024191;
    y[30] =  0.7071067811865476;
    x[31] =  0.9009688679024191;
    y[31] =  1.000000000000000;   
    x[32] =  1.000000000000000;
    y[32] = -0.9238795325112867;    
    x[33] =  1.000000000000000;
    y[33] = -0.3826834323650897;    
    x[34] =  1.000000000000000;
    y[34] =  0.3826834323650898;    
    x[35] =  1.000000000000000;
    y[35] =  0.9238795325112867; 
  }
  else if ( l == 8 )
  { 
    x[ 0] = -1.000000000000000;
    y[ 0] = -1.000000000000000;   
    x[ 1] = -1.000000000000000;
    y[ 1] = -0.7660444431189779;    
    x[ 2] = -1.000000000000000;
    y[ 2] = -0.1736481776669303;    
    x[ 3] = -1.000000000000000;
    y[ 3] =  0.5000000000000001;    
    x[ 4] = -1.000000000000000;
    y[ 4] =  0.9396926207859084;    
    x[ 5] = -0.9238795325112867;
    y[ 5] = -0.9396926207859083;    
    x[ 6] = -0.9238795325112867;
    y[ 6] = -0.4999999999999998;    
    x[ 7] = -0.9238795325112867;
    y[ 7] =  0.1736481776669304;    
    x[ 8] = -0.9238795325112867;
    y[ 8] =  0.7660444431189780;
    x[ 9] = -0.9238795325112867;
    y[ 9] =  1.000000000000000;
    x[10] = -0.7071067811865475;
    y[10] = -1.000000000000000;   
    x[11] = -0.7071067811865475;
    y[11] = -0.7660444431189779;    
    x[12] = -0.7071067811865475;
    y[12] = -0.1736481776669303;    
    x[13] = -0.7071067811865475;
    y[13] =  0.5000000000000001;    
    x[14] = -0.7071067811865475;
    y[14] =  0.9396926207859084;    
    x[15] = -0.3826834323650897;
    y[15] = -0.9396926207859083;    
    x[16] = -0.3826834323650897;
    y[16] = -0.4999999999999998;    
    x[17] = -0.3826834323650897;
    y[17] =  0.1736481776669304;    
    x[18] = -0.3826834323650897;
    y[18] =  0.7660444431189780;
    x[19] = -0.3826834323650897;
    y[19] =  1.000000000000000;    
    x[20] =  0.000000000000000;
    y[20] = -1.000000000000000;   
    x[21] =  0.000000000000000;
    y[21] = -0.7660444431189779;    
    x[22] =  0.000000000000000;
    y[22] = -0.1736481776669303;    
    x[23] =  0.000000000000000;
    y[23] =  0.5000000000000001;    
    x[24] =  0.000000000000000;
    y[24] =  0.9396926207859084;    
    x[25] =  0.3826834323650898;
    y[25] = -0.9396926207859083;    
    x[26] =  0.3826834323650898;
    y[26] = -0.4999999999999998;    
    x[27] =  0.3826834323650898;
    y[27] =  0.1736481776669304;    
    x[28] =  0.3826834323650898;
    y[28] =  0.7660444431189780;
    x[29] =  0.3826834323650898;
    y[29] =  1.000000000000000;
    x[30] =  0.7071067811865476;
    y[30] = -1.000000000000000;   
    x[31] =  0.7071067811865476;
    y[31] = -0.7660444431189779;    
    x[32] =  0.7071067811865476;
    y[32] = -0.1736481776669303;    
    x[33] =  0.7071067811865476;
    y[33] =  0.5000000000000001;    
    x[34] =  0.7071067811865476;
    y[34] =  0.9396926207859084;    
    x[35] =  0.9238795325112867;
    y[35] = -0.9396926207859083;    
    x[36] =  0.9238795325112867;
    y[36] = -0.4999999999999998;    
    x[37] =  0.9238795325112867;
    y[37] =  0.1736481776669304;    
    x[38] =  0.9238795325112867;
    y[38] =  0.7660444431189780;
    x[39] =  0.9238795325112867;
    y[39] =  1.000000000000000;    
    x[40] =  1.000000000000000;
    y[40] = -1.000000000000000;   
    x[41] =  1.000000000000000;
    y[41] = -0.7660444431189779;    
    x[42] =  1.000000000000000;
    y[42] = -0.1736481776669303;    
    x[43] =  1.000000000000000;
    y[43] =  0.5000000000000001;
    x[44] =  1.000000000000000;
    y[44] =  0.9396926207859084;   
  }
  else if ( l == 9 )
  { 
    x[ 0] = -1.000000000000000;
    y[ 0] = -1.000000000000000;   
    x[ 1] = -1.000000000000000;
    y[ 1] = -0.8090169943749473;    
    x[ 2] = -1.000000000000000;
    y[ 2] = -0.3090169943749473;    
    x[ 3] = -1.000000000000000;
    y[ 3] =  0.3090169943749475;
    x[ 4] = -1.000000000000000;
    y[ 4] =  0.8090169943749475;    
    x[ 5] = -1.000000000000000;
    y[ 5] =  1.000000000000000;   
    x[ 6] = -0.9396926207859083;
    y[ 6] = -0.9510565162951535;    
    x[ 7] = -0.9396926207859083;
    y[ 7] = -0.5877852522924730;
    x[ 8] = -0.9396926207859083;
    y[ 8] =  0.000000000000000;   
    x[ 9] = -0.9396926207859083;
    y[ 9] =  0.5877852522924731;    
    x[10] = -0.9396926207859083;
    y[10] =  0.9510565162951535;
    x[11] = -0.7660444431189779;
    y[11] = -1.000000000000000;   
    x[12] = -0.7660444431189779;
    y[12] = -0.8090169943749473;    
    x[13] = -0.7660444431189779;
    y[13] = -0.3090169943749473;    
    x[14] = -0.7660444431189779;
    y[14] =  0.3090169943749475;    
    x[15] = -0.7660444431189779;
    y[15] =  0.8090169943749475; 
    x[16] = -0.7660444431189779;
    y[16] =  1.000000000000000;   
    x[17] = -0.4999999999999998;
    y[17] = -0.9510565162951535;    
    x[18] = -0.4999999999999998;
    y[18] = -0.5877852522924730;
    x[19] = -0.4999999999999998;
    y[19] =  0.000000000000000;   
    x[20] = -0.4999999999999998;
    y[20] =  0.5877852522924731;    
    x[21] = -0.4999999999999998;
    y[21] =  0.9510565162951535;
    x[22] = -0.1736481776669303;
    y[22] = -1.000000000000000;   
    x[23] = -0.1736481776669303;
    y[23] = -0.8090169943749473;    
    x[24] = -0.1736481776669303;
    y[24] = -0.3090169943749473;    
    x[25] = -0.1736481776669303;
    y[25] =  0.3090169943749475;    
    x[26] = -0.1736481776669303;
    y[26] =  0.8090169943749475;
    x[27] = -0.1736481776669303;
    y[27] =  1.000000000000000;   
    x[28] =  0.1736481776669304;
    y[28] = -0.9510565162951535;    
    x[29] =  0.1736481776669304;
    y[29] = -0.5877852522924730;
    x[30] =  0.1736481776669304;
    y[30] =  0.000000000000000;   
    x[31] =  0.1736481776669304;
    y[31] =  0.5877852522924731;    
    x[32] =  0.1736481776669304;
    y[32] =  0.9510565162951535;
    x[33] =  0.5000000000000001;
    y[33] = -1.000000000000000;   
    x[34] =  0.5000000000000001;
    y[34] = -0.8090169943749473;    
    x[35] =  0.5000000000000001;
    y[35] = -0.3090169943749473;    
    x[36] =  0.5000000000000001;
    y[36] =  0.3090169943749475;    
    x[37] =  0.5000000000000001;
    y[37] =  0.8090169943749475;
    x[38] =  0.5000000000000001;
    y[38] =  1.000000000000000;   
    x[39] =  0.7660444431189780;
    y[39] = -0.9510565162951535;    
    x[40] =  0.7660444431189780;
    y[40] = -0.5877852522924730;
    x[41] =  0.7660444431189780;
    y[41] =  0.000000000000000;   
    x[42] =  0.7660444431189780;
    y[42] =  0.5877852522924731;    
    x[43] =  0.7660444431189780;
    y[43] =  0.9510565162951535;
    x[44] =  0.9396926207859084;
    y[44] = -1.000000000000000;   
    x[45] =  0.9396926207859084;
    y[45] = -0.8090169943749473;    
    x[46] =  0.9396926207859084;
    y[46] = -0.3090169943749473;    
    x[47] =  0.9396926207859084;
    y[47] =  0.3090169943749475;    
    x[48] =  0.9396926207859084;
    y[48] =  0.8090169943749475;
    x[49] =  0.9396926207859084;
    y[49] =  1.000000000000000;   
    x[50] =  1.000000000000000;
    y[50] = -0.9510565162951535;
    x[51] =  1.000000000000000;
    y[51] = -0.5877852522924730;    
    x[52] =  1.000000000000000;
    y[52] =  0.000000000000000;   
    x[53] =  1.000000000000000;
    y[53] =  0.5877852522924731;
    x[54] =  1.000000000000000;
    y[54] =  0.9510565162951535;    
  }
  else if ( l == 10 )
  {
    x[ 0] = -1.000000000000000;
    y[ 0] = -1.000000000000000;   
    x[ 1] = -1.000000000000000;
    y[ 1] = -0.8412535328311811;    
    x[ 2] = -1.000000000000000;
    y[ 2] = -0.4154150130018863;    
    x[ 3] = -1.000000000000000;
    y[ 3] =  0.1423148382732851;    
    x[ 4] = -1.000000000000000;
    y[ 4] =  0.6548607339452851;    
    x[ 5] = -1.000000000000000;
    y[ 5] =  0.9594929736144974;    
    x[ 6] = -0.9510565162951535;
    y[ 6] = -0.9594929736144974;    
    x[ 7] = -0.9510565162951535;
    y[ 7] = -0.6548607339452850;    
    x[ 8] = -0.9510565162951535;
    y[ 8] = -0.1423148382732850;    
    x[ 9] = -0.9510565162951535;
    y[ 9] =  0.4154150130018864;    
    x[10] = -0.9510565162951535;
    y[10] =  0.8412535328311812;
    x[11] = -0.9510565162951535;
    y[11] =  1.000000000000000;
    x[12] = -0.8090169943749473;
    y[12] = -1.000000000000000;   
    x[13] = -0.8090169943749473;
    y[13] = -0.8412535328311811;    
    x[14] = -0.8090169943749473;
    y[14] = -0.4154150130018863;    
    x[15] = -0.8090169943749473;
    y[15] =  0.1423148382732851;    
    x[16] = -0.8090169943749473;
    y[16] =  0.6548607339452851;    
    x[17] = -0.8090169943749473;
    y[17] =  0.9594929736144974;    
    x[18] = -0.5877852522924730;
    y[18] = -0.9594929736144974;    
    x[19] = -0.5877852522924730;
    y[19] = -0.6548607339452850;    
    x[20] = -0.5877852522924730;
    y[20] = -0.1423148382732850;    
    x[21] = -0.5877852522924730;
    y[21] =  0.4154150130018864;    
    x[22] = -0.5877852522924730;
    y[22] =  0.8412535328311812;
    x[23] = -0.5877852522924730;
    y[23] =  1.000000000000000;
    x[24] = -0.3090169943749473;
    y[24] = -1.000000000000000;   
    x[25] = -0.3090169943749473;
    y[25] = -0.8412535328311811;    
    x[26] = -0.3090169943749473;
    y[26] = -0.4154150130018863;    
    x[27] = -0.3090169943749473;
    y[27] =  0.1423148382732851;    
    x[28] = -0.3090169943749473;
    y[28] =  0.6548607339452851;    
    x[29] = -0.3090169943749473;
    y[29] =  0.9594929736144974;    
    x[30] =  0.000000000000000;
    y[30] = -0.9594929736144974;    
    x[31] =  0.000000000000000;
    y[31] = -0.6548607339452850;    
    x[32] =  0.000000000000000;
    y[32] = -0.1423148382732850;    
    x[33] =  0.000000000000000;
    y[33] =  0.4154150130018864;
    x[34] =  0.000000000000000;
    y[34] =  0.8412535328311812;    
    x[35] =  0.000000000000000;
    y[35] =  1.000000000000000; 
    x[36] =  0.3090169943749475;
    y[36] = -1.000000000000000;    
    x[37] =  0.3090169943749475;
    y[37] = -0.8412535328311811;    
    x[38] =  0.3090169943749475;
    y[38] = -0.4154150130018863;    
    x[39] =  0.3090169943749475;
    y[39] =  0.1423148382732851;    
    x[40] =  0.3090169943749475;
    y[40] =  0.6548607339452851;    
    x[41] =  0.3090169943749475;
    y[41] =  0.9594929736144974;    
    x[42] =  0.5877852522924731;
    y[42] = -0.9594929736144974;    
    x[43] =  0.5877852522924731;
    y[43] = -0.6548607339452850;    
    x[44] =  0.5877852522924731;
    y[44] = -0.1423148382732850;    
    x[45] =  0.5877852522924731;
    y[45] =  0.4154150130018864;    
    x[46] =  0.5877852522924731;
    y[46] =  0.8412535328311812;
    x[47] =  0.5877852522924731;
    y[47] =  1.000000000000000;
    x[48] =  0.8090169943749475;
    y[48] = -1.000000000000000;   
    x[49] =  0.8090169943749475;
    y[49] = -0.8412535328311811;    
    x[50] =  0.8090169943749475;
    y[50] = -0.4154150130018863;    
    x[51] =  0.8090169943749475;
    y[51] =  0.1423148382732851;    
    x[52] =  0.8090169943749475;
    y[52] =  0.6548607339452851;    
    x[53] =  0.8090169943749475;
    y[53] =  0.9594929736144974;    
    x[54] =  0.9510565162951535;
    y[54] = -0.9594929736144974;    
    x[55] =  0.9510565162951535;
    y[55] = -0.6548607339452850;    
    x[56] =  0.9510565162951535;
    y[56] = -0.1423148382732850;    
    x[57] =  0.9510565162951535;
    y[57] =  0.4154150130018864;    
    x[58] =  0.9510565162951535;
    y[58] =  0.8412535328311812;
    x[59] =  0.9510565162951535;
    y[59] =  1.000000000000000;    
    x[60] =  1.000000000000000;
    y[60] = -1.000000000000000;   
    x[61] =  1.000000000000000;
    y[61] = -0.8412535328311811;    
    x[62] =  1.000000000000000;
    y[62] = -0.4154150130018863;    
    x[63] =  1.000000000000000;
    y[63] =  0.1423148382732851;    
    x[64] =  1.000000000000000;
    y[64] =  0.6548607339452851; 
    x[65] =  1.000000000000000;
    y[65] =  0.9594929736144974;
  }
  else
  {
    std::cerr << "\n";
    std::cerr << "PADUA_POINT_SET - Fatal error!\n";
    std::cerr << "  Illegal value of L = " << l << "\n";
    std::cerr << "  Legal values are 1 through 10.\n";
    exit ( 1 );
  }

  return;
}
