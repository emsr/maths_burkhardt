void abwe1 ( int n, int m, double eps, double coef2, bool even, double b[], 
  double *x, double *w );
void abwe2 ( int n, int m, double eps, double coef2, bool even, double b[], 
  double *x, double *w1, double *w2 );
void kronrod ( int n, double eps, double x[], double w1[], double w2[] );
void kronrod_adjust ( double a, double b, int n, double x[], double w1[], double w2[] );
double r8_epsilon ( );
void timestamp ( );
