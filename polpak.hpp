
#include "export.h"

void agm_values ( int &n_data, double &a, double &b, double &fx );

BURKHARDT_EXPORT
double agud ( double gamma );
int align_enum ( int m, int n );

BURKHARDT_EXPORT
void bell ( int n, int b[] );
void bell_values ( int &n_data, int &n, int &c );
double benford ( int ival );

BURKHARDT_EXPORT
void bernoulli_number ( int n, double b[] );

BURKHARDT_EXPORT
void bernoulli_number2 ( int n, double b[] );

BURKHARDT_EXPORT
double bernoulli_number3 ( int n );

BURKHARDT_EXPORT
double bernoulli_poly ( int n, double x );

BURKHARDT_EXPORT
double bernoulli_poly2 ( int n, double x );
void bernoulli_number_values ( int &n_data, int &n, double &c );

BURKHARDT_EXPORT
void bernstein_poly ( int n, double x, double bern[] );
void bernstein_poly_values ( int &n_data, int &n, int &k, double &x, double &b );
void beta_values ( int *n_data, double *x, double *y, double *fxy );
void bpab ( int n, double x, double a, double b, double bern[] );
double *cardan_poly ( int n, double x, double s );
void cardan_poly_coef ( int n, double s, double c[] );
double *cardinal_cos ( int j, int m, int n, double t[] );
double *cardinal_sin ( int j, int m, int n, double t[] );
void catalan ( int n, int c[] );
void catalan_row_next ( bool next, int n, int irow[] );
void catalan_values ( int *n_data, int *n, int *c );
void charlier ( int n, double a, double x, double values[] );
double *cheby_t_poly ( int m, int n, double x[] );
double *cheby_t_poly_coef ( int n );
void cheby_t_poly_values ( int &n_data, int &n, double &x, double &fx );
double *cheby_t_poly_zero ( int n );
double *cheby_u_poly ( int m, int n, double x[] );
void cheby_u_poly_coef ( int n, double c[] );
void cheby_u_poly_values ( int &n_data, int &n, double &x, double &fx );
double *cheby_u_poly_zero ( int n );
void chebyshev_discrete ( int n, int m, double x, double v[] );
int collatz_count ( int n );
void collatz_count_max ( int n, int *i_max, int *j_max );
void collatz_count_values ( int *n_data, int *n, int *count );
void comb_row_next ( int n, int row[] );
int commul ( int n, int nfact, int iarray[] );
double complete_symmetric_poly ( int n, int r, double x[] );
double cos_power_int ( double a, double b, int n );
void cos_power_int_values ( int &n_data, double &a, double &b, int &n,
  double &fx );
void erf_values ( int &n_data, double &x, double &fx );

BURKHARDT_EXPORT
void euler_number ( int n, int e[] );

BURKHARDT_EXPORT
double euler_number2 ( int n );
void euler_number_values ( int *n_data, int *n, int *c );

BURKHARDT_EXPORT
double euler_poly ( int n, double x );
void eulerian ( int n, int e[] );
int f_hofstadter ( int n );
int fibonacci_direct ( int n );
void fibonacci_floor ( int n, int *f, int *i );
void fibonacci_recursive ( int n, int f[] );
int g_hofstadter ( int n );
void gamma_log_values ( int &n_data, double &x, double &fx );
void gamma_values ( int &n_data, double &x, double &fx );

BURKHARDT_EXPORT
void gegenbauer_poly ( int n, double alpha, double x, double cx[] );
void gegenbauer_poly_values ( int *n_data, int *n, double *a, double *x, double *fx );

BURKHARDT_EXPORT
void gen_hermite_poly ( int n, double x, double mu, double p[] );

BURKHARDT_EXPORT
void gen_laguerre_poly ( int n, double alpha, double x, double cx[] );

BURKHARDT_EXPORT
double gud ( double x );
void gud_values ( int *n_data, double *x, double *fx );
int h_hofstadter ( int n );
int hail ( int n );
void hermite_poly_phys ( int n, double x, double cx[] );
void hermite_poly_phys_coef ( int n, double c[] );
void hermite_poly_phys_values ( int *n_data, int *n, double *x, double *fx );
void hyper_2f1_values ( int *n_data, double *a, double *b, double *c, 
  double *x, double *fx );
int i4_choose ( int n, int k );
void i4_factor ( int n, int maxfactor, int &nfactor, int factor[], 
  int power[], int &nleft );
int i4_factorial ( int n );
void i4_factorial_values ( int *n_data, int *n, int *fn );
int i4_factorial2 ( int n );
void i4_factorial2_values ( int *n_data, int *n, int *fn );
int i4_gcd ( int i, int j );
bool i4_is_prime ( int n );
bool i4_is_triangular ( int i );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_partition_distinct_count ( int n );
int i4_sign ( int i );
void i4_swap ( int *i, int *j );
void i4_to_triangle ( int k, int *i, int *j );
int i4_uniform_ab ( int ilo, int ihi, int &seed );
void i4mat_print ( int m, int n, int a[], std::string title );
void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi,
  int jhi, std::string title );

BURKHARDT_EXPORT
double *jacobi_poly ( int n, double alpha, double beta, double x );
void jacobi_poly_values ( int &n_data, int &n, double &a, double &b, double &x,
  double &fx );
int jacobi_symbol ( int q, int p );
void krawtchouk ( int n, double p, double x, int m, double v[] );

BURKHARDT_EXPORT
void laguerre_associated ( int n, int m, double x, double cx[] );

BURKHARDT_EXPORT
void laguerre_poly ( int n, double x, double cx[] );
void laguerre_poly_coef ( int n, double c[] );
void laguerre_polynomial_values ( int *n_data, int *n, double *x, double *fx );

BURKHARDT_EXPORT
void legendre_associated ( int n, int m, double x, double cx[] );

BURKHARDT_EXPORT
void legendre_associated_normalized ( int n, int m, double x, double cx[] );

BURKHARDT_EXPORT
void legendre_associated_normalized_sphere_values ( int &n_data, int &n, int &m, 
  double &x, double &fx );
void legendre_associated_values ( int *n_data, int *n, int *m, double *x, 
  double *fx );

BURKHARDT_EXPORT
void legendre_function_q ( int n, double x, double cx[] );
void legendre_function_q_values ( int *n_data, int *n, double *x, double *fx );
void legendre_poly ( int n, double x, double cx[], double cpx[] );
void legendre_poly_coef ( int n, double c[] );
void legendre_poly_values ( int *n_data, int *n, double *x, double *fx );
int legendre_symbol ( int q, int p );

BURKHARDT_EXPORT
double lerch ( double z, int s, double a );
void lerch_values ( int *n_data, double *z, int *s, double *a, double *fx );
void lock ( int n, int a[] );

BURKHARDT_EXPORT
void meixner ( int n, double beta, double c, double x, double v[] );
int mertens ( int n );
void mertens_values ( int *n_data, int *n, int *c );
int moebius ( int n );
void moebius_values ( int *n_data, int *n, int *c );
void motzkin ( int n, int a[] );
double normal_01_cdf_inverse ( double p );
void normal_01_cdf_values ( int &n_data, double &x, double &fx );
int omega ( int n );
void omega_values ( int *n_data, int *n, int *c );
void partition_distinct_count_values ( int *n_data, int *n, int *c );
int pentagon_num ( int n );

BURKHARDT_EXPORT
int phi ( int n );
void phi_values ( int *n_data, int *n, int *c );
int plane_partition_num ( int n );
int poly_bernoulli ( int n, int k );
int poly_coef_count ( int dim, int degree );
int prime ( int n );
void psi_values ( int *n_data, double *x, double *fx );
int pyramid_num ( int n );
int pyramid_square_num ( int n );

BURKHARDT_EXPORT
double r8_agm ( double a, double b );

BURKHARDT_EXPORT
double r8_beta ( double x, double y );

BURKHARDT_EXPORT
double r8_choose ( int n, int k );
double r8_epsilon ( );

BURKHARDT_EXPORT
double r8_erf ( double x );

BURKHARDT_EXPORT
double r8_erf_inverse ( double y );

BURKHARDT_EXPORT
double r8_euler_constant ( );

BURKHARDT_EXPORT
double r8_factorial ( int n );

BURKHARDT_EXPORT
double r8_factorial_log ( int n );
void r8_factorial_log_values ( int *n_data, int *n, double *fn );
void r8_factorial_values ( int *n_data, int *n, double *fn );
double r8_gamma ( double x );
double r8_huge ( );
double r8_hyper_2f1 ( double a, double b, double c, double x );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_mop ( int i );
int r8_nint ( double x );
double r8_pi ( );

BURKHARDT_EXPORT
double r8_psi ( double xx );
double r8_uniform_01 ( int *seed );
int r8poly_degree ( int na, double a[] );
void r8poly_print ( int n, double a[], std::string title );
double r8poly_value_horner ( int n, double a[], double x );
double *r8vec_linspace_new ( int n, double a_first, double a_last );
void r8vec_print ( int n, double a[], std::string title );
void r8vec_zero ( int n, double a[] );
int s_len_trim ( char* s );
int sigma ( int n );
void sigma_values ( int *n_data, int *n, int *c );
int simplex_num ( int m, int n );
double sin_power_int ( double a, double b, int n );
void sin_power_int_values ( int *n_data, double *a, double *b, int *n, 
  double *fx );
int slice ( int dim_num, int slice_num );

BURKHARDT_EXPORT
void spherical_harmonic ( int l, int m, double theta, double phi, 
  double c[], double s[] );
void spherical_harmonic_values ( int *n_data, int *l, int *m, double *theta,
  double *phi, double *yr, double *yi ) ;

BURKHARDT_EXPORT
int *stirling1 ( int n, int m );

BURKHARDT_EXPORT
int *stirling2 ( int n, int m );
int tau ( int n );
void tau_values ( int *n_data, int *n, int *c );
int tetrahedron_num ( int n );
void timestamp ( );
int triangle_num ( int n );
int triangle_to_i4 ( int i, int j );
int trinomial ( int i, int j, int k );
int v_hofstadter ( int n );
void vibonacci ( int n, int &seed, int v[] );
void zeckendorf ( int n, int m_max, int *m, int i_list[], int f_list[] );

BURKHARDT_EXPORT
double zernike_poly ( int m, int n, double rho );
double *zernike_poly_coef ( int m, int n );

BURKHARDT_EXPORT
double zeta ( double p );
void zeta_values ( int *n_data, int *n, double *zeta );
