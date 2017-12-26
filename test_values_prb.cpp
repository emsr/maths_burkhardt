# include <cmath>
# include <complex>
# include <cstdlib>
# include <cstring>
# include <iomanip>
# include <iostream>

# include "test_values.hpp"

int main ( );

void abram0_values_test ( );
void abram1_values_test ( );
void abram2_values_test ( );
void agm_values_test ( );
void airy_ai_values_test ( );
void airy_ai_int_values_test ( );
void airy_ai_prime_values_test ( );
void airy_bi_values_test ( );
void airy_bi_int_values_test ( );
void airy_bi_prime_values_test ( );
void airy_cai_values_test ( );
void airy_cbi_values_test ( );
void airy_gi_values_test ( );
void airy_hi_values_test ( );
void arccos_values_test ( );
void arccosh_values_test ( );
void arcsin_values_test ( );
void arcsinh_values_test ( );
void arctan_values_test ( );
void arctan_int_values_test ( );
void arctanh_values_test ( );

void bei0_values_test ( );
void bei1_values_test ( );
void bell_values_test ( );
void ber0_values_test ( );
void ber1_values_test ( );
void bernoulli_number_values_test ( );
void bernoulli_poly_values_test ( );
void bernstein_poly_01_values_test ( );
void bessel_i0_values_test ( );
void bessel_i0_int_values_test ( );
void bessel_i0_spherical_values_test ( );
void bessel_i1_values_test ( );
void bessel_i1_spherical_values_test ( );
void bessel_in_values_test ( );
void bessel_ix_values_test ( );
void bessel_j0_values_test ( );
void bessel_j0_int_values_test ( );
void bessel_j0_spherical_values_test ( );
void bessel_j0_zero_values_test ( );
void bessel_j1_values_test ( );
void bessel_j1_spherical_values_test ( );
void bessel_jn_values_test ( );
void bessel_jx_values_test ( );
void bessel_k0_values_test ( );
void bessel_k0_int_values_test ( );
void bessel_k1_values_test ( );
void bessel_kn_values_test ( );
void bessel_kx_values_test ( );
void bessel_y0_values_test ( );
void bessel_y0_int_values_test ( );
void bessel_y0_spherical_values_test ( );
void bessel_y1_values_test ( );
void bessel_y1_spherical_values_test ( );
void bessel_yn_values_test ( );
void bessel_yx_values_test ( );
void beta_values_test ( );
void beta_cdf_values_test ( );
void beta_inc_values_test ( );
void beta_log_values_test ( );
void beta_noncentral_cdf_values_test ( );
void beta_pdf_values_test ( );
void binomial_values_test ( );
void binomial_cdf_values_test ( );
void binomial_pdf_values_test ( );
void bivariate_normal_cdf_values_test ( );

void catalan_values_test ( );
void cauchy_cdf_values_test ( );
void cbrt_values_test ( );
void cheby_t_poly_values_test ( );
void cheby_t01_poly_values_test ( );
void cheby_u_poly_values_test ( );
void cheby_u01_poly_values_test ( );
void cheby_v_poly_values_test ( );
void cheby_v01_poly_values_test ( );
void cheby_w_poly_values_test ( );
void cheby_w01_poly_values_test ( );
void chi_values_test ( );
void chi_square_cdf_values_test ( );
void chi_square_pdf_values_test ( );
void chi_square_noncentral_cdf_values_test ( );
void ci_values_test ( );
void cin_values_test ( );
void cinh_values_test ( );
void clausen_values_test ( );
void clebsch_gordan_values_test ( );
void collatz_count_values_test ( );
void cos_values_test ( );
void cos_degree_values_test ( );
void cos_power_int_values_test ( );
void cosh_values_test ( );
void cot_values_test ( );
void cp_values_test ( );

void dawson_values_test ( );
void debye1_values_test ( );
void debye2_values_test ( );
void debye3_values_test ( );
void debye4_values_test ( );
void dedekind_sum_values_test ( );
void dielectric_values_test ( );
void dilogarithm_values_test ( );

void e1_values_test ( );
void ei_values_test ( );
void elliptic_ea_values_test ( );
void elliptic_em_values_test ( );
void elliptic_ka_values_test ( );
void elliptic_km_values_test ( );
void erf_values_test ( );
void erfc_values_test ( );
void euler_number_values_test ( );
void euler_poly_values_test ( );
void exp_values_test ( );
void exp3_int_values_test ( );
void exponential_01_pdf_values_test ( );
void exponential_cdf_values_test ( );
void exponential_pdf_values_test ( );
void extreme_values_cdf_values_test ( );

void f_cdf_values_test ( );
void f_noncentral_cdf_values_test ( );
void fresnel_cos_values_test ( );
void fresnel_sin_values_test ( );
void frobenius_number_data_values_test ( );
void frobenius_number_order_values_test ( );
void frobenius_number_order2_values_test ( );

void gamma_values_test ( );
void gamma_01_pdf_values_test ( );
void gamma_cdf_values_test ( );
void gamma_inc_values_test ( );
void gamma_inc_p_values_test ( );
void gamma_inc_q_values_test ( );
void gamma_inc_tricomi_values_test ( );
void gamma_log_values_test ( );
void gamma_pdf_values_test ( );
void gegenbauer_poly_values_test ( );
void geometric_cdf_values_test ( );
void goodwin_values_test ( );
void gud_values_test ( );

void hermite_function_values_test ( );
void hermite_poly_phys_values_test ( );
void hermite_poly_prob_values_test ( );
void hyper_1f1_values_test ( );
void hyper_2f1_values_test ( );
void hypergeometric_cdf_values_test ( );
void hypergeometric_pdf_values_test ( );
void hypergeometric_u_values_test ( );

void i0ml0_values_test ( );
void i1ml1_values_test ( );
void i4_factorial_values_test ( );
void i4_factorial2_values_test ( );
void i4_fall_values_test ( );
void i4_rise_values_test ( );
void int_values_test ( );
void inverse_chi_square_pdf_values_test ( );
void inverse_gamma_pdf_values_test ( );

void jacobi_cn_values_test ( );
void jacobi_dn_values_test ( );
void jacobi_poly_values_test ( );
void jacobi_sn_values_test ( );
void jed_ce_values_test ( );
void jed_mjd_values_test ( );
void jed_rd_values_test ( );
void jed_weekday_values_test ( );

void kei0_values_test ( );
void kei1_values_test ( );
void ker0_values_test ( );
void ker1_values_test ( );

void laguerre_associated_values_test ( );
void laguerre_general_values_test ( );
void laguerre_polynomial_values_test ( );
void lambert_w_values_test ( );
void laplace_cdf_values_test ( );
void legendre_associated_values_test ( );
void legendre_associated_normalized_values_test ( );
void legendre_associated_normalized_sphere_values_test ( );
void legendre_poly_values_test ( );
void legendre_function_q_values_test ( );
void lerch_values_test ( );
void lobachevsky_values_test ( );
void lobatto_polynomial_values_test ( );
void lobatto_polynomial_derivatives_test ( );
void log_values_test ( );
void log_normal_cdf_values_test ( );
void log_series_cdf_values_test ( );
void log10_values_test ( );
void logarithmic_integral_values_test ( );
void logistic_cdf_values_test ( );

void mertens_values_test ( );
void moebius_values_test ( );
void multinomial_pdf_values_test ( );

void negative_binomial_cdf_values_test ( );
void nine_j_values_test ( );
void normal_01_cdf_values_test ( );
void normal_01_pdf_values_test ( );
void normal_cdf_values_test ( );
void normal_pdf_values_test ( );

void omega_values_test ( );
void owen_values_test ( );

void partition_count_values_test ( );
void partition_distinct_count_values_test ( );
void phi_values_test ( );
void pi_values_test ( );
void poisson_cdf_values_test ( );
void polylogarithm_values_test ( );
void prandtl_values_test ( );
void prime_values_test ( );
void psat_values_test ( );
void psi_values_test ( );

void r8_factorial_values_test ( );
void r8_factorial_log_values_test ( );
void r8_factorial2_values_test ( );
void r8_fall_values_test ( );
void r8_rise_values_test ( );
void rayleigh_cdf_values_test ( );

void scaled_inverse_chi_square_pdf_values_test ( );
void secvir_values_test ( );
void shi_values_test ( );
void si_values_test ( );
void sigma_values_test ( );
void sin_values_test ( );
void sin_degree_values_test ( );
void sin_power_int_values_test ( );
void sinh_values_test ( );
void six_j_values_test ( );
void sound_values_test ( );
void sphere_unit_area_values_test ( );
void sphere_unit_volume_values_test ( );
void spherical_harmonic_values_test ( );
void sqrt_values_test ( );
void stirling1_values_test ( );
void stirling2_values_test ( );
void stromgen_values_test ( );
void struve_h0_values_test ( );
void struve_h1_values_test ( );
void struve_l0_values_test ( );
void struve_l1_values_test ( );
void student_cdf_values_test ( );
void student_noncentral_cdf_values_test ( );
void subfactorial_values_test ( );
void surten_values_test ( );
void synch1_values_test ( );
void synch2_values_test ( );

void tan_values_test ( );
void tanh_values_test ( );
void tau_values_test ( );
void thercon_values_test ( );
void three_j_values_test ( );
void tran02_values_test ( );
void tran03_values_test ( );
void tran04_values_test ( );
void tran05_values_test ( );
void tran06_values_test ( );
void tran07_values_test ( );
void tran08_values_test ( );
void tran09_values_test ( );
void trigamma_values_test ( );
void truncated_normal_ab_cdf_test ( );
void truncated_normal_ab_pdf_test ( );
void truncated_normal_a_cdf_test ( );
void truncated_normal_a_pdf_test ( );
void truncated_normal_b_cdf_test ( );
void truncated_normal_b_pdf_test ( );
void tsat_values_test ( );

void van_der_corput_values_test ( );
void viscosity_values_test ( );
void von_mises_cdf_values_test ( );

void weekday_values_test ( );
void weibull_cdf_values_test ( );

void zeta_values_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TEST_VALUES_PRB.
//
//  Discussion:
//
//    TEST_VALUES_PRB tests the TEST_VALUE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 January 2016
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  std::cout << "\n";
  std::cout << "TEST_VALUES_PRB:\n";
  std::cout << "  C++ version,\n";
  std::cout << "  Test the TEST_VALUES library.\n";

  abram0_values_test ( );
  abram1_values_test ( );
  abram2_values_test ( );
  agm_values_test ( );
  airy_ai_values_test ( );
  airy_ai_int_values_test ( );
  airy_ai_prime_values_test ( );
  airy_bi_values_test ( );
  airy_bi_int_values_test ( );
  airy_bi_prime_values_test ( );
  airy_cai_values_test ( );
  airy_cbi_values_test ( );
  airy_gi_values_test ( );
  airy_hi_values_test ( );
  arccos_values_test ( );
  arccosh_values_test ( );
  arcsin_values_test ( );
  arcsinh_values_test ( );
  arctan_values_test ( );
  arctan_int_values_test ( );
  arctanh_values_test ( );

  bei0_values_test ( );
  bei1_values_test ( );
  bell_values_test ( );
  ber0_values_test ( );
  ber1_values_test ( );
  bernoulli_number_values_test ( );
  bernoulli_poly_values_test ( );
  bernstein_poly_01_values_test ( );
  bessel_i0_values_test ( );
  bessel_i0_int_values_test ( );
  bessel_i0_spherical_values_test ( );
  bessel_i1_values_test ( );
  bessel_i1_spherical_values_test ( );
  bessel_in_values_test ( );
  bessel_ix_values_test ( );
  bessel_j0_values_test ( );
  bessel_j0_int_values_test ( );
  bessel_j0_spherical_values_test ( );
  bessel_j0_zero_values_test ( );
  bessel_j1_values_test ( );
  bessel_j1_spherical_values_test ( );
  bessel_jn_values_test ( );
  bessel_jx_values_test ( );
  bessel_k0_values_test ( );
  bessel_k0_int_values_test ( );
  bessel_k1_values_test ( );
  bessel_kn_values_test ( );
  bessel_kx_values_test ( );
  bessel_y0_values_test ( );
  bessel_y0_int_values_test ( );
  bessel_y0_spherical_values_test ( );
  bessel_y1_values_test ( );
  bessel_y1_spherical_values_test ( );
  bessel_yn_values_test ( );
  bessel_yx_values_test ( );
  beta_values_test ( );
  beta_cdf_values_test ( );
  beta_inc_values_test ( );
  beta_log_values_test ( );
  beta_noncentral_cdf_values_test ( );
  beta_pdf_values_test ( );
  binomial_values_test ( );
  binomial_cdf_values_test ( );
  binomial_pdf_values_test ( );
  bivariate_normal_cdf_values_test ( );

  catalan_values_test ( );
  cauchy_cdf_values_test ( );
  cbrt_values_test ( );
  cheby_t_poly_values_test ( );
  cheby_t01_poly_values_test ( );
  cheby_u_poly_values_test ( );
  cheby_u01_poly_values_test ( );
  cheby_v_poly_values_test ( );
  cheby_v01_poly_values_test ( );
  cheby_w_poly_values_test ( );
  cheby_w01_poly_values_test ( );
  chi_values_test ( );
  chi_square_cdf_values_test ( );
  chi_square_pdf_values_test ( );
  chi_square_noncentral_cdf_values_test ( );
  ci_values_test ( );
  cin_values_test ( );
  cinh_values_test ( );
  clausen_values_test ( );
  clebsch_gordan_values_test ( );
  collatz_count_values_test ( );
  cos_values_test ( );
  cos_degree_values_test ( );
  cos_power_int_values_test ( );
  cosh_values_test ( );
  cot_values_test ( );
  cp_values_test ( );

  dawson_values_test ( );
  debye1_values_test ( );
  debye2_values_test ( );
  debye3_values_test ( );
  debye4_values_test ( );
  dedekind_sum_values_test ( );
  dielectric_values_test ( );
  dilogarithm_values_test ( );

  e1_values_test ( );
  ei_values_test ( );
  elliptic_ea_values_test ( );
  elliptic_em_values_test ( );
  elliptic_ka_values_test ( );
  elliptic_km_values_test ( );
  erf_values_test ( );
  erfc_values_test ( );
  euler_number_values_test ( );
  euler_poly_values_test ( );
  exp_values_test ( );
  exp3_int_values_test ( );
  exponential_01_pdf_values_test ( );
  exponential_cdf_values_test ( );
  exponential_pdf_values_test ( );
  extreme_values_cdf_values_test ( );

  f_cdf_values_test ( );
  f_noncentral_cdf_values_test ( );
  fresnel_cos_values_test ( );
  fresnel_sin_values_test ( );
  frobenius_number_data_values_test ( );
  frobenius_number_order_values_test ( );
  frobenius_number_order2_values_test ( );

  gamma_values_test ( );
  gamma_01_pdf_values_test ( );
  gamma_cdf_values_test ( );
  gamma_inc_values_test ( );
  gamma_inc_p_values_test ( );
  gamma_inc_q_values_test ( );
  gamma_inc_tricomi_values_test ( );
  gamma_log_values_test ( );
  gamma_pdf_values_test ( );
  gegenbauer_poly_values_test ( );
  geometric_cdf_values_test ( );
  goodwin_values_test ( );
  gud_values_test ( );

  hermite_function_values_test ( );
  hermite_poly_phys_values_test ( );
  hermite_poly_prob_values_test ( );
  hyper_1f1_values_test ( );
  hyper_2f1_values_test ( );
  hypergeometric_cdf_values_test ( );
  hypergeometric_pdf_values_test ( );
  hypergeometric_u_values_test ( );

  i0ml0_values_test ( );
  i1ml1_values_test ( );
  i4_factorial_values_test ( );
  i4_factorial2_values_test ( );
  i4_fall_values_test ( );
  i4_rise_values_test ( );
  int_values_test ( );
  inverse_chi_square_pdf_values_test ( );
  inverse_gamma_pdf_values_test ( );

  jacobi_cn_values_test ( );
  jacobi_dn_values_test ( );
  jacobi_poly_values_test ( );
  jacobi_sn_values_test ( );
  jed_ce_values_test ( );
  jed_mjd_values_test ( );
  jed_rd_values_test ( );
  jed_weekday_values_test ( );

  kei0_values_test ( );
  kei1_values_test ( );
  ker0_values_test ( );
  ker1_values_test ( );

  laguerre_associated_values_test ( );
  laguerre_general_values_test ( );
  laguerre_polynomial_values_test ( );
  lambert_w_values_test ( );
  laplace_cdf_values_test ( );
  legendre_associated_values_test ( );
  legendre_associated_normalized_values_test ( );
  legendre_associated_normalized_sphere_values_test ( );
  legendre_poly_values_test ( );
  legendre_function_q_values_test ( );
  lerch_values_test ( );
  lobachevsky_values_test ( );
  lobatto_polynomial_values_test ( );
  lobatto_polynomial_derivatives_test ( );
  log_values_test ( );
  log_normal_cdf_values_test ( );
  log_series_cdf_values_test ( );
  log10_values_test ( );
  logarithmic_integral_values_test ( );
  logistic_cdf_values_test ( );

  mertens_values_test ( );
  moebius_values_test ( );
  multinomial_pdf_values_test ( );

  negative_binomial_cdf_values_test ( );
  nine_j_values_test ( );
  normal_01_cdf_values_test ( );
  normal_01_pdf_values_test ( );
  normal_cdf_values_test ( );
  normal_pdf_values_test ( );

  omega_values_test ( );
  owen_values_test ( );

  partition_count_values_test ( );
  partition_distinct_count_values_test ( );
  phi_values_test ( );
  pi_values_test ( );
  poisson_cdf_values_test ( );
  polylogarithm_values_test ( );
  prandtl_values_test ( );
  prime_values_test ( );
  psat_values_test ( );
  psi_values_test ( );

  r8_factorial_values_test ( );
  r8_factorial_log_values_test ( );
  r8_factorial2_values_test ( );
  r8_fall_values_test ( );
  r8_rise_values_test ( );
  rayleigh_cdf_values_test ( );

  scaled_inverse_chi_square_pdf_values_test ( );
  secvir_values_test ( );
  shi_values_test ( );
  si_values_test ( );
  sigma_values_test ( );
  sin_values_test ( );
  sin_degree_values_test ( );
  sin_power_int_values_test ( );
  sinh_values_test ( );
  six_j_values_test ( );
  sound_values_test ( );
  sphere_unit_area_values_test ( );
  sphere_unit_volume_values_test ( );
  spherical_harmonic_values_test ( );
  sqrt_values_test ( );
  stirling1_values_test ( );
  stirling2_values_test ( );
  stromgen_values_test ( );
  struve_h0_values_test ( );
  struve_h1_values_test ( );
  struve_l0_values_test ( );
  struve_l1_values_test ( );
  student_cdf_values_test ( );
  student_noncentral_cdf_values_test ( );
  subfactorial_values_test ( );
  surten_values_test ( );
  synch1_values_test ( );
  synch2_values_test ( );

  tan_values_test ( );
  tanh_values_test ( );
  tau_values_test ( );
  thercon_values_test ( );
  three_j_values_test ( );
  tran02_values_test ( );
  tran03_values_test ( );
  tran04_values_test ( );
  tran05_values_test ( );
  tran06_values_test ( );
  tran07_values_test ( );
  tran08_values_test ( );
  tran09_values_test ( );
  trigamma_values_test ( );
  truncated_normal_ab_cdf_test ( );
  truncated_normal_ab_pdf_test ( );
  truncated_normal_a_cdf_test ( );
  truncated_normal_a_pdf_test ( );
  truncated_normal_b_cdf_test ( );
  truncated_normal_b_pdf_test ( );
  tsat_values_test ( );

  van_der_corput_values_test ( );
  viscosity_values_test ( );
  von_mises_cdf_values_test ( );

  weekday_values_test ( );
  weibull_cdf_values_test ( );

  zeta_values_test ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "TEST_VALUES_PRB:\n";
  std::cout << "  Normal end of execution.\n";
  std::cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void abram0_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ABRAM0_VALUES_TEST tests ABRAM0_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "ABRAM0_VALUES_TEST:\n";
  std::cout << "  ABRAM0_VALUES stores values of \n";
  std::cout << "  the Abramowitz function of order 0.\n";
  std::cout << "\n";
  std::cout << "                X                   ABRAM0(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    abram0_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void abram1_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ABRAM1_VALUES_TEST tests ABRAM1_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "ABRAM1_VALUES_TEST:\n";
  std::cout << "  ABRAM1_VALUES stores values of \n";
  std::cout << "  the Abramowitz function of order 1.\n";
  std::cout << "\n";
  std::cout << "                X                   ABRAM1(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    abram1_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void abram2_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ABRAM2_VALUES_TEST tests ABRAM2_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "ABRAM2_VALUES_TEST:\n";
  std::cout << "  ABRAM2_VALUES stores values of \n";
  std::cout << "  the Abramowitz function of order 2.\n";
  std::cout << "\n";
  std::cout << "                X                   ABRAM3(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    abram2_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void agm_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    AGM_VALUES_TEST tests AGM_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  int n_data;

  std::cout << "\n";
  std::cout << "AGM_VALUES_TEST:\n";
  std::cout << "  AGM_VALUES stores values of \n";
  std::cout << "  the arithmetic geometric mean function.\n";
  std::cout << "\n";
  std::cout << "           A          B              AGM(A,B)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    agm_values ( n_data, a, b, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(14) << std::setprecision (  6 ) << a  << "  "
         << std::setw(14) << std::setprecision (  6 ) << b  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void airy_ai_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    AIRY_AI_VALUES_TEST tests AIRY_AI_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double ai;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "AIRY_AI_VALUES_TEST:\n";
  std::cout << "  AIRY_AI_VALUES stores values of \n";
  std::cout << "  the Airy functions Ai(X).\n";
  std::cout << "\n";
  std::cout << "                X                     Ai(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    airy_ai_values ( n_data, x, ai );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << ai << "\n";
  }
  return;
}
//****************************************************************************80

void airy_ai_int_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    AIRY_AI_INT_VALUES_TEST tests AIRY_AI_INT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "AIRY_AI_INT_VALUES_TEST:\n";
  std::cout << "  AIRY_AI_INT_VALUES stores values of \n";
  std::cout << "  the integral of the Airy Ai function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    airy_ai_int_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void airy_ai_prime_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    AIRY_AI_PRIME_VALUES_TEST tests AIRY_AI_PRIME_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double aip;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "AIRY_AI_PRIME_VALUES_TEST:\n";
  std::cout << "  AIRY_AI_PRIME_VALUES stores values of \n";
  std::cout << "  the derivative of the Airy function Ai'(X).\n";
  std::cout << "\n";
  std::cout << "                X                    Ai'\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    airy_ai_prime_values ( n_data, x, aip );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                           << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x   << "  "
         << std::setw(24) << std::setprecision ( 16 ) << aip << "\n";
  }
  return;
}
//****************************************************************************80

void airy_bi_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    AIRY_BI_VALUES_TEST tests AIRY_BI_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bi;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "AIRY_BI_VALUES_TEST:\n";
  std::cout << "  AIRY_BI_VALUES stores values of \n";
  std::cout << "  the Airy function Bi.\n";
  std::cout << "\n";
  std::cout << "                X                     Bi\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    airy_bi_values ( n_data, x, bi );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << bi << "\n";
  }
  return;
}
//****************************************************************************80

void airy_bi_int_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    AIRY_BI_INT_VALUES_TEST tests AIRY_BI_INT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "AIRY_BI_INT_VALUES_TEST:\n";
  std::cout << "  AIRY_BI_INT_VALUES stores values of \n";
  std::cout << "  the integral of the Airy Bi function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    airy_bi_int_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void airy_bi_prime_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    AIRY_BI_PRIME_VALUES_TEST tests AIRY_BI_PRIME_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bip;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "AIRY_BI_PRIME_VALUES_TEST:\n";
  std::cout << "  AIRY_BI_PRIME_VALUES stores values of \n";
  std::cout << "  the derivative of Airy function Bi'(X).\n";
  std::cout << "\n";
  std::cout << "                X                     Bi'\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    airy_bi_prime_values ( n_data, x, bip );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                           << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x   << "  "
         << std::setw(24) << std::setprecision ( 16 ) << bip << "\n";
  }
  return;
}
//****************************************************************************80

void airy_cai_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    AIRY_CAI_VALUES_TEST tests AIRY_CAI_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  std::complex <double> cai;
  int n_data;
  std::complex <double> x;

  std::cout << "\n";
  std::cout << "AIRY_CAI_VALUES_TEST:\n";
  std::cout << "  AIRY_CAI_VALUES stores values of \n";
  std::cout << "  the Airy functions Ai(X) for std::complex argument.\n";
  std::cout << "\n";
  std::cout << "                X                     Ai\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    airy_cai_values ( n_data, x, cai );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                                   << "  "
         << std::setw(14) << std::setprecision ( 6 ) << std::real ( x )   << "  "
         << std::setw(14) << std::setprecision ( 6 ) << std::imag ( x )   << "  "
         << std::setw(24) << std::setprecision (16 ) << std::real ( cai ) << "  "
         << std::setw(24) << std::setprecision (16 ) << std::imag ( cai ) << "\n";
  }
  return;
}
//****************************************************************************80

void airy_cbi_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    AIRY_CBI_VALUES_TEST tests AIRY_CBI_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  std::complex <double> cbi;
  int n_data;
  std::complex <double> x;

  std::cout << "\n";
  std::cout << "AIRY_CBI_VALUES_TEST:\n";
  std::cout << "  AIRY_CBI_VALUES stores values of \n";
  std::cout << "  the Airy functions Bi(X) for std::complex argument.\n";
  std::cout << "\n";
  std::cout << "                X                     Bi\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    airy_cbi_values ( n_data, x, cbi );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                                   << "  "
         << std::setw(14) << std::setprecision ( 6 ) << std::real ( x )   << "  "
         << std::setw(14) << std::setprecision ( 6 ) << std::imag ( x )   << "  "
         << std::setw(24) << std::setprecision (16 ) << std::real ( cbi ) << "  "
         << std::setw(24) << std::setprecision (16 ) << std::imag ( cbi ) << "\n";
  }
  return;
}
//****************************************************************************80

void airy_gi_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    AIRY_GI_VALUES_TEST tests AIRY_GI_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bip;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "AIRY_GI_VALUES_TEST:\n";
  std::cout << "  AIRY_GI_VALUES stores values of \n";
  std::cout << "  the modified Airy function Gi(X).\n";
  std::cout << "\n";
  std::cout << "                X                     Gi\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    airy_gi_values ( n_data, x, bip );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                           << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x   << "  "
         << std::setw(24) << std::setprecision ( 16 ) << bip << "\n";
  }
  return;
}
//****************************************************************************80

void airy_hi_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    AIRY_HI_VALUES_TEST tests AIRY_HI_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bip;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "AIRY_HI_VALUES_TEST:\n";
  std::cout << "  AIRY_HI_VALUES stores values of \n";
  std::cout << "  the modified Airy function Hi(X).\n";
  std::cout << "\n";
  std::cout << "                X                     Hi\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    airy_hi_values ( n_data, x, bip );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                           << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x   << "  "
         << std::setw(24) << std::setprecision ( 16 ) << bip << "\n";
  }
  return;
}
//****************************************************************************80

void arccos_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ARCCOS_VALUES_TEST tests ARCCOS_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "ARCCOS_VALUES_TEST:\n";
  std::cout << "  ARCCOS_VALUES stores values of the arc cosine function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    arccos_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void arccosh_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ARCCOSH_VALUES_TEST tests ARCCOSH_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "ARCCOSH_VALUES_TEST:\n";
  std::cout << 
    "  ARCCOSH_VALUES stores values of the hyperbolic arc cosine function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    arccosh_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void arcsin_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ARCSIN_VALUES_TEST tests ARCSIN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "ARCSIN_VALUES_TEST:\n";
  std::cout << "  ARCSIN_VALUES stores values of the arc sine function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    arcsin_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void arcsinh_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ARCSINH_VALUES_TEST tests ARCSINH_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "ARCSINH_VALUES_TEST:\n";
  std::cout << 
    "  ARCSINH_VALUES stores values of the hyperbolic arc sine function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    arcsinh_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void arctan_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ARCTAN_VALUES_TEST tests ARCTAN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "ARCTAN_VALUES_TEST:\n";
  std::cout << "  ARCTAN_VALUES stores values of the arc tangent function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    arctan_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void arctan_int_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ARCTAN_INT_VALUES_TEST tests ARCTAN_INT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bip;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "ARCTAN_INT_VALUES_TEST:\n";
  std::cout << "  ARCTAN_INT_VALUES stores values of \n";
  std::cout << "  the arctangent integral.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    arctan_int_values ( n_data, x, bip );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                           << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x   << "  "
         << std::setw(24) << std::setprecision ( 16 ) << bip << "\n";
  }
  return;
}
//****************************************************************************80

void arctanh_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ARCTANH_VALUES_TEST tests ARCTANH_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "ARCTANH_VALUES_TEST:\n";
  std::cout << 
    "  ARCTANH_VALUES stores values of the hyperbolic arc tangent function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    arctanh_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bei0_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BEI0_VALUES_TEST tests BEI0_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BEI0_VALUES_TEST:\n";
  std::cout << "  BEI0_VALUES stores values of \n";
  std::cout << "  the Kelvin function BEI of order 0.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bei0_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bei1_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BEI1_VALUES_TEST tests BEI1_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BEI1_VALUES_TEST:\n";
  std::cout << "  BEI1_VALUES stores values of \n";
  std::cout << "  the Kelvin function BEI of order 1.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bei1_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bell_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BELL_VALUES_TEST tests BELL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "BELL_VALUES_TEST:\n";
  std::cout << "  BELL_VALUES returns values of \n";
  std::cout << "  the Bell numbers.\n";
  std::cout << "\n";
  std::cout << "     N        BELL(N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bell_values ( n_data, n, c );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                  << "  "
         << std::setw(6)  << n << "  "
         << std::setw(10) << c << "\n";
  }
  return;
}
//****************************************************************************80

void ber0_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BER0_VALUES_TEST tests BER0_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BER0_VALUES_TEST:\n";
  std::cout << "  BER0_VALUES stores values of \n";
  std::cout << "  the Kelvin function BER of order 0.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    ber0_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void ber1_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BER1_VALUES_TEST tests BER1_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BER1_VALUES_TEST:\n";
  std::cout << "  BER1_VALUES stores values of \n";
  std::cout << "  the Kelvin function BER of order 1.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    ber1_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bernoulli_number_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_NUMBER_VALUES_TEST tests BERNOULLI_NUMBER_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double c;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "BERNOULLI_NUMBER_VALUES_TEST:\n";
  std::cout << "  BERNOULLI_NUMBER_VALUES returns values of \n";
  std::cout << "  the Bernoulli numbers.\n";
  std::cout << "\n";
  std::cout << "     N              B(N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bernoulli_number_values ( n_data, n, c );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                  << "  "
         << std::setw(6)  << n << "  "
         << std::setw(12) << c << "\n";
  }
  return;
}
//****************************************************************************80

void bernoulli_poly_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_POLY_VALUES_TEST tests BERNOULLI_POLY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double b;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BERNOULLI_POLY_VALUES_TEST:\n";
  std::cout << "  BERNOULLI_POLY_VALUES returns values of \n";
  std::cout << "  the Bernoulli Polynomials.\n";
  std::cout << "\n";
  std::cout << "     N     X      BERNOULLI(N)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bernoulli_poly_values ( n_data, n, x, b );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                  << "  "
         << std::setw(6)  << n << "  "
         << std::setw(12) << x << "  "
         << std::setw(12) << b << "\n";
  }
  return;
}
//****************************************************************************80

void bernstein_poly_01_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BERNSTEIN_POLY_01_VALUES_TEST tests BERNSTEIN_POLY_01_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double b;
  int k;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BERNSTEIN_POLY_01_VALUES_TEST:\n";
  std::cout << "  BERNSTEIN_POLY_01_VALUES returns values of \n";
  std::cout << "  the Bernstein Polynomials.\n";
  std::cout << "\n";
  std::cout << "     N     K       X      BERNSTEIN(N,K)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bernstein_poly_01_values ( n_data, n, k, x, b );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                  << "  "
         << std::setw(6)  << n << "  "
         << std::setw(6)  << k << "  "
         << std::setw(12) << x << "  "
         << std::setw(12) << b << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_i0_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_I0_VALUES_TEST tests BESSEL_I0_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_I0_VALUES_TEST:\n";
  std::cout << "  BESSEL_I0_VALUES stores values of \n";
  std::cout << "  the Bessel I0 function.\n";
  std::cout << "\n";
  std::cout << "      X         I0(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_i0_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_i0_int_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_I0_INT_VALUES_TEST tests BESSEL_I0_INT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bip;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_I0_INT_VALUES_TEST:\n";
  std::cout << "  BESSEL_I0_INT_VALUES stores values of \n";
  std::cout << "  the integral of the Bessel I0 function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_i0_int_values ( n_data, x, bip );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                           << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x   << "  "
         << std::setw(24) << std::setprecision ( 16 ) << bip << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_i0_spherical_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_I0_SPHERICAL_VALUES_TEST tests BESSEL_I0_SPHERICAL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_I0_SPHERICAL_VALUES_TEST:\n";
  std::cout << "  BESSEL_I0_SPHERICAL_VALUES stores values of\n";
  std::cout << "  the spherical Bessel i0 function.\n";
  std::cout << "\n";
  std::cout << "      X            i0(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   bessel_i0_spherical_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_i1_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_I1_VALUES_TEST tests BESSEL_I1_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_I1_VALUES_TEST:\n";
  std::cout << "  BESSEL_I1_VALUES stores values of \n";
  std::cout << "  the Bessel I1 function.\n";
  std::cout << "\n";
  std::cout << "      X         I1(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_i1_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_i1_spherical_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_I1_SPHERICAL_VALUES_TEST tests BESSEL_I1_SPHERICAL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_I1_SPHERICAL_VALUES_TEST:\n";
  std::cout << "  BESSEL_I1_SPHERICAL_VALUES stores values of\n";
  std::cout << "  the spherical Bessel i1 function.\n";
  std::cout << "\n";
  std::cout << "      X            i1(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   bessel_i1_spherical_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_in_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_IN_VALUES_TEST tests BESSEL_IN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_IN_VALUES_TEST:\n";
  std::cout << "  BESSEL_IN_VALUES stores values of \n";
  std::cout << "  the Bessel In function.\n";
  std::cout << "\n";
  std::cout << "      N     X         IN(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_in_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_ix_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_IX_VALUES_TEST tests BESSEL_IX_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double nu;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_IX_VALUES_TEST:\n";
  std::cout << "  BESSEL_IX_VALUES stores values of \n";
  std::cout << "  the Bessel In function for NONINTEGER order.\n";
  std::cout << "\n";
  std::cout << "      NU    X         IN(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_ix_values ( n_data, nu, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << nu << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_j0_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_J0_VALUES_TEST tests BESSEL_J0_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_J0_VALUES_TEST:\n";
  std::cout << "  BESSEL_J0_VALUES stores values of \n";
  std::cout << "  the Bessel J0 function.\n";
  std::cout << "\n";
  std::cout << "      X         J0(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_j0_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_j0_int_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_J0_INT_VALUES_TEST tests BESSEL_J0_INT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bip;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_J0_INT_VALUES_TEST:\n";
  std::cout << "  BESSEL_J0_INT_VALUES stores values of \n";
  std::cout << "  the integral of the Bessel J0 function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_j0_int_values ( n_data, x, bip );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                           << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x   << "  "
         << std::setw(24) << std::setprecision ( 16 ) << bip << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_j0_spherical_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_J0_SPHERICAL_VALUES_TEST tests BESSEL_J0_SPHERICAL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_J0_SPHERICAL_VALUES_TEST:\n";
  std::cout << "  BESSEL_J0_SPHERICAL_VALUES stores values of\n";
  std::cout << "  the spherical Bessel j0 function.\n";
  std::cout << "\n";
  std::cout << "      X            j0(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   bessel_j0_spherical_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_j0_zero_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_J0_ZERO_VALUES_TEST tests BESSEL_J0_ZERO_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 January 2016
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int k;
  int n_data;

  std::cout << "\n";
  std::cout << "BESSEL_J0_ZERO_VALUES_TEST:\n";
  std::cout << "  BESSEL_J0_ZERO_VALUES stores values of zeros of\n";
  std::cout << "  the Bessel J0 function.\n";
  std::cout << "\n";
  std::cout << "       K         X(K)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_j0_zero_values ( n_data, k, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout << "  " << std::setw(6) << k
         << "  " << std::setw(24) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_j1_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_J1_VALUES_TEST tests BESSEL_J1_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_J1_VALUES_TEST:\n";
  std::cout << "  BESSEL_J1_VALUES stores values of \n";
  std::cout << "  the Bessel J1 function.\n";
  std::cout << "\n";
  std::cout << "      X         J1(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_j1_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_j1_spherical_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_J1_SPHERICAL_VALUES_TEST tests BESSEL_J1_SPHERICAL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_J1_SPHERICAL_VALUES_TEST:\n";
  std::cout << "  BESSEL_J1_SPHERICAL_VALUES stores values of\n";
  std::cout << "  the spherical Bessel j1 function.\n";
  std::cout << "\n";
  std::cout << "      X            j1(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   bessel_j1_spherical_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_jn_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_JN_VALUES_TEST tests BESSEL_JN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_JN_VALUES_TEST:\n";
  std::cout << "  BESSEL_JN_VALUES stores values of \n";
  std::cout << "  the Bessel Jn function.\n";
  std::cout << "\n";
  std::cout << "      N     X         JN(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_jn_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_jx_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_JX_VALUES_TEST tests BESSEL_JX_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double nu;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_JX_VALUES_TEST:\n";
  std::cout << "  BESSEL_JX_VALUES stores values of \n";
  std::cout << "  the Bessel Jn function for NONINTEGER order.\n";
  std::cout << "\n";
  std::cout << "      NU      X         JN(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_jx_values ( n_data, nu, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << nu << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_k0_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_K0_VALUES_TEST tests BESSEL_K0_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_K0_VALUES_TEST:\n";
  std::cout << "  BESSEL_K0_VALUES stores values of \n";
  std::cout << "  the Bessel K0 function.\n";
  std::cout << "\n";
  std::cout << "      X         K0(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_k0_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_k0_int_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_K0_INT_VALUES_TEST tests BESSEL_K0_INT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bip;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_K0_INT_VALUES_TEST:\n";
  std::cout << "  BESSEL_K0_INT_VALUES stores values of \n";
  std::cout << "  the integral of the Bessel K0 function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_k0_int_values ( n_data, x, bip );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                           << "  "
         << std: std::setw(24) << std::setprecision ( 16 ) << x   << "  "
         << std::setw(24) << std::setprecision ( 16 ) << bip << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_k1_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_K1_VALUES_TEST tests BESSEL_K1_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_K1_VALUES_TEST:\n";
  std::cout << "  BESSEL_K1_VALUES stores values of \n";
  std::cout << "  the Bessel K1 function.\n";
  std::cout << "\n";
  std::cout << "      X         K1(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_k1_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_kn_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_KN_VALUES_TEST tests BESSEL_KN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_KN_VALUES_TEST:\n";
  std::cout << "  BESSEL_KN_VALUES stores values of \n";
  std::cout << "  the Bessel Kn function.\n";
  std::cout << "\n";
  std::cout << "      N      X         KN(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_kn_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_kx_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_KX_VALUES_TEST tests BESSEL_KX_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double nu;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_KX_VALUES_TEST:\n";
  std::cout << "  BESSEL_KX_VALUES stores values of \n";
  std::cout << "  the Bessel Kn function for NONINTEGER order.\n";
  std::cout << "\n";
  std::cout << "      NU     X         KN(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_kx_values ( n_data, nu, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << nu << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_y0_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_Y0_VALUES_TEST tests BESSEL_Y0_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_Y0_VALUES_TEST:\n";
  std::cout << "  BESSEL_Y0_VALUES stores values of \n";
  std::cout << "  the Bessel Y0 function.\n";
  std::cout << "\n";
  std::cout << "      X         Y0(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_y0_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_y0_int_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_Y0_INT_VALUES_TEST tests BESSEL_Y0_INT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_Y0_INT_VALUES_TEST:\n";
  std::cout << "  BESSEL_Y0_INT_VALUES stores values of \n";
  std::cout << "  the integral of the Bessel Y0 function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_y0_int_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_y0_spherical_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_Y0_SPHERICAL_VALUES_TEST tests BESSEL_Y0_SPHERICAL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_Y0_SPHERICAL_VALUES_TEST:\n";
  std::cout << "  BESSEL_Y0_SPHERICAL_VALUES stores values of\n";
  std::cout << "  the spherical Bessel y0 function.\n";
  std::cout << "\n";
  std::cout << "                X                      y0(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   bessel_y0_spherical_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_y1_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_Y1_VALUES_TEST tests BESSEL_Y1_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_Y1_VALUES_TEST:\n";
  std::cout << "  BESSEL_Y1_VALUES stores values of \n";
  std::cout << "  the Bessel Y1 function.\n";
  std::cout << "\n";
  std::cout << "                X                   Y1(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_y1_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_y1_spherical_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_Y1_SPHERICAL_VALUES_TEST tests BESSEL_Y1_SPHERICAL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_Y1_SPHERICAL_VALUES_TEST:\n";
  std::cout << "  BESSEL_Y1_SPHERICAL_VALUES stores values of\n";
  std::cout << "  the spherical Bessel y1 function.\n";
  std::cout << "\n";
  std::cout << "                X                      y1(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   bessel_y1_spherical_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_yn_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_YN_VALUES_TEST tests BESSEL_YN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_YN_VALUES_TEST:\n";
  std::cout << "  BESSEL_YN_VALUES stores values of \n";
  std::cout << "  the Bessel Yn function.\n";
  std::cout << "\n";
  std::cout << "      N     X         YN(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_yn_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bessel_yx_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_YX_VALUES_TEST tests BESSEL_YX_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double nu;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BESSEL_YX_VALUES_TEST:\n";
  std::cout << "  BESSEL_YX_VALUES stores values of \n";
  std::cout << "  the Bessel Yn function for NONINTEGER order.\n";
  std::cout << "\n";
  std::cout << "      NU    X         YN(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_yx_values ( n_data, nu, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << nu << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void beta_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_CDF_VALUES_TEST tests BETA_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BETA_CDF_VALUES_TEST:\n";
  std::cout << "  BETA_CDF_VALUES stores values of\n";
  std::cout << "  the Beta CDF.\n";
  std::cout << "\n";
  std::cout << "      A            B            X            CDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_cdf_values ( n_data, a, b, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                 << "  "
         << std::setw(12)                     << a  << "  "
         << std::setw(12)                     << b  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void beta_inc_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_INC_VALUES_TEST tests BETA_INC_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BETA_INC_VALUES_TEST:\n";
  std::cout << "  BETA_INC_VALUES stores values of\n";
  std::cout << "  the incomplete Beta function.\n";
  std::cout << "\n";
  std::cout << "      A            B            X            BETA_INC(A,B)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( n_data, a, b, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                 << "  "
         << std::setw(12)                     << a  << "  "
         << std::setw(12)                     << b  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void beta_log_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_LOG_VALUES_TEST tests BETA_LOG_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fxy;
  int n_data;
  double x;
  double y;

  std::cout << "\n";
  std::cout << "BETA_LOG_VALUES_TEST:\n";
  std::cout << "  BETA_LOG_VALUES stores values of\n";
  std::cout << "  the logarithm of the Beta function.\n";
  std::cout << "\n";
  std::cout << "      X              Y         BETA_LOG(X,Y)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_log_values ( n_data, x, y, fxy );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                 << "  "
         << std::setw(12)                     << x   << "  "
         << std::setw(12)                     << y   << "  "
         << std::setw(24) << std::setprecision(16) << fxy << "\n";
  }
  return;
}
//****************************************************************************80

void beta_noncentral_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_NONCENTRAL_CDF_VALUES_TEST tests BETA_NONCENTRAL_CDF_VALUES.
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
  double a;
  double b;
  double fx;
  double lambda;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BETA_NONCENTRAL_CDF_VALUES_TEST:\n";
  std::cout << "  BETA_NONCENTRAL_CDF_VALUES stores values of\n";
  std::cout << "  the noncentral Beta CDF.\n";
  std::cout << "\n";
  std::cout << "      A            B       LAMBDA             X            CDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_noncentral_cdf_values ( n_data, a, b, lambda, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                 << "  "
         << std::setw(12)                     << a  << "  "
         << std::setw(12)                     << b  << "  "
         << std: std::setw(12)                     << lambda << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void beta_pdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_PDF_VALUES_TEST tests BETA_PDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  double beta;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "BETA_PDF_VALUES_TEST:\n";
  std::cout << "  BETA_PDF_VALUES stores values of\n";
  std::cout << "  the Beta PDF.\n";
  std::cout << "\n";
  std::cout << "      ALPHA        BETA         X            PDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_pdf_values ( n_data, alpha, beta, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                 << "  "
         << std::setw(12)                     << alpha  << "  "
         << std::setw(12)                     << beta  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void beta_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_VALUES_TEST tests BETA_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fxy;
  int n_data;
  double x;
  double y;

  std::cout << "\n";
  std::cout << "BETA_VALUES_TEST:\n";
  std::cout << "  BETA_VALUES stores values of\n";
  std::cout << "  the Beta function.\n";
  std::cout << "\n";
  std::cout << "      X              Y         BETA(X,Y)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_values ( n_data, x, y, fxy );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                 << "  "
         << std::setw(12)                     << x   << "  "
         << std::setw(12)                     << y   << "  "
         << std::setw(24) << std::setprecision(16) << fxy << "\n";
  }
  return;
}
//****************************************************************************80

void binomial_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL_VALUES_TEST tests BINOMIAL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int b;
  int c;
  int n_data;

  std::cout << "\n";
  std::cout << "BINOMIAL_VALUES_TEST:\n";
  std::cout << "  BINOMIAL_VALUES returns values of\n";
  std::cout << "  the binomial numbers.\n";
  std::cout << "\n";
  std::cout << "     A     B        C(A,B)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    binomial_values ( n_data, a, b, c );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                  << "  "
         << std::setw(6)  << a << "  "
         << std::setw(6)  << b << "  "
         << std::setw(12) << c << "\n";
  }
  return;
}
//****************************************************************************80

void binomial_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL_CDF_VALUES_TEST tests BINOMIAL_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  double b;
  double fx;
  int n_data;
  int x;

  std::cout << "\n";
  std::cout << "BINOMIAL_CDF_VALUES_TEST:\n";
  std::cout << "  BINOMIAL_CDF_VALUES returns values of \n";
  std::cout << "  the Binomial Cumulative Density Function.\n";
  std::cout << "\n";
  std::cout << "     A      B        X   CDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    binomial_cdf_values ( n_data, a, b, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                  << "  "
         << std::setw(6)                      << a  << "  "
         << std::setw(8)                      << b  << "  "
         << std::setw(4)                      << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void binomial_pdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL_PDF_VALUES_TEST tests BINOMIAL_PDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  double b;
  double fx;
  int n_data;
  int x;

  std::cout << "\n";
  std::cout << "BINOMIAL_PDF_VALUES_TEST:\n";
  std::cout << "  BINOMIAL_PDF_VALUES returns values of \n";
  std::cout << "  the Binomial Probability Density Function.\n";
  std::cout << "\n";
  std::cout << "       A          B                    X       PDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    binomial_pdf_values ( n_data, a, b, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                  << "  "
         << std::setw(6)                      << a  << "  "
         << std::setw(24) << std::setprecision(16) << b  << "  "
         << std::setw(4)                      << x  << "  "
         << std: std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void bivariate_normal_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BIVARIATE_NORMAL_CDF_VALUES_TEST tests BIVARIATE_NORMAL_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2009
//
//  Author:
//
//    John Burkardt
//
{
  double fxy;
  int n_data;
  double r;
  double x;
  double y;

  std::cout << "\n";
  std::cout << "BIVARIATE_NORMAL_CDF_VALUES_TEST:\n";
  std::cout << "  BIVARIATE_NORMAL_CDF_VALUES stores values of\n";
  std::cout << "  the bivariate normal CDF.\n";
  std::cout << "\n";
  std::cout << "      X            Y            R            F(R)(X,Y)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bivariate_normal_cdf_values ( n_data, x, y, r, fxy );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                 << "  "
         << std::setw(12)                     << x   << "  "
         << std::setw(12)                     << y   << "  "
         << std::setw(12)                     << r   << "  "
         << std::setw(24) << std::setprecision(16) << fxy << "\n";
  }
  return;
}
//****************************************************************************80

void catalan_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CATALAN_VALUES_TEST tests CATALAN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "CATALAN_VALUES_TEST:\n";
  std::cout << "  CATALAN_VALUES returns values of \n";
  std::cout << "  the Catalan numbers.\n";
  std::cout << "\n";
  std::cout << "     N        C(N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    catalan_values ( n_data, n, c );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                  << "  "
         << std::setw(6)  << n << "  "
         << std::setw(10) << c << "\n";
  }
  return;
}
//****************************************************************************80

void cauchy_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CAUCHY_CDF_VALUES_TEST tests CAUCHY_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  std::cout << "\n";
  std::cout << "CAUCHY_CDF_VALUES_TEST:\n";
  std::cout << "  CAUCHY_CDF_VALUES returns values of \n";
  std::cout << "  the Cauchy Cumulative Density Function.\n";
  std::cout << "\n";
  std::cout << "     Mu      Sigma        X   CDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cauchy_cdf_values ( n_data, mu, sigma, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                 << "  "
         << std::setw(8)                      << mu    << "  "
         << std::setw(8)                      << sigma << "  "
         << std::setw(8)                      << x     << "  "
         << std::setw(24) << std::setprecision(16) << fx    << "\n";
  }
  return;
}
//****************************************************************************80

void cbrt_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CBRT_VALUES_TEST tests CBRT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "CBRT_VALUES_TEST:\n";
  std::cout << "  CBRT_VALUES stores values of the cube root function.\n";
  std::cout << "\n";
  std::cout << "      X            CBRT(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cbrt_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void cheby_t_poly_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_T_POLY_VALUES_TEST tests CHEBY_T_POLY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "CHEBY_T_POLY_VALUES_TEST:\n";
  std::cout << "  CHEBY_T_POLY_VALUES returns values of\n";
  std::cout << "  the Chebyshev T polynomials.\n";
  std::cout << "\n";
  std::cout << "     N       X      T(N)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cheby_t_poly_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(8)  << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void cheby_t01_poly_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_T01_POLY_VALUES_TEST tests CHEBY_T01_POLY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "CHEBY_T01_POLY_VALUES_TEST:\n";
  std::cout << "  CHEBY_T01_POLY_VALUES returns values of\n";
  std::cout << "  the shifted Chebyshev T polynomials.\n";
  std::cout << "\n";
  std::cout << "     N       X      T01(N)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cheby_t01_poly_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(8)  << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void cheby_u_poly_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_U_POLY_VALUES_TEST tests CHEBY_U_POLY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "CHEBY_U_POLY_VALUES_TEST:\n";
  std::cout << "  CHEBY_U_POLY_VALUES returns values of\n";
  std::cout << "  the Chebyshev U polynomials.\n";
  std::cout << "\n";
  std::cout << "     N       X      U(N)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cheby_u_poly_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(8)  << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void cheby_u01_poly_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_U01_POLY_VALUES_TEST tests CHEBY_U01_POLY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "CHEBY_U01_POLY_VALUES_TEST:\n";
  std::cout << "  CHEBY_U01_POLY_VALUES returns values of\n";
  std::cout << "  the shifted Chebyshev U polynomials.\n";
  std::cout << "\n";
  std::cout << "     N       X      U01(N)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cheby_u01_poly_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(8)  << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void cheby_v_poly_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_V_POLY_VALUES_TEST tests CHEBY_V_POLY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "CHEBY_V_POLY_VALUES_TEST:\n";
  std::cout << "  CHEBY_V_POLY_VALUES returns values of\n";
  std::cout << "  the Chebyshev V polynomials.\n";
  std::cout << "\n";
  std::cout << "     N       X      V(N)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cheby_v_poly_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(8)  << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void cheby_v01_poly_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_V01_POLY_VALUES_TEST tests CHEBY_V01_POLY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "CHEBY_V01_POLY_VALUES_TEST:\n";
  std::cout << "  CHEBY_V01_POLY_VALUES returns values of\n";
  std::cout << "  the shifted Chebyshev V polynomials.\n";
  std::cout << "\n";
  std::cout << "     N       X      V01(N)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cheby_v01_poly_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(8)  << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void cheby_w_poly_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_W_POLY_VALUES_TEST tests CHEBY_W_POLY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "CHEBY_W_POLY_VALUES_TEST:\n";
  std::cout << "  CHEBY_W_POLY_VALUES returns values of\n";
  std::cout << "  the Chebyshev W polynomials.\n";
  std::cout << "\n";
  std::cout << "     N       X      W(N)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cheby_w_poly_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(8)  << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void cheby_w01_poly_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_W01_POLY_VALUES_TEST tests CHEBY_W01_POLY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "CHEBY_W01_POLY_VALUES_TEST:\n";
  std::cout << "  CHEBY_W01_POLY_VALUES returns values of\n";
  std::cout << "  the shifted Chebyshev W polynomials.\n";
  std::cout << "\n";
  std::cout << "     N       X      W01(N)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cheby_w01_poly_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(8)  << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void chi_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_VALUES_TEST tests CHI_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "CHI_VALUES_TEST:\n";
  std::cout << "  CHI_VALUES stores values of\n";
  std::cout << "  the Hyperbolic Cosine Integral function CHI(X).\n";
  std::cout << "\n";
  std::cout << "      X            CHI(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    chi_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void chi_square_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_CDF_VALUES_TEST tests CHI_SQUARE_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "CHI_SQUARE_CDF_VALUES_TEST:\n";
  std::cout << "  CHI_SQUARE_CDF_VALUES returns values of \n";
  std::cout << "  the Chi-Squared Cumulative Density Function.\n";
  std::cout << "\n";
  std::cout << "     N       X    CDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    chi_square_cdf_values ( n_data, a, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << a  << "  "
         << std::setw(8)  << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void chi_square_pdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_PDF_VALUES_TEST tests CHI_SQUARE_PDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2015
//
//  Author:
//
//    John Burkardt
//
{
  double df;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "CHI_SQUARE_PDF_VALUES_TEST:\n";
  std::cout << "  CHI_SQUARE_PDF_VALUES returns values of \n";
  std::cout << "  the Chi-Squared Probability Density Function.\n";
  std::cout << "\n";
  std::cout << "     DF         X    PDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    chi_square_pdf_values ( n_data, df, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(8)  << df  << "  "
         << std::setw(8)  << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void chi_square_noncentral_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_NONCENTRAL_CDF_VALUES_TEST tests CHI_SQUARE_NONCENTRAL_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int df;
  double fx;
  double lambda;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "CHI_SQUARE_NONCENTRAL_CDF_VALUES_TEST:\n";
  std::cout << "  CHI_SQUARE_NONCENTRAL_CDF_VALUES returns values of\n";
  std::cout << "  the noncentral Chi-Squared Cumulative Density Function.\n";
  std::cout << "\n";
  std::cout << "      X      LAMBDA     DF     CDF\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    chi_square_noncentral_cdf_values ( n_data, df, lambda, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                       << "  "
         << std::setw(10) << x      << "  "
         << std: std::setw(8)  << lambda << "  "
         << std::setw(4)  << df     << "  "
         << std::setw(12) << fx     << "\n";
  }
  return;
}
//****************************************************************************80

void ci_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CI_VALUES_TEST tests CI_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "CI_VALUES_TEST:\n";
  std::cout << "  CI_VALUES stores values of\n";
  std::cout << "  the Cosine Integral function CI(X).\n";
  std::cout << "\n";
  std::cout << "      X            CI(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    ci_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void cin_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CIN_VALUES_TEST tests CIN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "CIN_VALUES_TEST:\n";
  std::cout << "  CIN_VALUES stores values of\n";
  std::cout << "  the Cosine Integral function CIN(X).\n";
  std::cout << "\n";
  std::cout << "      X            CIN(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cin_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void cinh_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CINH_VALUES_TEST tests CINH_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "CINH_VALUES_TEST:\n";
  std::cout << "  CINH_VALUES stores values of\n";
  std::cout << "  the Hyperbolic Cosine Integral function CINH(X).\n";
  std::cout << "\n";
  std::cout << "      X            CINH(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cinh_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void clausen_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CLAUSEN_VALUES_TEST tests CLAUSEN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "CLAUSEN_VALUES_TEST:\n";
  std::cout << "  CLAUSEN_VALUES stores values of \n";
  std::cout << "  Clausen's integral function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    clausen_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void clebsch_gordan_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CLEBSCH_GORDAN_VALUES_TEST tests CLEBSCH_GORDAN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double j1;
  double j2;
  double j3;
  double m1;
  double m2;
  double m3;
  int n_data;

  std::cout << "\n";
  std::cout << "CLEBSCH_GORDAN_VALUES_TEST:\n";
  std::cout << "  CLEBSCH_GORDAN_VALUES returns values of\n";
  std::cout << "  the Clebsch Gordan coefficient.\n";
  std::cout << "\n";
  std::cout << "      J1      J2      J3      M1      M2      M3        CG\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    clebsch_gordan_values ( n_data, j1, j2, j3, m1, m2, m3, fx );

    if ( n_data == 0 )
    {
      break;
    }

    std::cout << "  " << std::setw(6) << j1
         << "  " << std::setw(6) << j2
         << "  " << std::setw(6) << j3
         << "  " << std::setw(6) << m1
         << "  " << std: std::setw(6) << m2
         << "  " << std::setw(6) << m3
         << "  " << std::setprecision(16) << std::setw(24) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void collatz_count_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    COLLATZ_COUNT_VALUES_TEST tests COLLATZ_COUNT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 March 2006
//
//  Author:
//
//    John Burkardt
//
{
  int count;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "COLLATZ_COUNT_VALUES_TEST:\n";
  std::cout << "  COLLATZ_COUNT_VALUES returns values of\n";
  std::cout << "  the length of the Collatz sequence that\n";
  std::cout << "  starts at N.\n";
  std::cout << "\n";
  std::cout << "         N      COLLATZ_COUNT(N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    collatz_count_values ( n_data, n, count );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout << "  " << std::setw(8)  << n
         << "  " << std::setw(12) << count << "\n";
  }

  return;
}
//****************************************************************************80

void cos_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    COS_VALUES_TEST tests COS_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "COS_VALUES_TEST:\n";
  std::cout << "   COS_VALUES stores values of the cosine function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cos_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void cos_degree_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    COS_DEGREE_VALUES_TEST tests COS_DEGREE_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "COS_DEGREE_VALUES_TEST:\n";
  std::cout << "   COS_DEGREE_VALUES stores values of the cosine function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cos_degree_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void cos_power_int_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    COS_POWER_INT_VALUES_TEST tests COS_POWER_INT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "COS_POWER_INT_VALUES_TEST:\n";
  std::cout << "  COS_POWER_INT_VALUES returns values of\n";
  std::cout << "  the integral of the N-th power of the cosine function.\n";
  std::cout << "\n";
  std::cout << "         A         B       N        FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   cos_power_int_values ( n_data, a, b, n, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                 << "  "
         << std::setw(8)                      << a  << "  "
         << std::setw(8)                      << b  << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void cosh_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    COSH_VALUES_TEST tests COSH_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "COSH_VALUES_TEST:\n";
  std::cout << "   COSH_VALUES stores values of the hyperbolic cosine function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cosh_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void cot_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    COT_VALUES_TEST tests COT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "COT_VALUES_TEST:\n";
  std::cout << "   COT_VALUES stores values of the cotangent function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cot_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void cp_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CP_VALUES_TEST tests CP_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double cp;
  int n_data;
  double p;
  double tc;

  std::cout << "\n";
  std::cout << "CP_VALUES_TEST:\n";
  std::cout << "  CP_VALUES stores values of\n";
  std::cout << "  the specific heat CP\n";
  std::cout << "  as a function of temperature and pressure.\n";
  std::cout << "\n";
  std::cout << "      T            P            CP(T,P)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cp_values ( n_data, tc, p, cp );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << tc << "  "
         << std::setw(12) << p  << "  "
         << std::setw(12) << cp << "\n";
  }
  return;
}
//****************************************************************************80

void dawson_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    DAWSON_VALUES_TEST tests DAWSON_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "DAWSON_VALUES_TEST:\n";
  std::cout << "  DAWSON_VALUES stores values of\n";
  std::cout << "  Dawson's integral function.\n";
  std::cout << "\n";
  std::cout << "      X          DAWSON(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    dawson_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void debye1_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    DEBYE1_VALUES_TEST tests DEBYE1_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "DEBYE1_VALUES_TEST:\n";
  std::cout << "  DEBYE1_VALUES stores values of \n";
  std::cout << "  the Debye function of order 1.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    debye1_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void debye2_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    DEBYE2_VALUES_TEST tests DEBYE2_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "DEBYE2_VALUES_TEST:\n";
  std::cout << "  DEBYE2_VALUES stores values of \n";
  std::cout << "  the Debye function of order 2.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    debye2_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void debye3_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    DEBYE3_VALUES_TEST tests DEBYE3_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "DEBYE3_VALUES_TEST:\n";
  std::cout << "  DEBYE3_VALUES stores values of \n";
  std::cout << "  the Debye function of order 3.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    debye3_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void debye4_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    DEBYE4_VALUES_TEST tests DEBYE4_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "DEBYE4_VALUES_TEST:\n";
  std::cout << "  DEBYE4_VALUES stores values of \n";
  std::cout << "  the Debye function of order 4.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    debye4_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void dedekind_sum_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    DEDEKIND_SUM_VALUES_TEST tests DEDEKIND_SUM_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 July 2009
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  int n;
  int n_data;
  int p;
  int q;

  std::cout << "\n";
  std::cout << "DEDEKIND_SUM_VALUES_TEST:\n";
  std::cout << "  DEDEKIND_SUM_VALUES stores values of the Dedekind sum\n";
  std::cout << "  (N/D) = Dedekind_Sum(P,Q).\n";
  std::cout << "\n";
  std::cout << "       P       Q       N       D\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    dedekind_sum_values ( n_data, p, q, n, d );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                  << "  "
         << std::setw(6) << p  << "  "
         << std::setw(6) << q  << "  "
         << std::setw(6) << n  << "  "
         << std::setw(6) << d  << "\n";
  }
  return;
}
//****************************************************************************80

void dielectric_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    DIELECTRIC_VALUES_TEST tests DIELECTRIC_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double eps;
  int n_data;
  double p;
  double tc;

  std::cout << "\n";
  std::cout << "DIELECTRIC_VALUES_TEST:\n";
  std::cout << "  DIELECTRIC_VALUES stores values of\n";
  std::cout << "  the dielectric function.\n";
  std::cout << "\n";
  std::cout << "      T           P            EPS(T,P)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    dielectric_values ( n_data, tc, p, eps );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                    << "  "
         << std::setw(12) << tc  << "  "
         << std::setw(12) << p   << "  "
         << std::setw(12) << eps << "\n";
  }
  return;
}
//****************************************************************************80

void dilogarithm_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    DILOGARITHM_VALUES_TEST tests DILOGARITHM_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "DILOGARITHM_VALUES_TEST:\n";
  std::cout << "  DILOGARITHM_VALUES stores values of\n";
  std::cout << "  the dilogarithm function.\n";
  std::cout << "\n";
  std::cout << "      X          DILOGARITHM(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    dilogarithm_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void e1_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    E1_VALUES_TEST tests E1_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "E1_VALUES_TEST:\n";
  std::cout << "  E1_VALUES stores values of\n";
  std::cout << "  the exponential integral function E1(X).\n";
  std::cout << "\n";
  std::cout << "      X          E1(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    e1_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void ei_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    EI_VALUES_TEST tests EI_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "EI_VALUES_TEST:\n";
  std::cout << "  EI_VALUES stores values of\n";
  std::cout << "  the exponential integral function EI(X).\n";
  std::cout << "\n";
  std::cout << "      X          EI(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    ei_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void elliptic_ea_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_EA_VALUES_TEST tests ELLIPTIC_EA_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "ELLIPTIC_EA_VALUES_TEST:\n";
  std::cout << "  ELLIPTIC_EA_VALUES stores values of\n";
  std::cout << "  the complete elliptic integral of the second\n";
  std::cout << "  kind, with parameter angle ALPHA in degrees.\n";
  std::cout << "\n";
  std::cout << "    ALPHA        EA(ALPHA)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    elliptic_ea_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void elliptic_em_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_EM_VALUES_TEST tests ELLIPTIC_EM_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "ELLIPTIC_EM_VALUES_TEST:\n";
  std::cout << "  ELLIPTIC_EM_VALUES stores values of\n";
  std::cout << "  the complete elliptic integral of the second\n";
  std::cout << "  kind, with parameter modulus M.\n";
  std::cout << "\n";
  std::cout << "      M            EM(M)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    elliptic_em_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void elliptic_ka_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_KA_VALUES_TEST tests ELLIPTIC_KA_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "ELLIPTIC_KA_VALUES_TEST:\n";
  std::cout << "  ELLIPTIC_KA_VALUES stores values of\n";
  std::cout << "  the complete elliptic integral of the first\n";
  std::cout << "  kind, with parameter angle ALPHA in degrees.\n";
  std::cout << "\n";
  std::cout << "    ALPHA        KA(ALPHA)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    elliptic_ka_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void elliptic_km_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_KM_VALUES_TEST tests ELLIPTIC_KM_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "ELLIPTIC_KM_VALUES_TEST:\n";
  std::cout << "  ELLIPTIC_KM_VALUES stores values of\n";
  std::cout << "  the complete elliptic integral of the first\n";
  std::cout << "  kind, with parameter modulus M.\n";
  std::cout << "\n";
  std::cout << "      M            KM(M)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    elliptic_km_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void erf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ERF_VALUES_TEST tests ERF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "ERF_VALUES_TEST:\n";
  std::cout << "  ERF_VALUES stores values of\n";
  std::cout << "  the error function ERF(X).\n";
  std::cout << "\n";
  std::cout << "      X          ERF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    erf_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void erfc_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ERFC_VALUES_TEST tests ERFC_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "ERFC_VALUES_TEST:\n";
  std::cout << "  ERFC_VALUES stores values of\n";
  std::cout << "  the complementary error function ERFC(X).\n";
  std::cout << "\n";
  std::cout << "      X          ERFC(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    erfc_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void euler_number_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    EULER_NUMBER_VALUES_TEST tests EULER_NUMBER_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "EULER_NUMBER_VALUES_TEST:\n";
  std::cout << "  EULER_NUMBER_VALUES returns values of\n";
  std::cout << "  the Euler numbers.\n";
  std::cout << "\n";
  std::cout << "     N        EULER_NUMBER(N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    euler_number_values ( n_data, n, c );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                  << "  "
         << std::setw(6)  << n << "  "
         << std::setw(10) << c << "\n";
  }
  return;
}
//****************************************************************************80

void euler_poly_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    EULER_POLY_VALUES_TEST tests EULER_POLY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "EULER_POLY_VALUES_TEST:\n";
  std::cout << "  EULER_POLY_VALUES returns values of\n";
  std::cout << "  the Euler numbers.\n";
  std::cout << "\n";
  std::cout << "     N     X       EULER_POLY(N)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    euler_poly_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(8)  << x  << "  "
         << std::setw(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void exp_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    EXP_VALUES_TEST tests EXP_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "EXP_VALUES_TEST:\n";
  std::cout << "   EXP_VALUES stores values of the exponential function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    exp_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void exp3_int_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    EXP3_INT_VALUES_TEST tests EXP3_INT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "EXP3_INT_VALUES_TEST:\n";
  std::cout << "  EXP3_INT_VALUES stores values of \n";
  std::cout << "  the exponential integral function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    exp3_int_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void exponential_01_pdf_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    EXPONENTIAL_01_PDF_VALUES_TEST tests EXPONENTIAL_01_PDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "EXPONENTIAL_01_PDF_VALUES_TEST:\n";
  std::cout << "  EXPONENTIAL_01_PDF_VALUES stores values of\n";
  std::cout << "  the standard exponential Probability Density Function.\n";
  std::cout << "\n";
  std::cout << "            X                   PDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    exponential_01_pdf_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void exponential_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_CDF_VALUES_TEST tests EXPONENTIAL_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double lambda;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "EXPONENTIAL_CDF_VALUES_TEST:\n";
  std::cout << "  EXPONENTIAL_CDF_VALUES stores values of \n";
  std::cout << "  the exponential CDF.\n";
  std::cout << "\n";
  std::cout << "       LAMBDA         X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    exponential_cdf_values ( n_data, lambda, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                              << "  "
         << std::setw(24) << std::setprecision ( 8 )  << lambda << "  "
         << std::setw(24) << std::setprecision ( 8 )  << x      << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx     << "\n";
  }
  return;
}
//****************************************************************************80

void exponential_pdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_PDF_VALUES_TEST tests EXPONENTIAL_PDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double lambda;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "EXPONENTIAL_PDF_VALUES_TEST:\n";
  std::cout << "  EXPONENTIAL_PDF_VALUES stores values of \n";
  std::cout << "  the exponential PDF.\n";
  std::cout << "\n";
  std::cout << "       LAMBDA         X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    exponential_pdf_values ( n_data, lambda, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                              << "  "
         << std::setw(24) << std::setprecision ( 8 )  << lambda << "  "
         << std::setw(24) << std::setprecision ( 8 )  << x      << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx     << "\n";
  }
  return;
}
//****************************************************************************80

void extreme_values_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    EXTREME_VALUES_CDF_VALUES_TEST tests EXTREME_VALUES_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  double beta;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "EXTREME_VALUES_CDF_VALUES_TEST:\n";
  std::cout << "  EXTREME_VALUES_CDF_VALUES stores values of \n";
  std::cout << "  the extreme values CDF.\n";
  std::cout << "\n";
  std::cout << "        Alpha    Beta        X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    extreme_values_cdf_values ( n_data, alpha, beta, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                              << "  "
         << std::setw(12)                        << alpha  << "  "
         << std::setw(12)                        << beta   << "  "
         << std::setw(12)                        << x      << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx     << "\n";
  }
  return;
}
//****************************************************************************80

void f_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    F_CDF_VALUES_TEST tests F_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int b;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "F_CDF_VALUES_TEST:\n";
  std::cout << "  F_CDF_VALUES stores values of\n";
  std::cout << "  the F cumulative density function.\n";
  std::cout << "\n";
  std::cout << "     A       B            X            CDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    f_cdf_values ( n_data, a, b, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << a  << "  "
         << std::setw(6)  << b  << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void f_noncentral_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    F_NONCENTRAL_CDF_VALUES_TEST tests F_NONCENTRAL_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int b;
  double fx;
  double lambda;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "F_NONCENTRAL_CDF_VALUES_TEST:\n";
  std::cout << "  F_NONCENTRAL_CDF_VALUES stores values of\n";
  std::cout << "  the F cumulative density function.\n";
  std::cout << "\n";
  std::cout << "     A       B            LAMBDA    X            CDF\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    f_noncentral_cdf_values ( n_data, a, b, lambda, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                       << "  "
         << std::setw(6)  << a      << "  "
         << std::setw(6)  << b      << "  "
         << std::setw(8)  << lambda << "  "
         << std::setw(12) << x      << "  "
         << std::setw(12) << fx     << "\n";
  }
  return;
}
//****************************************************************************80

void fresnel_cos_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    FRESNEL_COS_VALUES_TEST tests FRESNEL_COS_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "FRESNEL_COS_VALUES_TEST:\n";
  std::cout << "  FRESNEL_COS_VALUES stores values of\n";
  std::cout << "  the Fresnel cosine integral C(X).\n";
  std::cout << "\n";
  std::cout << "      X           C(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    fresnel_cos_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void fresnel_sin_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    FRESNEL_SIN_VALUES_TEST tests FRESNEL_SIN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "FRESNEL_SIN_VALUES_TEST:\n";
  std::cout << "  FRESNEL_SIN_VALUES stores values of\n";
  std::cout << "  the Fresnel sine integral S(X).\n";
  std::cout << "\n";
  std::cout << "      X           S(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    fresnel_sin_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void frobenius_number_data_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    FROBENIUS_NUMBER_DATA_VALUES_TEST tests FROBENIUS_NUMBER_DATA_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  int *c;
  int f;
  int i;
  int n_data;
  int order;

  std::cout << "\n";
  std::cout << "FROBENIUS_NUMBER_DATA_VALUES_TEST:\n";
  std::cout << "  FROBENIUS_NUMBER_DATA_VALUES returns the corresponding\n";
  std::cout << "  coin denominations.\n";

  n_data = 0;

  for ( ; ; )
  {
    frobenius_number_order_values ( n_data, order );

    if ( n_data == 0 )
    {
      break;
    }

    c = new int[order];

    frobenius_number_data_values ( n_data, order, c, f );

    std::cout << "\n";
    std::cout << "  Order = " << order << "\n";
    for ( i = 0; i < order; i++ )
    {
      std::cout << "  " << std::setw(8) << c[i];
    }
    std::cout << "\n";
    std::cout << "  Frobenius number = " << f << "\n";

    delete [] c;
  }
  return;
}
//****************************************************************************80

void frobenius_number_order_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    FROBENIUS_NUMBER_ORDER_VALUES tests FROBENIUS_NUMBER_ORDER_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n_data;
  int order;

  std::cout << "\n";
  std::cout << "FROBENIUS_NUMBER_ORDER_VALUES_TEST:\n";
  std::cout << "  FROBENIUS_NUMBER_ORDER_VALUES returns the order for\n";
  std::cout << "  a Frobenius problem;\n";
  std::cout << "\n";
  std::cout << "   Problem   ORDER\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    frobenius_number_order_values ( n_data, order );

    if ( n_data == 0 )
    {
      break;
    }

    std::cout << "\n";
    std::cout << "  " << std::setw(4) << n_data
         << "  " << std::setw(4) << order << "\n";
  }
  return;
}
//****************************************************************************80

void frobenius_number_order2_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    FROBENIUS_NUMBER_ORDER2_VALUES_TEST tests FROBENIUS_NUMBER_ORDER2_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c1;
  int c2;
  int f;
  int n_data;

  std::cout << "\n";
  std::cout << "FROBENIUS_NUMBER_ORDER2_VALUES_TEST:\n";
  std::cout << "  FROBENIUS_NUMBER_ORDER2_VALUES returns values of \n";
  std::cout << "  the Frobenius number of order 2.\n";
  std::cout << "\n";
  std::cout << "         C1        C2          F(C1,C2)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    frobenius_number_order2_values ( n_data, c1, c2, f );

    if ( n_data == 0 )
    {
      break;
    }

    std::cout << "  " << std::setw(8) << c1
         << "  " << std::setw(8) << c2
         << "  " << std::setw(8) << f << "\n";
  }

  return;
}
//****************************************************************************80

void gamma_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_VALUES_TEST tests GAMMA_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "GAMMA_VALUES_TEST:\n";
  std::cout << "  GAMMA_VALUES stores values of the Gamma function.\n";
  std::cout << "\n";
  std::cout << "      X            GAMMA(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void gamma_01_pdf_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    GAMMA_01_PDF_VALUES_TEST tests GAMMA_01_PDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "GAMMA_01_PDF_VALUES_TEST:\n";
  std::cout << "  GAMMA_01_PDF_VALUES stores values of\n";
  std::cout << "  the standard Gamma Probability Density Function.\n";
  std::cout << "\n";
  std::cout << "       ALPHA        X                   PDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_01_pdf_values ( n_data, alpha, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                             << "  "
         << std::setw(12)                        << alpha << "  "
         << std::setw(12)                        << x     << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx    << "\n";
  }
  return;
}
//****************************************************************************80

void gamma_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_CDF_VALUES_TEST tests GAMMA_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  std::cout << "\n";
  std::cout << "GAMMA_CDF_VALUES_TEST:\n";
  std::cout << "  GAMMA_CDF_VALUES stores values of\n";
  std::cout << "  the Gamma CDF.\n";
  std::cout << "\n";
  std::cout << "      M    Sigma      X            CDF((X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_cdf_values ( n_data, mu, sigma, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                           << "  "
         << std::setw(12)                     << mu     << "  "
         << std::setw(12)                     << sigma  << "  "
         << std::setw(12)                     << x      << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void gamma_inc_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_INC_VALUES_TEST tests GAMMA_INC_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "GAMMA_INC_VALUES_TEST:\n";
  std::cout << "   GAMMA_INC_VALUES stores values of\n";
  std::cout << "   the incomplete Gamma function.\n";
  std::cout << "\n";
  std::cout << "      A            X            GAMMA_INC(A)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_values ( n_data, a, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(12)                     << a  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void gamma_inc_p_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_INC_P_VALUES_TEST tests GAMMA_INC_P_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "GAMMA_INC_P_VALUES_TEST:\n";
  std::cout << "   GAMMA_INC_P_VALUES stores values of\n";
  std::cout << "   the incomplete Gamma P function.\n";
  std::cout << "\n";
  std::cout << "      A            X            F(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_p_values ( n_data, a, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(12)                     << a  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void gamma_inc_q_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_INC_Q_VALUES_TEST tests GAMMA_INC_Q_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "GAMMA_INC_Q_VALUES_TEST:\n";
  std::cout << "   GAMMA_INC_Q_VALUES stores values of\n";
  std::cout << "   the incomplete Gamma Q function.\n";
  std::cout << "\n";
  std::cout << "      A            X            F(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_q_values ( n_data, a, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(12)                     << a  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void gamma_inc_tricomi_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_INC_TRICOMI_VALUES_TEST tests GAMMA_INC_TRICOMI_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "GAMMA_INC_TRICOMI_VALUES_TEST:\n";
  std::cout << "   GAMMA_INC_TRICOMI_VALUES stores values of\n";
  std::cout << "   the incomplete Tricomi Gamma function.\n";
  std::cout << "\n";
  std::cout << "      A            X            F(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_tricomi_values ( n_data, a, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(12)                     << a  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void gamma_log_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_LOG_VALUES_TEST tests GAMMA_LOG_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "GAMMA_LOG_VALUES_TEST:\n";
  std::cout << "  GAMMA_LOG_VALUES stores values of\n";
  std::cout << "  the logarithm of the Gamma function.\n";
  std::cout << "\n";
  std::cout << "      X            GAMMA_LOG(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_log_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void gamma_pdf_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    GAMMA_PDF_VALUES_TEST tests GAMMA_PDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  double beta;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "GAMMA_PDF_VALUES_TEST:\n";
  std::cout << "  GAMMA_PDF_VALUES stores values of\n";
  std::cout << "  a Gamma Probability Density Function.\n";
  std::cout << "\n";
  std::cout << "       BETA          ALPHA        X                   PDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_pdf_values ( n_data, beta, alpha, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                             << "  "
         << std::setw(12)                        << beta  << "  "
         << std::setw(12)                        << alpha << "  "
         << std::setw(12)                        << x     << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx    << "\n";
  }
  return;
}
//****************************************************************************80

void gegenbauer_poly_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GEGENBAUER_POLY_VALUES_TEST tests GEGENBAUER_POLY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "GEGENBAUER_POLY_VALUES_TEST:\n";
  std::cout << "  GEGENBAUER_POLY_VALUES returns values of\n";
  std::cout << "  the Gegenbauer polynomials.\n";
  std::cout << "\n";
  std::cout << "       N       A       X       G(N,A)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gegenbauer_poly_values ( n_data, n, a, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(10)                     << a  << "  "
         << std::setw(10)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }

  return;
}
//****************************************************************************80

void geometric_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GEOMETRIC_CDF_VALUES_TEST tests GEOMETRIC_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double cdf;
  int n_data;
  double p;
  int x;

  std::cout << "\n";
  std::cout << "GEOMETRIC_CDF_VALUES_TEST:\n";
  std::cout << "  GEOMETRIC_CDF_VALUES stores values of\n";
  std::cout << "  the Geometric Probability Cumulative Density Function.\n";
  std::cout << "\n";
  std::cout << "      X      P       CDF\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    geometric_cdf_values ( n_data, x, p, cdf );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                        << "  "
         << std::setw(6)                      << x   << "  "
         << std::setw(8)                      << p   << "  "
         << std::setw(24) << std::setprecision(16) << cdf << "\n";
  }
  return;
}
//****************************************************************************80

void goodwin_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GOODWIN_VALUES_TEST tests GOODWIN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "GOODWIN_VALUES_TEST:\n";
  std::cout << "  GOODWIN_VALUES stores values of \n";
  std::cout << "  the Goodwin function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    goodwin_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void gud_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GUD_VALUES_TEST tests GUD_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "GUD_VALUES_TEST:\n";
  std::cout << "  GUD_VALUES stores values of\n";
  std::cout << "  the Gudermannian function.\n";
  std::cout << "\n";
  std::cout << "      X            GUD(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gud_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void hermite_function_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_FUNCTION_VALUES_TEST tests HERMITE_FUNCTION_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "HERMITE_FUNCTION_VALUES_TEST\n";
  std::cout << "  HERMITE_FUNCTION_VALUES stores values of\n";
  std::cout << "  the Hermite function.\n";
  std::cout << "\n";
  std::cout << "     N      X            Hf(N,X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    hermite_function_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void hermite_poly_phys_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLY_PHYS_VALUES_TEST tests HERMITE_POLY_PHYS_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "HERMITE_POLY_PHYS_VALUES_TEST\n";
  std::cout << "  HERMITE_POLY_PHYS_VALUES stores values of\n";
  std::cout << "  the physicist's Hermite polynomials.\n";
  std::cout << "\n";
  std::cout << "     N      X            H(N,X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    hermite_poly_phys_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void hermite_poly_prob_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLY_PROB_VALUES_TEST tests HERMITE_POLY_PROB_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "HERMITE_POLY_PROB_VALUES_TEST\n";
  std::cout << "  HERMITE_POLY_PROB_VALUES stores values of\n";
  std::cout << "  the probabilist's Hermite polynomials.\n";
  std::cout << "\n";
  std::cout << "     N      X            He(N,X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    hermite_poly_prob_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void hyper_1f1_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HYPER_1F1_VALUES_TEST tests HYPER_1F1_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "HYPER_1F1_VALUES_TEST:\n";
  std::cout << "  HYPER_1F1_VALUES stores values of\n";
  std::cout << "  the hypergeometric function 1F1.\n";
  std::cout << "\n";
  std::cout << "      A      B      X   Hyper_1F1(A,B,X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    hyper_1f1_values ( n_data, a, b, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(8)                      << a  << "  "
         << std::setw(8)                      << b  << "  "
         << std::setw(8)                      << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void hyper_2f1_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HYPER_2F1_VALUES_TEST tests HYPER_2F1_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 September 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "HYPER_2F1_VALUES_TEST:\n";
  std::cout << "  HYPER_2F1_VALUES stores values of\n";
  std::cout << "  the hypergeometric function 2F1.\n";
  std::cout << "\n";
  std::cout << "      A      B     C      X   Hyper_2F1(A,B,C,X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    hyper_2f1_values ( n_data, a, b, c, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(8)                      << a  << "  "
         << std::setw(8)                      << b  << "  "
         << std::setw(8)                      << c  << "  "
         << std::setw(8)                      << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void hypergeometric_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_CDF_VALUES_TEST tests HYPERGEOMETRIC_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  int pop;
  int sam;
  int succ;
  int x;

  std::cout << "\n";
  std::cout << "HYPERGEOMETRIC_CDF_VALUES_TEST:\n";
  std::cout << "  HYPERGEOMETRIC_CDF_VALUES stores values of\n";
  std::cout << "  the Hypergeometric CDF.\n";
  std::cout << "\n";
  std::cout << "     SAM    SUC   POP     X   HyperCDF(S,S,P)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    hypergeometric_cdf_values ( n_data, sam, succ, pop, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                         << "  "
         << std::setw(8)                      << sam  << "  "
         << std::setw(8)                      << succ << "  "
         << std::setw(8)                      << pop  << "  "
         << std::setw(8)                      << x    << "  "
         << std::setw(24) << std::setprecision(16) << fx   << "\n";
  }
  return;
}
//****************************************************************************80

void hypergeometric_pdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_PDF_VALUES_TEST tests HYPERGEOMETRIC_PDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 January 2008
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  int pop;
  int sam;
  int succ;
  int x;

  std::cout << "\n";
  std::cout << "HYPERGEOMETRIC_PDF_VALUES_TEST:\n";
  std::cout << "  HYPERGEOMETRIC_PDF_VALUES stores values of\n";
  std::cout << "  the Hypergeometric PDF.\n";
  std::cout << "\n";
  std::cout << "     SAM    SUC   POP     X   HyperPDF(S,S,P)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    hypergeometric_pdf_values ( n_data, sam, succ, pop, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                         << "  "
         << std::setw(8)                      << sam  << "  "
         << std::setw(8)                      << succ << "  "
         << std::setw(8)                      << pop  << "  "
         << std::setw(8)                      << x    << "  "
         << std::setw(24) << std::setprecision(16) << fx   << "\n";
  }
  return;
}
//****************************************************************************80

void hypergeometric_u_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_U_VALUES_TEST tests HYPERGEOMETRIC_U_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "HYPERGEOMETRIC_U_VALUES_TEST:\n";
  std::cout << "  HYPERGEOMETRIC_U_VALUES stores values of\n";
  std::cout << "  the hypergeometric function U.\n";
  std::cout << "\n";
  std::cout << "      A      B      X   HyperU(A,B,X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    hypergeometric_u_values ( n_data, a, b, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(8)                      << a  << "  "
         << std::setw(8)                      << b  << "  "
         << std::setw(8)                      << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void i0ml0_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I0ML0_VALUES_TEST tests I0ML0_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "I0ML0_VALUES_TEST:\n";
  std::cout << "  I0ML0_VALUES stores values of \n";
  std::cout << "  the I0-L0 function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    i0ml0_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void i1ml1_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I1ML1_VALUES_TEST tests I1ML1_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "I1ML1_VALUES_TEST:\n";
  std::cout << "  I1ML1_VALUES stores values of \n";
  std::cout << "  the I1-L1 function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    i1ml1_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void i4_factorial_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_FACTORIAL_TEST tests I4_FACTORIAL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  int fn;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "I4_FACTORIAL_TEST:\n";
  std::cout << "   I4_FACTORIAL_VALUES returns values of\n";
  std::cout << "   the factorial function.\n";
  std::cout << "\n";
  std::cout << "      N        Factorial(N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    i4_factorial_values ( n_data, n, fn );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(12) << fn << "\n";
  }
  return;
}
//****************************************************************************80

void i4_factorial2_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_FACTORIAL2_TEST tests I4_FACTORIAL2_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int fn;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "I4_FACTORIAL2_TEST:\n";
  std::cout << "   I4_FACTORIAL2_VALUES return;s values of\n";
  std::cout << "   the double factorial function.\n";
  std::cout << "\n";
  std::cout << "      N         DoubleFactorial(N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    i4_factorial2_values ( n_data, n, fn );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(12) << fn << "\n";
  }
  return;
}
//****************************************************************************80

void i4_fall_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_FALL_VALUES_TEST tests I4_FALL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  int fmn;
  int m;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "I4_FALL_VALUES_TEST:\n";
  std::cout << "  I4_FALL_VALUES returns some exact values\n";
  std::cout << "  of the integer falling factorial function:\n";
  std::cout << "\n";
  std::cout << "     M     N      I4_FALL(M,N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    i4_fall_values ( n_data, m, n, fmn );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << m   << "  "
         << std::setw(6)  << n   << "  "
         << std::setw(12) << fmn << "\n";
  }
  return;
}
//****************************************************************************80

void i4_rise_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_RISE_VALUES_TEST tests I4_RISE_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  int fmn;
  int m;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "I4_RISE_VALUES_TEST:\n";
  std::cout << "  I4_RISE_VALUES returns some exact values\n";
  std::cout << "  of the integer rising factorial function:\n";
  std::cout << "\n";
  std::cout << "     M     N      I4_RISE(M,N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    i4_rise_values ( n_data, m, n, fmn );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << m   << "  "
         << std::setw(6)  << n   << "  "
         << std::setw(12) << fmn << "\n";
  }
  return;
}
//****************************************************************************80

void int_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    INT_VALUES_TEST tests INT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "INT_VALUES_TEST:\n";
  std::cout << "  INT_VALUES stores values of the integer part of a real number.\n";
  std::cout << "\n";
  std::cout << "      X            INT(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    int_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void inverse_chi_square_pdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    INVERSE_CHI_SQUARE_PDF_VALUES_TEST tests INVERSE_CHI_SQUARE_PDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2015
//
//  Author:
//
//    John Burkardt
//
{
  double df;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "INVERSE_CHI_SQUARE_PDF_VALUES_TEST:\n";
  std::cout << "  INVERSE_CHI_SQUARE_PDF_VALUES returns values of \n";
  std::cout << "  the inverse Chi-Square Probability Density Function.\n";
  std::cout << "\n";
  std::cout << "     DF        X    PDF\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    inverse_chi_square_pdf_values ( n_data, df, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(8)  << df  << "  "
         << std::setw(8)  << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void inverse_gamma_pdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    INVERSE_GAMMA_PDF_VALUES_TEST tests INVERSE_GAMMA_PDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2015
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  double beta;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "INVERSE_GAMMA_PDF_VALUES_TEST:\n";
  std::cout << "  INVERSE_GAMMA_PDF_VALUES returns values of \n";
  std::cout << "  the inverse gamma Probability Density Function.\n";
  std::cout << "\n";
  std::cout << "     ALPHA        BETA      X    PDF\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    inverse_gamma_pdf_values ( n_data, alpha, beta, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(8) << alpha << "  "
         << std::setw(8) << beta << "  "
         << std::setw(8) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void jacobi_cn_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_CN_VALUES_TEST tests JACOBI_CN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "JACOBI_CN_VALUES_TEST:\n";
  std::cout << "  JACOBI_CN_VALUES returns values of \n";
  std::cout << "  the Jacobi elliptic CN function.\n";
  std::cout << "\n";
  std::cout << "      A         X       CN(A,X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    jacobi_cn_values ( n_data, a, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(10)                     << a  << "  "
         << std::setw(10)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void jacobi_dn_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_DN_VALUES_TEST tests JACOBI_DN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "JACOBI_DN_VALUES_TEST:\n";
  std::cout << "  JACOBI_DN_VALUES returns values of \n";
  std::cout << "  the Jacobi elliptic DN function.\n";
  std::cout << "\n";
  std::cout << "      A         X       DN(A,X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    jacobi_dn_values ( n_data, a, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                 << "  "
         << std::setw(10)                     << a  << "  "
         << std::setw(10)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void jacobi_poly_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_POLY_VALUES_TEST tests JACOBI_POLY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "JACOBI_POLY_VALUES_TEST:\n";
  std::cout << "  JACOBI_POLY_VALUES returns values of\n";
  std::cout << "  the Jacobi polynomial.\n";
  std::cout << "\n";
  std::cout << "       N         A         B      X       J(N,A,B)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    jacobi_poly_values ( n_data, n, a, b, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(8)                      << a  << "  "
         << std::setw(8)                      << b  << "  "
         << std::setw(10)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }

  return;
}
//****************************************************************************80

void jacobi_sn_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_SN_VALUES_TEST tests JACOBI_SN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "JACOBI_SN_VALUES_TEST:\n";
  std::cout << "  JACOBI_SN_VALUES returns values of \n";
  std::cout << "  the Jacobi elliptic SN function.\n";
  std::cout << "\n";
  std::cout << "      A         X       SN(A,X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    jacobi_sn_values ( n_data, a, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(10)                     << a  << "  "
         << std::setw(10)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void jed_ce_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    JED_CE_VALUES_TEST tests JED_CE_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 May 2006
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  double f;
  double jed;
  int n_data;
  int m;
  int y;

  std::cout << "\n";
  std::cout << "JED_CE_VALUES_TEST:\n";
  std::cout << "  JED_CE_VALUES returns:\n";
  std::cout << "  JED, a Julian Ephemeris Date, and\n";
  std::cout << "  YMDF, the corresponding year, month, day, fraction.\n";
  std::cout << "\n";
  std::cout << "        JED          Y   M   D    F\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    jed_ce_values ( n_data, jed, y, m, d, f );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout << "  " << std::setw(12) << jed
         << "  " << std::setw(6)  << y
         << "  " << std::setw(2)  << m
         << "  " << std::setw(2)  << d
         << "  " << std::setw(6)  << f << "\n";
  }
  return;
}
//****************************************************************************80

void jed_mjd_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    JED_MJD_VALUES_TEST tests JED_MJD_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double jed;
  int n_data;
  double mjd;

  std::cout << "\n";
  std::cout << "JED_MJD_VALUES_TEST:\n";
  std::cout << "  JED_MJD_VALUES returns:\n";
  std::cout << "  JED, a Julian Ephemeris Date, and\n";
  std::cout << "  MJD, the corresponding Modified Julian Day count.\n";
  std::cout << "\n";
  std::cout << "   JED      MJD\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    jed_mjd_values ( n_data, jed, mjd );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                    << "  "
         << std::setw(12) << jed << "  "
         << std::setw(12) << mjd << "\n";
  }
  return;
}
//****************************************************************************80

void jed_rd_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    JED_RD_VALUES_TEST tests JED_RD_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double jed;
  int n_data;
  double rd;

  std::cout << "\n";
  std::cout << "JED_RD_VALUES_TEST:\n";
  std::cout << "  JED_RD_VALUES returns:\n";
  std::cout << "  JED, a Julian Ephemeris Date, and\n";
  std::cout << "  RD, the corresponding Reingold Dershowitz Day count.\n";
  std::cout << "\n";
  std::cout << "   JED      RD\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    jed_rd_values ( n_data, jed, rd );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                    << "  "
         << std::setw(12) << jed << "  "
         << std::setw(12) << rd  << "\n";
  }
  return;
}
//****************************************************************************80

void jed_weekday_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    JED_WEEKDAY_VALUES_TEST tests JED_WEEKDAY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double jed;
  int n_data;
  int weekday;
  std::string weekday_name[7] = {
    "Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", 
    "Friday", "Saturday" };

  std::cout << "\n";
  std::cout << "JED_WEEKDAY_VALUES_TEST:\n";
  std::cout << "  JED_WEEKDAY_VALUES returns Julian Ephemeris Dates \n";
  std::cout << "  (JED) and the corresponding weekday\n";
  std::cout << "\n";
  std::cout << "   JED      #  Weekday\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    jed_weekday_values ( n_data, jed, weekday );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                            << "  "
         << std::setw(12) << jed         << "  "
         << std::setw(1)  << weekday     << "  "
         << weekday_name[weekday-1] << "\n";
  }
  return;
}
//****************************************************************************80

void kei0_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    KEI0_VALUES_TEST tests KEI0_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "KEI0_VALUES_TEST:\n";
  std::cout << "  KEI0_VALUES stores values of \n";
  std::cout << "  the Kelvin function KEI of order 0.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    kei0_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void kei1_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    KEI1_VALUES_TEST tests KEI1_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "KEI1_VALUES_TEST:\n";
  std::cout << "  KEI1_VALUES stores values of \n";
  std::cout << "  the Kelvin function KEI of order 1.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    kei1_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void ker0_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    KER0_VALUES_TEST tests KER0_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "KER0_VALUES_TEST:\n";
  std::cout << "  KER0_VALUES stores values of \n";
  std::cout << "  the Kelvin function KER of order 0.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    ker0_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void ker1_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    KER1_VALUES_TEST tests KER1_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 June 2006
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "KER1_VALUES_TEST:\n";
  std::cout << "  KER1_VALUES stores values of \n";
  std::cout << "  the Kelvin function KER of order 1.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    ker1_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void laguerre_associated_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_ASSOCIATED_VALUES_TEST tests LAGUERRE_ASSOCIATED_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int m;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "LAGUERRE_ASSOCIATED_VALUES_TEST:\n";
  std::cout << "  LAGUERRE_ASSOCIATED_VALUES stores values of\n";
  std::cout << "  the associated Laguerre polynomials.\n";
  std::cout << "\n";
  std::cout << "     N     M    X             L(N,M)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    laguerre_associated_values ( n_data, n, m, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                 << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(6)                      << m  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void laguerre_general_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_GENERAL_VALUES_TEST tests LAGUERRE_GENERAL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "LAGUERRE_GENERAL_VALUES_TEST:\n";
  std::cout << "  LAGUERRE_GENERAL_VALUES stores values of\n";
  std::cout << "  the generalized Laguerre function.\n";
  std::cout << "\n";
  std::cout << "     N     A    X             L(N,A)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    laguerre_general_values ( n_data, n, a, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                 << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(12)                      << a  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void laguerre_polynomial_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_POLYNOMIAL_VALUES_TEST tests LAGUERRE_POLYNOMIAL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "LAGUERRE_POLYNOMIAL_VALUES_TEST:\n";
  std::cout << "  LAGUERRE_POLYNOMIAL_VALUES stores values of \n";
  std::cout << "  the Laguerre polynomials.\n";
  std::cout << "\n";
  std::cout << "     N     X            L(N)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    laguerre_polynomial_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void lambert_w_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LAMBERT_W_VALUES_TEST tests LAMBERT_W_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "LAMBERT_W_VALUES_TEST:\n";
  std::cout << "  LAMBERT_W_VALUES stores values of \n";
  std::cout << "  the Lambert W function.\n";
  std::cout << "\n";
  std::cout << "                X                     W(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    lambert_w_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void laplace_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LAPLACE_CDF_VALUES_TEST tests LAPLACE_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double beta;
  double fx;
  double mu;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "LAPLACE_CDF_VALUES_TEST:\n";
  std::cout << "  LAPLACE_CDF_VALUES returns values of \n";
  std::cout << "  the Laplace Cumulative Density Function.\n";
  std::cout << "\n";
  std::cout << "     Mu      Beta         X   CDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    laplace_cdf_values ( n_data, mu, beta, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                         << "  "
         << std::setw(8)                      << mu   << "  "
         << std::setw(8)                      << beta << "  "
         << std::setw(8)                      << x    << "  "
         << std::setw(24) << std::setprecision(16) << fx   << "\n";
  }
  return;
}
//****************************************************************************80

void legendre_associated_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_ASSOCIATED_VALUES_TEST tests LEGENDRE_ASSOCIATED_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int m;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "LEGENDRE_ASSOCIATED_VALUES_TEST:\n";
  std::cout << "  LEGENDRE_ASSOCIATED_VALUES stores values of\n";
  std::cout << "  the associated Legendre polynomials.\n";
  std::cout << "\n";
  std::cout << "     N     M    X             P(N,M)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    legendre_associated_values ( n_data, n, m, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(6)                      << m  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void legendre_associated_normalized_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_ASSOCIATED_NORMALIZED_VALUES_TEST tests LEGENDRE_ASSOCIATED_NORMALIZED_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2010
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int m;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "LEGENDRE_ASSOCIATED_NORMALIZED_VALUES_TEST:\n";
  std::cout << "  LEGENDRE_ASSOCIATED_NORMALIZED_VALUES stores values of\n";
  std::cout << "  the normalized associated Legendre polynomials.\n";
  std::cout << "\n";
  std::cout << "     N     M    X             P(N,M)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    legendre_associated_normalized_values ( n_data, n, m, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(6)                      << m  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void legendre_associated_normalized_sphere_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES_TEST tests LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int m;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES_TEST:\n";
  std::cout << "  LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES stores values of\n";
  std::cout << "  the associated Legendre polynomials, normalized for the unit sphere.\n";
  std::cout << "\n";
  std::cout << "     N     M    X             P(N,M)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    legendre_associated_normalized_sphere_values ( n_data, n, m, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(6)                      << m  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void legendre_poly_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_POLY_VALUES_TEST tests LEGENDRE_POLY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "LEGENDRE_POLY_VALUES_TEST:\n";
  std::cout << "  LEGENDRE_POLY_VALUES stores values of \n";
  std::cout << "  the Legendre polynomials.\n";
  std::cout << "\n";
  std::cout << "     N    X             P(N)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    legendre_poly_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void legendre_function_q_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_FUNCTION_Q_VALUES_TEST tests LEGENDRE_FUNCTION_Q_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "LEGENDRE_FUNCTION_Q_VALUES_TEST:\n";
  std::cout << "  LEGENDRE_FUNCTION_Q_VALUES stores values of\n";
  std::cout << "  the Legendre Q function.\n";
  std::cout << "\n";
  std::cout << "     N    X             Q(N)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    legendre_function_q_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                 << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void lerch_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LERCH_VALUES_TEST tests LERCH_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  int n_data;
  int s;
  double z;

  std::cout << "\n";
  std::cout << "LERCH_VALUES_TEST:\n";
  std::cout << "  LERCH_VALUES returns values of\n";
  std::cout << "  the Lerch transcendent function.\n";
  std::cout << "\n";
  std::cout << "      Z      S      A      Fx\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    lerch_values ( n_data, z, s, a, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(24) << std::setprecision(16) << z  << "  "
         << std::setw(6)                      << s  << "  "
         << std::setw(12)                     << a  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void lobachevsky_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LOBACHEVSKY_VALUES_TEST tests LOBACHEVSKY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "LOBACHEVSKY_VALUES_TEST:\n";
  std::cout << "  LOBACHEVSKY_VALUES stores values of \n";
  std::cout << "  the Lobachevsky function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    lobachevsky_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(12)                        << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void lobatto_polynomial_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LOBATTO_POLYNOMIAL_VALUES_TEST tests LOBATTO_POLYNOMIAL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2013
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "LOBATTO_POLYNOMIAL_VALUES_TEST:\n";
  std::cout << "  LOBATTO_POLYNOMIAL_VALUES stores values of \n";
  std::cout << "  the completed Lobatto polynomials.\n";
  std::cout << "\n";
  std::cout << "     N    X             Lo(N)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    lobatto_polynomial_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void lobatto_polynomial_derivatives_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LOBATTO_POLYNOMIAL_DERIVATIVES_TEST tests LOBATTO_POLYNOMIAL_DERIVATIVES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 November 2014
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "LOBATTO_POLYNOMIAL_DERIVATIVES_TEST:\n";
  std::cout << "  LOBATTO_POLYNOMIAL_VALUES stores derivatives of \n";
  std::cout << "  the completed Lobatto polynomials.\n";
  std::cout << "\n";
  std::cout << "     N    X             Lo'(N)(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    lobatto_polynomial_derivatives ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void log_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_VALUES_TEST tests LOG_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "LOG_VALUES_TEST:\n";
  std::cout << "   LOG_VALUES stores values of the natural logarithm function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    log_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void log_normal_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_CDF_VALUES_TEST tests LOG_NORMAL_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  std::cout << "\n";
  std::cout << "LOG_NORMAL_CDF_VALUES_TEST:\n";
  std::cout << "  LOG_NORMAL_CDF_VALUES returns values of \n";
  std::cout << "  the Log Normal Cumulative Density Function.\n";
  std::cout << "\n";
  std::cout << "     Mu      Sigma        X   CDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    log_normal_cdf_values ( n_data, mu, sigma, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                 << "  "
         << std::setw(8)                      << mu    << "  "
         << std::setw(8)                      << sigma << "  "
         << std::setw(8)                      << x     << "  "
         << std::setw(24) << std::setprecision(16) << fx    << "\n";
  }
  return;
}
//****************************************************************************80

void log_series_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_SERIES_CDF_VALUES_TEST tests LOG_SERIES_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double t;

  std::cout << "\n";
  std::cout << "LOG_SERIES_CDF_VALUES_TEST:\n";
  std::cout << "  LOG_SERIES_CDF_VALUES returns values of \n";
  std::cout << "  the Log Series Cumulative Density Function.\n";
  std::cout << "\n";
  std::cout << "     T      N   CDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    log_series_cdf_values ( n_data, t, n, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(24) << std::setprecision(16) << t  << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void log10_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LOG10_VALUES_TEST tests LOG10_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "LOG10_VALUES_TEST:\n";
  std::cout << "   LOG10_VALUES stores values of the base 10 logarithm function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    log10_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void logarithmic_integral_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LOGARITHMIC_INTEGRAL_VALUES_TEST tests LOGARITHMIC_INTEGRAL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "LOGARITHMIC_INTEGRAL_VALUES_TEST:\n";
  std::cout << "  LOGARITHMIC_INTEGAL_VALUES stores values of\n";
  std::cout << "  the logarithmic integral function.\n";
  std::cout << "\n";
  std::cout << "      X            LI(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    logarithmic_integral_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                 << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void logistic_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LOGISTIC_CDF_VALUES_TEST tests LOGISTIC_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double beta;
  double fx;
  double mu;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "LOGISTIC_CDF_VALUES_TEST:\n";
  std::cout << "  LOGISTIC_CDF_VALUES returns values of \n";
  std::cout << "  the Logistic Cumulative Density Function.\n";
  std::cout << "\n";
  std::cout << "     Mu      Beta         X   CDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    logistic_cdf_values ( n_data, mu, beta, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                 << "  "
         << std::setw(8)                      << mu   << "  "
         << std::setw(8)                      << beta << "  "
         << std::setw(8)                      << x    << "  "
         << std::setw(24) << std::setprecision(16) << fx   << "\n";
  }
  return;
}
//****************************************************************************80

void mertens_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    MERTENS_VALUES_TEST tests MERTENS_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 October 2007
//
//  Author:
//
//    John Burkardt
//
{
  int fn;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "MERTENS_VALUES_TEST:\n";
  std::cout << "  MERTENS_VALUES returns values of\n";
  std::cout << "  the Mertens function.\n";
  std::cout << "\n";
  std::cout << "     N         MERTENS(N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    mertens_values ( n_data, n, fn );

    if ( n_data == 0 )
    {
      break;
    } 
    std::cout                   << "  "
         << std::setw(8)  << n  << "  "
         << std::setw(12) << fn << "\n";
  }
  return;
}
//****************************************************************************80

void moebius_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    MOEBIUS_VALUES_TEST tests MOEBIUS_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int fn;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "MOEBIUS_VALUES_TEST:\n";
  std::cout << "  MOEBIUS_VALUES returns values of\n";
  std::cout << "  the Moebius function.\n";
  std::cout << "\n";
  std::cout << "     N         MU(N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    moebius_values ( n_data, n, fn );

    if ( n_data == 0 )
    {
      break;
    } 
    std::cout                   << "  "
         << std::setw(8)  << n  << "  "
         << std::setw(12) << fn << "\n";
  }
  return;
}
//****************************************************************************80

void multinomial_pdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MULTINOMIAL_PDF_VALUES_TEST tests MULTINOMIAL_PDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int m;
  int n;
  int n_data1;
  int n_data2;
  double *p;
  double pdf;
  int *x;

  std::cout << "\n";
  std::cout << "MULTINOMIAL_PDF_VALUES_TEST:\n";
  std::cout << "  MULTINOMIAL_PDF_VALUES stores values of the Multinomial PDF.\n";
  std::cout << "  Given M possible outcomes on a single trial,\n";
  std::cout << "  with each outcome having probability P,\n";
  std::cout << "  PDF is the probability that after N trials,\n";
  std::cout << "  outcome I occurred X(I) times.\n";
  std::cout << "\n";
  std::cout << "     N     M     I      P        X        PDF()\n";

  n_data1 = 0;
  n_data2 = 0;

  for ( ; ; )
  {
    multinomial_pdf_sizes ( n_data1, m );

    if ( n_data1 == 0 )
    {
      break;
    }

    p = new double[m];
    x = new int[m];
    multinomial_pdf_values ( n_data2, m, n, p, x, pdf );

    std::cout << "\n";
    for ( i = 0; i < m; i++ )
    {
      std::cout << "              " << std::setw(4) << i
           << "  " << std::setw(8) << p[i]
           << "  " << std::setw(4) << x[i] << "\n";
    }
    std::cout << "  " << std::setw(4) << n
         << "  " << std::setw(4) << m
         << "                        " << std::setw(14) << pdf << "\n";

    delete [] p;
    delete [] x;
  }
  return;
}
//****************************************************************************80

void negative_binomial_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NEGATIVE_BINOMIAL_CDF_VALUES_TEST tests NEGATIVE_BINOMIAL_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double cdf;
  int f;
  int n_data;
  double p;
  int s;

  std::cout << "\n";
  std::cout << "NEGATIVE_BINOMIAL_CDF_VALUES_TEST:\n";
  std::cout << "  NEGATIVE_BINOMIAL_CDF_VALUES stores values of\n";
  std::cout << "  the Negative Binomial Cumulative Density Function.\n";
  std::cout << "\n";
  std::cout << "     F     S         P         CDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    negative_binomial_cdf_values ( n_data, f, s, p, cdf );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                        << "  "
         << std::setw(6)                      << f   << "  "
         << std::setw(6)                      << s   << "  "
         << std::setw(12)                     << p   << "  "
         << std::setw(24) << std::setprecision(16) << cdf << "\n";
  }
  return;
}
//****************************************************************************80

void nine_j_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NINE_J_VALUES_TEST demonstrates NINE_J_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double j1;
  double j2;
  double j3;
  double j4;
  double j5;
  double j6;
  double j7;
  double j8;
  double j9;
  int n_data;

  std::cout << "\n";
  std::cout << "NINE_J_VALUES_TEST:\n";
  std::cout << "  NINE_J_VALUES returns values of\n";
  std::cout << "  the Wigner 9J coefficient.\n";
  std::cout << "\n";
  std::cout << "      J1      J2      J3      J4      J5      J6"
       << "      J7      J8      J9        NINE_J\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    nine_j_values ( n_data, j1, j2, j3, j4, j5, j6, j7, j8, j9, fx );

    if ( n_data == 0 )
    {
      break;
    }

    std::cout << "  " << std::setw(6) << j1
         << "  " << std::setw(6) << j2
         << "  " << std::setw(6) << j3
         << "  " << std::setw(6) << j4
         << "  " << std::setw(6) << j5
         << "  " << std::setw(6) << j6
         << "  " << std::setw(6) << j7
         << "  " << std::setw(6) << j8
         << "  " << std::setw(6) << j9
         << "  " << std::setprecision(16) << std::setw(24) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void normal_01_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    NORMAL_01_CDF_VALUES_TEST tests NORMAL_01_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "NORMAL_01_CDF_VALUES_TEST:\n";
  std::cout << "  NORMAL_01_CDF_VALUES stores values of\n";
  std::cout << "  the Normal 01 Cumulative Density Function.\n";
  std::cout << "\n";
  std::cout << "            X                   CDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void normal_01_pdf_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    NORMAL_01_PDF_VALUES_TEST tests NORMAL_01_PDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "NORMAL_01_PDF_VALUES_TEST:\n";
  std::cout << "  NORMAL_01_PDF_VALUES stores values of\n";
  std::cout << "  the Normal 01 Probability Density Function.\n";
  std::cout << "\n";
  std::cout << "            X                   PDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    normal_01_pdf_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void normal_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    NORMAL_CDF_VALUES_TEST tests NORMAL_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  std::cout << "\n";
  std::cout << "NORMAL_CDF_VALUES_TEST:\n";
  std::cout << "  NORMAL_CDF_VALUES stores values of\n";
  std::cout << "  the Normal Cumulative Density Function.\n";
  std::cout << "\n";
  std::cout << "            X                   CDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    normal_cdf_values ( n_data, mu, sigma, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                             << "  "
         << std::setw(12)                        << mu    << "  "
         << std::setw(12)                        << sigma << "  "
         << std::setw(12)                        << x     << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx    << "\n";
  }
  return;
}
//****************************************************************************80

void normal_pdf_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    NORMAL_PDF_VALUES_TEST tests NORMAL_PDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 July 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  std::cout << "\n";
  std::cout << "NORMAL_PDF_VALUES_TEST:\n";
  std::cout << "  NORMAL_PDF_VALUES stores values of\n";
  std::cout << "  the Normal Probability Density Function.\n";
  std::cout << "\n";
  std::cout << "            X                   PDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    normal_pdf_values ( n_data, mu, sigma, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                             << "  "
         << std::setw(12)                        << mu    << "  "
         << std::setw(12)                        << sigma << "  "
         << std::setw(12)                        << x     << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx    << "\n";
  }
  return;
}
//****************************************************************************80

void omega_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    OMEGA_VALUES_TEST tests OMEGA_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int fn;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "OMEGA_VALUES_TEST:\n";
  std::cout << "  OMEGA_VALUES returns values of\n";
  std::cout << "  the Omega function.\n";
  std::cout << "\n";
  std::cout << "     N           OMEGA(N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    omega_values ( n_data, n, fn );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << n  << "  "
         << std::setw(10) << fn << "\n";
  }
  return;
}
//****************************************************************************80

void owen_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    OWEN_VALUES_TEST tests OWEN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double h;
  int n_data;
  double t;

  std::cout << "\n";
  std::cout << "OWEN_VALUES_TEST\n";
  std::cout << "  OWEN_VALUES stores values of\n";
  std::cout << "  Owen's T function.\n";
  std::cout << "\n";
  std::cout << "          H            A            T\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    owen_values ( n_data, h, a, t );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                      << "  "
         << std::setw(12)                     << h << "  "
         << std::setw(12)                     << a << "  "
         << std::setw(24) << std::setprecision(16) << t << "\n";
  }
  return;
}
//****************************************************************************80

void partition_count_values_test ( )
 
//****************************************************************************80
//
//  Purpose:
//
//    PARTITION_COUNT_VALUES_TEST tests PARTITION_COUNT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int fn;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "PARTITION_COUNT_VALUES_TEST:\n";
  std::cout << "  PARTITION_COUNT_VALUES returns values of \n";
  std::cout << "  the integer partition count function.\n";
  std::cout << "\n";
  std::cout << "     N         P(N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    partition_count_values ( n_data, n, fn );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(12) << fn << "\n";
  }
  return;
}
//****************************************************************************80

void partition_distinct_count_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    PARTITION_DISTINCT_COUNT_VALUES_TEST tests PARTITION_DISTINCT_COUNT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int fn;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "PARTITION_DISTINCT_COUNT_VALUES_TEST:\n";
  std::cout << "  PARTITION_DISTINCT_COUNT_VALUES returns values of \n";
  std::cout << "  the integer distinct partition count function.\n";
  std::cout << "\n";
  std::cout << "     N         Q(N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    partition_distinct_count_values ( n_data, n, fn );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(12) << fn << "\n";
  }
  return;
}
//****************************************************************************80

void phi_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    PHI_VALUES_TEST tests PHI_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int fn;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "PHI_VALUES_TEST:\n";
  std::cout << "  PHI_VALUES returns values of\n";
  std::cout << "  the PHI function.\n";
  std::cout << "\n";
  std::cout << "     N         PHI(N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    phi_values ( n_data, n, fn );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(10) << fn << "\n";
  }
  return;
}
//****************************************************************************80

void pi_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    PI_VALUES_TEST tests PI_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int fn;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "PI_VALUES_TEST:\n";
  std::cout << "  PI_VALUES returns values of\n";
  std::cout << "  the PI function.\n";
  std::cout << "\n";
  std::cout << "     N         PI(N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    pi_values ( n_data, n, fn );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << n  << "  "
         << std::setw(10) << fn << "\n";
  }
  return;
}
//****************************************************************************80

void poisson_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    POISSON_CDF_VALUES_TEST tests POISSON_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  int n_data;
  int x;

  std::cout << "\n";
  std::cout << "POISSON_CDF_VALUES_TEST:\n";
  std::cout << "  POISSON_CDF_VALUES returns values of\n";
  std::cout << "  the Poisson Cumulative Density Function.\n";
  std::cout << "\n";
  std::cout << "      A     X       CDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    poisson_cdf_values ( n_data, a, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(8)                      << a  << "  "
         << std::setw(4)                      << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void polylogarithm_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLYLOGARITHM_VALUES_TEST tests POLYLOGARITHM_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n;
  int n_data;
  double z;

  std::cout << "\n";
  std::cout << "POLYLOGARITHM_VALUES_TEST:\n";
  std::cout << "  POLYLOGARITHM_VALUES returns values of \n";
  std::cout << "  the polylogarithm function.\n";
  std::cout << "\n";
  std::cout << "     N      Z          Fx\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    polylogarithm_values ( n_data, n, z, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(24) << std::setprecision(16) << z  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void prandtl_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    PRANDTL_VALUES_TEST tests PRANDTL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int n_data;
  double p;
  double pr;
  double tc;

  std::cout << "\n";
  std::cout << "PRANDTL_VALUES_TEST:\n";
  std::cout << "  PRANDTL_VALUES stores values of\n";
  std::cout << "  the Prandtl number of water\n";
  std::cout << "  as a function of temperature and pressure.\n";
  std::cout << "\n";
  std::cout << "      T            P            Pr(T,P)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    prandtl_values ( n_data, tc, p, pr );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << tc << "  "
         << std::setw(12) << p  << "  "
         << std::setw(12) << pr << "\n";
  }
  return;
}
//****************************************************************************80

void prime_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    PRIME_VALUES_TEST tests PRIME_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int n_data;
  int p;

  std::cout << "\n";
  std::cout << "PRIME_VALUES_TEST:\n";
  std::cout << "  PRIME_VALUES returns values of\n";
  std::cout << "  the prime function.\n";
  std::cout << "\n";
  std::cout << "           N          P[N]\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    prime_values ( n_data, n, p );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                  << "  "
         << std::setw(12) << n << "  "
         << std::setw(12) << p << "\n";
  }

  return;
}
//****************************************************************************80

void psat_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    PSAT_VALUES_TEST tests PSAT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int n_data;
  double psat;
  double tc;

  std::cout << "\n";
  std::cout << "PSAT_VALUES_TEST:\n";
  std::cout << "  PSAT_VALUES stores values of\n";
  std::cout << "  the saturation pressure of water\n";
  std::cout << "  as a function of temperature.\n";
  std::cout << "\n";
  std::cout << "      T            PSAT(T)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    psat_values ( n_data, tc, psat );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                         << "  "
         << std::setw(24) << std::setprecision(16) << tc   << "  "
         << std::setw(24) << std::setprecision(16) << psat << "\n";
  }
  return;
}
//****************************************************************************80

void psi_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    PSI_VALUES_TEST tests PSI_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "PSI_VALUES_TEST\n";
  std::cout << "  PSI_VALUES stores values of\n";
  std::cout << "  the PSI function.\n";
  std::cout << "\n";
  std::cout << "      X            PSI(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    psi_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void r8_factorial_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    R8_FACTORIAL_VALUES_TEST tests R8_FACTORIAL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fn;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "R8_FACTORIAL_VALUES_TEST:\n";
  std::cout << "  R8_FACTORIAL_VALUES stores values of\n";
  std::cout << "  the factorial function (using double arithmetic).\n";
  std::cout << "\n";
  std::cout << "      N       Factorial(N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial_values ( n_data, n, fn );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(24) << std::setprecision(16) << fn << "\n";
  }
  return;
}
//****************************************************************************80

void r8_factorial_log_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL_LOG_VALUES_TEST tests R8_FACTORIAL_LOG_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fn;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "R8_FACTORIAL_LOG_VALUES_TEST:\n";
  std::cout << "  R8_FACTORIAL_LOG_VALUES stores values of\n";
  std::cout << "  the logarithm of the factorial function\n";
  std::cout << "  (using real arithmetic).\n";
  std::cout << "\n";
  std::cout << "      N       Log(Factorial(N))\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial_log_values ( n_data, n, fn );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(24) << std::setprecision(16) << fn << "\n";
  }
  return;
}
//****************************************************************************80

void r8_factorial2_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    R8_FACTORIAL2_VALUES_TEST tests R8_FACTORIAL2_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  double f;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "R8_FACTORIAL2_VALUES_TEST:\n";
  std::cout << "  R8_FACTORIAL2_VALUES stores values of\n";
  std::cout << "  the double factorial function (using double arithmetic).\n";
  std::cout << "\n";
  std::cout << "      N       F\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial2_values ( n_data, n, f );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(24) << std::setprecision(16) << f << "\n";
  }
  return;
}
//****************************************************************************80

void r8_fall_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FALL_VALUES_TEST tests R8_FALL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  double f;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "R8_FALL_VALUES_TEST:\n";
  std::cout << "  R8_FALL_VALUES returns some exact values\n";
  std::cout << "  of the falling factorial function:\n";
  std::cout << "\n";
  std::cout << "     X     N      R8_FALL(X,N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    r8_fall_values ( n_data, x, n, f );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(8)  << x << "  "
         << std::setw(8)  << n << "  "
         << std::setw(12) << f << "\n";
  }
  return;
}
//****************************************************************************80

void r8_rise_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_RISE_VALUES_TEST tests R8_RISE_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  double f;
  int n;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "R8_RISE_VALUES_TEST:\n";
  std::cout << "  R8_RISE_VALUES returns some exact values\n";
  std::cout << "  of the  rising factorial function:\n";
  std::cout << "\n";
  std::cout << "     X     N      R8_RISE(X,N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    r8_rise_values ( n_data, x, n, f );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(8)  << x << "  "
         << std::setw(8)  << n << "  "
         << std::setw(12) << f << "\n";
  }
  return;
}
//****************************************************************************80

void rayleigh_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    RAYLEIGH_CDF_VALUES_TEST tests RAYLEIGH_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double sigma;
  double x;

  std::cout << "\n";
  std::cout << "RAYLEIGH_CDF_VALUES_TEST:\n";
  std::cout << "  RAYLEIGH_CDF_VALUES stores values of\n";
  std::cout << "  the Rayleigh CDF.\n";
  std::cout << "\n";
  std::cout << "      SIGMA        X            CDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    rayleigh_cdf_values ( n_data, sigma, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                 << "  "
         << std::setw(12)                     << sigma  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void scaled_inverse_chi_square_pdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SCALED_INVERSE_CHI_SQUARE_PDF_VALUES_TEST tests SCALED_INVERSE_CHI_SQUARE_PDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2015
//
//  Author:
//
//    John Burkardt
//
{
  double df;
  double fx;
  int n_data;
  double x;
  double xi;

  std::cout << "\n";
  std::cout << "SCALED_INVERSE_CHI_SQUARE_PDF_VALUES_TEST:\n";
  std::cout << "  SCALED_INVERSE_CHI_SQUARE_PDF_VALUES returns values of \n";
  std::cout << "  the scaled inverse Chi-Square Probability Density Function.\n";
  std::cout << "\n";
  std::cout << "     DF        XI        X    PDF\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    scaled_inverse_chi_square_pdf_values ( n_data, df, xi, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(8) << df << "  "
         << std::setw(8) << xi << "  "
         << std::setw(8) << x  << "  "
         << std::setw(12) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void secvir_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SECVIR_VALUES_TEST tests SECVIR_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int n_data;
  double tc;
  double vir;

  std::cout << "\n";
  std::cout << "SECVIR_VALUES_TEST:\n";
  std::cout << "  SECVIR_VALUES stores values of\n";
  std::cout << "  the second virial coefficient of water\n";
  std::cout << "  as a function of temperature.\n";
  std::cout << "\n";
  std::cout << "      T            VIR(T)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   secvir_values ( n_data, tc, vir );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                    << "  "
         << std::setw(12) << tc  << "  "
         << std::setw(12) << vir << "\n";
  }
  return;
}
//****************************************************************************80

void shi_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SHI_VALUES_TEST tests SHI_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "SHI_VALUES_TEST:\n";
  std::cout << "  SHI_VALUES stores values of\n";
  std::cout << "  the hyperbolic sine integral function.\n";
  std::cout << "\n";
  std::cout << "      X            SHI(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   shi_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void si_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SI_VALUES_TEST tests SI_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "SI_VALUES_TEST:\n";
  std::cout << "  SI_VALUES stores values of\n";
  std::cout << "  the sine integral function.\n";
  std::cout << "\n";
  std::cout << "      X            SI(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   si_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void sigma_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SIGMA_VALUES_TEST tests SIGMA_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int fn;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "SIGMA_VALUES_TEST:\n";
  std::cout << "  SIGMA_VALUES returns values of\n";
  std::cout << "  the SIGMA function.\n";
  std::cout << "\n";
  std::cout << "       N         SIGMA(N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   sigma_values ( n_data, n, fn );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(12) << fn << "\n";
  }
  return;
}
//****************************************************************************80

void sin_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SIN_VALUES_TEST tests SIN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "SIN_VALUES_TEST:\n";
  std::cout << "   SIN_VALUES stores values of the sine function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    sin_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void sin_degree_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SIN_DEGREE_VALUES_TEST tests SIN_DEGREE_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "SIN_DEGREE_VALUES_TEST:\n";
  std::cout << "   SIN_DEGREE_VALUES stores values of the sine function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    sin_degree_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void sin_power_int_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SIN_POWER_INT_VALUES_TEST tests SIN_POWER_INT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "SIN_POWER_INT_VALUES_TEST:\n";
  std::cout << "  SIN_POWER_INT_VALUES returns values of\n";
  std::cout << "  the integral of the N-th power of the sine function.\n";
  std::cout << "\n";
  std::cout << "         A         B       N        FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   sin_power_int_values ( n_data, a, b, n, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                 << "  "
         << std::setw(8)                      << a  << "  "
         << std::setw(8)                      << b  << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void sinh_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SINH_VALUES_TEST tests SINH_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "SINH_VALUES_TEST:\n";
  std::cout << "   SINH_VALUES stores values of the hyperbolic sine function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    sinh_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void six_j_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SIX_J_VALUES_TEST tests SIX_J_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double j1;
  double j2;
  double j3;
  double j4;
  double j5;
  double j6;
  int n_data;

  std::cout << "\n";
  std::cout << "SIX_J_VALUES_TEST:\n";
  std::cout << "  SIX_J_VALUES returns values of \n";
  std::cout << "  the Wigner 6J coefficient.\n";
  std::cout << "\n";
  std::cout << 
    "      J1      J2      J3      J4      J5      J6        SIX_J\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    six_j_values ( n_data, j1, j2, j3, j4, j5, j6, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout << "  " << std::setw(6) << j1
         << "  " << std::setw(6) << j2
         << "  " << std::setw(6) << j3
         << "  " << std::setw(6) << j4
         << "  " << std::setw(6) << j5
         << "  " << std::setw(6) << j6
         << "  " << std::setprecision(16) << std::setw(24) << fx << "\n";
  }

  return;
}
//****************************************************************************80

void sound_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SOUND_VALUES_TEST tests SOUND_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double c;
  int n_data;
  double p;
  double tc;

  std::cout << "\n";
  std::cout << "SOUND_VALUES_TEST:\n";
  std::cout << "  SOUND_VALUES stores values of\n";
  std::cout << "  the spead of sound in water\n";
  std::cout << "  as a function of temperature and pressure.\n";
  std::cout << "\n";
  std::cout << "      T            P            C(T,P)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   sound_values ( n_data, tc, p, c );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << tc << "  "
         << std::setw(12) << p  << "  "
         << std::setw(12) << c  << "\n";
  }
  return;
}
//****************************************************************************80

void sphere_unit_area_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_AREA_VALUES_TEST tests SPHERE_UNIT_AREA_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2007
//
//  Author:
//
//    John Burkardt
//
{  
  double fx;
  int n_data;
  int n;

  std::cout << "\n";
  std::cout << "SPHERE_UNIT_AREA_VALUES_TEST:\n";
  std::cout << "  SPHERE_UNIT_AREA_VALUES stores values of\n";
  std::cout << "  the area of the unit sphere in various dimensions.\n";
  std::cout << "\n";
  std::cout << "      N           AREA\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    sphere_unit_area_values ( n_data, n, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void sphere_unit_volume_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_VOLUME_VALUES_TEST tests SPHERE_UNIT_VOLUME_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2007
//
//  Author:
//
//    John Burkardt
//
{  
  double fx;
  int n_data;
  int n;

  std::cout << "\n";
  std::cout << "SPHERE_UNIT_VOLUME_VALUES_TEST:\n";
  std::cout << "  SPHERE_UNIT_VOLUME_VALUES stores values of\n";
  std::cout << "  the volume of the unit sphere in various dimensions.\n";
  std::cout << "\n";
  std::cout << "      N           VOLUME\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    sphere_unit_volume_values ( n_data, n, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(6)                      << n  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void spherical_harmonic_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERICAL_HARMONIC_VALUES_TEST tests SPHERICAL_HARMONIC_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2007
//
//  Author:
//
//    John Burkardt
//
{

  int l;
  int m;
  int n_data;
  double phi;
  double theta;
  double yi;
  double yr;

  std::cout << "\n";
  std::cout << "SPHERICAL_HARMONIC_VALUES_TEST:\n";
  std::cout << "  SPHERICAL_HARMONIC_VALUES stores values of\n";
  std::cout << "  the spherical harmonic function.\n";
  std::cout << "\n";
  std::cout << "   L   M    THETA       PHI           Yr                    Yi\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    spherical_harmonic_values ( n_data, l, m, theta, phi, yr, yi );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(2)                      << l     << "  "
         << std::setw(2)                      << m     << "  "
         << std::setw(8)  << std::setprecision(4)  << theta << "  "
         << std::setw(8)  << std::setprecision(4)  << phi   << "  "
         << std::setw(24) << std::setprecision(16) << yr    << "  "
         << std::setw(24) << std::setprecision(16) << yi    << "\n";
  }
  return;
}
//****************************************************************************80

void sqrt_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SQRT_VALUES_TEST tests SQRT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "SQRT_VALUES_TEST:\n";
  std::cout << "  SQRT_VALUES returns some exact values.\n";
  std::cout << "\n";
  std::cout << "     X       Fx\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   sqrt_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    std::cout                                       << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void stirling1_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    STIRLING1_VALUES_TEST tests STIRLING1_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int m;
  int n;
  int n_data;
  int s1;

  std::cout << "\n";
  std::cout << "STIRLING1_VALUES_TEST:\n";
  std::cout << "  STIRLING1_VALUES returns values of\n";
  std::cout << "  the Stirling numbers of the first kind.\n";
  std::cout << "\n";
  std::cout << "     N     N        S1\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    stirling1_values ( n_data, n, m, s1 );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(6)  << m  << "  "
         << std::setw(12) << s1 << "\n";
  }
  return;
}
//****************************************************************************80

void stirling2_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    STIRLING2_VALUES_TEST tests STIRLING2_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int m;
  int n;
  int n_data;
  int s2;

  std::cout << "\n";
  std::cout << "STIRLING2_VALUES_TEST:\n";
  std::cout << "  STIRLING2_VALUES returns values of\n";
  std::cout << "  the Stirling numbers of the second kind.\n";
  std::cout << "\n";
  std::cout << "     N     N        S2\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    stirling1_values ( n_data, n, m, s2 );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(6)  << m  << "  "
         << std::setw(12) << s2 << "\n";
  }
  return;
}
//****************************************************************************80

void stromgen_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    STROMGEN_VALUES_TEST tests STROMGEN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "STROMGEN_VALUES_TEST:\n";
  std::cout << "  STROMGEN_VALUES stores values of \n";
  std::cout << "  the Stromgen function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    stromgen_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void struve_h0_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    STRUVE_H0_VALUES_TEST tests STRUVE_H0_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2007
//
//  Author:
//
//    John Burkardt
//
{  
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "STRUVE_H0_VALUES_TEST:\n";
  std::cout << "  STRUVE_H0_VALUES stores values of\n";
  std::cout << "  the Struve H0 function.\n";
  std::cout << "\n";
  std::cout << "      X            H0(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    struve_h0_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void struve_h1_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    STRUVE_H1_VALUES_TEST tests STRUVE_H1_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "STRUVE_H1_VALUES_TEST:\n";
  std::cout << "  STRUVE_H1_VALUES stores values of\n";
  std::cout << "  the Struve H1 function.\n";
  std::cout << "\n";
  std::cout << "      X            H1(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   struve_h1_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void struve_l0_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    STRUVE_L0_VALUES_TEST tests STRUVE_L0_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2007
//
//  Author:
//
//    John Burkardt
//
{  
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "STRUVE_L0_VALUES_TEST:\n";
  std::cout << "  STRUVE_L0_VALUES stores values of\n";
  std::cout << "  the Struve L0 function.\n";
  std::cout << "\n";
  std::cout << "      X            L0(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    struve_l0_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void struve_l1_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    STRUVE_L1_VALUES_TEST tests STRUVE_L1_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "STRUVE_L1_VALUES_TEST:\n";
  std::cout << "  STRUVE_L1_VALUES stores values of\n";
  std::cout << "  the Struve L1 function.\n";
  std::cout << "\n";
  std::cout << "      X            L1(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   struve_l1_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void student_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    STUDENT_CDF_VALUES_TEST tests STUDENT_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2005
//
//  Author:
//
//    John Burkardt
//
{
  double c;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "STUDENT_CDF_VALUES_TEST:\n";
  std::cout << "  STUDENT_CDF_VALUES returns values of\n";
  std::cout << "  the Student T Cumulative Density Function.\n";
  std::cout << "\n";
  std::cout << "      C     X       CDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   student_cdf_values ( n_data, c, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                       << "  "
         << std::setw(16)                     << c  << "  "
         << std::setw(16)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void student_noncentral_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    STUDENT_NONCENTRAL_CDF_VALUES_TEST tests STUDENT_NONCENTRAL_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int df;
  double fx;
  double lambda;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "STUDENT_NONCENTRAL_CDF_VALUES_TEST:\n";
  std::cout << "  STUDENT_NONCENTRAL_CDF_VALUES returns values of\n";
  std::cout << "  the noncentral Student T Cumulative Density Function.\n";
  std::cout << "\n";
  std::cout << "    DF     LAMBDA        X        CDF\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    student_noncentral_cdf_values ( n_data, df, lambda, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                           << "  "
         << std::setw(6)                      << df     << "  "
         << std::setw(8)                      << lambda << "  "
         << std::setw(8)                      << x      << "  "
         << std::setw(24) << std::setprecision(16) << fx     << "\n";
  }
  return;
}
//****************************************************************************80

void subfactorial_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SUBFACTORIAL_VALUES_TEST tests SUBFACTORIAL_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  int fn;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "SUBFACTORIAL_VALUES_TEST:\n";
  std::cout << "  SUBFACTORIAL_VALUES returns values of\n";
  std::cout << "  the subfactorial function.\n";
  std::cout << "\n";
  std::cout << "      N       Subfactorial[N]\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    subfactorial_values ( n_data, n, fn );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(12) << fn << "\n";
  }
  return;
}
//****************************************************************************80

void surten_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    SURTEN_VALUES_TEST tests SURTEN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int n_data;
  double sigma;
  double tc;

  std::cout << "\n";
  std::cout << "SURTEN_VALUES_TEST:\n";
  std::cout << "  SURTEN_VALUES stores values of\n";
  std::cout << "  the surface tension of water\n";
  std::cout << "  as a function of temperature.\n";
  std::cout << "\n";
  std::cout << "      T            SIGMA(T)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   surten_values ( n_data, tc, sigma );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                      << "  "
         << std::setw(12) << tc    << "  "
         << std::setw(12) << sigma << "\n";
  }
  return;
}
//****************************************************************************80

void synch1_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SYNCH1_VALUES_TEST tests SYNCH1_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "SYNCH1_VALUES_TEST:\n";
  std::cout << "  SYNCH1_VALUES stores values of \n";
  std::cout << "  the Synchrotron function of order 1.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    synch1_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void synch2_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SYNCH2_VALUES_TEST tests SYNCH2_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "SYNCH2_VALUES_TEST:\n";
  std::cout << "  SYNCH2_VALUES stores values of \n";
  std::cout << "  the Synchrotron function of order 2.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    synch2_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void tan_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TAN_VALUES_TEST tests TAN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TAN_VALUES_TEST:\n";
  std::cout << "   TAN_VALUES stores values of the tangent function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    tan_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void tanh_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TANH_VALUES_TEST tests TANH_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TANH_VALUES_TEST:\n";
  std::cout << "   TANH_VALUES stores values of the hyperbolic tangent function.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    tanh_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void tau_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    TAU_VALUES_TEST tests TAU_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int fn;
  int n;
  int n_data;

  std::cout << "\n";
  std::cout << "TAU_VALUES_TEST:\n";
  std::cout << "  TAU_VALUES returns values of\n";
  std::cout << "  the TAU function.\n";
  std::cout << "\n";
  std::cout << "     N         TAU(N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   tau_values ( n_data, n, fn );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(6)  << n  << "  "
         << std::setw(12) << fn << "\n";
  }
  return;
}
//****************************************************************************80

void thercon_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    THERCON_VALUES_TEST tests THERCON_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double lambda;
  int n_data;
  double p;
  double tc;

  std::cout << "\n";
  std::cout << "THERCON_VALUES_TEST:\n";
  std::cout << "  THERCON_VALUES stores values of\n";
  std::cout << "  the thermal conductivity of water\n";
  std::cout << "  as a function of temperature and pressure.\n";
  std::cout << "\n";
  std::cout << "      T            P            LAMBDA(T,P)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   thercon_values ( n_data, tc, p, lambda );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                       << "  "
         << std::setw(12) << tc     << "  "
         << std::setw(12) << p      << "  "
         << std::setw(12) << lambda << "\n";
  }
  return;
}
//****************************************************************************80

void three_j_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    THREE_J_VALUES_TEST tests THREE_J_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double j1;
  double j2;
  double j3;
  double m1;
  double m2;
  double m3;
  int n_data;

  std::cout << "\n";
  std::cout << "THREE_J_VALUES_TEST:\n";
  std::cout << "  THREE_J_VALUES returns values of\n";
  std::cout << "  the Wigner 3J coefficient.\n";
  std::cout << "\n";
  std::cout << "      J1      J2      J3      M1      M2      M3        THREE_J\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    three_j_values ( n_data, j1, j2, j3, m1, m2, m3, fx );

    if ( n_data == 0 )
    {
      break;
    }

    std::cout << "  " << std::setw(6) << j1
         << "  " << std::setw(6) << j2
         << "  " << std::setw(6) << j3
         << "  " << std::setw(6) << m1
         << "  " << std::setw(6) << m2
         << "  " << std::setw(6) << m3
         << "  " << std::setprecision(16) << std::setw(24) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void tran02_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRAN02_VALUES_TEST tests TRAN02_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TRAN02_VALUES_TEST:\n";
  std::cout << "  TRAN02_VALUES stores values of \n";
  std::cout << "  the Transport function of order 2.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    tran02_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void tran03_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRAN03_VALUES_TEST tests TRAN03_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TRAN03_VALUES_TEST:\n";
  std::cout << "  TRAN03_VALUES stores values of \n";
  std::cout << "  the Transport function of order 3.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    tran03_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void tran04_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRAN04_VALUES_TEST tests TRAN04_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TRAN04_VALUES_TEST:\n";
  std::cout << "  TRAN04_VALUES stores values of \n";
  std::cout << "  the Transport function of order 4.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    tran04_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void tran05_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRAN05_VALUES_TEST tests TRAN05_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TRAN05_VALUES_TEST:\n";
  std::cout << "  TRAN05_VALUES stores values of \n";
  std::cout << "  the Transport function of order 5.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    tran05_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void tran06_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRAN06_VALUES_TEST tests TRAN06_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TRAN06_VALUES_TEST:\n";
  std::cout << "  TRAN06_VALUES stores values of \n";
  std::cout << "  the Transport function of order 6.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    tran06_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void tran07_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRAN07_VALUES_TEST tests TRAN07_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TRAN07_VALUES_TEST:\n";
  std::cout << "  TRAN07_VALUES stores values of \n";
  std::cout << "  the Transport function of order 7.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    tran07_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void tran08_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRAN08_VALUES_TEST tests TRAN08_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TRAN08_VALUES_TEST:\n";
  std::cout << "  TRAN08_VALUES stores values of \n";
  std::cout << "  the Transport function of order 8.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    tran08_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout << "  " << std::setw(24) << std::setprecision ( 16 ) << x  
         << "  " << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void tran09_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRAN09_VALUES_TEST tests TRAN09_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TRAN09_VALUES_TEST:\n";
  std::cout << "  TRAN09_VALUES stores values of \n";
  std::cout << "  the Transport function of order 9.\n";
  std::cout << "\n";
  std::cout << "                X                     FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    tran09_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                          << "  "
         << std::setw(24) << std::setprecision ( 16 ) << x  << "  "
         << std::setw(24) << std::setprecision ( 16 ) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void trigamma_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    TRIGAMMA_VALUES_TEST tests TRIGAMMA_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "TRIGAMMA_VALUES_TEST\n";
  std::cout << "  TRIGAMMA_VALUES stores values of\n";
  std::cout << "  the TriGamma function.\n";
  std::cout << "\n";
  std::cout << "      X            FX\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    trigamma_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout << "  " << std::setw(12)                     << x  
         << "  " << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void truncated_normal_ab_cdf_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    TRUNCATED_NORMAL_AB_CDF_VALUES_TEST tests TRUNCATED_NORMAL_AB_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  std::cout << "\n";
  std::cout << "TRUNCATED_NORMAL_AB_CDF_VALUES_TEST:\n";
  std::cout << "  TRUNCATED_NORMAL_AB_CDF_VALUES stores values of\n";
  std::cout << "  the Truncated Normal Cumulative Density Function.\n";
  std::cout << "\n";
  std::cout << "        MU     SIGMA       A         B         X        CDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_ab_cdf_values ( n_data, mu, sigma, a, b, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout << "  " << std::setw(8) << mu
         << "  " << std::setw(8) << sigma
         << "  " << std::setw(8) << a
         << "  " << std::setw(8) << b
         << "  " << std::setw(8) << x
         << "  " << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void truncated_normal_ab_pdf_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    TRUNCATED_NORMAL_AB_PDF_VALUES_TEST tests TRUNCATED_NORMAL_AB_PDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  std::cout << "\n";
  std::cout << "TEST1547:\n";
  std::cout << "  TRUNCATED_NORMAL_AB_PDF_VALUES stores values of\n";
  std::cout << "  the Truncated Normal Probability Density Function.\n";
  std::cout << "\n";
  std::cout << "        MU     SIGMA       A         B         X        PDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_ab_pdf_values ( n_data, mu, sigma, a, b, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout << "  " << std::setw(8) << mu
         << "  " << std::setw(8) << sigma
         << "  " << std::setw(8) << a
         << "  " << std::setw(8) << b
         << "  " << std::setw(8) << x
         << "  " << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void truncated_normal_a_cdf_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    TRUNCATED_NORMAL_A_CDF_VALUES_TEST tests TRUNCATED_NORMAL_A_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  std::cout << "\n";
  std::cout << "TEST1548:\n";
  std::cout << "  TRUNCATED_NORMAL_A_CDF_VALUES stores values of\n";
  std::cout << "  the lower Truncated Normal Cumulative Density Function.\n";
  std::cout << "\n";
  std::cout << "        MU     SIGMA       A         X        CDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_a_cdf_values ( n_data, mu, sigma, a, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout << "  " << std::setw(8) << mu
         << "  " << std::setw(8) << sigma
         << "  " << std::setw(8) << a
         << "  " << std::setw(8) << x
         << "  " << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void truncated_normal_a_pdf_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    TRUNCATED_NORMAL_A_PDF_VALUES_TEST tests TRUNCATED_NORMAL_A_PDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  std::cout << "\n";
  std::cout << "TRUNCATED_NORMAL_A_PDF_VALUES_TEST:\n";
  std::cout << "  TRUNCATED_NORMAL_A_PDF_VALUES stores values of\n";
  std::cout << "  the lower Truncated Normal Probability Density Function.\n";
  std::cout << "\n";
  std::cout << "        MU     SIGMA       A         X        PDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_a_pdf_values ( n_data, mu, sigma, a, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout << "  " << std::setw(8) << mu
         << "  " << std::setw(8) << sigma
         << "  " << std::setw(8) << a
         << "  " << std::setw(8) << x
         << "  " << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void truncated_normal_b_cdf_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    TRUNCATED_NORMAL_B_CDF_VALUES_TEST tests TRUNCATED_NORMAL_B_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double b;
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  std::cout << "\n";
  std::cout << "TRUNCATED_NORMAL_B_CDF_VALUES_TEST:\n";
  std::cout << "  TRUNCATED_NORMAL_B_CDF_VALUES stores values of\n";
  std::cout << "  the upper Truncated Normal Cumulative Density Function.\n";
  std::cout << "\n";
  std::cout << "        MU     SIGMA       B         X        CDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_b_cdf_values ( n_data, mu, sigma, b, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout << "  " << std::setw(8) << mu
         << "  " << std::setw(8) << sigma
         << "  " << std::setw(8) << b
         << "  " << std::setw(8) << x
         << "  " << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void truncated_normal_b_pdf_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    TRUNCATED_NORMAL_B_PDF_VALUES_TEST tests TRUNCATED_NORMAL_B_PDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  std::cout << "\n";
  std::cout << "TRUNCATED_NORMAL_B_PDF_VALUES_TEST:\n";
  std::cout << "  TRUNCATED_NORMAL_B_PDF_VALUES stores values of\n";
  std::cout << "  the upper Truncated Normal Probability Density Function.\n";
  std::cout << "\n";
  std::cout << "        MU     SIGMA       B         X        PDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_b_pdf_values ( n_data, mu, sigma, b, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout << "  " << std::setw(8) << mu
         << "  " << std::setw(8) << sigma
         << "  " << std::setw(8) << b
         << "  " << std::setw(8) << x
         << "  " << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void tsat_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    TSAT_VALUES_TEST tests TSAT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int n_data;
  double p;
  double tc;

  std::cout << "\n";
  std::cout << "TSAT_VALUES_TEST:\n";
  std::cout << "  TSAT_VALUES stores values of\n";
  std::cout << "  the saturation temperature\n";
  std::cout << "  as a function of pressure.\n";
  std::cout << "\n";
  std::cout << "      P           Tsat(P)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   tsat_values ( n_data, p, tc );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << p  << "  "
         << std::setw(12) << tc << "\n";
  }
  return;
}
//****************************************************************************80

void van_der_corput_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    VAN_DER_CORPUT_VALUES_TEST tests VAN_DER_CORPUT_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int base;
  int n_data;
  int seed;
  double value;

  std::cout << "\n";
  std::cout << "VAN_DER_CORPUT_VALUES_TEST:\n";
  std::cout << "  VAN_DER_CORPUT_VALUES stores values of\n";
  std::cout << "  the van der Corput sequence in a given base.\n";
  std::cout << "\n";
  std::cout << "      BASE      SEED    VDC(BASE,SEED)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    van_der_corput_values ( n_data, base, seed, value );

    if ( n_data == 0 )
    {
      break;
    }

    std::cout                      << "  "
         << std::setw(8)  << base  << "  "
         << std::setw(8)  << seed  << "  "
         << std::setw(14) << value << "\n";

  }

  return;
}
//****************************************************************************80

void viscosity_values_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    VISCOSITY_VALUES_TEST tests VISCOSITY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double eta;
  int n_data;
  double p;
  double tc;

  std::cout << "\n";
  std::cout << "VISCOSITY_VALUES_TEST:\n";
  std::cout << "  VISCOSITY_VALUES stores values of\n";
  std::cout << "  the viscosity of water\n";
  std::cout << "  as a function of temperature and pressure.\n";
  std::cout << "\n";
  std::cout << "      T            P            ETA(T,P)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   viscosity_values ( n_data, tc, p, eta );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                   << "  "
         << std::setw(12) << tc  << "  "
         << std::setw(12) << p   << "  "
         << std::setw(12) << eta << "\n";
  }
  return;
}
//****************************************************************************80

void von_mises_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    VON_MISES_CDF_VALUES_TEST tests VON_MISES_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "VON_MISES_CDF_VALUES_TEST:\n";
  std::cout << "  VON_MISES_CDF_VALUES stores values of\n";
  std::cout << "  the von Mises CDF.\n";
  std::cout << "\n";
  std::cout << "      A            B            X            CDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    von_mises_cdf_values ( n_data, a, b, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                 << "  "
         << std::setw(12)                     << a  << "  "
         << std::setw(12)                     << b  << "  "
         << std::setw(12)                     << x  << "  "
         << std::setw(24) << std::setprecision(16) << fx << "\n";
  }
  return;
}
//****************************************************************************80

void weekday_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    WEEKDAY_VALUES_TEST tests WEEKDAY_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  int m;
  int n_data;
  int w;
  int y;

  std::cout << "\n";
  std::cout << "WEEKDAY_VALUES_TEST:\n";
  std::cout << "  WEEKDAY_VALUES returns values of \n";
  std::cout << "  the weekday for a given Y/M/D date.\n";
  std::cout << "\n";
  std::cout << "     Y     M     D     W\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    weekday_values ( n_data, y, m, d, w );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout << "  " << std::setw(4) << y
         << "  " << std::setw(4) << m
         << "  " << std::setw(4) << d 
         << "  " << std::setw(4) << w << "\n";
  }
  return;
}
//****************************************************************************80

void weibull_cdf_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_CDF_VALUES_TEST tests WEIBULL_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  double beta;
  double fx;
  int n_data;
  double x;

  std::cout << "\n";
  std::cout << "WEIBULL_CDF_VALUES_TEST:\n";
  std::cout << "  WEIBULL_CDF_VALUES returns values of \n";
  std::cout << "  the Weibull Cumulative Density Function.\n";
  std::cout << "\n";
  std::cout << "     Alpha   Beta        X   CDF(X)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    weibull_cdf_values ( n_data, alpha, beta, x, fx );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                 << "  "
         << std::setw(8)                      << alpha << "  "
         << std::setw(8)                      << beta  << "  "
         << std::setw(8)                      << x     << "  "
         << std::setw(24) << std::setprecision(16) << fx    << "\n";
  }
  return;
}
//****************************************************************************80

void zeta_values_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ZETA_VALUES_TEST tests ZETA_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int n_data;
  double zeta;

  std::cout << "\n";
  std::cout << "ZETA_VALUES_TEST:\n";
  std::cout << "  ZETA_VALUES returns values of \n";
  std::cout << "  the Riemann Zeta function.\n";
  std::cout << "\n";
  std::cout << "     N        ZETA(N)\n";
  std::cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    zeta_values ( n_data, n, zeta );

    if ( n_data == 0 )
    {
      break;
    }
    std::cout                                         << "  "
         << std::setw(6)                      << n    << "  "
         << std::setw(24) << std::setprecision(16) << zeta << "\n";
  }
  return;
}
