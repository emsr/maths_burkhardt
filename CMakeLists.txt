cmake_minimum_required (VERSION 3.12)

project(maths_burkhardt LANGUAGES CXX C Fortran)

if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++11")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif(MSVC)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")

add_library(maths_burkhardt SHARED
	    special_functions.f90
	    special_functions_module.f90
	    prob.cpp
	    asa109.cpp
	    asa152.cpp
	    asa226.cpp
	    asa241.cpp
	    asa243.cpp
	    asa310.cpp
	    bernstein_polynomial.cpp
	    cdflib.cpp
	    chebyshev_polynomial.cpp
	    divdif.cpp
	    filon.cpp
	    gegenbauer_polynomial.cpp
	    hermite_polynomial.cpp
	    jacobi_polynomial.cpp
	    kronrod.cpp
	    lobatto_polynomial.cpp
	    minmax.cpp
	    polpak.cpp
	    timestamp.cpp
	    toms443.cpp
	    toms462.cpp
	    toms917.cpp
	    quadrule.cpp
	    line_monte_carlo.cpp
	    line_ncc_rule.cpp
	    line_nco_rule.cpp
	    tetrahedron_monte_carlo.cpp
	    triangle_fekete_rule.cpp
	    triangle_felippa_rule.cpp
	    triangle_monte_carlo.cpp
	    blas0.cpp
	    blas1_c.cpp
	    blas1_d.cpp
	    blas1_s.cpp
	    blas1_z.cpp
	    blas2_d.cpp
	    blas3_d.cpp
)
target_include_directories(maths_burkhardt PUBLIC .)
target_link_libraries(maths_burkhardt gfortran)
