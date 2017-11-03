cmake_minimum_required (VERSION 3.2)

project(maths_burkhardt LANGUAGES CXX C Fortran)

if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++11")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif(MSVC)

add_library(maths_burkhardt SHARED
            maths_burkhardt.h
	    special_functions.f90
	    asa109.cpp
	    asa243.cpp
	    asa310.cpp
	    bernstein_polynomial.cpp
	    cdflib.cpp
	    chebyshev_polynomial.cpp
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
           )
target_link_libraries(maths_burkhardt gfortran)