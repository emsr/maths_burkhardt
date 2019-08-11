module special_functions

contains

  ! airy_ai airy_bi
  subroutine airyab(x, ai, bi, ad, bd) bind(c)
    use iso_c_binding
    implicit none
    real(c_double), intent(in), value :: x
    real(c_double), intent(out) :: ai, bi, ad, bd
    call airya(x, ai, bi, ad, bd)
  end subroutine

  ! struve_h
  subroutine struveh(nu, x, sh) bind(c)
    use iso_c_binding
    implicit none
    real(c_double), intent(in), value :: nu, x
    real(c_double), intent(out) :: sh
    call stvhv(nu, x, sh)
  end subroutine

  ! struve_l
  subroutine struvel(nu, x, sl) bind(c)
    use iso_c_binding
    implicit none
    real(c_double), intent(in), value :: nu, x
    real(c_double), intent(out) :: sl
    call stvlv(nu, x, sl)
  end subroutine

  ! conf_hyperg
  subroutine confhyp(a, b, x, chg) bind(c)
    use iso_c_binding
    implicit none
    real(c_double), intent(in), value :: a, b, x
    real(c_double), intent(out) :: chg
    call cchg(a, b, x, chg)
  end subroutine

  ! tricomi_u
  subroutine tricomiu(a, b, x, hu, md) bind(c)
    use iso_c_binding
    implicit none
    real(c_double), intent(in), value :: a, b, x
    real(c_double), intent(out) :: hu
    integer(c_int), intent(out) :: md
    call chgu(a, b, x, hu, md)
  end subroutine

  ! hyperg
  subroutine hyperg2f1(a, b, c, x, hf) bind(c)
    use iso_c_binding
    implicit none
    real(c_double), intent(in), value :: a, b, c, x
    real(c_double), intent(out) :: hf
    call hygfx(a, b, c, x, hf)
  end subroutine

end module
