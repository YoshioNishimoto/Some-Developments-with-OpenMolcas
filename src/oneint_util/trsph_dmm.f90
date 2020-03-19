Subroutine dmm_tranform(l,Fi2, dmm)

  Implicit None

  Integer, intent(in) :: l
  Real, intent(in) :: Fi2
  Real*8, intent(out) :: dmm(dim,dim)

  Integer :: k_min, k_max, k, m, mp, a1, a2
  Real :: fact1, fact2, fact3, fact4, coeff
  Integer :: dim
  dim = 2*l + 1
!====================================================================

  Do m = -l, l
     Do mp = -l, l
        fact1 = factorial(l+mp)
        fact2 = factorial(l-mp)
        fact3 = factorial(l+m)
        fact4 = factorial(l-m)
        coeff = (-1)**(mp -m) * sqrt(fact1 * fact2 * fact3 * fact4)

        Call k_interval(l, m, mp, k_min, k_max)

        Do k = k_min, k_max
           a1 = choose(l+m, k)
           a2 = choose(l-m, l-mp-k)

           dmm = dmm + (-1)**k * a1 * a2 * cos((Fi2)/2)**(2*l-mp+m-2*k)*  &
                sin((Fi2/2))**(2*k-m+mp)
        End Do
     End Do
  End Do

!====================================================================
Contains

  Function factorial (n) result (res)

    Implicit None
    Integer, intent (in) :: n
    Integer :: res, i

    res = product ((/(i, i = 1, n)/))

  End Function factorial

!------------------------------------------------------

  Function choose (n, k) result (res)

    Implicit None
    Integer, intent (in) :: n
    Integer, intent (in) :: k
    Integer :: res

    res = factorial (n) / (factorial (k) * factorial (n - k))

  End Function choose

!-----------------------------------------------------

  Subroutine k_interval(l, m, mp, k_min, k_max)

    Integer, intent(in) :: l, m, mp
    Integer, intent(out) :: k_min, k_max

    k_min = max(0, m-mp)
    k_max = min(l-mp, l+m)

  End Subroutine k_interval

!======================================================
End Subroutine dmm_tranform
