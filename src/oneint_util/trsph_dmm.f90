Subroutine dmm_transform(l,Fi1,Fi2,Fi3,dmm,lda)

  Implicit None

  Integer, intent(in) :: l, lda
  Real, intent(in) :: Fi1, Fi2, Fi3
  Real*8, intent(out) :: dmm(lda,-l:l)  !Dmm

  Integer :: k_min, k_max, k, m, mp, a1, a2
  Integer :: i
  Real :: fact1, fact2, fact3, fact4, coeff, tmp
  Real :: re1, re3
!====================================================================

  Do m = -l, l
     Do mp = -l, l
        fact1 = factorial(l+mp)
        fact2 = factorial(l-mp)
        fact3 = factorial(l+m)
        fact4 = factorial(l-m)
        coeff = (-1)**(mp -m) * sqrt(fact1 * fact2 * fact3 * fact4)

        Call k_interval(l, m, mp, k_min, k_max)

        tmp=0.0D0
        Do k = k_min, k_max
           a1 = choose(l+m, k)
           a2 = choose(l-m, l-mp-k)

           tmp = tmp + (-1)**k * a1 * a2 * cos((Fi2)/2)**(2*l-mp+m-2*k)*  &
                sin((Fi2/2))**(2*k-m+mp)
        End Do

        If (mp.gt.0) Then
           re1 = cos(-mp*Fi1)/Sqrt(2.0D0)
        Else If (mp.lt.0) Then
           re1 = sin(-mp*Fi1)/Sqrt(2.0D0)
        Else
           re1 = 1.0D0
        End If

        If (m .gt.0) Then
           re3 = cos(-m *Fi3)/Sqrt(2.0D0)
        Else If (m .lt.0) Then
           re3 = sin(-m *Fi3)/Sqrt(2.0D0)
        Else
           re3 = 1.0D0
        End If

        dmm(mp+l+1,m) = Real(re1*tmp*re3)

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
End Subroutine dmm_transform
