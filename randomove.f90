       MODULE randomove
       !
       !
       IMPLICIT NONE
       !
       contains
       !
       SUBROUTINE VECTDIST(v1, v2, r12)
       implicit none
       REAL(KIND=8),intent(in) :: v1(3), v2(3)
       REAL(KIND=8),intent(out) :: r12
       REAL(KIND=8) :: x21, y21, z21
       x21 = v2(1)-v1(1)
       y21 = v2(2)-v1(2)
       z21 = v2(3)-v1(3)
       r12 = DSQRT(x21**2 + y21**2 + z21**2)
       END SUBROUTINE VECTDIST
       !
       !
       SUBROUTINE VECTDISTV2(v1, v2, r12, L)
       implicit none
       INTEGER(KIND=4) :: I
       REAL(KIND=8), intent(in) :: L, v1(3), v2(3)
       REAL(KIND=8), intent(out) :: r12
       REAL(KIND=8) :: dr21(3)
       dr21 = v2 - v1
       DO I = 1,3
       IF (DABS(dr21(I)) .GT. L/2) then
          IF (dr21(I) .LT. 0) then
             dr21(I) = dr21(I) + L
          ELSE IF (dr21(I) .GT. 0) then
             dr21(I) = dr21(I) - L
          ENDIF
       ENDIF
       ENDDO
       r12 = DSQRT(dr21(1)**2 + dr21(2)**2 + dr21(3)**2)
       ENDSUBROUTINE VECTDISTV2
       !
       !
       SUBROUTINE ranmove(N, L, mcoor, fcoor, k)
       implicit none
       INTEGER(KIND=4) :: I, J, N
       INTEGER(KIND=4), intent(out) :: k
       REAL(KIND=8), intent(in) :: L
       REAL(KIND=8) :: wx, wy, wz, rk, movex, movey, movez, drmax, M
       REAL(KIND=8), intent(in) :: mcoor(3,N)
       REAL(KIND=8), intent(out) :: fcoor(3,N)
       REAL(KIND=8) :: wi(3)
       !
       !
       drmax = 0.10
       !
       !
       call random_number(rk)
       !
       M = rk * N + 1
       k = int(M)
       !
       call random_number(wx)
       call random_number(wy)
       call random_number(wz)
       !
       !
       !
       movex = 2 * drmax * (wx - 0.5)
       movey = 2 * drmax * (wy - 0.5)
       movez = 2 * drmax * (wz - 0.5)
       ! 
       !
       wi = (/movex, movey, movez /)
       !
       !
       !
       !
       DO I = 1, 3
          fcoor(I, k) = mcoor(I, k) + wi(I)
          IF (fcoor(I, k) .GE. L) then
             WRITE(6,*) "atom off limits"
             fcoor(I, k) = fcoor(I, k) - L
          ELSE IF (fcoor(I, k) .LT. 0.0) then
             WRITE(6,*) "atom off limits"
             fcoor(I, k) = fcoor(I, k) + L
          ENDIF
       ENDDO
       !
       WRITE(6,*) "MOVED ATOM:", k
       WRITE(6,*) "MOVED IT TO:", (fcoor(I,k),I = 1,3)
       !
       ENDSUBROUTINE ranmove
       !
       !
       !
       !
       SUBROUTINE IKvec(N, L, mcoor, k, vrMI)
       implicit none
       INTEGER(KIND=4), intent(in) :: k, N
       INTEGER(KIND=4) :: X, I
       REAL(KIND=8), intent(in) :: L
       REAL(KIND=8), intent(in) :: mcoor(3,N)
       REAL(KIND=8), allocatable, intent(out) :: vrMI(:)
       REAL(KIND=8) :: vi(3), vk(3)
       REAL(KIND=8), allocatable :: vrIK(:)
       REAL(KIND=8) :: rik, L2
       !
       Allocate (vrMI(N-1))
       Allocate (vrIK(N-1))
       !
       !WRITE(6,*) "Aqui no es"
       DO I = 1, 3
          vk = mcoor(I,k)
       ENDDO
       !
       X = 1
       !
       DO WHILE (X .LT. N) 
          IF (X .EQ. k) then
             X = X + 1
             cycle
          ELSE
             DO I = 1, 3
                vi(I) = mcoor(I,X)
                call VECTDIST(vi, vk, rik)
                vrik(X) = rik
             ENDDO
          ENDIF
          X = X + 1
       ENDDO
       !
       !
       !Minimum image
       !
       L2 = L/2
       !
       DO I = 1, N-1
          IF (vrIK(I) .GT. L2) then
             vrMI(I) = DABS(vrIK(I) - L)
          ELSE 
             vrMI(I) = vrIK(I)
          ENDIF
       ENDDO
       !
       !
       ENDSUBROUTINE IKvec
       !
       !
       !
       SUBROUTINE deltV(rmk, rfk, fm)
       implicit none
       REAL(KIND=8), intent(in) :: rmk, rfk 
       REAL(KIND=8), intent(out) :: fm
       REAL(KIND=8) :: m, f
       !
       !
       !
       m = (1/rmk)**12 - (1/rmk)**6
       f = (1/rfk)**12 - (1/rfk)**6
       !
       fm = f - m
       !
       ENDSUBROUTINE deltV
       !
       !
       SUBROUTINE alpha(dV, T, a)
       implicit none
       REAL(KIND=8), intent(in) :: dV, T
       REAL(KIND=8), intent(out) :: a
       REAL(KIND=8), parameter :: kb = 8.314462618D-3
       REAL(KIND=8) :: expon
       !
       expon = dV / T
       a = dexp(-expon)
       !
       ENDSUBROUTINE
       !
       !
       !
       SUBROUTINE acceptornot(a, accept)
       implicit none
       REAL(KIND=8), intent(in) :: a
       REAL(KIND=8) :: w
       LOGICAL, intent(out) :: accept
       !
       call random_number(w)
       !
       IF (a .GE. w) then
          accept = .TRUE.
       ELSE
          accept = .FALSE.
       ENDIF
       !
       !
       ENDSUBROUTINE
       !
       !
       SUBROUTINE RADdist(N, dens, Nacc, dr, Ni, rvec, U, g)
       implicit none
       INTEGER(KIND=4), intent(in) :: N, U, Nacc
       REAL(KIND=8), intent(in) ::  rvec(U), dens, dr
       INTEGER(KIND=4), intent(in) :: Ni(U)
       INTEGER(KIND=4) :: I
       REAL(KIND=8) :: rNacc, rU, rN, cte, rpart
       REAL(KIND=8), PARAMETER :: pi = 4 * ATAN(1.0)
       REAL(KIND=8), intent(out) :: g(U)
       !
       rNacc = DBLE(Nacc)
       rU = DBLE(U) 
       rN = DBLE(N)
       !
       cte = rNacc * (rN / 2.0) * (4.0 / 3.0) * pi * dens
       !
       !WRITE(6,*) cte
       DO I = 1, U
          rpart = rvec(I)**3 - (rvec(I) - dr)**3
     !     WRITE(6,*) "rpart"
     !     WRITE(6,*) rpart
          g(I) = Ni(I) / (cte * rpart)
       ENDDO
       !
       ENDSUBROUTINE
       !
       !
       END MODULE randomove
