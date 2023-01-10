      module srandom1
      !
      use randomove
      !
      implicit none
      !
      contains
      !
      !
      SUBROUTINE random1(N, L, coor)
      IMPLICIT NONE
      INTEGER(KIND=4), intent(in) :: N
      INTEGER(KIND=4) ::  M, I, J
      REAL(KIND=8), intent(in) :: L
      REAL(KIND=8) :: rij, mx, my, mz, rc
      REAL(KIND=8) :: ra(3), rb(3)
      REAL(KIND=8), intent(out) :: coor(N,3)
      !
      !
      !Box size in reduced units of length
      !
      !WRITE(6,*) "Length of box"
      !READ(5,*) L
      !
      !WRITE(6,*)"How many atoms"
      !READ(5,*) N
      !
      !L = real(LI)
      !
      !Allocate (coor(N,3))
      !
      call random_number(mx)
      call random_number(my)
      call random_number(mz)
      ra(1) = mx*L
      ra(2) = my*L
      ra(3) = mz*L
      !
      coor(1,1) = ra(1)
      coor(1,2) = ra(2)
      coor(1,3) = ra(3)
      !
      WRITE(6,*) ra(1), ra(2), ra(3)
      !
      !minimum threshold proposed for the distance
      rc = 0.9
      !
      OPEN (unit=9, file="incoord.xyz", action="write")
      !
      M = 2
      DO WHILE (M .LE. N)
         call random_number(mx)
         call random_number(my)
         call random_number(mz)
         rb(1) = mx*L
         rb(2) = my*L
         rb(3) = mz*L
         call VECTDIST(ra, rb, rij)
         IF (rij .GE. rc) then
            DO J = 1,3
               coor(M,J) = rb(J)
            ENDDO
            M = M + 1
         ELSE
            WRITE(6,"(A4,I2.1,A26)") "The ", M, "th atom had to be &
            repeated"
            CYCLE
         END IF
      ENDDO
      !
      WRITE(9,*) N
      WRITE(9,*) "INITIAL SYSTEM "
      DO I = 1,N
         WRITE(9,*) "Ar", (coor(I,J),J=1,3)
      ENDDO
      !
      !
      !CLOSE(9)
      !
      END SUBROUTINE random1
      !
      !
      !
      END MODULE srandom1
