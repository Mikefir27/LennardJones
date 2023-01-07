      MODULE srandom2
      !
      IMPlICIT NONE
      !
      !
      contains
      !
      !
      SUBROUTINE random2(N, L, cv)
      implicit none
      INTEGER(KIND=4), intent(in) :: N
      REAL(KIND=8), intent(in) :: L
      REAL(KIND=8), allocatable, intent(out) :: cv(:,:)
      INTEGER(KIND=4) :: I, J, V, RA, ix, Z, Y, II, JJ
      REAL(KIND=8) :: K, A, B, C, M, LM, x
      REAL(KIND=8),allocatable :: coor(:,:), cd(:,:)
      !
      !
      !Calculate a value of M:
      K = N**(1.00/3.00)
      !
      M = 1
      DO WHILE (M .LT. K)
         M = M + 1
      ENDDO
      !WRITE(6,*) M
      !
      !Create the initial state of L*L*L:
      !
      V = M**3
      Allocate (coor(V, 3))
      !
      LM = L/M
      !WRITE(6,*) LM
      A = 1.0
      B = 1.0
      C = 1.0
      I = 1
      M = real(M)
      DO WHILE (A .LE. M)
         coor(I,1) = A*LM
         coor(I,2) = B*LM
         coor(I,3) = C*LM
        ! B = B + 1
        ! I = I + 1
         DO WHILE (B .LE. M)
            I = I + 1
            C = C + 1
            coor(I,1) = A*LM
            coor(I,2) = B*LM
            coor(I,3) = C*LM
          !  C = C + 1
            I = I + 1
            DO WHILE (C .LE. M)
               C = C + 1
               IF (C .GT. M) then
                  cycle
               Else
                  coor(I,1) = A*LM
                  coor(I,2) = B*LM
                  coor(I,3) = C*LM
           !    C = C + 1
                  I = I + 1
               ENDIF
            END DO
            B = B + 1
            C = 1
            IF (B .GT. M) then
               cycle
            Else
               coor(I,1) = A*LM
               coor(I,2) = B*LM
               coor(I,3) = C*LM
            ENDIF
         END DO
         B = 1
         A = A + 1
      END DO
      !WRITE(6,*) I
      !
      !WRITE(6,*) coor
      OPEN (unit=9, file="incoordfull.xyz", action="write")
      !
      WRITE(9,*) V 
      WRITE(9,*) "INITIAL SYSTEM "
      DO I = 1,V
         WRITE(9,*) "Ar", (coor(I,J),J=1,3)
      ENDDO
      !
      !
      !REMOVE RANDOM ATOMS
      !
      !
      RA = V - N
      Allocate (cv(V, 3))
      Allocate (cd(V, 3))
      !
      DO I = 1, V
         DO J = 1, 3
            cv(I, J) = coor(I,J)
      !      WRITE(6,*) cv(I,J)
         ENDDO
      ENDDO
      !
      !
      !
      DO I = 1, RA
         JJ = 1
         call random_number(x)
         Z = V + 1 - I
      !   WRITE(6,*) "Z = ", Z
         x = x * Z
         ix = nint(x)
      !   WRITE(6,*) "RN = ", ix
         DO WHILE (ix .EQ. 0)
      !      WRITE(6,*) ix
            call random_number(x)
            x = x * Z
            ix = nint(x)
         ENDDO
            DO J = 1, Z
               IF (J .NE. ix) then
                  cd(JJ,1) = cv(J,1)
                  cd(JJ,2) = cv(J,2)
                  cd(JJ,3) = cv(J,3)
               ELSE
       !           WRITE(6,*) "The random number is:", J
                  JJ = JJ - 1
                  CONTINUE
               ENDIF
               JJ = JJ + 1
            ENDDO
            WRITE(6,*) "Loop", I
            DO II = 1, V-I
       !        WRITE(6,*) "Atom", II
               DO J = 1,3
       !           WRITE(6,*) "cd", cd(II,J)
               ENDDO
            ENDDO
            cv = 0.0
        !    WRITE(6,*) "Convertir la cd a cv:"
            DO J = 1, V-I
        !       WRITE(6,*) "Atom", J
               DO Y = 1, 3
                  cv(J, Y) = cd(J, Y)
         !         WRITE(6,*) "cv", cv(J,Y)
               ENDDO
            ENDDO
            cd = 0.0
      ENDDO
      !
      !
      !
      !
      OPEN (unit=10, file="incoord2.xyz", action="write")
      !
      WRITE(10,*) N
      WRITE(10,*) "INITIAL SYSTEM "
      DO I = 1,N
         WRITE(10,*) "Ar", (cv(I,J),J=1,3)
      ENDDO
      !
      !
      !
      END SUBROUTINE random2

      END MODULE srandom2
