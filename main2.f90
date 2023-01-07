       PROGRAM main
       !
       Use srandom1
       Use srandom2
       Use randomove
       !
       !
       !
       IMPLICIT NONE
       !
       !
       INTEGER(KIND=4) :: A, N, k, I, J, B, M, X, modx, U, Nacc, UU,&
       unitn
       REAL(KIND=8) :: dV, Temp, alp, L, L2, rfk,rmk, SUMdV, mm, f, fm,&
       NOTACC, ACC, alphaprop, dr, rmax, ri, r, dens
       LOGICAL :: accept
       REAL(KIND=8), allocatable :: coor(:,:), mcoor(:,:), fcoor(:,:), &
       vrMIm(:), vrMIf(:), vecr(:),  g(:), vecfkdist(:)
       INTEGER(KIND=4), allocatable :: Ni(:)
       REAL(KIND=8) :: point(3), pointkf(3), pointkm(3)
       CHARACTER(14) :: filename, cNacc
       !
       WRITE(6,*) "This is the main program. The first step is to gene&
       rate an initial state of Argon molecules. Select the method in &
       which this state is going to be generated:"
       WRITE(6,*) "(1) Randomly generating N atoms in a L*L*L cube."
       WRITE(6,*) "(2) Periodically filling a L*L*L cube and eliminatin&
       g atoms until N atoms are left."
       !
       !
       !
       READ(5,*) A
       DO WHILE ((A .GT. 3) .or. (A .LE. 0))
          WRITE(6,*) "Please select a correct method"
          READ(5,*) B
          A = B
       END DO
       !
       !
       !
       WRITE(6,*) "Select number of atoms"
       READ(5,*) N
       !
       WRITE(6,*) "Select length of box"
       READ(5,*) L
       !
       !
       Allocate (coor(N,3))
       Allocate (mcoor(3,N))
       !
       !
       IF (A .EQ. 1) then
          call random1(N, L, coor)
       ELSE
          call random2(N, L, coor)
       ENDIF
       !
       !
       Allocate (fcoor(3,N))
       mcoor = TRANSPOSE(coor)
       !
       !
       !
       !EQUILIBRATION PART
       !
       !
       !
       Temp = 1.2
       !
       OPEN (unit=11, file="LJcoor.xyz", action="write")
       !
       Allocate (vrMIm(N-1), vrMIf(N-1))
       !
       M = N - 1
       !
       X = 1
       !
       SUMdv = 0
       !
       L2 = L/2
       !
       !TO CALCULATE ALPHA
       ACC = 0.0
       NOTACC = 0.0
       !

       !
       WRITE(6,*) "EQUILIBRATION PART"
       !DO X = 1,5
       DO WHILE (X .LE. 2000)
!          modx = mod(X, 5)
          dV = 0
          SUMdV = 0
          call ranmove(N, L, mcoor, fcoor, k)
          DO I = 1, N
             IF (I .NE. k) then
                DO J = 1, 3
                   pointkm(J) = mcoor(J,k)
                   pointkf(J) = fcoor(J,k)
                   point(J) = mcoor(J,I)
                ENDDO
                call VECTDISTV2(point, pointkf, rfk, L)
                call VECTDISTV2(point, pointkm, rmk, L)
                call deltV(rmk, rfk, fm)
                SUMdv = SUMdv + fm
             ELSE
                cycle
             ENDIF
          ENDDO
          dV = 4 * SUMdv
          IF (dV .LE. 0) then
             X = X + 1
          ELSE
             call alpha(dV, Temp, alp)
             call acceptornot(alp, accept)
             IF (accept .EQV. .TRUE.) then
                continue
             ELSE
                cycle
             ENDIF
             X = X + 1
          ENDIF
          DO J = 1, 3
             mcoor(J,k) = fcoor(J,k)
          ENDDO
       ENDDO
       !
       !
       OPEN (unit=16, file="Equilibrcoord.xyz", action="write")
       WRITE(16,*) N
       WRITE(16,*) "Equilibration Coordinates"
       DO I = 1, N
          WRITE(16,*) "Ar", (fcoor(:,I))
       ENDDO
       !
       !PRODUCTION PART
       !!
       Allocate (vecfkdist(N-1))
       WRITE(6,*) "PRODUCTION PART"
       !TO CALCULATE g(ri)
       Nacc = 0
       rmax = 4.0
       dr = 0.01
       !
       U = INT(rmax / dr)
       Allocate (Ni(U), vecr(U))
       vecr = 0.0
       Ni = 0
       !CREATE vecr
       !
       ri = 0
       DO I = 1, U
          vecr(I) = ri + dr
          ri = vecr(I)
       ENDDO
       !
       unitn = 15
       !
       DO X = 1, 500*N
          modx = mod(X, 5)
          ri = 0.0
          dV = 0
          SUMdV = 0
          !WRITE(6,*) "----------------------------------------------"
          WRITE(6,*) "Cycle ", X
          call ranmove(N, L, mcoor, fcoor, k)
          !WRITE(6,*) "k = ", k
          DO I = 1, N
             IF (I .NE. k) then
                DO J = 1, 3
                   pointkm(J) = mcoor(J,k)
                   pointkf(J) = fcoor(J,k)
                   point(J) = mcoor(J,I)
                ENDDO
                call VECTDISTV2(point, pointkf, rfk, L)
                vecfkdist(I) = rfk
                call VECTDISTV2(point, pointkm, rmk, L)
                call deltV(rmk, rfk, fm)
                SUMdv = SUMdv + fm
             ELSE
                cycle
             ENDIF
          ENDDO
          dV = 4 * SUMdv
         ! WRITE(6,*) "The difference of energy is: ", dV, "kj/mol)"
!
          IF (dV .LE. 0) then
             accept = .TRUE.
             ACC = ACC + 1
             mcoor = fcoor
          ELSE
        !     WRITE(6,*) "This dV is positive"
             call alpha(dV, Temp, alp)
        !     WRITE(6,*) "alpha is :", alp
             call acceptornot(alp, accept)
             IF (accept .EQV. .TRUE.) then
                ACC = ACC + 1
                mcoor = fcoor
             ELSE
                NOTACC = NOTACC + 1
             ENDIF
          ENDIF
!
       !      WRITE(6,*) "Change Accepted!"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!This part creates the histogram by adding 1 to every Ni in the interval of bond lengths defined in vecr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
          IF (modx .EQ. 0) then
             DO I = 1, N-1
        !        WRITE(6,*) "distance between", k, "and", I, "=",&
      !vecfkdist(I)
                DO J = 1, U
            !       WRITE(6,*) "histogram",vecr(J), "---", (vecr(J) + dr)
                   IF ((vecfkdist(I) .GE. vecr(J)).and.&
      (vecfkdist(I) .LT.(vecr(J)+dr))) then
         !             WRITE(6,*) "rfk lies between",vecr(J),&
      !"and", vecr(J)+dr
       !            Nacc = Nacc + 1
                      Ni(J) = Ni(J) + 1
                   ELSE
                      continue
                   ENDIF
                ENDDO
             ENDDO
             Nacc = Nacc + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
             OPEN (unit=14, file="RaddistcdeltVoord.xyz", action=&
      "write")
             WRITE(14,*) N
             WRITE(14,*) "Coordinates for radial dist.", Nacc
             DO J = 1, N
                WRITE(14,*) "Ar", (fcoor(I,J),I = 1,3)
             ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
          ELSE
             continue
          ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
          WRITE(6,*) "----------------------------------------------"
       ENDDO
       !
       !
       !
       alphaprop = 100 * ACC / (ACC + NOTACC)
       WRITE(6,*) "alpha proportion = ", alphaprop, "%"
       !
       !
       !GET g(ri)
       !
     !   DO I = 1, U
     !      WRITE(6,*) "U = ", I
     !      WRITE(6,*) Ni(I)
     !   ENDDO
       WRITE(6,*) "Nacc", Nacc
       !
       dens = 0.7
       Allocate (g(U))
!
       call RADdist(N, dens, Nacc, dr, Ni, vecr, U, g)
       !WRITE(6,*) "HOLA"
      ! WRITE(6,*) "r", "g(ri)"
       DO I = 1, U
          WRITE(6, *) Ni(I)
       ENDDO
      !
       OPEN (unit=13, file="radialdist.xyz", action="write")
       WRITE(13,'(A, 6X, A)') "r", "g(ri)"
       DO I = 1, U
          WRITE(13, '(F1.2, 6X, F4.4)') vecr(I), g(I)
       ENDDO
       !
       !WRITE IN AN XYZ DOCUMENT
       !
       WRITE(11,*) N
       WRITE(11,*) "FINAL STATE"
       DO I = 1, N
          WRITE(11,*) "Ar", (fcoor(J,I),J = 1,3)
       ENDDO
       !
       !
       !
       !
       END PROGRAM main
