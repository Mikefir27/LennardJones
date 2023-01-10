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
       unitn, W, Q, Tr, Y, indice
       REAL(KIND=8) :: dV, Temp, alp, L, L2, rfk,rmk, SUMdV, mm, f, fm,&
       NOTACC, ACC, alphaprop, dr, rmax, ri, r, dens, rij
       LOGICAL :: accept
       REAL(KIND=8), allocatable :: coor(:,:), mcoor(:,:), fcoor(:,:), &
       vrMIm(:), vrMIf(:), vecr(:),  g(:), vecfkdist(:), alldist(:)
       INTEGER(KIND=4), allocatable :: Ni(:)
       REAL(KIND=8) :: point(3), pointkf(3), pointkm(3), ivec(3), jvec(3)
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
       WRITE(6,*) "Select number of atoms!"
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
       WRITE(9,*) N
       WRITE(9,*) "Transpose "
       DO I = 1,N
         WRITE(9,*) "Ar", (mcoor(J,I),J=1,3)
       ENDDO
!
       WRITE(9,*) N
       WRITE(9,*) "INITIAL SYSTEM again"
       DO I = 1,N
          WRITE(9,*) "Ar", (coor(I,J),J=1,3)
       ENDDO
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
       DO Y = 1, 2000*N
       !DO WHILE (X .LE. 2000)
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
       !         WRITE(6,*) "Atom", I
                call VECTDIST(point, pointkf, rfk)
      !          WRITE(6,*) "rfk before : ", rfk
                call VECTDISTV2(point, pointkf, rfk, L)
       !         WRITE(6,*) "rfk after: ", rfk
                call VECTDISTV2(point, pointkm, rmk, L)
       !         WRITE(6,*) "______________________________________"
                call deltV(rmk, rfk, fm)
                SUMdv = SUMdv + fm
             ELSE
                cycle
             ENDIF
          ENDDO
          dV = 4 * SUMdv
          IF (dV .LE. 0) then
            accept = .TRUE.
          ELSE
             call alpha(dV, Temp, alp)
             call acceptornot(alp, accept)
          ENDIF
          IF (accept .EQV. .TRUE.) then
            mcoor = fcoor
          ENDIF
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
       !
       Tr = ((N - 1)*N)/2
       U = INT(rmax / dr)
       Allocate (Ni(U), vecr(U))
       Allocate (alldist(Tr))
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
          modx = mod(X, N)
          ri = 0.0
          dV = 0
          SUMdV = 0
          !WRITE(6,*) "----------------------------------------------"
          !WRITE(6,*) "Cycle ", X
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
                call VECTDISTV2(point, pointkm, rmk, L)
               ! IF (isnan(rfk)) print*, "atom ", I, "for k = ", k, &
               ! "coord f =", point, "coord k =", pointkf
                call deltV(rmk, rfk, fm)
                WRITE(6,*) "Atom", I
                WRITE(6,*) "rfk : ", rfk, "rmk : ", rmk
                WRITE(6,*) "______________________________________"
                SUMdv = SUMdv + fm
             ELSE
                cycle
             ENDIF
          ENDDO
          dV = 4 * SUMdv
          WRITE(6,*) "The difference of energy is: ", dV, "kj/mol)"

          IF (dV .LE. 0) then
             accept = .TRUE.
          ELSE
             call alpha(dV, Temp, alp)
        !     WRITE(6,*) "This dV is positive"
             call alpha(dV, Temp, alp)
        !     WRITE(6,*) "alpha is :", alp
             call acceptornot(alp, accept)
          ENDIF
          !
          IF (accept .EQV. .TRUE.) then
             ACC = ACC + 1
             mcoor = fcoor
          ELSE
             NOTACC = NOTACC + 1
          ENDIF
       !      WRITE(6,*) "Change Accepted!"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!This part creates the histogram by adding 1 to every Ni in the interval of bond lengths defined in vecr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
          IF (modx .EQ. 0) then
             !WRITE(6,*) "In this cycle g is calculated"
             W = 0
             DO I = 1, N-1
             DO J = I+1, N
                W = W + 1
                DO Q = 1, 3 
                   ivec(Q) = mcoor(Q, I)
                   jvec(Q) = mcoor(Q, J)
                ENDDO
                call VECTDISTV2(ivec, jvec, rij, L)
                   IF ((rij .GE. 0.0) .and.&
   (rij .LT. 4.0)) then
                      indice = INT((rij/dr)+1)
                      Ni(indice) = Ni(indice) + 1     
                     ! WRITE(6,*) "rij = ", rij, "INDEX", indice           
                   ENDIF
              !WRITE(6,*) "rij is", rij
 !               DO Y = 1, U
 !                 IF ((rij .GE. vecr(Y)) .and.&
 !  (rij .LT. (vecr(Y)+dr))) then
 !                    Ni(Y) = Ni(Y) + 1
 !                  ENDIF
 !               ENDDO
             ENDDO 
             ENDDO
             Nacc = Nacc + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
             OPEN (unit=14, file="Raddistcoord.xyz", action=&
    "write")
             WRITE(14,*) N
             WRITE(14,*) "Coordinates for radial dist.", Nacc
             DO J = 1, N
                WRITE(14,*) "Ar", (fcoor(I,J),I = 1,3)
             ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
          ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
         ! WRITE(6,*) "----------------------------------------------"
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
          WRITE(6, *) "r = ", vecr(I), "Ni = ", Ni(I), "g = ", g(I)
       ENDDO
      !
       OPEN (unit=13, file="radialdist.xyz", action="write")
       WRITE(13,'(A, 6X, A)') "r", "g(ri)"
       DO I = 1, U
          WRITE(13, '(F4.2, 6X, F8.4)') vecr(I), g(I)
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
