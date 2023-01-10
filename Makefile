FO=gfortran
ARC=randomove.f90 sinit1.f90 sinit2.f90 main2.f90
OBJ=${ARC:.f90=.o}

%.o: %.f90
	$(FO) -o $@ -c $<

executable.out: $(OBJ) 
	$(FO) -o $@ $(OBJ)

clean:
	rm *.mod *.o *.out