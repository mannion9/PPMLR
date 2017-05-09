CF	  = gfortran
FLAGS = -fdefault-real-8
debugFLAGS  = -O0 -g -Wall  -Wtabs -Wuninitialized -Wextra -Wconversion -fbacktrace -fbounds-check
objects = Main.o io.o Inputs.o Init.o Boundary.o PPM.o Timestep.o EffectiveStates.o Remap.o RiemannSolvers.o EOS.o

exe : $(objects)
	@$(CF) $(FLAGS) $(debugFLAGS) -o exe $(objects)
Inputs.mod : Inputs.o Inputs.f90
	@$(CF) $(FLAGS) $(debugFLAGS) -c Inputs.f90
Inputs.o : Inputs.f90
	@$(CF) $(FLAGS) $(debugFLAGS) -c Inputs.f90
%.o: %.f90 Inputs.mod
	@$(CF) $(FLAGS) $(debugFLAGS) -c $<
clean :
	@$(RM) *.o *.mod Inputs.f90
