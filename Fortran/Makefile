#MAKE

# Variables
OUTFILE=solveur_gfortran solveur_ifort

# Usaual compilation flags
CFLAGS1=-Wall -Wextra -O3
CFLAGS2=-O3

# Special rules and targets
.PHONY: all clean cleanall help $(OUTFILE)

# Rules and targets
all: $(OUTFILE)

solveur_gfortran: 
	gfortran $(CFLAGS1) -c methode_resolution.f90
	gfortran $(CFLAGS1) -c main.f90 methode_resolution.f90
	gfortran $(CFLAGS1) -o $@ main.o methode_resolution.o

solveur_ifort: 
	ifort $(CFLAGS2) -c methode_resolution.f90
	ifort $(CFLAGS2) -c main.f90 methode_resolution.f90
	ifort $(CFLAGS2) -o $@ main.o methode_resolution.o

clean:
	@(rm -f *.gcov *.gcda *.gcno *.dat *mod *.o *~ *# $(OUTFILE))

help :
	@(echo -e "Usage:")
	@(echo -e "  make [all]\t\tBuild executable ")
	@(echo -e "  make clean\t\tRemove all object and temporary files")
	@(echo -e "  make help\t\tDisplay this help")
