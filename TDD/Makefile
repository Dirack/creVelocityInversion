# Makefile
# 
# Purpose: Run all tests with 'make'.
# Or one test as example 'make test_diffmig'.
# 
# Site: https://dirack.github.io
# 
# Version 1.0
# 
# Programer: Rodolfo A C Neves (Dirack) 05/10/2020
# 
# Email: rodolfo_profissional@hotmail.com
# 
# License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

SRC = diffmig diffsimul
DEP = $(addprefix test_, $(SRC)) 
SRCNIPTOMO= ../Mnipmodsetup.c ../raytrace.c ../grid2.c ../atela.c
SRCNIPMODEL= ../Mniptimecurve.c ../raytrace.c ../grid2.c ../atela.c
SRCNIPMODEL2= ../Mniptimecurve2.c ../raytrace.c ../grid2.c ../atela.c
SRCGETPARAM= ../Mgetparameter.c
SRCSTEREO = ../Mstereoniptomo.c ../tomography.c ../raytrace.c ../grid2.c ../atela.c

all:	$(DEP) integ

test_%:
	@echo $@
	@python -m unittest $@

integ:
	cd ../experiments/modelingAndStack/ && scons
	cd ../experiments/algorithmVelocityInversion/ && scons

# Version of the sfnipmodsetup program for GDB (GNU debugger)
nipmodsetup.x:	$(SRCNIPTOMO)
	gcc -g $^ -I$$RSFSRC/include -I.. -L$$RSFSRC/lib -o $@ -lrsf -lm
	@echo "Program nipmodsetup.x for GDB built"
	gdb -q --command=gdb_script1 --args $@

# Version of the sfniptimecurve program for GDB (GNU debugger)
niptimecurve.x:	$(SRCNIPMODEL)
	gcc -g $^ -I$$RSFSRC/include -I.. -L$$RSFSRC/lib -o $@ -lrsf -lm
	@echo "Program niptimecurve.x for GDB built"
	gdb -q --command=gdb_script2 --args $@

# Version of the sfniptimecurve2 program for GDB (GNU debugger)
niptimecurve2.x:	$(SRCNIPMODEL2)
	gcc -g $^ -I$$RSFSRC/include -I.. -L$$RSFSRC/lib -o $@ -lrsf -lm
	@echo "Program niptimecurve2.x for GDB built"
	gdb -q --command=gdb_script3 --args $@

# Version of the sfgetparameter program for GDB (GNU debugger)
getparameter.x:	setup1 $(SRCGETPARAM)
	gcc -g $(SRCGETPARAM) -I$$RSFSRC/include -I.. -L$$RSFSRC/lib -o $@ -lrsf -lm
	@echo "Program getparameter.x for GDB built"
	gdb -q --command=gdb_script4 --args $@

setup1:
	@echo "Setup input files for getparameter.x test..."
	cd stereoniptomo && scons -f SConscript1

setup2:	setup1
	@echo "Setup input files for NIP sources plot..."
	cd stereoniptomo && scons nipsources.rsf && scons -f SConscript2

# Plot NIP sources using matplotlib python library
plotnip:	setup2
	cd stereoniptomo && ./viewNipSources.py

# Version of the sfstereoniptomo program for GDB (GNU debugger)
stereoniptomo.x:	setup1 $(SRCSTEREO)
	gcc -DGDB_DEBUG -g $(SRCSTEREO) -I$$RSFSRC/include -I.. -L$$RSFSRC/lib -o $@ -lrsf -lm
	@echo "Program getparameter.x for GDB built"
	gdb -q --command=gdb_script5 --args $@

clean:
	rm *.x

.PHONY: clean nipmodsetup.x niptimecurve.x niptimecurve2.x getparameter.x stereoniptomo.x
