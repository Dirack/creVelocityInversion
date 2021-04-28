# Makefile
# 
# Purpose: Run all tests with 'make'.
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

SRCNIPTOMO= ../Mnipmodsetup.c ../raytrace.c ../grid2.c ../atela.c
SRCGETPARAM= ../Mgetparameter.c
SRCSTEREO = ../Mstereoniptomo.c ../tomography.c ../raytrace.c ../grid2.c ../atela.c ../vfsacrsnh_lib.c

all:
	@echo "No unit and integration tests found, run 'make help' for more details" 

# Version of the sfnipmodsetup program for GDB (GNU debugger)
nipmodsetup.x:	$(SRCNIPTOMO)
	gcc -g $^ -I$$RSFSRC/include -I.. -L$$RSFSRC/lib -o $@ -lrsf -lm
	@echo "Program nipmodsetup.x for GDB built"
	gdb -q --command=gdb_script1 --args $@

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

.PHONY: clean nipmodsetup.x getparameter.x stereoniptomo.x

help:
	@echo "Use 'make' to run unit tests and intergration test"
	@echo "Use 'make progname.x' to generate a GDB version of progname for debug"
	@echo "Programs available:"
	@echo "\tgetparameter.x\n\tnipmodsetup.x\n\tstereoniptomo.x"
