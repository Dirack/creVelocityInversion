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

DEP= ../raytrace.c ../grid2.c ../atela.c
SRCNIPTOMO= ../Mnipmodsetup.c $(DEP)
SRCGETPARAM= ../Mgetparameter.c
SRCSTEREO = ../Mstereoniptomo.c ../tomography.c $(DEP) ../vfsacrsnh_lib.c

# model files: Stacked section and parameters cube from CRE stack
PARCUBE=model2/parametersCube.rsf
STKSEC=model2/stackedSection.rsf

# Extract parameters from stacked section
NP=$(shell (<stereoniptomo/stackedSection.asc sfget n1 n2; echo "n1*n2") | bc)
N1=$(shell <stereoniptomo/stackedSection.asc sfget n1 parform=n)
D1=$(shell <stereoniptomo/stackedSection.asc sfget d1 parform=n)
O1=$(shell <stereoniptomo/stackedSection.asc sfget o1 parform=n)
N2=$(shell <stereoniptomo/stackedSection.asc sfget n2 parform=n)
D2=$(shell <stereoniptomo/stackedSection.asc sfget d2 parform=n)
O2=$(shell <stereoniptomo/stackedSection.asc sfget o2 parform=n)
PARAM=np=$(NP) n1=$(N1) d1=$(D1) o1=$(O1) n2=$(N2) d2=$(D2) o2=$(O2)

# Function to call gdb scripts
# 1 is the number of the script
# 2 are the parameters
# 3 is the file to debug
define GDB
	gdb -q --command=gdb_scripts/gdb_script$(1) $(2) --args $(3)
endef

all:	nipmodsetup.x getparameter.x stereoniptomo.x
	@echo "Running integration tests, run 'make help' for more details"
	$(call GDB,1,-ex "r" -ex "q",nipmodsetup.x)
	$(call GDB,2,-ex "r" -ex "q",getparameter.x)
	$(call GDB,3,-ex "r" -ex "q",stereoniptomo.x)
	

# Version of the sfnipmodsetup program for GDB (GNU debugger)
nipmodsetup:	nipmodsetup.x GDB1

nipmodsetup.x:	$(SRCNIPTOMO)
	@echo "Program nipmodsetup.x for GDB built"
	gcc -g $^ -I$$RSFSRC/include -I.. -L$$RSFSRC/lib -o $@ -lrsf -lm

GDB1:
	@echo "Running GDB script for nipmodsetup.x"
	$(call GDB,1,,nipmodsetup.x)

# Version of the sfgetparameter program for GDB (GNU debugger)
getparameter:	getparameter.x GDB2

getparameter.x:	setup1 $(SRCGETPARAM)
	@echo "Program getparameter.x for GDB built"
	gcc -g $(SRCGETPARAM) -I$$RSFSRC/include -I.. -L$$RSFSRC/lib -o $@ -lrsf -lm

GDB2:
	@echo "Running GDB script for getparameter.x"
	$(call GDB,2,,getparameter.x)

setup1:
	@echo "Setup input files for getparameter.x test..."
	cd stereoniptomo && scons -f SConscript1 $(PARAM) && scons model_0.rsf && scons sz.rsf && scons vspline_0.rsf

# Version of the sfstereoniptomo program for GDB (GNU debugger)
stereoniptomo:	stereoniptomo.x GDB3

stereoniptomo.x:	setup1 $(SRCSTEREO)
	@echo "Program stereoniptomo.x for GDB built"
	gcc -DGDB_DEBUG -g $(SRCSTEREO) -I$$RSFSRC/include -I.. -L$$RSFSRC/lib -o $@ -lrsf -lm

GDB3:
	@echo "Running GDB script for stereoniptomo.x"
	$(call GDB,3,-ex "b 117" -ex "r",stereoniptomo.x)

clean:
	rm *.x

.PHONY: clean nipmodsetup.x getparameter.x stereoniptomo.x

help:
	@echo "Use 'make' to run unit tests and intergration test"
	@echo "Use 'make progname' to generate a GDB version of progname for debug"
	@echo "Programs available:"
	@echo "\tgetparameter\n\tnipmodsetup\n\tstereoniptomo"
