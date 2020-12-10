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
SRCNIPTOMO= ../Mniptomo.c ../raytrace.c ../grid2.c ../atela.c

all:	$(DEP)

test_%:
	@echo $@
	@python -m unittest $@

integ:
	cd ../experiments/modelingAndStack/ && scons
	cd ../experiments/algorithmVelocityInversion/ && scons

# Version of the sfniptomo program for GDB (GNU debugger)
niptomo.x:	$(SRCNIPTOMO)
	gcc -g $^ -I$$RSFSRC/include -I.. -L$$RSFSRC/lib -o $@ -lrsf -lm
	@echo "Program niptomo.x for GDB built"

clean:
	rm *.x
