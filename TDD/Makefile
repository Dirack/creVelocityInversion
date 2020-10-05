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

all:	$(DEP)

test_%:
	@echo $@
	@python -m unittest $@
