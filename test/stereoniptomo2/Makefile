# Makefile
# 
# Purpose: Run and control SConscripts and do model backup.
# 
# Site: https://dirack.github.io
# 
# Version 1.0
# 
# Programer: Rodolfo A C Neves (Dirack) 23/03/2021
# 
# Email: rodolfo_profissional@hotmail.com
# 
# License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

# Number of iterations in inversion
NI=10

# model files: Stacked section and parameters cube from CRE stack
PARCUBE=model/parametersCube.rsf
STKSEC=model/stackedSection.rsf

# Extract parameters from stacked section
NP=$(shell (<stackedSection.asc sfget n1 n2; echo "n1*n2") | bc)
N1=$(shell <stackedSection.asc sfget n1 parform=n)
D1=$(shell <stackedSection.asc sfget d1 parform=n)
O1=$(shell <stackedSection.asc sfget o1 parform=n)
N2=$(shell <stackedSection.asc sfget n2 parform=n)
D2=$(shell <stackedSection.asc sfget d2 parform=n)
O2=$(shell <stackedSection.asc sfget o2 parform=n)
PARAM=np=$(NP) n1=$(N1) d1=$(D1) o1=$(O1) n2=$(N2) d2=$(D2) o2=$(O2)

# Run the experiment
all:	stackedSection.asc parametersCube.asc
	scons -f SConscript1 $(PARAM)
	scons ns=$(shell wc -l pick.txt | cut -d" " -f1) ni=$(NI)
	scons -f SConscript2
	./viewNipSources.py

# Get help
help:
	@echo "Run:"
	@echo "'make picking' to do picking in stacked section"
	@echo "'make' to run stereoniptomo experiment"
	@echo "'make backup' to convert stacked section and parameters cube to ascii"

# Build an ascii file of the stacked section and parameters cube for backup
backup:
	<$(PARCUBE) sfdisfil format="%f " number=n > $(PARCUBE:.rsf=.asc)
	echo "n1=8 d1=1 o1=0 n2=$(NP) d2=1 o2=0 in=$(notdir $(PARCUBE:.rsf=.asc)) data_format=ascii_float" >> $(PARCUBE:.rsf=.asc)
	<$(STKSEC) sfdisfil format="%f " number=n > $(STKSEC:.rsf=.asc)
	echo "n1=$(N1) d1=$(D1) o1=$(O1) n2=$(N2) d2=$(D2) o2=$(O2) in=$(notdir $(STKSEC:.rsf=.asc)) data_format=ascii_float" >> $(STKSEC:.rsf=.asc)

picking:	stackedSection.asc
	scons -f SConscript1 $(PARAM) stackedSection.rsf
	<stackedSection.rsf sfipick > pick.txt

stackedSection.asc:	$(STKSEC:.rsf=.asc)
	cp $< $@

parametersCube.asc:	$(PARCUBE:.rsf=.asc)
	cp $< $@

clean:
	scons -c
	sfrm *.rsf
	rm *.asc *.txt

.PHONY:	clean
