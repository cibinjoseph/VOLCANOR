##################################
#  USE CMAKE GENERATED MAKEFILE  #
##################################

resultspath=./Results

all:
	@echo
	@<&2 echo "  Use cmake generated Makefile in build directory"
	@echo

fileclean:
	-rm -f $(resultspath)/*.plt
	-rm -f $(resultspath)/*.curve
	-rm -f $(resultspath)/*.dat
	-rm -f $(resultspath)/*.txt
	-rm -f status.txt
