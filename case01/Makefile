##################################
#  USE CMAKE GENERATED MAKEFILE  #
##################################

resultspath=./Results

all:
	make warning

warning:
	@echo
	@<&2 echo "  Use cmake generated Makefile in build/"
	@echo

fileclean:
	-rm -f $(resultspath)/*.plt
	-rm -f $(resultspath)/*.curve
	-rm -f $(resultspath)/*.dat
	-rm -f $(resultspath)/*.txt
	-rm -f status.txt

run:
	mkdir -p $(resultspath)
	make fileclean
	./volcanor
