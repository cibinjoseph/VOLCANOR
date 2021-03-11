##################################
#  USE CMAKE GENERATED MAKEFILE  #
#    FOR BUILDING EXECUTABLE     #
##################################

resultspath=./Results

all:
	make run

fileclean:
	-rm -f $(resultspath)/*.plt
	-rm -f $(resultspath)/*.curve
	-rm -f $(resultspath)/*.dat
	-rm -f $(resultspath)/*.txt
	-rm -f volcanor.log

run:
	mkdir -p $(resultspath)
	make fileclean
	./volcanor 2>&1 | tee volcanor.log
