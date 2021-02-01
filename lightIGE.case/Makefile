##################################
#  USE CMAKE GENERATED MAKEFILE  #
#    FOR BUILDING EXECUTABLE     #
##################################

resultspath=./Results
restartpath=./Restart

all:
	make run

fileclean:
	-rm -f $(resultspath)/*.plt
	-rm -f $(resultspath)/*.curve
	-rm -f $(resultspath)/*.dat
	-rm -f $(resultspath)/*.txt
	-rm -f status.txt

run:
	mkdir -p $(resultspath) ${restartpath}
	make fileclean
	./volcanor
