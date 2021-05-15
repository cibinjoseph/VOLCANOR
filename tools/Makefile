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
	-rm -f $(resultspath)/*.csv
	-rm -f $(resultspath)/*.json
	-rm -f $(resultspath)/*.txt
	-rm -f volcanor.log

run:
	mkdir -p $(resultspath)
	mkdir -p $(restartpath)
	make fileclean
	time -p ./bin/volcanor 2>&1 | tee volcanor.log
