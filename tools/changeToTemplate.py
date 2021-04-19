#!/usr/bin/python3
""" Converts input file comments to format given in template file """
""" The template file is assumed to be the last in the arguments """
""" All inputs of the input file is retained """

import sys
import time

def changeFile(inputFile, templateFile):
    with open(inputFile, 'r') as inptF:
        inptLines = inptF.readlines()
    # Remove comments from input file
    onlyInputs = []
    for line in inptLines:
        if line[0] != '#':
            onlyInputs.append(line)

    with open(templateFile, 'r') as tmplF:
        tmplLines = tmplF.readlines()

    with open(inputFile, 'w') as inptF:
        for line in tmplLines:
            if line[0] == '#':
                inptF.write(line)
            else:
                inptF.write(onlyInputs[0])
                del(onlyInputs[0])
    if len(onlyInputs) != 0:
        print('Warning: The following input lines have not been used')
        print(onlyInputs)
    return

args = sys.argv

if len(args) == 1:
    print('Usage: changeFile file1 [file2 ...] templateFile')
    print('Converts input file comments to format given in template file')
    print('The template file is assumed to be the last in the arguments')
    print('All inputs of the input file is retained')
    quit()

inputFileList = args[1:-1]
templateFile = args[-1]

print('Template file is: ' + templateFile)
print('Changing files in 5 seconds. Ctrl+C to stop')
time.sleep(5)

for inputFile in inputFileList:
    print(inputFile)
    changeFile(inputFile, templateFile)
