#!/usr/bin/python3
""" Converts input file comments to format given in template file """
""" The template file is assumed to be the last in the arguments """
""" All inputs of the input file is retained """

import sys
import time
from pathlib import Path


def getInputLines(inputFile):
    """ Returns lines ignoring comments """
    with open(inputFile, 'r') as fh:
        lines = fh.readlines()
    inputLines = []
    for line in lines:
        if line[0] != '#':
            inputLines.append(line)

    return inputLines

def changeFile(inputFile, templateFile, retainTemplateLines=[]):
    inptLines = getInputLines(inputFile)
    tmplLines = getInputLines(templateFile)

    # Ensure no. of lines match
    if retainTemplateLines == []:
        if len(inptLines) != len(tmplLines):
            raise ValueError('No. of lines do not match template in: ' \
                             + inputFile)

    with open(templateFile, 'r') as tmplF:
        tmplLines = tmplF.readlines()

    i = 0
    with open(inputFile, 'w') as inptF:
        for line in tmplLines:
            # Write comments from template file
            if line[0] == '#':
                inptF.write(line)
            else:
                # If input line,
                i = i + 1
                # check if specified to be retained
                if i in retainTemplateLines:
                    inptF.write(line)
                # else copy line from input file
                else:
                    inptF.write(inptLines[0])
                    del(inptLines[0])
    return

def main():
    scriptPath = str(Path(__file__).parent)
    args = sys.argv

    if len(args) == 1:
        print('Usage: changeFile [-l 1,12,...] file1 [file2 ...]')
        print('Converts input file comments to format given in template file')
        print('The template file is guessed from input filename.')
        print('All inputs of the input file is retained.')
        print('Input line no.s given after -l are copied from the template.')
        quit()

    if args[1] == '-l':
        retainTemplateLines = [int(i) for i in args[2].split(',')]
        inputFileList = args[3:]
    else:
        retainTemplateLines = []
        inputFileList = args[1:]

    for inputFile in inputFileList:
        print(inputFile)
        if 'config.in' in inputFile:
            templateFile = scriptPath + '/template.case/config.in'
        elif 'geom' in inputFile:
            # Assume geom file
            templateFile = scriptPath + '/template.case/geom01.in'
        else:
            raise ValueError('Invalid input filename')

        print('Template: ' + templateFile)
        changeFile(inputFile, templateFile, retainTemplateLines)

if __name__ == '__main__':
    main()
