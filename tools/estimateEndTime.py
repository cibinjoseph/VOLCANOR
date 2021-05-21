#!/usr/bin/python3
""" Pass log file that contains "iter MaxIter" being recorded """
""" to estimate time taken for MaxIter iterations to be completed """

import os
import sys
import subprocess
import datetime
import numpy as np

polyDeg = 2

def getdata(filename):
    line = subprocess.check_output(['tail', '-1', filename]).split()
    if len(line) == 2:
        itr, tot = line
    else:
        itr, tot = 0, 0
    mtime = os.path.getmtime(filename)
    return int(itr), int(tot), mtime

if __name__ == "__main__":
    try:
        filename = sys.argv[1]
    except:
        raise ValueError('No file supplied as argument')

    mList = (polyDeg+1)*[0]
    iList = (polyDeg+1)*[0]
    counter = 0
    while True:
        itr, tot, mtime = getdata(filename)
        if itr > iList[-1]:
            counter += 1
            _ = mList.pop(0)
            mList.append(mtime)
            _ = iList.pop(0)
            iList.append(itr)

            if counter > polyDeg:
                timeAt = np.polyval( \
                            np.polyfit(iList, mList, deg=polyDeg), tot)
                dt = timeAt - mtime
                estTime = datetime.datetime.fromtimestamp(timeAt)
                print(estTime)
