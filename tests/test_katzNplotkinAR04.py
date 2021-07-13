import unittest
import subprocess
import numpy as np
import parseResults as pr
import json
import os
import pandas as pd


class TestAR04(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """ Run case and store reference results """
        # setupClass is used instead of setUp so that code is only run once
        scriptDir = os.path.dirname(os.path.realpath(__file__)) + '/'
        cls.caseDir = scriptDir + 'katzNplotkin-AR04.case/'
        cls.ResultsDir = cls.caseDir + 'Results/'
        cls.refResultsDir = cls.caseDir + 'referenceResults/'
        p = subprocess.run(['make'], cwd=cls.caseDir)
        cls.returncode = p.returncode

        pr.ResultsDir = cls.ResultsDir
        # Read reference params file
        with open(cls.refResultsDir + 'r01Params.json' + '.ref') as fh:
            cls.params_ref = json.load(fh)

        # Read reference results
        dataHist = np.loadtxt(cls.refResultsDir + \
                              'r01ForceNonDim.csv' + '.ref', \
                              skiprows=1)
        CT = dataHist[:, 1]
        cls.CTlast_ref = CT[-1]
        cls.CTavg_ref = np.average(CT[-10:])  # Over last 10 iterations

        dataPD = pd.read_table(cls.refResultsDir + \
                          'r01b01ForceDist00160.csv' + '.ref', \
                             delim_whitespace=True)
        data = dataPD.to_dict(orient='list')

        cls.secSpan_ref = data['secSpan']
        cls.secCL_ref = data['secCL']
        cls.secArea_ref = data['secArea']
        cls.secVel_ref = data['secVel']
        cls.secAlpha_ref = data['secAlpha']

    def testProcess(self):
        self.assertEqual(self.returncode, 0, 'Exit code is non-zero')

    def testParams(self):
        params = pr.getParams()
        for k in self.params_ref.keys():
            if isinstance(self.params_ref[k], float): 
                self.assertAlmostEqual(params[k], self.params_ref[k], \
                                places=6, \
                                msg='Parameter ' + k + ' does not match')
            else:
                self.assertEqual(params[k], self.params_ref[k], \
                                'Parameter ' + k + ' does not match')

    def testCT(self):
        """ Check CT history values """
        dataHist = np.loadtxt(self.ResultsDir + 'r01ForceNonDim.csv', \
                              skiprows=1)
        it = dataHist[:, 0]
        CT = dataHist[:, 1]

        # Check recorded number of iterations
        if it[0] == 0:
            self.assertEqual(CT.size, 161)
        if it[0] == 1:
            self.assertEqual(CT.size, 160)

        # Average CT value
        CTavg = np.average(CT[-10:])  # Over last 10 iterations
        self.assertAlmostEqual(CTavg, self.CTavg_ref, places=6, \
                               msg='Average CT does not match')

        # Final CT value
        CTlast = CT[-1]
        self.assertAlmostEqual(CTlast, self.CTlast_ref, places=6, \
                               msg='Last CT does not match')

    def testSect(self):
        """ Check sectional quantities """
        data, filename_ = pr.getForceDist()

        for val, val_ref in zip(data['secSpan'], self.secSpan_ref):
            self.assertAlmostEqual(val, val_ref, places=6, \
                                  msg='Sectional span does not match')
        for val, val_ref in zip(data['secCL'], self.secCL_ref):
            self.assertAlmostEqual(val, val_ref, places=6, \
                                  msg='Sectional CL does not match')
        for val, val_ref in zip(data['secArea'], self.secArea_ref):
            self.assertAlmostEqual(val, val_ref, places=6, \
                                  msg='Sectional area does not match')
        for val, val_ref in zip(data['secVel'], self.secVel_ref):
            self.assertAlmostEqual(val, val_ref, places=6, \
                                  msg='Sectional vel does not match')
        for val, val_ref in zip(data['secAlpha'], self.secAlpha_ref):
            self.assertAlmostEqual(val, val_ref, places=6, \
                                  msg='Sectional alpha does not match')



if __name__ == "__main__":
    unittest.main()
