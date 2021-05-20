import unittest
import subprocess
import numpy as np
import parseResults as pr
import json
import os


class TestElevate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """ Run case and store reference results """
        # setupClass is used instead of setUp so that code is only run once
        scriptDir = os.path.dirname(os.path.realpath(__file__)) + '/'
        cls.caseDir = scriptDir + 'elevateTest.case/'
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
        cls.CTavg_ref = np.average(CT[-15:])  # Over last 15 iterations

        data = np.loadtxt(cls.refResultsDir + \
                          'r01b01ForceDist00300.csv' + '.ref', \
                          skiprows=1)
        cls.secSpan_ref = data[:, 0]
        cls.secCL_ref = data[:, 1]
        cls.secArea_ref = data[:, 5]
        cls.secVel_ref = data[:, 6]
        cls.secAlpha_ref = data[:, 8]

    def testProcess(self):
        self.assertEqual(self.returncode, 0, 'Exit code is non-zero')

    def testParams(self):
        params = pr.getParams()
        for k in self.params_ref.keys():
            self.assertEqual(params[k], self.params_ref[k], \
                            'A parameter does not match')

    def testCT(self):
        """ Check CT history values """
        dataHist = np.loadtxt(self.ResultsDir + 'r01ForceNonDim.csv', \
                              skiprows=1)
        it = dataHist[:, 0]
        CT = dataHist[:, 1]
        niter = CT.size

        # Check recorded number of iterations
        if it[0] == 0:
            self.assertEqual(niter, 301)
        if it[0] == 1:
            self.assertEqual(niter, 300)

        # Average CT value
        CTavg = np.average(CT[-15:])  # Over last 15 iterations
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
