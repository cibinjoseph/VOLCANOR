import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

vlm = pd.read_csv('vlmResults.csv', \
                  delim_whitespace=True).to_dict(orient='list')
vsp = pd.read_csv('vspResults.csv', \
                  delim_whitespace=True).to_dict(orient='list')

denom = 0.2054*np.pi*74.18**2
vsp['secLift'] = np.array(vsp['CT_H'])*denom
vsp['secSpan'] = np.array(vsp['S'])*0.7+0.3
vsp['dx'] = np.array(vsp['Area'])/np.array(vsp['Chord'])
vsp['dLdx'] = vsp['secLift']/vsp['dx']

vsp_old = pd.read_csv('vspResults_old.csv', \
                  delim_whitespace=True).to_dict(orient='list')

denom = 0.2054*np.pi*74.18**2
vsp_old['secLift'] = np.array(vsp_old['CT_H'])*denom
vsp_old['secSpan'] = np.array(vsp_old['S'])*0.7+0.3
vsp_old['dx'] = np.array(vsp_old['Area'])/np.array(vsp_old['Chord'])
vsp_old['dLdx'] = vsp_old['secLift']/vsp_old['dx']

vlm['dx'] = np.array(vlm['secArea'])/np.array(vlm['secChord'])
vlm['dLdx'] = vlm['secLift']/vlm['dx']

plt.plot(vlm['secSpan'], vlm['dLdx'], 'r', label='VLM')
plt.plot(vsp['secSpan'], vsp['dLdx'], 'b', label='OpenVSP')
plt.plot(vsp_old['secSpan'], vsp_old['dLdx'], 'k', label='OpenVSP older version')
plt.legend()
plt.xlabel('r')
plt.ylabel('$dT/dr$')
plt.title('Venus rotor')

plt.show()
