#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from nonstdlib import *

nucleotide = input('Which nucleotide? ')
df = pd.read_csv('merged_spectra.txt', sep='\t', header=0)
pn = df[df['sample'] == 'protein'].Abs.data
p = np.array(df[df['sample'] == 'unfolded_protein_ref'].Abs.data)
n = np.array(df[df['sample'] == nucleotide].Abs.data)
nm = np.array(df[df['sample'] == nucleotide].nm.data)

def combine_spectra(spectra, *concentrations):
    return sum(k * s for s, k in zip(spectra, concentrations))

spectra = np.vstack([p,n])
conc = curve_fit(combine_spectra, spectra, pn, p0=[1,1])[0]
print(conc)


nm280 = np.nonzero(nm == 280.0)[0][0]
nm260 = np.nonzero(nm == 260.0)[0][0]
protein_corr = conc[0] * p
nuc_corr = conc[1] * n

plt.plot(nm, pn, 'k-', label='Protein + Nuc')
plt.plot(nm, protein_corr, 'b-', label='Protein (k={})'.format(sci(conc[0])))
plt.plot(nm, nuc_corr, 'r-', label='Nuc (k={})'.format(sci(conc[1])))
plt.plot(nm, combine_spectra(spectra, *conc), 'k--', label='Best Fit')
plt.legend(loc='best')
plt.xlim(min(nm), max(nm))
plt.show()

print('Protein conc/uM:', 1e6*(protein_corr[nm280]/(29910*1)))
print('260/280 ratio: ', pn[nm260]/pn[nm280])

