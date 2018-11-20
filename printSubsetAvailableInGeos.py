#!/usr/bin/env python3

import h5py, sys, os  
import numpy as np

def main(whichSpectrum, idx):
    thisScriptPath = os.path.split(os.path.realpath(__file__))[0]
    coefPath = os.path.join(thisScriptPath,'etc')
    if(whichSpectrum == 'CRISFSR'):
        h5 = h5py.File( os.path.join(coefPath,'rtcoef_jpss_0_cris-fsr_so2.H5') ) 
    elif(whichSpectrum == 'CRISNPP'):
        h5 = h5py.File( os.path.join(coefPath, 'rtcoef_jpss_0_cris_so2.H5') )
    elif(whichSpectrum == 'IASI'):
        h5 = h5py.File( os.path.join(coefPath, 'rtcoef_metop_2_iasi_so2.H5') )
    elif(whichSpectrum == 'AIRS'): 
        h5 = h5py.File( os.path.join(coefPath, 'rtcoef_eos_2_airs_so2.H5') )
    print("Printing Channel, Wavenumber, Wavelength for "+whichSpectrum)   
    if('CRIS' in whichSpectrum or 'AIRS' in whichSpectrum): 
        for l in list(h5['INSTRUMENT']['FILTERS'].keys()):
            #print(list(h5['INSTRUMENT']['FILTERS'][l].keys()))
            if(int(l) in idx):
                print("{:4d} {:10.4f} {:10.4f}".format(np.asarray(h5['INSTRUMENT']['FILTERS'][l]['ICHAN'])[0], np.asarray(h5['INSTRUMENT']['FILTERS'][l]['CWV'])[0], 10000.0/np.asarray(h5['INSTRUMENT']['FILTERS'][l]['CWV'])[0]) )
    elif( 'IASI' in whichSpectrum):
        for i,wvno in enumerate(h5['INSTRUMENT']['INTERF_WV']):
            if (int(i+1) in idx):
                print("{:4d} {:10.4f} {:10.4f}".format(int(i+1), wvno, 10000.0/wvno) ) 
if __name__ == "__main__":
    fn = sys.argv[1]
    idx = []
    with open(fn) as f:
        lines = f.readlines()
    for l in lines:
        idx.append( int(l.split()[1]) )
    if('airs' in fn): whichSpectrum = "AIRS"
    elif('cris' in fn and 'fsr' in fn): whichSpectrum = "CRISFSR"
    elif('iasi' in fn): whichSpectrum = "IASI"
    elif('cris' in fn and not 'fsr' in fn): whichSpectrum = "CRISNPP"
    main(whichSpectrum, idx)
