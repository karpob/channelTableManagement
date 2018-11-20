#!/usr/bin/env python3

import h5py, sys, os  
import numpy as np
def nucapsIdxToSpectralIdx(whichSpectrum, idxNucaps):
    if(whichSpectrum == 'CRISFSR'):
        band1Idx = idxNucaps[ np.where( idxNucaps <= 715 ) ] - 2
        band2Idx = idxNucaps[ np.where( (idxNucaps > 715 ) & (idxNucaps<=1584))] - 6
        band3Idx = idxNucaps[ np.where( idxNucaps > 1584 ) ] - 10
        b1, b2, b3 = band1Idx.tolist(), band2Idx.tolist(), band3Idx.tolist()
        b1.extend(b2)
        b1.extend(b3)
    elif(whichSpectrum == 'CRISNPP'):
        band1Idx = idxNucaps[ np.where( idxNucaps <= 715 ) ] - 2
        band2Idx = idxNucaps[ np.where( (idxNucaps > 715 ) & (idxNucaps<=864))] - 6
        band3Idx = idxNucaps[ np.where( idxNucaps > 864 ) ] - 10
        b1, b2, b3 = band1Idx.tolist(), band2Idx.tolist(), band3Idx.tolist()
        b1.extend(b2)
        b1.extend(b3)
    else: b1 = idxNucaps
    return np.asarray(b1)

def spectralIdxToNucapsIdx(whichSpectrum, idx):
    if(whichSpectrum == 'CRISFSR'):
        band1Idx = idx[ np.where( idx<=713 )] + 2
        band2Idx = idx[ np.where( (idx>713) & (idx<=1578) )] + 6
        band3Idx = idx[ np.where( idx > 1578)] + 10
        b1, b2, b3 = band1Idx.tolist(), band2Idx.tolist(), band3Idx.tolist()
        b1.extend(b2)
        b1.extend(b3)
    elif(whichSpectrum == 'CRISNPP'):
        band1Idx = idx[ np.where( idx<=713 )] + 2
        band2Idx = idx[ np.where( (idx>713) & (idx<=864) )] + 6
        band3Idx = idx[ np.where( idx > 864)] + 10
        b1, b2, b3 = band1Idx.tolist(), band2Idx.tolist(), band3Idx.tolist()
        b1.extend(b2)
        b1.extend(b3)
    else:b1 = idx
    
    return np.asarray(b1)


def main(whichSpectrum):
    thisScriptPath = os.path.split(os.path.realpath(__file__))[0]
    coefPath = os.path.join(thisScriptPath,'../etc')
    if(whichSpectrum == 'CRISFSR'):
        h5 = h5py.File( os.path.join(coefPath,'rtcoef_jpss_0_cris-fsr_so2.H5') ) 
    elif(whichSpectrum == 'CRISNPP'):
        h5 = h5py.File( os.path.join(coefPath, 'rtcoef_jpss_0_cris_so2.H5') )
    elif(whichSpectrum == 'IASI'):
        h5 = h5py.File( os.path.join(coefPath, 'rtcoef_metop_2_iasi_so2.H5') )
    elif(whichSpectrum == 'AIRS'): 
        h5 = h5py.File( os.path.join(coefPath, 'rtcoef_eos_2_airs_so2.H5') )
    print("Printing Channel, NucapsIdx, Wavenumber, Wavelength for "+whichSpectrum)   
    if('CRIS' in whichSpectrum or 'AIRS' in whichSpectrum): 
        spectralIdx = np.arange(1,len(h5['INSTRUMENT']['FILTERS'])+1)
        nucapsIdx = spectralIdxToNucapsIdx(whichSpectrum, spectralIdx)
        for i,l in enumerate(list(h5['INSTRUMENT']['FILTERS'].keys())):
            print("{:4d} {:4d} {:10.4f} {:10.4f}".format(np.asarray(h5['INSTRUMENT']['FILTERS'][l]['ICHAN'])[0],nucapsIdx[i], np.asarray(h5['INSTRUMENT']['FILTERS'][l]['CWV'])[0], 10000.0/np.asarray(h5['INSTRUMENT']['FILTERS'][l]['CWV'])[0]) )
    elif( 'IASI' in whichSpectrum):
        for i,wvno in enumerate(h5['INSTRUMENT']['INTERF_WV']):
            print("{:4d} {:4d} {:10.4f} {:10.4f}".format(int(i+1), int(i+1), wvno, 10000.0/wvno) ) 
if __name__ == "__main__":
    whichSpectrum = 'CRISFSR'
    if(len(sys.argv)>1):
        whichSpectrum = sys.argv[1].upper()
    main(whichSpectrum)
