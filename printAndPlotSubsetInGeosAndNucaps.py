#!/usr/bin/env python3

import h5py, sys, os  
import numpy as np
from matplotlib import pyplot as plt
from collections import OrderedDict
ecmwfOzoneChannels = {}
ecmwfOzoneChannels['cris'] =[590, 594, 598,  602,  606,  610, 614, 618, 622, 626, 645, 649, 653, 657, 661, 665]
# assume they didn't do anything new for FSR.
ecmwfOzoneChannels['cris-fsr'] =[590, 594, 598,  602,  606,  610, 614, 618, 622, 626, 645, 649, 653, 657, 661, 665]
ecmwfOzoneChannels['airs'] =[1012, 1019, 1024, 1030, 1038, 1048, 1069, 1079, 1082, 1088, 1090, 1092, 1104, 1111, 1115, 1116, 1119, 1120, 1123]
ecmwfOzoneChannels['iasi'] = [1479, 1509, 1513, 1521, 1536, 1574, 1578, 1579, 1585, 1587, 1626, 1639, 1643, 1652, 1658, 1671] 
def main_plot(freqMeta, nucapsTable, geosAssimilated, instrument, platform ):
    """
    Go through and generate some plots.
    """
    print("number in geos",len(geosAssimilated))
    #plot everything in nucaps table

    plotNucapsRetrievalBands(freqMeta, nucapsTable, instrument, platform)

    # zoom in and plot near the ozone 9.6 um band.
    plotNucapsRetrievalBands(freqMeta, nucapsTable, instrument, platform, zoom =np.asarray([980.0,1080]))
    plotNucapsRatios(freqMeta, nucapsTable, instrument, platform)
    plotNucapsRatios(freqMeta, nucapsTable, instrument, platform, zoom = np.asarray([980.0,1080]))
    plotGeosAssimilated(nucapsTable, geosAssimilated, instrument, platform)   
    plotEcmwfOzone(nucapsTable, instrument, platform, zoom = np.asarray([980.0,1080]))

def plotEcmwfOzone(nucapsTable, instrument, platform, zoom = []):
    Tb = np.asarray(nucapsTable['BT'])
    Tb[ np.where(Tb <= 50) ] = np.nan
     
    fig = plt.figure(figsize =(10,5))
    ax = plt.subplot(111)
    ecmwfAssimilatedTransformed = instrumentIdxToNucapsIdx(instrument, np.asarray(ecmwfOzoneChannels[instrument]))
    for i in ecmwfAssimilatedTransformed: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='orange',label='ECMWF Assimilated Ozone')
    ax.plot(nucapsTable['freq'], Tb,'k')
    ax.set_xlabel('Wavenumber [cm$^{\mathrm{-1}}$]')
    ax.set_ylabel('Brightness Temperature [K]')
    #shrink axes to fit legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # get just the first instance of each axvline
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    if(len(zoom)>0): ax.set_xlim(zoom)
    # drop in legend.
    ax.legend(by_label.values(), by_label.keys(),loc='center left', bbox_to_anchor=(1, 0.5) )

    instrument = os.path.join('plots',instrument)
    if(len(zoom)>0):plt.savefig(instrument+'_'+platform+'_ecmwf_bands_{}_{}.png'.format(zoom[0],zoom[1]))
    else: plt.savefig(instrument+'_'+platform+'_ecmwf_bands.png')
    
def plotGeosAssimilated(nucapsTable, geosAssimilated, instrument, platform):
    Tb = np.asarray(nucapsTable['BT'])
    Tb[ np.where(Tb <= 50) ] = np.nan
     
    fig = plt.figure(figsize =(10,5))
    ax = plt.subplot(111)
    geosAssimilatedTransformed = instrumentIdxToNucapsIdx(instrument, np.asarray(geosAssimilated))
    for i in geosAssimilatedTransformed: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='red',label='Geos Assimilated')
    ax.plot(nucapsTable['freq'], Tb,'k')
    ax.set_xlabel('Wavenumber [cm$^{\mathrm{-1}}$]')
    ax.set_ylabel('Brightness Temperature [K]')
    #shrink axes to fit legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # get just the first instance of each axvline
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    # drop in legend.
    ax.legend(by_label.values(), by_label.keys(),loc='center left', bbox_to_anchor=(1, 0.5) )
    instrument = os.path.join('plots',instrument)
    plt.savefig(instrument+'_'+platform+'_nucaps_bands_geos.png')
    
def plotNucapsRatios(freqMeta, nucapsTable, instrument, platform, zoom = [], zoomy=[]):
    fig = plt.figure(figsize =(10,5))
    ax = plt.subplot(111)
    ax.plot(nucapsTable['freq'], nucapsTable['Dhno3'], color='purple', linestyle=':', label='D$_{\mathrm{HNO}_3}$')
    ax.plot(nucapsTable['freq'], nucapsTable['Dn2o'], color='magenta', linestyle=':', label='D$_{\mathrm{N}_2\mathrm{O}}$')
    #ax.plot(nucapsTable['freq'],nucapsTable['Dso2'], color='gray' )
    for i in freqMeta['ivozon']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='orange',label='O$_3$')
    for i in freqMeta['ivHNO3']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='purple',label='HNO$_3$')
    for i in freqMeta['ivch4']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='goldenrod',label='CH$_4$')
    for i in freqMeta['ivwatr']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='blue', label='H$_2$O')
    for i in freqMeta['ivN2O']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='magenta',label='N$_2$O')
    for i in freqMeta['ivSO2']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='gray',label='SO$_2$')
    for i in freqMeta['ivco']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='yellow', label = 'CO')
    for i in freqMeta['ivch4']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='brown', label='CH$_4$')
    for i in freqMeta['ivtemp']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='red', label = 'Temperature')
    for i in freqMeta['ivtemp2']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='red', label = 'Temperature')
    for i in freqMeta['ivstrat']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='maroon',label = 'Stratospheric')

    ax.plot(nucapsTable['freq'], nucapsTable['Dhno3'], color='purple', linestyle=':', label='D$_{\mathrm{HNO}_3}$')
    ax.plot(nucapsTable['freq'], nucapsTable['Dn2o'], color='magenta', linestyle=':', label='D$_{\mathrm{N}_2\mathrm{O}}$')
    ax.plot(nucapsTable['freq'], nucapsTable['wat'], color='blue', linestyle=':', label='wat')
    ax.plot(nucapsTable['freq'], nucapsTable['ozo'], color='orange', linestyle=':', label='ozo')
    ax.plot(nucapsTable['freq'], nucapsTable['ch4'], color='brown', linestyle=':', label='ch4')
    ax.plot(nucapsTable['freq'], nucapsTable['co'], color='yellow', linestyle=':', label='co')
    'idx    freq mod CTTUSWOdMmO B    NEDT TUNING RTAERR     BT    O-C |   fix   wat   ozo   ch4    co Dhno3  Dn2o  Dso2 | P_tot P_fix P_wat P_ozo P_ch4  P_co'
    ax.set_xlabel('Wavenumber [cm$^{\mathrm{-1}}$]')
    ax.set_ylabel('Ratio')
    #shrink axes to fit legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # get just the first instance of each axvline
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    if(len(zoom)>0): ax.set_xlim(zoom)
    if(len(zoomy)>0): ax.set_ylim(zoomy)
    # drop in legend.
    ax.legend(by_label.values(), by_label.keys(),loc='center left', bbox_to_anchor=(1, 0.5) )

    instrument = os.path.join('plots',instrument)
    if(len(zoom)>0):plt.savefig(instrument+'_'+platform+'_nucaps_ratios_{}_{}.png'.format(zoom[0],zoom[1]))
    else: plt.savefig(instrument+'_'+platform+'_nucaps_ratios.png')

def plotNucapsRetrievalBands(freqMeta, nucapsTable, instrument, platform, zoom = []):
    'idx    freq mod CTTUSWOdMmO B    NEDT TUNING RTAERR     BT    O-C |   fix   wat   ozo   ch4    co Dhno3  Dn2o  Dso2 | P_tot P_fix P_wat P_ozo P_ch4  P_co'
    Tb = np.asarray(nucapsTable['BT'])
    Tb[ np.where(Tb <= 50) ] = np.nan
     
    fig = plt.figure(figsize =(10,5))
    ax = plt.subplot(111)
    for i in freqMeta['ivozon']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='orange',label='O$_3$')
    for i in freqMeta['ivHNO3']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='purple',label='HNO$_3$')
    for i in freqMeta['ivch4']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='goldenrod',label='CH$_4$')
    for i in freqMeta['ivwatr']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='blue', label='H$_2$O')
    for i in freqMeta['ivN2O']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='magenta',label='N$_2$O')
    for i in freqMeta['ivSO2']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='gray',label='SO$_2$')
    for i in freqMeta['ivco']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='yellow', label = 'CO')
    for i in freqMeta['ivch4']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='brown', label='CH$_4$')
    for i in freqMeta['ivtemp']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='red', label = 'Temperature')
    for i in freqMeta['ivtemp2']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='red', label = 'Temperature')
    for i in freqMeta['ivstrat']: ax.axvline(x = nucapsTable['freq'][int(i)-1],color='maroon',label = 'Stratospheric')
    ax.plot(nucapsTable['freq'], Tb,'k')
    ax.set_xlabel('Wavenumber [cm$^{\mathrm{-1}}$]')
    ax.set_ylabel('Brightness Temperature [K]')
    #shrink axes to fit legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # get just the first instance of each axvline
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    if(len(zoom)>0): ax.set_xlim(zoom)
    # drop in legend.
    ax.legend(by_label.values(), by_label.keys(),loc='center left', bbox_to_anchor=(1, 0.5) )

    instrument = os.path.join('plots',instrument)
    if(len(zoom)>0):plt.savefig(instrument+'_'+platform+'_nucaps_bands_{}_{}.png'.format(zoom[0],zoom[1]))
    else: plt.savefig(instrument+'_'+platform+'_nucaps_bands.png')

def getInstrumentWavenumbersFromRttov(instrument, idxBufrSubset, geosAssimilated, nucapsTableFreq, idxNucapsOzone):    
    thisScriptPath = os.path.split(os.path.realpath(__file__))[0]
    coefPath = os.path.join(thisScriptPath,'etc')
 
    if(instrument == 'cris-fsr'):
        h5 = h5py.File( os.path.join(coefPath,'rtcoef_jpss_0_cris-fsr_so2.H5') ) 
    elif(instrument == 'cris'):
        h5 = h5py.File( os.path.join(coefPath, 'rtcoef_jpss_0_cris_so2.H5') )
    elif(instrument  == 'iasi'):
        h5 = h5py.File( os.path.join(coefPath, 'rtcoef_metop_2_iasi_so2.H5') )
    elif(instrument == 'airs'): 
        h5 = h5py.File( os.path.join(coefPath, 'rtcoef_eos_2_airs_so2.H5') )
    else: sys.exit("I don't know this instrument:{}".format(instrument))
    wavenumbers = [] 
    if('cris' in instrument or 'airs' in instrument):
        instrumentIdx = np.arange(1,len( list(h5['INSTRUMENT']['FILTERS'].keys()) ) + 1)
        channelNames = list(h5['INSTRUMENT']['FILTERS'].keys())
        nucapsIdx = instrumentIdxToNucapsIdx(instrument, instrumentIdx)
        for l in instrumentIdx:
            wavenumbers.append(h5['INSTRUMENT']['FILTERS'][channelNames[l-1]]['CWV'][0])
    elif( 'iasi' in instrument):
        for wvno in list(h5['INSTRUMENT']['INTERF_WV']):
            wavenumbers.append(wvno)

    outpath = os.path.join(coefPath, instrument+'_wavenumbers.h5')
    idxBufrSubsetInNucaps = instrumentIdxToNucapsIdx(instrument, idxBufrSubset)
    geosAssimilatedInNucaps = instrumentIdxToNucapsIdx(instrument, geosAssimilated)
    idxNucapsOzoneInInstrument = nucapsIdxToInstrumentIdx(instrument, idxNucapsOzone)
    idxEcmwfOzoneInNucaps = instrumentIdxToNucapsIdx(instrument, np.asarray(ecmwfOzoneChannels[instrument]))
    with h5py.File( outpath, "w" ) as f:
        dset = f.create_dataset("wavenumbers",data = np.asarray(wavenumbers))
        dset = f.create_dataset("wavenumbersNucaps",data = np.asarray(nucapsTableFreq) )
        dset = f.create_dataset("idxNucapsOzoneInNucaps",data = np.asarray(idxNucapsOzone) )
        dset = f.create_dataset("idxNucapsOzoneInInstrument",data = np.asarray(idxNucapsOzoneInInstrument) )
        dset = f.create_dataset("idxBufrSubset",data = np.asarray(idxBufrSubset) )
        dset = f.create_dataset("idxBufrSubsetInNucaps",data = np.asarray(idxBufrSubsetInNucaps) )
        dset = f.create_dataset("geosAssimilated",data = np.asarray(geosAssimilated) )
        dset = f.create_dataset("geosAssimilatedInNucaps",data = np.asarray(geosAssimilatedInNucaps) )
        dest = f.create_dataset("idxEcmwfOzoneInInstrument", data = np.asarray(ecmwfOzoneChannels[instrument]) )
        dest = f.create_dataset("idxEcmwfOzoneInNucaps", data = np.asarray(idxEcmwfOzoneInNucaps)) 
    return np.asarray(wavenumbers)
      
def main_print(instrument, nucapsTable, idxBufrSubset, idxNucapsInInstrumentIdx, freqMeta, instrumentWavenumbersFromRttov):
    nchans = 0
    thisScriptPath = os.path.split(os.path.realpath(__file__))[0]
    coefPath = os.path.join(thisScriptPath,'etc')
    print("Printing Common Ozone Subset.")
    print("Printing Channel, Nucaps idx, Wavenumber (from RTTOV), Wavenumber (from NUCAPS), Wavelength for "+instrument)   
    if('cris' in instrument or 'airs' in instrument):
        instrumentWavenumberIndex = np.arange(1,len( instrumentWavenumbersFromRttov ) + 1)
        nucapsIdx = instrumentIdxToNucapsIdx(instrument, instrumentWavenumberIndex)
        for l in instrumentWavenumberIndex:
            if(l in idxBufrSubset and l in idxNucapsInInstrumentIdx ):
                wavenumber = instrumentWavenumbersFromRttov[l-1]
                print("{:4d} {:4d} {:10.4f} {:10.4f}  {:10.4f}".format(l, nucapsIdx[l-1],wavenumber, freqMeta['freqozon'][ np.where(l == idxNucapsInInstrumentIdx)[0][0] ], 10000.0/wavenumber) )
                nchans+=1
    elif( 'iasi' in instrument):
        for i,wvno in enumerate(instrumentWavenumbersFromRttov):
            if (int(i+1) in idxBufrSubset and i+1 in idxNucapsInInstrumentIdx):
                print("{:4d} {:4d} {:10.4f} {:10.4f} {:10.4f}".format(int(i+1), int(i+1), wvno, freqMeta['freqozon'][ np.where( i+1 == idxNucapsInInstrumentIdx)[0][0] ], 10000.0/wvno) ) 
                nchans+=1
    print ( "Number of ozone channels in common subset between NUCAPS and BUFR subset: {:d}".format(nchans) )
    """
    print ("Sanity Check of Wavenumber and index transforms.")
    print ("idx from BUFR subset, NUCAPS idx from BUFR subset passed through transform, wavenumber from rttov, wavenumber from NUCAPS") 

    nucapsIdxFromBufrSubset = instrumentIdxToNucapsIdx(instrument, np.asarray(idxBufrSubset))
    for i,ibuf in enumerate(idxBufrSubset):
        print("{:d} {:d} {:10.4f} {:10.4f}".format(ibuf, nucapsIdxFromBufrSubset[i], instrumentWavenumbersFromRttov[ibuf-1], nucapsTable['freq'][ nucapsIdxFromBufrSubset[i]-1 ] ) ) 
    """
def nucapsIdxToInstrumentIdx(instrument, idxNucaps):
    if(instrument == 'cris-fsr'):
        band1Idx = idxNucaps[ np.where( idxNucaps <= 715 ) ] - 2
        band2Idx = idxNucaps[ np.where( (idxNucaps > 715 ) & (idxNucaps<=1584))] - 6
        band3Idx = idxNucaps[ np.where( idxNucaps > 1584 ) ] - 10
        b1, b2, b3 = band1Idx.tolist(), band2Idx.tolist(), band3Idx.tolist()
        b1.extend(b2)
        b1.extend(b3)
    elif(instrument == 'cris'):
        band1Idx = idxNucaps[ np.where( idxNucaps <= 715 ) ] - 2
        band2Idx = idxNucaps[ np.where( (idxNucaps > 715 ) & (idxNucaps<=864))] - 6
        band3Idx = idxNucaps[ np.where( idxNucaps > 864 ) ] - 10
        b1, b2, b3 = band1Idx.tolist(), band2Idx.tolist(), band3Idx.tolist()
        b1.extend(b2)
        b1.extend(b3)
    else: b1 = idxNucaps
    return np.asarray(b1)

def parseGeosActive(instrument, platform, path2tbl):
    """
    input: instrument/platform (e.g., iasi/metop-a), and path to the GEOS-DAS gmao_satinfo.db/active_channels.tbl 
    output: assimilated instrument channels.
    """
    lines = []
    with open(path2tbl) as f:
        lines = f.readlines()
    splitLines = []
    assimilatedChannels =[]
    for l in lines:
        splitLine = l.strip().split()
        if(len(splitLine)>5):
            if(instrument in splitLine[5] and platform in splitLine[0]):
                if (int(splitLine[3][0:4]) >= 2100):
                    if(int(splitLine[6])>0):
                        for i in range(int(splitLine[6])):
                            assimilatedChannels.append(int(splitLine[6+i+1]))
    return assimilatedChannels

def parseF89(fpath):
    """
    Parse data from NUCAPS F89 files. Return two dictionaries. The first is the channel/frequency information, 
    the second is flags/spectra for what looks to be a sample profile, with some peturbations for HNO3, NO2, SO2, etc...
    """
    # first bit is grabbing the cryptic metadata that tells you channels used in parts of retrieval
    # We'll call them "keys"    
    keys = ['nchsort','ivsort','freqsort',\
            'nchstrat','ivstrat','freqstrat',\
            'nchsurf','ivsurf','freqsurf',\
            'nchtemp','ivtemp','freqtemp',\
            'nchtemp2','ivtemp2','freqtemp2',\
            'nchwatr','ivwatr','freqwatr',\
            'nchozon','ivozon','freqozon',\
            'nchcld','ivcldccr','freqcldccr','cldhgtidx','ivcldhgt','freqcldhgt',\
            'nchco2','ivco2','freqco2',\
            'nchsfovco2','ivsfovco2','freqsfovco2','masksfovco2',\
            'nchch4', 'ivch4', 'freqch4',\
            'nchco', 'ivco', 'freqco',\
            'nchHNO3', 'ivHNO3','freqHNO3',\
            'nchN2O', 'ivN2O','freqN2O',\
            'nchSO2', 'ivSO2','freqSO2',\
            'nchdustscore','ivdustscore','freqdustscore']
    # go through the file and read the lines.
    with open(fpath) as f:
        lines = f.readlines()

    # mark the lines associated with one of the keys above.
    keyLines = []
    for i,l in enumerate(lines):
        for k in keys:
            if k in l and i not in keyLines:
                keyLines.append(i)
            if k == 'freqdustscore':
                if '# of temperature.1' in l:
                    keyLines.append(i)
    # go through and make chunks associated with each key.
    dataLines = {}
    for i,k in enumerate(keys):
        start = keyLines[i]
        end = keyLines[i+1]
        dataLines[k] = lines[start:end]
    # pass through again, this time putting data associated with the key...super messy.
    # don't ask me what I did here...it works.
    data = {}
    for k in list(keys):
        buf = dataLines[k]
        bufOut = []
        for l in buf:
            line = l.strip('\n').replace('=','').replace(k,'')
            bufOut.append(line)
        data[k] = []
        for l in bufOut:
            array = l.split(',')
            for item in array:
                if not item == '': 
                    if 'mask' not in k and not item.isspace() and k[0] !='n' : data[k].append(float(item))
                    elif('mask' in k): data[k].append(item)
                    elif(k[0] =='n'): data[k] = int(item)
    # next part is to get the table of stuff, which I think might be useful? Unless it's extra stuff associated with the microwave sounder, in which case...less useful.     
    channelData = data
      
    tableMarker = 'idx    freq mod CTTUSWOdMmO B    NEDT TUNING RTAERR     BT    O-C |   fix   wat   ozo   ch4    co Dhno3  Dn2o  Dso2 | P_tot P_fix P_wat P_ozo P_ch4  P_co'
    tableStarts = []
    
    for i,l in enumerate(lines):
        if (tableMarker[0:27] in l):
            tableStarts.append(i)
            # Stop looking after we hit microwave sounder (it won't find the full marker because the microwave header is slightly different).
            # we only want to read one table. Getting this far for one table was painful enough!
            if(not tableMarker in l): break
    tableBuf = []
    for idx,start in enumerate(tableStarts):
        if(idx+1 < len(tableStarts)):
            tableBuf.extend(lines[start+1:tableStarts[idx+1]-1])
        # otherwise it's the microwave sounder, which we don't want here.
        #else:
        #    tableBuf.append(lines[start+1::])
    tableData = {}            
    tableDataKeys = tableMarker.replace('|','').replace('mod','').split()
    for k in tableDataKeys:
        tableData[k] = []
    tableData['flagCloudClearing'] = []
    tableData['flagTemperaturePass1'] = []
    tableData['flagTemperaturePass2'] = []
    tableData['flagUpper'] = []
    tableData['flagH2O'] = []
    tableData['flagO3'] = []
    tableData['flagCO2'] = []
    tableData['flagCH4'] = []
    tableData['flagCO'] = []
    tableData['flagHNO3'] = []
    tableData['flagN2O'] = []
    tableData['flagSO2'] = []
    tableData['flagUsed'] = []
    for l in tableBuf:
        tableLine = l.strip().replace('|','').split()
        if( len(tableLine) == 24):
            # we actually have mod data, drop it! Not relevant to what I'm doing (I think).
            del tableLine[2]
        for i,k in enumerate(tableLine):
            if tableDataKeys[i] == 'idx':
                tableData[ tableDataKeys[i] ].append(int(k))
            elif tableDataKeys[i] == 'B':
                if(k =='.'): tableData[ tableDataKeys[i] ].append(False)
                else: tableData[ tableDataKeys[i] ].append(True)
            elif tableDataKeys[i] == 'CTTUSWOdMmO':
                if('C' in k): tableData['flagCloudClearing'].append(True)
                else: tableData['flagCloudClearing'].append(False)

                if(k[1] == 'T'): tableData['flagTemperaturePass1'].append(True)
                else: tableData['flagTemperaturePass1'].append(False)

                if(k[2] == 'T'): tableData['flagTemperaturePass2'].append(True)
                else: tableData['flagTemperaturePass2'].append(False)

                if('U' in k ): tableData['flagUpper'].append(True)
                else: tableData['flagUpper'].append(False)

                if('W' in k): tableData['flagH2O'].append(True)
                else: tableData['flagH2O'].append(False)

                if('O' in k ): tableData['flagO3'].append(True)
                else: tableData['flagO3'].append(False)

                if('d' in k ): tableData['flagCO2'].append(True)
                else: tableData['flagCO2'].append(False)
                
                if('M' in k ): tableData['flagCH4'].append(True)
                else: tableData['flagCH4'].append(False)

                if('m' in k): tableData['flagCO'].append(True)
                else: tableData['flagCO'].append(False)

                if('h' in k): tableData['flagHNO3'].append(True)
                else: tableData['flagHNO3'].append(False)

                if('n' in k): tableData['flagN2O'].append(True)
                else: tableData['flagN2O'].append(False)

                if('s' in k): tableData['flagSO2'].append(True)
                else: tableData['flagSO2'].append(False)
                
                if('N' in k): tableData['flagUsed'].append(True)
                else: tableData['flagUsed'].append(False)

                tableData[ tableDataKeys[i] ].append(k)
            else:
                if(k != '.' and k != 'BAD'):
                    tableData[ tableDataKeys[i] ].append(float(k))
                else:
                    tableData[ tableDataKeys[i] ].append(np.nan)

         
    return channelData,tableData     

def instrumentIdxToNucapsIdx(instrument, idx):
    if(instrument == 'cris-fsr'):
        band1Idx = idx[ np.where( idx<=713 )] + 2
        band2Idx = idx[ np.where( (idx>713) & (idx<=1578) )] + 6
        band3Idx = idx[ np.where( idx > 1578)] + 10
        b1, b2, b3 = band1Idx.tolist(), band2Idx.tolist(), band3Idx.tolist()
        b1.extend(b2)
        b1.extend(b3)
    elif(instrument == 'cris'):
        band1Idx = idx[ np.where( idx<=713 )] + 2
        band2Idx = idx[ np.where( (idx>713) & (idx<=864) )] + 6
        band3Idx = idx[ np.where( idx > 864)] + 10
        b1, b2, b3 = band1Idx.tolist(), band2Idx.tolist(), band3Idx.tolist()
        b1.extend(b2)
        b1.extend(b3)
    else:b1 = idx
    
    return np.asarray(b1)


if __name__ == "__main__":

    # input file contains index (start with 1) 
    fn = sys.argv[1]
    idxBufrSubset = []
    with open(fn) as f: lines = f.readlines()
    for l in lines: idxBufrSubset.append( int(l.split()[1]) )

    if('airs' in fn): 
        instrument, platform = 'airs','aqua'
        f89File = 'etc/airs_v5p9.f89'
    elif('cris' in fn and 'fsr' in fn): 
        instrument, platform = 'cris-fsr','n20'
        f89File = 'etc/cris_npp_fsr.f89'
    elif('cris' in fn ): 
        instrument, platform = 'cris','npp'
        f89File = 'etc/cris_npp_nsr.f89'
 
    elif('iasi' in fn): 
        instrument, platform = 'iasi','metop-a'
        f89File = 'etc/d07bin_oc_day.f89'
    else:
        sys.exit("I don't recognize this input file. Exiting.")

    geosAssimilated = parseGeosActive(instrument, platform, os.path.join('etc','active_channels.tbl') )
    freqMeta, nucapsTable = parseF89(f89File)
    potentiallyCachedFile = os.path.join( os.path.split(os.path.realpath(__file__))[0], 'etc',instrument+'_wavenumbers.h5')

    idxNucapsOzone = freqMeta['ivozon']
    freqNucapsOzone = freqMeta['freqozon']
    if(not os.path.exists( potentiallyCachedFile ) ):
        print("Reading from RTTOV coefficients.")
        instrumentWavenumbersFromRttov = getInstrumentWavenumbersFromRttov(instrument, np.asarray(idxBufrSubset), np.asarray(geosAssimilated), nucapsTable['freq'], np.asarray(idxNucapsOzone ) )
    else:
        print("Reading from Cached files.")
        h5 = h5py.File( potentiallyCachedFile )
        instrumentWavenumbersFromRttov = h5['wavenumbers']



    idxNucapsInInstrumentOzone = nucapsIdxToInstrumentIdx(instrument, np.asarray(idxNucapsOzone) )
    main_print(instrument, nucapsTable, idxBufrSubset, idxNucapsInInstrumentOzone, freqMeta, instrumentWavenumbersFromRttov)
    main_plot(freqMeta, nucapsTable, geosAssimilated, instrument, platform)
