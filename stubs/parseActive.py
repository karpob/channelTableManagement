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
                    print("read count",int(splitLine[6]))
                    if(int(splitLine[6])>0):
                        for i in range(int(splitLine[6])):
                            assimilatedChannels.append(int(splitLine[6+i+1]))
    return assimilatedChannels
if __name__ == "__main__":
    assimilatedChannels = parseGeosActive('iasi','metop-a','active_channels.tbl')
    print(assimilatedChannels)
    
    
