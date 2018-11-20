import struct
f = open('cris-fsr431_n20.SpcCoeff.bin','rb')
buf1 = f.read(44)
buf2 = f.read(37)
bb = f.read(3)
fsrArray = []
for i in range(431):
    buf = f.read(4)
    val, = struct.unpack('i',buf)
    print(i+1, val)
    fsrArray.append(val)
print(fsrArray)
print(len(fsrArray))
"""
buf3 = f.read(4)
print(buf3)
val1 = struct.unpack('i',buf3)
bb = f.read(4)
print(bb)
c = struct.unpack('i',bb)
print(val1,c)
"""
#print(val1,val2)
#for n,aa in enumerate(a):
#    print(n,aa)

