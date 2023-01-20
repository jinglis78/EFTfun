import os
coeffs='cW','cll1'
path='/tmp/test_example12/c_scan_{coeff}.txt'

def parse(fn):
    if not os.path.exists(fn):
        print 'WARNING:',fn,'not found'
        return None
    f=open(fn)
    txt=f.read().strip()    
    return [[float(x.strip()) for x in t.split(',')] for t in txt.split()]

for c in coeffs:
    fn=path.format(coeff=c)
    if fn is None:
        continue
    central,(low68,high68),(low95,high95)=parse(fn)
    print c,'central: {}, 68% CL: [{}:{}], 95% CL: [{}:{}]'.format(central,low68,high68,low95,high95)
