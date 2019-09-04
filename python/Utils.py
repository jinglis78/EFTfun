import math

#TODO: logger

def isFile(data):
    if data.split(':')[0].lower().endswith('.root') or data.split(':')[0].lower().endswith('.yoda'):
        return True
    return False

def saveToTree(results,outfile):
    from array import array
    from ROOT import TTree
    outfile.cd()
    t=TTree('tree','tree')
    arrays={}
    for r in results:
        arrays[r]=array( 'f', [ 0. ] )
        t.Branch( r, arrays[r], r+'/F' )
    for i in range(len(results[r])):
        for r in results:
            arrays[r][0]=results[r][i]
        t.Fill()
    t.Write()
    


def round2(number,digits=3):
    if number==0:
        return 0
    if math.log10(abs(number))>digits-1:
        return int(round(number))   
    return round(number, digits-int(math.floor(math.log10(abs(number))))-1)

def getOptionsAndData(data):
    if '{' in data:
        options=set(o.lower() for o in data[data.find('{')+1:data.find('}')].split())
    else:
        options=set()
    data=data[data.find('}')+1:].strip()
    return options,data

