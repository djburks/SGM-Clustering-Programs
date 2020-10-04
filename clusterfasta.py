# Takes in MCL clusters and outputs fasta files for each segment post clustering.

import sys
import os

inabc = sys.argv[1]
clus = str(sys.argv[2])
clusext = str(sys.argv[2]).replace('.','')
if float(clus) >= 10:
    clusext = clusext + '0'

insfna = inabc.split('.abc')[0] + '.network'
root = inabc.split('.abc')[0]

com1 = 'mcxload -abc ' + inabc + ' --stream-mirror -write-tab ' + root + '.tab -o ' + root + '.mci'
com2 = 'mcl ' + root + '.mci -I ' + clus
com3 = 'mcxdump -icl out.' + root + '.mci.I' + clusext + ' -tabr ' + root + '.tab -o dump.' + root + '.mci.I' + clusext

os.system(com1)
os.system(com2)
os.system(com3)



segments = {}
clusters = {}

with open(insfna) as infile:
    for lines in infile:
        lines = lines.rstrip()
        if lines == '###...###':
            break
        if lines.startswith('>'):
            name = lines[1:].split('/')[-1]
            segments[name] = ''
            curseg = name
        else:
            segments[curseg] += lines
            
c = 0
with open('dump.' + root + '.mci.I' + clusext) as infile:
    for lines in infile:
        clusters[c] = []
        lines = lines.rstrip()
        values = lines.split('\t')
        for v in values:
            clusters[c].append(v)
        c += 1
outfasta = inabc.split('.abc')[0] + '.sfnab'
out = open(outfasta,'w')
out.close()

for c in clusters:
    seglist = sorted(clusters[c])
    outfasta = inabc.split('.abc')[0] + '.sfnab'
    out = open(outfasta,'a')
    out.write('>' + sys.argv[1].split('.')[0] + '.fna.' + str(c) + '\n')
    fastastring = ''
    for s in seglist:
        fastastring += segments[s]
    out.write(fastastring + '\n')

os.system('rm *' + root + '*mci*')
os.system('rm *' + root + '*tab')
