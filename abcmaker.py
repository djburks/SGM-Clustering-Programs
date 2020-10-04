import sys

segments = []
switch = False
thresh = 1
outfile = sys.argv[1].split('.network')[0] + '.abc'
o = open(outfile,'w')
numz = []

with open(sys.argv[1]) as infile:
    for lines in infile:
        lines = lines.rstrip()
        if (lines == '###...###'):
            switch = True
            continue
        if (switch == True):
            numz.append(float(lines))
            
            
        
maxn = max(numz)
minn = min(numz)

switch = False


with open(sys.argv[1]) as infile:
    for lines in infile:
        lines = lines.rstrip()
        if lines.startswith('>'):
            segments.append((lines[1:]).split('/')[-1])
        elif (lines == '###...###'):
            switch = True
            x = 0
            y = x + 1
            continue
        if (switch == True):
            div = float(lines)
            if (div < thresh):
                #div = 1
                div = 1 - ((div - minn) / (maxn - minn))
                #div = thresh - div
                o.write(segments[x] + '\t' + segments[y] + '\t' + str(div) + '\n')
                y += 1
            else:
                y += 1

            if (y == len(segments)):
                x += 1
                y = x + 1
        
