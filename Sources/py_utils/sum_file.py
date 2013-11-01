


import sys


if len(sys.argv)<2:
    print sys.argv[0], '[file to sum] [file to sum2]'
    exit()


sum_v=[]
for l in open(sys.argv[1]):
    s=l.split()
    if len(s)>0:
        sum_v.append(float(s[0]))


sum_v2=[]

for l in open(sys.argv[2]):
    s=l.split()
    if len(s)>0:
        sum_v2.append(float(s[0]))

for i in range(len(sum_v)):
    print sum_v[i], sum_v2[i]
    
print '#', sum(sum_v), sum(sum_v2)
