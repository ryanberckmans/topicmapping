

label_file='parallel_0_0_3/word_wn_count.txt'
sig_word_file='sig_words.edges'

wn_word={}


for l in open(label_file):
    s=l.split()
    wn_word[int(s[1])]=s[0]
    
#print wn_word


wn=-1
neighs=[]
for l in open(sig_word_file):
    s=l.split()
    wnew=int(s[0])
    
    # if new wn
    if wnew!=wn:
        if wn>=0:
            sorted_neighs=sorted(neighs, reverse=True)
            for i in range(min(len(sorted_neighs), 10)):
                print wn_word[wn], wn_word[sorted_neighs[i][1]], sorted_neighs[i]
            print '----------', len(sorted_neighs), wn
        neighs=[]
        wn=wnew
    neighs.append((float(s[2]), int(s[1])))

