

def hist(l):
    
    h={}
    for x in l:
        if x not in h:
            h[x]=0
        h[x]+=1
    return h

def get_word_wn_count(corpus):

    '''
        returns {word : [new_label, occurrences]}
    '''
    
    word_wn_count={}
        
    for s in corpus:
        for w in s:
            if w not in word_wn_count:
                new_entry=[len(word_wn_count), 0]
                word_wn_count[w]=new_entry
            word_wn_count[w][1]+=1
    
    return word_wn_count


def arrange_list(b, word_wn_count):
    
    newb=[0]*len(word_wn_count)
    for k in range(len(b)):
        if k in word_wn_count:
            newb[word_wn_count[k][0]]=b[k]
    
    sum_newb=sum(newb)
    if sum_newb>0:
        newb=[v/sum_newb for v in newb]
    else:
        print 'WARNING: for some topics, no words have drawn at all. You should have used too many topics or too few documents'
    return newb


def write_corpus_format(corpus, filename, betas=[], overwrite_corpus=True):

    '''
        corpus is supposed to be a list of list of words - strings or ints
        [['word1','word2','word3'], ['word1', 'word1'], [...], ...]
        BUT int format should be used  if betas is passed
        
        betas[t][w] is p(w|t) (not the logarithm)
        which is renamed according to the labels used in the corpus file (non existing words are skipped)
        
        filename is where to write the corpus in the LDA-corpus format
        
        overwrite_corpus simply renames the word labels so that also in the syn_corpus.txt
        labels match
        
    '''

    word_wn_count=get_word_wn_count(corpus)
    
    f=open(filename, 'w')
        
    for d in corpus:
        
        d_renamed=[]
        for w in d:
            d_renamed.append(word_wn_count[w][0])
        
        h=hist(d_renamed)
        
        f.write(str(len(h))+' ')
        for k in sorted(h.keys()):
            f.write('%d:%d '%(k, h[k]))
        f.write('\n')
    
    f.close()
    
    if overwrite_corpus:
        for i in range(len(corpus)):
            for j in range(len(corpus[i])):
                corpus[i][j]=word_wn_count[corpus[i][j]][0]
    
    if len(betas)>0:
        newbetas=[]
        for b in betas:
            newbetas.append(arrange_list(b, word_wn_count))
        
        return newbetas


def get_pw_and_pd_from_corpus_file(corpus_text_file):


    '''
        takes a corpus file in format text, and return pw and pd
        pw is using the usual labelling convention of naming the words
        in order of appearance
    '''
    
    corpus=[v.split() for v in open(corpus_text_file, 'r').readlines()]
    word_wn_count=get_word_wn_count(corpus)
    
    pw=[0.]*len(word_wn_count)
    for k,v in word_wn_count.iteritems():
        pw[v[0]]=float(v[1])

    pd=[]
    for d in corpus:
        if len(d)>0:
            pd.append(float(len(d)))
    
    pw=[v/sum(pw) for v in pw]
    pd=[v/sum(pd) for v in pd]

    return pw, pd
    



if __name__=='__main__':
    
    print get_pw_and_pd_from_corpus_file('syn_corpus.txt')

    
    
  