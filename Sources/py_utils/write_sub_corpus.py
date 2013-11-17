


'''
    takes a file in format
    'wn word_str occ topic'
    and a topic number
    and select all words coming 
    from that topic 
    to write a new corpus
    also, the list of selected documents
    is given
'''

import sys


if __name__=='__main__':
    
    
    if len(sys.argv)<3:
        print 'python', sys.argv[0], '[word_assignments] [topic_number]'
        exit()
    
    infile= sys.argv[1]
    print 'word_assignments', infile
    topic_number= int(sys.argv[2])
    print 'topic_number', topic_number
    
    sub_corpus=open('sub_corpus.txt', 'w')
    doc_list=open('doc_list.txt', 'w')
    
    doc_num=0
    for l in open(infile):
        s=l.split()
        word=''
        occ=0
        # (word, occ)
        words=[]
        for i in range(len(s)):
            if i%4==1:
                word=s[i]
            if i%4==2:
                occ=int(s[i])
            if i%4==3:
                topic=int(s[i])
                if topic==topic_number:
                    words.append((word,occ))
        
        if len(words)>1:
            # writing doc
            doc_list.write(str(doc_num)+'\n')
            for w,occ in words:
                for o in range(occ):
                    sub_corpus.write(w+' ')
            sub_corpus.write('\n')
        else:
            print 'doc skipped', doc_num
        
        
        doc_num+=1
    
    doc_list.close()
    sub_corpus.close()
    
