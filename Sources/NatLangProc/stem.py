


'''
    stems and removes stop words
'''

import sys
import os
from stemming import porter2
import re
import string


def get_names_from_list(listfile):
    
    
    x=[]
    f=open(listfile, 'r')
    for l in f.readlines():
        s=l.split()
        x.append(s[0].lower())
    f.close()
    
    return x


def stem(file_obj):
    raw_data = file_obj.readlines()
    
    corpus = []
    syms_rm = ['.', ',', '\\', '\n', '/', '{', '}', '[', ']', '\'', '\"',
               '(', ')', '$', '^', '<', '>', '_', '-', '~', '=', '+', '*', '%', '&', '#']
    for doc in raw_data:
        re.sub(r'[^\w]', '', doc)
        for sym in syms_rm:
            doc = doc.replace(sym, '')
        doc = str.lower(doc)
        doc_stem = ' '.join([porter2.stem(word) for word in doc.split(' ')])
        corpus.append(doc_stem)
    return corpus



if __name__=='__main__':

    if len(sys.argv)<3:
        print sys.argv[0], '[file to clean] [output file]'
        exit()
    
    f=open(sys.argv[1], "r")
    c=stem(f)
    f.close()
    
    black_list_file = os.path.dirname(os.path.realpath(__file__)) + '/blacklist129.txt'        
    black_list=get_names_from_list(black_list_file)
    black_list=set(black_list)
    
    f2=open(sys.argv[2], "w")
    for d in c:
        l= d.split()
        
        for s in l:
            news=s.translate(string.maketrans("",""), string.punctuation)
            if news.lower() not in black_list:
                f2.write(news.lower()+" ")
        f2.write("\n")

        
    f2.close()
    