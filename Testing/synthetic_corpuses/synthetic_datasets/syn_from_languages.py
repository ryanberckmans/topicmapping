#! /usr/bin/env python


import os
import math
import sys
import math
import random
import json

from syn_utilities import cumulative_from_list
from syn_utilities import print_matrix
from syn_utilities import get_thetas
from syn_utilities import draw
from syn_utilities import get_log_matrix
from write_corpus_format import write_corpus_format



current_folder = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, current_folder+'/../')
language_file=current_folder+'/languages.json'



def write_syn_corpus_languages(number_of_documents, number_of_words_per_doc, rule):
    
    corpus_file='syn_corpus.corpus'
    beta_file='syn_betas.txt'
    theta_file='syn_thetas.txt'
    
    languages=json.load(open(language_file))    
    for l in languages:
        number_of_words_per_language=len(languages[l])
    
    number_of_topics=len(languages)
    top_topics = number_of_topics/5
    other_topics = number_of_topics - top_topics
    pt=[(1.-rule)/other_topics]*other_topics + [rule/top_topics]*top_topics
    print 'rule', rule, 'number_of_words_per_language', number_of_words_per_language
    pt=[float(v)/sum(pt) for v in pt]
    print '*** planted p(t):', pt
    pt_cum=cumulative_from_list(pt)
    
    # setting betas from languages
    betas=[]
    for l in languages:
        language_index=len(betas)
        b=  [0]*number_of_words_per_language * language_index \
            + languages[l] \
            + [0]* number_of_words_per_language * (number_of_topics - language_index-1)
        sumb=sum(b)
        print 'words:', len(b), 'norm:', sumb
        if sumb>0:
            b=[float(v)/sumb for v in b]
        betas.append(b)
    
    # beta cumulative
    pw_givent_cum=[]
    for v in betas:
        pw_givent_cum.append(cumulative_from_list(v))
    
    thetas=[]
    corpus=[]
    for doc in range(number_of_documents):
        # decide the language
        t=draw(pt_cum)
        theta=[0]*number_of_topics
        theta[t]=1
        thetas.append(theta)
        doc=[]
        # decide the words given the language
        for w in range(number_of_words_per_doc):
            doc.append(draw(pw_givent_cum[t]))
        corpus.append(doc)
    
    new_betas=write_corpus_format(corpus, corpus_file, betas=betas)
    print_matrix(corpus, corpus_file.split('.corpus')[0]+'.txt')
    print_matrix(thetas, theta_file)
    print_matrix(get_log_matrix(new_betas), beta_file)



if __name__=='__main__':
    
    
    if len(sys.argv)<2:
        print sys.argv[0], ' [number_of_documents] [rule]'
        exit()
    
    number_of_documents=int(sys.argv[1])
    # rule=x means that the top 20% topics gets x probability
    # 0.2 => equiprobable topics
    rule=float(sys.argv[2])
    number_of_words_per_doc=100
    write_syn_corpus_languages(number_of_documents, number_of_words_per_doc, rule)








