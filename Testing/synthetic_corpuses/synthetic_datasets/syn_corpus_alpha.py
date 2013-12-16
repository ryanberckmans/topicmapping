#! /usr/bin/env python


import os
import math
import sys
import math
import random

from syn_utilities import cumulative_from_list
from syn_utilities import print_matrix
from syn_utilities import get_thetas
from syn_utilities import draw
from syn_utilities import get_log_matrix

current_folder = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, current_folder+'/../')

from write_corpus_format import write_corpus_format


def generate_betas_and_gammas(alpha_words, pt, pw, alpha_doc, N, w_per_doc):

    """
        alpha_words[w] are the concentration parameters for word w
        pt[t] is the topic size of topic t
        pw[w] is the probability pf word w
        alpha_doc is the concentration parameter for the docs
        N is the number of docs,
    """
    
    # normalizing pt and pw
    pt=[float(v)/sum(pt) for v in pt]
    
    print '*** planted p(t):', pt
    pt_cum=cumulative_from_list(pt)
    pw=[float(v)/sum(pw) for v in pw]
    
    
    # getting pt_given_w from alpha_words
    pt_given_w=[]
    for alpha in alpha_words:
        params=[alpha * p * len(pt) for p in pt]
        pt_gw=get_thetas(params, pt_cum)
        pt_given_w.append(pt_gw)
    
    
    # getting pw_given_t inverting pt_given_w    
    betas=[]
    for t in range(len(pt)):
        b=[]
        for w in range(len(pt_given_w)):
            b.append(pt_given_w[w][t]*pw[w])
        sumb=sum(b)
        if sumb>0:
            b=[float(v)/sumb for v in b]
        betas.append(b)
    

    
    pw_givent_cum=[]
    for v in betas:
        pw_givent_cum.append(cumulative_from_list(v))
    
    
    doc_parameters=[alpha_doc * p * len(pt) for p in pt]
    print 'alphas:', doc_parameters
    thetas=[]
    corpus=[]
    
    for doc_index in range(N):
        
        thetas_d=get_thetas(doc_parameters, pt_cum)

        thetas.append(thetas_d)
        theta_cum=cumulative_from_list(thetas_d)
        doc=[]
        for wj in range(w_per_doc):
            t=draw(theta_cum)
            doc.append(draw(pw_givent_cum[t]))
        corpus.append(doc)

    return betas,thetas,corpus


def write_syn_corpus(alpha_doc, fraction_of_generic_words, \
                           pt, pw, \
                           number_of_unique_words, number_of_topics, \
                           number_of_documents, number_of_words_per_doc, \
                           equally_probable_topics=False, equally_probable_words=False, \
                           alpha_for_generic_words=1., \
                           corpus_file='syn_corpus.corpus', \
                           beta_file='syn_betas.txt', \
                           theta_file='syn_thetas.txt'):
    '''
    
        INPUT:
        * alpha_doc is a float number which is used to compute the memberships of the documents as 
            well as non-generic words
        * fraction_of_generic_words
        * pt[t] is the probability of using topic t 
        * pw[w] is the probability of using word w
        * equally_probable_words & equally_probable_topics are flags to set pw and pt to equal probabilities
        
        The function writes down three files in the standard LDA format.
        
    '''
    
    
    if equally_probable_words:
        pw=number_of_unique_words*[1.]
    if equally_probable_topics:
        pt=number_of_topics*[1.]
    
    assert len(pt)==number_of_topics, 'pt doesnt match the number of topics'
    assert len(pw)==number_of_unique_words, 'pw doesnt match the number of unique words'
            
    generic_words=int(fraction_of_generic_words*number_of_unique_words)
    alpha_words=[alpha_for_generic_words]*generic_words+[alpha_doc]*(number_of_unique_words-generic_words)
    
    print '*** Number of documents ', number_of_documents
    print '*** Number of topics ', number_of_topics
    print '*** Number of unique words ', number_of_unique_words
    print '*** Number of generic words ', generic_words

    
    p_word_given_topic, p_topic_given_document, corpus = generate_betas_and_gammas(alpha_words, \
                                                                        pt, pw, \
                                                                        alpha_doc, \
                                                                        number_of_documents, \
                                                                        number_of_words_per_doc)
    
    print 'results can be found in:', corpus_file, corpus_file.split('.corpus')[0]+'.txt',theta_file, beta_file 
    
    new_betas=write_corpus_format(corpus, corpus_file, betas=p_word_given_topic)
    print_matrix(corpus, corpus_file.split('.corpus')[0]+'.txt')
    print_matrix(p_topic_given_document, theta_file)
    print_matrix(get_log_matrix(new_betas), beta_file)

    
if __name__=='__main__':
    
    '''
        Most parameters are set in this file (see below).
        You are just required to input alpha and the fraction of generic words.
        See ReadMe.txt for details
    '''
    
    
    if len(sys.argv)<4:
        print sys.argv[0], '[alpha] [fraction of generic words] [rule]'
        exit()
    
    
    
    alpha_doc=float(sys.argv[1])
    fraction_of_generic_words=float(sys.argv[2])
    # rule=x means that the top 20% topics gets x probability
    rule=float(sys.argv[3])

    
    ############# parameters ##################
    number_of_documents=1000
    number_of_words_per_doc=50
    number_of_topics=20
    # number of words in the vocabulary
    number_of_unique_words=2000
    ############# parameters ##################

    # top topics are the biggest topic
    # dividing by 5, means we consider 20% top topics
    top_topics = number_of_topics/5
    other_topics = number_of_topics - top_topics
    # p(t) i.e. topic probability
    pt=[(1.-rule)/other_topics]*other_topics + [rule/top_topics]*top_topics
    print rule
    
    
    assert alpha_doc>0, 'alpha need to be positive'
    assert 0<=fraction_of_generic_words<=1, 'fraction_of_generic_words should be between 0 and 1'
    assert 0<=rule<=1, 'rule should be between 0 and 1'
    
    
    write_syn_corpus(alpha_doc= alpha_doc, \
                     fraction_of_generic_words= fraction_of_generic_words, \
                     pt = pt, \
                     pw=[], \
                     number_of_unique_words=number_of_unique_words,\
                     number_of_topics= len(pt), \
                     number_of_documents=number_of_documents, \
                     number_of_words_per_doc=number_of_words_per_doc, \
                     equally_probable_words=True)
    