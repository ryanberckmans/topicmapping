


import math
import random
import bisect
import numpy as np


def get_thetas(params, pt_cum):
    
    
    '''
        given params as alpha hyperparameters, 
        the function returns a samples from 
        a dirichlet distribution
        
    '''
    
    sample = [random.gammavariate(a,1) for a in params]
    sum_sample=sum(sample)
    
    if sum_sample>0:
        sample = [v/sum_sample for v in sample]
    else:
        select_topic=draw(pt_cum)
        sample[select_topic]=1
        sum_sample=sum(sample)
        sample = [v/sum_sample for v in sample]
    
    return sample



def draw(l):
    if abs(l[-1]-1)>1e-6:
        print l
    assert abs(l[-1]-1)<1e-6, '\Some topics have no words! Please reduce the number of topics or increase the number of words'
    return bisect.bisect(l, random.uniform(0,1))



def cumulative_from_list(thetas):
    
    thetas_cum=[]
    sum=0
    for v in thetas:
        sum+=v
        thetas_cum.append(sum)
    
    return thetas_cum


def print_matrix(l, outfile):
    f=open(outfile, 'w')
    for v in l:
        for k in range(len(v)):
            f.write(str(v[k])+' ')
        f.write('\n')
    f.close()


def get_log_matrix(p_word_given_topic):

    log_betas=[]
    for v in p_word_given_topic:
        log_p=[]
        for v2 in v:
            if v2>0:
                log_p.append(math.log(v2))
            else:
                log_p.append(-1000)
        log_betas.append(log_p)
    
    return log_betas

