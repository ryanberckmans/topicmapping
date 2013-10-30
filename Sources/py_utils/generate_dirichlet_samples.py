import math
import random
import bisect
import numpy as np


def get_thetas(params):
    
    
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
        print 'error'
        exit()
    
    return sample

params=[0.003, 0.003, 0.003, 0.003]
params=[0.004]*20+[0.1]


f=open('sample.txt', 'w')
for i in range(20000):
    s=get_thetas(params)
    for k in range(len(s)):
        f.write(str(s[k]+params[k])+' ')
    f.write('\n')
f.close()
        