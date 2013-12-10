
'''
    this script is to start LDA with random initial conditions
    writes a random model with num_topics and num_terms
'''
import random
import sys

if __name__=='__main__':
    
    if len(sys.argv)<2:
        print 'python', sys.argv[0], '[num_topics] [num_terms]'
        exit()
        
    num_topics=int(sys.argv[1])
    num_terms=int(sys.argv[2])
    
    f=open('random_model.txt', 'w')
    for i in range(num_topics):
        f.write(str(i)+' ')
        for wn in range(num_terms):
            f.write(str(wn)+' '+str(1./num_terms+random.uniform(0,1))+' ')
        f.write('\n')
    f.close()