

'''
    this script is meant to run topicmapping
    for building the network of significant words
    for very large datasets
    '''

import sys
import os
import inspect
import time


# no. jobs needed is jobs_square_root**2
jobs_square_root = 3
pvalue = 0.05     

if __name__=='__main__':
    
    if len(sys.argv)<2:
        print 'python', sys.argv[0], '[corpus_file]'
        print '# jobs is set to:', jobs_square_root**2
        print 'open this file to change the #jobs or pvalue'
        exit()
    
    file_dir= os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    topicmap_path=file_dir+'/../../bin/topicmap'
    
    working_dir=os.getcwd()
    print 'working directory:', working_dir
    
    for i in range(jobs_square_root):
        for j in range(jobs_square_root):
            print 'job::', i, j
            new_folder='parallel_'+str(i)+'_'+str(j)+'_'+str(jobs_square_root)
            os.system('mkdir '+new_folder)
            os.chdir(new_folder)
            os.system(topicmap_path+' -p '+str(pvalue)+' -f ../'+sys.argv[1]+' -parall '+str(i)+\
                      ':'+str(j)+':'+str(jobs_square_root)+' > parall.log &\n')
            
            # waiting a few minutes before running another job
            time.sleep(30*60)
            os.chdir(working_dir)
    
    print 'when jobs will be done, you can collect the networks with:'
    print 'cat parallel_*/sig_words*.edges > sig_words.edges'
    print 'Write the network in pajek format with:'
    print 'mkdir sig_words_results/'
    print 'Run Infomap with:'
    print './bin/Infomap pajek.net sig_words_results/ --two-level --undirected --num-trials 10 --seed 10 > infomap.log &'
