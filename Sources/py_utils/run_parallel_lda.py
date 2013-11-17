
import sys
import os
import inspect
import time
import glob
import math
import time
from time import gmtime, strftime


def update_activity(activity, first=False):
    
    if first:
        ftim_act=open('time_activity.log', 'w')
    else:
        ftim_act=open('time_activity.log', 'app')
    ftim_act.write(activity+' '+strftime("%Y-%m-%d %H:%M:%S", gmtime())+'\n')
    ftim_act.close()



def slice_docs(no_docs, no_jobs):
    '''
        returns a list with the No. of docs in each folder
    '''
    l=[0]*no_jobs
    for i in xrange(no_docs):
        l[i%no_jobs]+=1
    

    prev=0
    cuml=[]
    for i, v in enumerate(l):
        cuml.append(v+prev)
        prev+=v

    assert cuml[-1]==no_docs
    return cuml
    

def split_docs_in_folders(corpus_file, folders, no_docs, no_jobs):
    '''
        writes smaller corpora
        in folders
    '''
    ranges=slice_docs(no_docs, no_jobs)
    print 'docs in folders:', ranges
    doc_counter=0
    folder_counter=0
    for l in open(corpus_file):
        if doc_counter==0:
            outfile=open(folders[folder_counter]+'/sliced_corpus.txt', 'w')
    
        if doc_counter>ranges[folder_counter]:
            print 'done with copying all docs in ', folders[folder_counter]
            print 'docs copied::', doc_counter
            outfile.close()
            folder_counter+=1
            assert folder_counter<len(folders)
            print 'moving to folder', folders[folder_counter]
            outfile=open(folders[folder_counter]+'/sliced_corpus.txt', 'w')
        
        outfile.write(l)
        doc_counter+=1
    
    print 'done with copying docs. Docs copied::', doc_counter


def get_likelihoods(folders, filename='lda_log_likelihood.txt'):
    
    tot_lik=0.
    for f in folders:
        #print ('cat '+f+'/'+filename)
        #os.system('cat '+f+'/'+filename)
        #print ('cat '+f+'/'+filename)
        #os.system('cat '+f+'/'+filename)
        tot_lik+=float(open(f+'/'+filename).readlines()[-1].split()[-1])
    return tot_lik
    


def all_jobs_are_done(files, waiting_time, filename='log.log'):
    
    '''
        waits until 
        '--done!--'
        appears in all 'files/filename'
        waiting_time is simply 
        how often to check
    '''
    done=False
    while done==False:
        done=True
        
        for f in files:
        
            last_word='not done!'
            try:
                last_word=open(f+'/'+filename).readlines()[-1].split()[0]
            except:
                last_word='not done!'
            # we are done yet *somewhere*
            if last_word!='--done!--':
                done=False
                print 'not done yet in', f

        
        if done==False:
            update_activity('checking again in '+str(waiting_time)+' secs')
            time.sleep(waiting_time)
    

def get_lines(infile):
    os.system('wc '+infile+' > word_count.tmp')
    no_lines = int(open('word_count.tmp').readlines()[-1].split()[0])
    return no_lines



if __name__=='__main__':
    
    
    
    if len(sys.argv)<4:
        print 'python', sys.argv[0], '[corpus_file] [initial_model] [no_jobs] [word_wn_count]'
        exit()
    
    
    corpus_file= sys.argv[1]
    print 'corpus_file', corpus_file
    model_file= sys.argv[2]
    print 'model_file', model_file
    no_jobs= int(sys.argv[3])
    print 'no_jobs', no_jobs
    word_wn_file= sys.argv[4]
    
    update_activity('starting', True)
    #========= default parameters ===============
    initial_alpha=0.01
    option_t='-t 100'
    #secs between checks until all_jobs_are_done
    waiting_time=3
    # waits these seconds for files to be closed
    writing_time=2
    #========= default parameters ===============

    # deleting before starting
    os.system('rm -r parallel_lda_*')

    # creating folders
    for i in range(no_jobs):
        print 'job::', i
        new_folder='parallel_lda_'+str(i)
        os.system('mkdir '+new_folder)
    
    # how many docs?
    no_docs=get_lines(corpus_file)
    print 'No. docs', no_docs
    # how many topics?
    no_topics=get_lines(model_file)
    # writing alpha
    alpha_file=open('alphas.txt', 'w')
    for i in xrange(no_topics):
        alpha_file.write('0.01 ')
    alpha_file.write('\n')
    alpha_file.close()
    # getting folders
    folders=glob.glob('parallel_lda_*')
    print 'folders::', folders
    # slicing corpus to folders
    split_docs_in_folders(corpus_file, folders, no_docs, no_jobs)
    
    # No. iterations
    iter=0
    MAX_iter=100
    old_lik=-1e100
    
    likout=open('par_likelihood.txt', 'w')
    while iter < MAX_iter:
        
        for f in folders:
            # move to folder
            # E step
            os.chdir(f)
            command_line='nohup ../bin/topicmap -f sliced_corpus.txt -infer -model ../'+\
                          model_file+' -alpha_file ../alphas.txt -word_wn ../'+\
                          word_wn_file+' '+option_t+' > log.log &'
            #print 'running::', command_line
            os.system(command_line)
            os.chdir('../')

        
        # this lets you wait until all jobs are done
        all_jobs_are_done(folders, waiting_time)
        time.sleep(writing_time)
        
        # updating models and alphas
        lik=get_likelihoods(folders)
        
        # collecting
        os.system('cat parallel_lda_*/lda_class_words.txt > model.txt')
        os.system('cat parallel_lda_*/lda_gammas.txt > all_gammas.txt')
        time.sleep(writing_time)
        
        # optimizing alpha
        os.system('./bin/opt_alpha all_gammas.txt')
        time.sleep(writing_time)
    
        # updates file name and counter
        model_file='model.txt'
        iter+=1
        
        # are we done?
        print 'likelihood::', lik, 'prev::', old_lik, 'iter::', iter
        update_activity('likelihood:: '+str(lik)+' prev:: '+\
                        str(old_lik)+' iter:: '+str(iter))
        likout.write(str(iter)+' '+str(lik)+'\n')
        # at least 9 iterations
        if iter>9:
            if math.fabs((lik-old_lik)/old_lik)<1e-5 or lik<old_lik:
                break
        old_lik=lik
        
    print 'done'
    likout.close()
    
    
    
    