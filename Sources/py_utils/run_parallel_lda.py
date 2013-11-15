
import sys
import os
import inspect
import time
import glob


def slice_docs(no_docs, no_jobs):
    
    div


def all_jobs_are_done(files, filename='log.log', waiting_time=10):
    
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
            print 'checking again in ', waiting_time, ' secs'
            time.sleep(waiting_time)
    
        



if __name__=='__main__':
    
    if len(sys.argv)<4:
        print 'python', sys.argv[0], '[corpus_file] [initial_model] [no_jobs]'
        exit()
    
    
    corpus_file=sys.argv[1]
    model_file=sys.argv[2]
    no_jobs= int(sys.argv[3])

    for i in range(no_jobs):
        print 'job::', i
        new_folder='parallel_lda_'+str(i)
        os.system('mkdir '+new_folder)
    
    # slicing corpus to folders
    # how many docs?
    os.system('wc '+sys.argv[1]+' > word_count.tmp')
    no_docs = int(open('word_count.tmp').readlines()[-1].split()[0])
    print 'no_docs', no_docs
    
    dones=glob.glob('parallel_lda_*')
    
    # splits docs in folder
    
    while condition:
        
        
        for s in folders:
            # move to folder
            # E step
            #os.system('./bin/topicmap -f 'datasets'  -infer -model 'initialmodel' -alpha_file 'initialalpha' -word_wn word_wn_count.txt &')
            pass
        # this lets you wait until all jobs are done
        all_jobs_are_done(dones)
        # cat folders*/lda_class_words.txt > model.txt
        # cat folders*/lda_gammas.txt > all_gammas.txt
        # sum folders*/lda_log_likelihood.txt
        os.system('./bin/opt_alpha gammafile')
        
    
    
    print 'done!!'
    
    
    
    