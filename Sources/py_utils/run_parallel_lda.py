
import sys
import os
import inspect
import time
import glob


def all_jobs_are_done(files, filename='done.txt', waiting_time=10):
    
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
            if last_word!='done!!':
                done=False
                print 'not done yet in', f

        
        if done==False:
            print 'checking again in ', waiting_time, ' secs'
            time.sleep(waiting_time)
    
        



if __name__=='__main__':
    
    if len(sys.argv)<4:
        print 'python', sys.argv[0], '[corpus_file] [beta_file] [alpha_file]'
        print '# jobs is set to:', jobs_square_root
        print 'open this file to change the #jobs'
        exit()
        
    # no. jobs needed is jobs_square_root**2
    jobs_square_root = 3

    for i in range(jobs_square_root):
        print 'job::', i
        new_folder='parallel_lda_'+str(i)
        os.system('mkdir '+new_folder)
    
    
    dones=glob.glob('parallel_lda_*')
    
    # this lests you wait until all jobs are done
    all_jobs_are_done(dones)
    
    print 'done!!'
    
    
    
    