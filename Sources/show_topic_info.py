
import sys

def parse_file(filename):
    
    
    topic_dict={}
    subtopic_dict={}
    
    next_is_a_topic=False
    next_is_a_subtopic=False
    num_docs=0
    topic_id=0
    subtopic_id=-1
    
    for l in open(filename):
        s=l.split()
        if len(s)>0:
            if s[0]=='#docs':
                num_docs=int(s[1])
            if next_is_a_topic:
                # setting topic name
                topic_dict[topic_id][1]=l
                next_is_a_topic=False
            if next_is_a_subtopic:
                # setting subtopic name
                subtopic_dict[subtopic_id][1]=l
                topic_dict[topic_id][2].append(subtopic_id)
                next_is_a_subtopic=False

            if s[0]=='topic:' or s[0]=='****topic:':
                next_is_a_topic=True
                pt=float(s[-1])
                topic_id=int(s[1])
                # setting new topic
                topic_dict[topic_id]=[pt, "", []]
            if s[0]=='->':
                next_is_a_subtopic=True
                pt=float(s[-1])
                subtopic_id=int(s[2])
                # setting new subtopic
                subtopic_dict[subtopic_id]=[pt, ""]
            if s[0]=='no' and s[1]=='subtopics':
                # if there are no subtopics 
                pt=topic_dict[topic_id][0]
                subtopic_id+=1
                subtopic_dict[subtopic_id]=[pt, topic_dict[topic_id][1]]
                
    return num_docs, topic_dict, subtopic_dict



if __name__=='__main__':
    
    if len(sys.argv)<6:
        print sys.argv[0]+' [(sub)topic_summary.txt] [min_docs] [num_words_per_topic] [corpus file] [doc_topic_file]'
        exit()
    
    # topic_dict[topic_id] = [pt, name, subtopic_list]
    # subtopic_dict[topic_id] = [pt, name]
    min_docs=float(sys.argv[2])
    num_words_per_topic=int(sys.argv[3])
    
    # parsing 
    num_docs, topic_dict, subtopic_dict= parse_file(sys.argv[1])
    print 'num_docs', num_docs
    #print 'Number of subtopics: ', len(subtopic_dict)
    topic_flag=len(subtopic_dict)==0



    if topic_flag:
        f=open('short_summary.txt', 'w')
    else:
        f=open('short_summary_sub.txt', 'w')
        
    
    # writing summary
    for topic_id, topic_stuff in topic_dict.iteritems():
        if topic_stuff[0]*num_docs> min_docs:
            f.write('*** topic: '+str(topic_id)+' pt: ' +str(topic_stuff[0])+'\n')
            f.write(' '.join(topic_stuff[1].split()[:num_words_per_topic])+'\n')
            for subtopic_id in topic_stuff[2]:
                subtopic_stuff=subtopic_dict[subtopic_id]
                if subtopic_stuff[0] * num_docs > min_docs:
                    f.write('   -> subtopic: '+str(subtopic_id)+' pt: ' +str(subtopic_stuff[0])+'\n')
                    f.write('     '+' '.join(subtopic_stuff[1].split()[:num_words_per_topic])+'\n')
    f.close()


    # writing back
    if topic_flag:
        f=open('doc_topics_text.txt', 'w')
    else:
        f=open('doc_subtopics_text.txt', 'w')
        
    corpus_file=sys.argv[4]
    thetas=open(sys.argv[5]).readlines()
    
    print "corpus file: ", corpus_file
    print "topic file: ", sys.argv[5]
    counter=0
    
    for l in open(corpus_file):
        if len(l.split())>0:
            #print counter, len(thetas)
            if counter>=len(thetas):
                print 'error: topic file lenght does not match corpus file. Are you sure ', sys.argv[5], 'is the right file?'
                exit()
            
            theta=thetas[counter]
            f.write(l)
            f.write(theta)
            topic_prs=[v.split(':') for v in theta.split()]
            for t in topic_prs:
                if topic_flag:
                    topic_name=topic_dict[int(t[0])][1]
                    topic_name=' '.join(topic_name.split()[:num_words_per_topic])
                else:
                    #print t
                    topic_name=subtopic_dict[int(t[0])][1]
                    topic_name=' '.join(topic_name.split()[:num_words_per_topic])

                f.write('\''+topic_name+'\':'+t[1]+'\t')
            f.write('\n\n')
            counter+=1
        else:
            print 'empty!'
    
    f.close()
    assert counter==len(thetas), 'topics and corpus have different lengths!'
    
    
    
    