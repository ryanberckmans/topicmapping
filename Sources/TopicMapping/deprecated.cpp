

void set_corpuses_from_topics(deque<mapii> & doc_assignments, \
                              map<int, word_corpus> & old_topic_corpus, \
                              map<int, mapii> & topic_old_docs, word_corpus & oC) {
    
    
    
    // doc_topic_newid[doc][topic] is the new id of the doc in that corpus
    deque<mapii> doc_topic_newid;
    
    RANGE_loop(i, doc_assignments) {
        
        mapii topic_newid;
        IT_loop(mapii, itm, doc_assignments[i]) {
            
            // itm->second is the topic
            // itm->first is the word
            // i is the doc
            
            // adding this topic to the corpuses
            if(old_topic_corpus.count(itm->second)==0) {
                word_corpus word_corpus_;
                old_topic_corpus[itm->second]=word_corpus_;
            }
            
            // adding this doc to the corpus
            // it is the first time, we see this topic for this doc
            if(topic_newid.count(itm->second)==0) {
                topic_newid[itm->second]=old_topic_corpus[itm->second].docs.size();
                doc doc_;
                old_topic_corpus[itm->second].docs.push_back(doc_);
            }
            
            old_topic_corpus[itm->second].docs[topic_newid[itm->second]].wn_occurences.insert(make_pair(itm->first, oC.docs[i].wn_occurences[itm->first]));
        }
        // storing the doc ids
        doc_topic_newid.push_back(topic_newid);
    }
    
    invert_doc_topic_newid(doc_topic_newid, topic_old_docs);
    
    // asserting old_topic_corpus have labels in order
    for(UI i=0; i<old_topic_corpus.size(); i++) {
        if(old_topic_corpus.count(i)==0) {
            cerr<<"error in second level"<<endl;
            exit(-1);
        }
    } 
    
}


