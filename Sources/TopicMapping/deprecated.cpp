

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
                topic_newid[itm->second]=old_topic_corpus[itm->second].docs_.size();
                doc doc_;
                old_topic_corpus[itm->second].docs_.push_back(doc_);
            }
            
            old_topic_corpus[itm->second].docs_[topic_newid[itm->second]].wn_occurences_.insert(make_pair(itm->first, oC.docs_[i].wn_occurences_[itm->first]));
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



//P.set_int("-subt", 0., false, "[int]: minimum number of documents per subtopic. Default is 0, but 10 is recommended for big corpuses.");    
//P.set_double("-conv", 1e-8, false, "[double]: if infomap relative gain is smaller than this, Infomap stops. Default: 1e-8.");
//P.set_bool("-nos", false, false, ": no subtopics are provided.");
//P.set_int("-subdocs", 10, false, "[int]: minimum size of each subtopic (#docs_). Default: 10.");
//P.set_int("-subwords", 10, false, "[int]: minimum size of each subtopic (#words). Default: 10.");
//P.set_bool("-fullout", false, false, ": writes thetas and betas file as well. Not recommended for big corpuses (lots of zeros)");

