


void word_corpus::set_class_words_to_zeros() {
    
    
    // initializing class_word_ldav_ and class_total_ldav_
    // class_word_ldav_[topic][wn] 
    // class_total_ldav_[topic]
    

    // ============ this should be made more efficient =======
    class_word_ldav_.clear();
    class_total_ldav_.clear();

    
    DD void_dd_class;
    void_dd_class.assign(word_occurrences_.size(), 0.);
    for(int k=0; k<num_topics_ldav_; k++) {
        class_word_ldav_.push_back(void_dd_class);
    }
    class_total_ldav_.assign(num_topics_ldav_, 0.);
    
}




void word_corpus::initialize_lda_data(deque<mapid> & doc_topic, 
                                      map<int, mapid> & topic_word) {
    
    // this function is getting all data structure for lda topics ready
    // some asserts are also done and should be removed later
    
    // ---------- TEMPORARY ----------
    topic_word.clear();
    ifstream gin("run1/final.beta");
    string gins;
    int count_line=0;
    while(getline(gin, gins)) {
        DD vv;
        cast_string_to_doubles(gins, vv);
        mapid topic_word_mapid;
        double check_sum=0.;
        RANGE_loop(i, vv) {
            topic_word_mapid[i]=exp(vv[i]);
            check_sum+=topic_word_mapid[i];
        }
        if(fabs(check_sum-1)>1e-4) {
            cerr<<"error in check_sum "<<check_sum-1<<endl; exit(-1);
        }
        topic_word[count_line]=topic_word_mapid;
        count_line+=1;
    }
    // TEMPORARY ---------------------
    
    
    // asserting everything is starting from zero and being consecutive
    DI all_topics;
    for (map<int, mapid>::iterator topic_itm= topic_word.begin(); 
         topic_itm!=topic_word.end(); topic_itm++) {
        all_topics.push_back(topic_itm->first);
    }
    assert_consecutive(all_topics);
    // assertation done
    
    //lda data structures
    
    num_topics_ldav_=all_topics.size();
    cout<<"number of topics for lda em: "<<num_topics_ldav_<<endl;
    
    phis_ldav_.clear();
    betas_ldav_.clear();
    alphas_ldav_.clear();
    class_word_ldav_.clear();
    class_total_ldav_.clear();
    gammas_ldav_.clear();
    
    // betas_ldav_ is inverted respect with the original code
    // I think it's better this way (one long vector of short vectors rather 
    // than few long vectors)
    // should I do the same thing on class_word_ldav_?
    
    DD void_dd_numtops;
    void_dd_numtops.assign(num_topics_ldav_, 0.);
    
    RANGE_loop(wn, word_occurrences_) {
        betas_ldav_.push_back(void_dd_numtops);
        // this could be made more efficient (making it depend on single documents)
        phis_ldav_.push_back(void_dd_numtops);
    }
    
    
    // TODO 
    // the initialization of alphas_ldav_ is quite arbitrary 
    // and should be done
    // using doc_topic
    alphas_ldav_.assign(num_topics_ldav_, 0.01);
    
    // copying topic_word in betas_ldav_
    // this should be avoided and pass betas_ldav_ directly
    // to the lda model
    for (map<int, mapid>::iterator topic_itm= topic_word.begin(); 
         topic_itm!=topic_word.end(); topic_itm++) {
        
        int word_wn=0;
        IT_loop(mapid, itm2, topic_itm->second) { 
            assert_ints(itm2->first, word_wn);
            betas_ldav_[word_wn][topic_itm->first]=log(itm2->second);
            ++word_wn;
        }
    }
    
    
    set_class_words_to_zeros();
    
    // initializing gammas
    RANGE_loop(doc_number, docs_) {
        gammas_ldav_.push_back(void_dd_numtops);
    }
    
    
    cout<<"============== DIMENSIONS =============="<<endl;
    cout<<"number of topics for lda em: "<<num_topics_ldav_<<endl;
    
    cout<<"phis_ldav_ "<<phis_ldav_.size()<<" x "<<phis_ldav_[0].size()<<endl;
    cout<<"betas_ldav_ "<<betas_ldav_.size()<<" x "<<betas_ldav_[0].size()<<endl;
    
    cout<<"alphas_ldav_ "<<alphas_ldav_.size()<<endl;
    
    cout<<"class_word_ldav_ "<<class_word_ldav_.size()<<" x "<<class_word_ldav_[0].size()<<endl;
    cout<<"class_total_ldav_ "<<class_total_ldav_.size()<<endl;
    
    cout<<"gammas_ldav_ "<<gammas_ldav_.size()<<" x "<<gammas_ldav_[0].size()<<endl;

}




