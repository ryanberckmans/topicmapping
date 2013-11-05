

// values below this threshold are 
// not included in betas
# define SPARSE_limit 1e-12



void word_corpus::set_class_words_to_zeros_map() {
    
    // initializing class_word_ldav_ and class_total_ldav_
    // class_word_ldav_map_[wn][topic] 
    // class_total_ldav_map_[topic]
    
    // class_word_ldav_map_ should be already initialized
    assert_ints(word_occurrences_.size(), class_word_ldav_map_.size());
    RANGE_loop(wn, class_word_ldav_map_) {
        class_word_ldav_map_[wn].clear();
    }
    class_total_ldav_map_.clear();
    
}





void word_corpus::initialize_lda_data(map<int, mapid> & topic_word, double alphas_init) {
    
    // this function is getting all data structure for lda topics ready
    // some asserts are also done
    
    
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
    
    // initializing alphas_ldav_
    alphas_ldav_.clear();
    cout<<"alpha initialized with:: "<<alphas_init<<endl;
    alphas_ldav_.assign(num_topics_ldav_, alphas_init);


    cout<<"number of topics for LDA:: "<<num_topics_ldav_<<endl;    
        
    // =============== sparse data structure initialization ======================= //
    
    phis_ldav_map_.clear();
    betas_ldav_map_.clear();
    class_word_ldav_map_.clear();
    class_total_ldav_map_.clear();
    gammas_ldav_map_.clear();

    
    // betas_ldav_map_, phis_ldav_map_ and class_word_ldav_map_ have all the same structure
    // betas_ldav_map_[wn]= { topic : value }
    mapid void_mapid;
    deqid void_deqid;
    RANGE_loop(wn, word_occurrences_) {        
        betas_ldav_map_.push_back(void_mapid);
        phis_ldav_map_.push_back(void_deqid);
        class_word_ldav_map_.push_back(void_mapid);
    }
    
    // copying topic_word in betas_ldav_
    // this is basically topic word
    // but stored as [wn][topic]
    // rather than [topic][wn]
    // also I take log
    // and I remove words below SPARSE_limit
    // (which should be none actually)
    
    for (map<int, mapid>::iterator topic_itm= topic_word.begin(); 
         topic_itm!=topic_word.end(); topic_itm++) {
        IT_loop(mapid, itm2, topic_itm->second) { 
            int word_wn=itm2->first;
            assert_ints(itm2->first, word_wn);
            // words which are below this threshold should never 
            // be relevant in this topic
            if(itm2->second>SPARSE_limit)
                betas_ldav_map_[word_wn][topic_itm->first] = log(itm2->second);
        }
    }

    
    
    // ========== check betas is normalized ===============
    cout<<"check betas normalization"<<endl;
    mapid topic_norm;
    RANGE_loop(wn, word_occurrences_) {
        mapid & betas_ldav_wn = betas_ldav_map_.at(wn);
        IT_loop(mapid, topic_pr, betas_ldav_wn) {
            int_histogram(topic_pr->first, topic_norm, exp(topic_pr->second));
        }
    }
    IT_loop(mapid, itm, topic_norm) {
        assert_floats(itm->second, 1., \
                      "error in betas normalization", 1e-5);
    }
    // ========== check betas is normalized ===============
    
    
    set_class_words_to_zeros_map();
    
    RANGE_loop(doc_number, docs_) {
        // we only keep track of topics
        // which actually have words in doc_number
        // for the rest:
        // gammas_ldav_map_[missing_topic]=alphas_ldav_[missing_topic]
        // but this is not recorded
        gammas_ldav_map_.push_back(void_mapid);
        IT_loop(deqii, wn_occ, docs_[doc_number].wn_occs_) {
            IT_loop(mapid, topic_prob, betas_ldav_map_[wn_occ->first]) {
                // different normalizations could be possible here
                int_histogram(topic_prob->first, gammas_ldav_map_[doc_number], 0.);
            }
        }
    }
    
    
    

    
    
    
    
}




