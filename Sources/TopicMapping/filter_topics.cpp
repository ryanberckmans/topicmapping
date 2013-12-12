


void word_corpus::get_rid_of_non_prevalent_topics(mapii & hard_mems, DI & doc_prevalent_topics) {
    
    // this function is deleting from hard_mems words which appear in topics
    // which are not prevalent in any document
    
    
    //cout<<"total words in partition: "<<hard_mems.size()<<endl;
    
    SI prevalent_topics;
    RANGE_loop(i, docs_) {
        prevalent_topics.insert(doc_prevalent_topics[i]);
    }
    
    
    //cout<<"prevalent topics: "<<prevalent_topics.size()<<endl;
    mapii new_hard_mems;
    IT_loop(mapii, itm, hard_mems) {
        // if the topic is prevalent in at least one doc
        if(prevalent_topics.find(itm->second)!=prevalent_topics.end()) {
            new_hard_mems.insert(*itm);
        }
    }
    
    hard_mems=new_hard_mems;
    //cout<<"words in partition after removing small topics: "<<hard_mems.size()<<endl;
    
    
    
}


void word_corpus::initial_ptopic(deque<mapid> & doc_topic, map<int, mapii> & word_topic, \
                                 const mapii & hard_mems, const DI & doc_prevalent_topics) {
    
    // we set doc_topic[doc] as p(t|doc)
    // and word_topic[word] as n(w,t)
    // pt[topic] is n(t)
    
    doc_topic.clear();
    word_topic.clear();
    
    RANGE_loop(i, docs_) {
        
        mapid topic_distr;
        
        IT_loop(deqii, itm, docs_[i].wn_occs_) {
            
            if(hard_mems.count(itm->first)>0) {
                // if the word is in the partition
                int topic_num=hard_mems.at(itm->first);
                int_histogram(topic_num, topic_distr, double(itm->second)/docs_[i].num_words_);
                if(word_topic.count(itm->first)==0) {
                    mapii newmapii;
                    word_topic[itm->first]=newmapii;
                }
                int_histogram(topic_num, word_topic[itm->first], itm->second);
            }            
        }
        
        int prevalent_topic = doc_prevalent_topics[i];
        
        IT_loop(deqii, itm, docs_[i].wn_occs_) {
            
            if(hard_mems.count(itm->first)==0) {
                
                // if the word is not here, I use the most significant topic in this doc
                int_histogram(prevalent_topic, topic_distr, double(itm->second)/docs_[i].num_words_);
                if(word_topic.count(itm->first)==0) {
                    mapii newmapii;
                    word_topic[itm->first]=newmapii;
                }
                int_histogram(prevalent_topic, word_topic[itm->first], itm->second);  
            }            
        }
        
        doc_topic.push_back(topic_distr);
        
    }
    
}


void word_corpus::get_betas(map<int, mapii> & word_topic, map<int, mapid> & topic_word, mapid & pt) {
    
    pt.clear();
    topic_word.clear();
    
    for(map<int, mapii>::iterator itm = word_topic.begin(); itm!=word_topic.end(); itm++) {
        
        IT_loop(mapii, itm2, itm->second) if(itm2->second>1e-10) {
            if(topic_word.count(itm2->first)==0) {
                mapid newmapid;
                topic_word[itm2->first]=newmapid;
            }
            int_histogram(itm->first, topic_word[itm2->first], double(itm2->second));
        }
    }
    
    // 
    for(map<int, mapid>::iterator itm = topic_word.begin(); itm!=topic_word.end(); itm++) {
        pt[itm->first]=normalize_mapid(itm->second);
    }
    normalize_mapid(pt);
}


double word_corpus::compute_likelihood(map<int, mapid> & topic_word, map<int, mapii> & word_topic, deque<mapid> & doc_topic) {
    
    double loglik=0.;
    RANGE_loop(i, docs_) {
        
        
        IT_loop(deqii, itm, docs_[i].wn_occs_) {
            double pr=0;
            mapii & possible_topics = word_topic[itm->first];
            
            // intersection of possible_topics & doc_topic[i]
            mapid intermap;
            mapid::iterator out_itr( intermap.begin() );
            set_intersection(doc_topic[i].begin(), doc_topic[i].end(), 
                             possible_topics.begin(), possible_topics.end(), 
                             inserter( intermap, out_itr ), doc_topic[i].value_comp());
            
            IT_loop(mapid, itm2, intermap) {
                if(topic_word.count(itm2->first)>0) {
                    pr+= itm2->second * get_from_mapid(topic_word.at(itm2->first), itm->first);
                }
            }
            
            if(pr==0) {
                
                cout<<"doc "<<i<<endl;
                cout<<"word "<<itm->first<<endl;
                prints(doc_topic[i]);
                prints(word_topic[itm->first]);                 
                cerr<<"error in pr"<<endl;
                exit(-1);
            }
            loglik+= itm->second * log(pr);
        }
    }
    return loglik;
}



double word_corpus::likelihood_filter(map<int, mapii> & word_topic, deque<mapid> & doc_topic, \
                                      map<int, mapid> & topic_word, mapid & pt, \
                                      const mapii & hard_mems, \
                                      double filtering_par, const DI & doc_prevalent_topics,\
                                      deque<mapii> & doc_assignments) {
    
    topic_word.clear();
    pt.clear();
    
    // we set doc_topic[doc] as p(t|doc)
    // and word_topic[word] as n(w,t)
    
    //doc_assignments[doc][wn] is the topic to which the word has been assigned
    doc_assignments.clear();
    RANGE_loop(i, docs_) {
        mapii void_mapii_;
        doc_assignments.push_back(void_mapii_);
    }
    
    
    RANGE_loop(i, docs_) {
        
        int prevalent_topic = doc_prevalent_topics[i];
        
        IT_loop(deqii, itm, docs_[i].wn_occs_) {
            
            // topic num is either the prevalent_topic or the topic
            // the word belongs to originally 
            int topic_num = prevalent_topic;
            if(hard_mems.count(itm->first)>0) {
                topic_num=hard_mems.at(itm->first);
            }
            // if the topic is small, we move words to prevalent_topic
            if(doc_topic[i].at(topic_num)<filtering_par and topic_num!=prevalent_topic) {
                word_topic[itm->first][topic_num]-=itm->second;
                int_histogram(prevalent_topic, word_topic[itm->first], itm->second);
            }
            // now topic_num is prevalent_topic, if it is small
            if(doc_topic[i].at(topic_num)<filtering_par) {
                topic_num=prevalent_topic;
            }
            // we insert the word assignement
            doc_assignments[i][itm->first]=topic_num;
        }
        
        IT_loop(mapid, itm, doc_topic[i]) if(itm->second<filtering_par and itm->first!=prevalent_topic) {
            doc_topic[i][prevalent_topic]+=itm->second;
            itm->second=0;
        }
    }
    
    get_betas(word_topic, topic_word, pt);
    
    /**
     // checking word_topic is correct
     for(map<int, mapii>::iterator itm = word_topic.begin(); itm!=word_topic.end(); itm++) {
     int occ=0; IT_loop(mapii, itm2, itm->second) occ+=itm2->second;
     if(occ!=word_occurrences_[itm->first]) {
     cerr<<"occ is wrong for word  "<<itm->first<<endl;
     exit(-1);
     }
     }*/
    
    
    if(false) {
        
        for(map<int, mapii>::iterator itm = word_topic.begin(); itm!=word_topic.end(); itm++) {
            cout<<"word:: "<<itm->first<<endl;
            prints(itm->second);
        }
        RANGE_loop(i, docs_) {
            cout<<"doc:: "<<i<<endl;
            prints(doc_topic[i]);
        }
        for(map<int, mapid>::iterator itm = topic_word.begin(); itm!=topic_word.end(); itm++) {
            cout<<"topic:: "<<itm->first<<endl;
            prints(itm->second);
        }
        cout<<"pt"<<endl;
        prints(pt);
    }
    
    return compute_likelihood(topic_word, word_topic, doc_topic);
    
}








void word_corpus::get_prevalent_topics(DI & doc_prevalent_topics, mapii & hard_mems) {
    
    /* first I need to compute the probability of using a topic
     in the null model, where we draw words at random according to 
     simply their probability */
    
    
    log_fact_table log_fact_table_;
    int max_doc_size=0;
    RANGE_loop(i, docs_) max_doc_size=max(max_doc_size, int(docs_[i].num_words_));
    log_fact_table_._set_(max_doc_size);
    // for each topic, the probability we use it in the null model
    mapid topic_pr;
    
    RANGE_loop(wn, word_occurrences_) {
        if (hard_mems.count(wn)>0) {
            int topic_num=hard_mems.at(wn);
            int_histogram(topic_num, topic_pr, double(word_occurrences_[wn]));
        }
    }
    
    // so far, for each topic I am counting the occurrences of 
    // the words in that topic
    normalize_mapid(topic_pr);
    
    double min_number_of_docs=double(min_docs_);
    // topics which are used in a very few documents
    SI banned_topics;
    IT_loop(mapid, itm, topic_pr) if( itm->second < min_number_of_docs / docs_.size() ) {
        banned_topics.insert(itm->first);
    }
    
    
    RANGE_loop(i, docs_) {
        
        mapii topic_usage;
        // getting topic_usage for this doc
        IT_loop(deqii, itm, docs_[i].wn_occs_) {
            if(hard_mems.count(itm->first)>0) {
                int topic_num=hard_mems.at(itm->first);
                int_histogram(topic_num, topic_usage, itm->second);
            }
        }
        
        double smallest_binpvalue=2;
        double smallest_binpvalue_non_banned=2;
        int prevalent_topic=-1;
        int prevalent_topic_non_banned=-1;
        
        IT_loop(mapii, itm, topic_usage) {            
            // topic_id is itm->first
            // usage is itm->second            
            
            double bin_pvalue= log_fact_table_.cum_binomial_right(itm->second, docs_[i].num_words_, topic_pr.at(itm->first));
            if(bin_pvalue<smallest_binpvalue_non_banned and banned_topics.find(itm->first) == banned_topics.end()) {
                prevalent_topic_non_banned=itm->first;
                smallest_binpvalue_non_banned=bin_pvalue;
            }
            if(bin_pvalue<smallest_binpvalue) {
                prevalent_topic=itm->first;
                smallest_binpvalue=bin_pvalue;
            }
        }
        // prevalent_topic==-1 means that the document had just filtered out words
        if(prevalent_topic_non_banned!=-1) {
            doc_prevalent_topics.push_back(prevalent_topic_non_banned);
        } else {
            doc_prevalent_topics.push_back(prevalent_topic);
        }
    }
    
    
    get_rid_of_non_prevalent_topics(hard_mems, doc_prevalent_topics);
    int number_of_topics=0; IT_loop(mapii, itm, hard_mems) number_of_topics=max(number_of_topics, itm->second);
    
    RANGE_loop(i, doc_prevalent_topics) {
        int prevalent_topic=doc_prevalent_topics[i];
        if(prevalent_topic==-1) {
            ++number_of_topics;
            prevalent_topic=number_of_topics;
            doc_prevalent_topics[i]=prevalent_topic;
        }
    }
}



double word_corpus::optimal_filtering(mapii & hard_mems, \
                                      mapid & pt_best, \
                                      deque<mapid> & doc_topic_best, \
                                      map<int, mapid> & topic_word_best, \
                                      bool verbose, double step, string out_dir) {
    
    
    // for each document, the topic which we believe
    // the document mostly belongs to
    DI doc_prevalent_topics;
    get_prevalent_topics(doc_prevalent_topics, hard_mems);    
    
    // this is the starting point (the same for all values of filtering)
    deque<mapid> doc_topic_initial;
    map<int, mapii> word_topic_initial;
    if(verbose) cout<<"initial step"<<endl;
    initial_ptopic(doc_topic_initial, word_topic_initial, hard_mems, doc_prevalent_topics);
    if(verbose) cout<<"initial step done"<<endl;
    
    // these are the parameters of the best model
    double max_ll=-1e300;
    double optimal_par=0.;
    
    double eff_ntopics=0.;

    //doc_assignments[doc][wn] is the topic to which the word has been assigned
    deque<mapii> doc_assignments_best;

    for(double filtering_par=min_filter_; filtering_par<max_filter_+step; filtering_par+=step) {
        
        map<int, mapii> word_topic; // for each word, {topic:usage}
        deque<mapid> doc_topic;     // for each doc, {topic:probability}
        map<int, mapid> topic_word;
        mapid pt;
        deque<mapii> doc_assignments;
        
        word_topic= word_topic_initial;
        doc_topic= doc_topic_initial;
        
        // filtering
        double loglikelihood=likelihood_filter(word_topic, doc_topic, topic_word, pt, hard_mems, filtering_par, doc_prevalent_topics, doc_assignments);
        
        if(verbose) {
            cout<<"filtering: "<<filtering_par<<" loglikelihood: "<<loglikelihood<<" #topics: "<<topic_word.size()<<endl;
        }
        
        if(loglikelihood>max_ll) {
            max_ll=loglikelihood;
            optimal_par=filtering_par;
            doc_topic_best=doc_topic;
            topic_word_best=topic_word;
            doc_assignments_best=doc_assignments;
            pt_best=pt;
            eff_ntopics=compute_eff_num_topics(pt);
            if(verbose) cout<<"best filtering so far: "<<optimal_par<<endl;
        }
        
    }
    
    // sorting topic names so that they start from zero and there are no gaps
    make_topic_names_consecutive(doc_topic_best, topic_word_best, doc_assignments_best, pt_best);
    
    if(verbose) cout<<"optimal filtering: "<<optimal_par<<endl;
    
    ofstream asgout((out_dir+"/plsa_word_assignments.txt").c_str());
    RANGE_loop(doc_number, doc_assignments_best) {

        IT_loop(deqii, wn_occ, docs_[doc_number].wn_occs_) {
            int best_topic= doc_assignments_best[doc_number].at(wn_occ->first);
            asgout<<wn_occ->first<<" "<<word_strings_.at(wn_occ->first)<<" ";
            asgout<<wn_occ->second<<" "<<best_topic<<" ";
        }
        asgout<<endl;
    }
    asgout.close();
    
    return eff_ntopics;
    
}




double word_corpus::dimap(int Nruns, \
                          double step,
                          bool print_sig_words,
                          mapid & pt, \
                          deque<mapid> & doc_topic_best, \
                          map<int, mapid> & topic_word_best, \
                          string out_dir) {
    
    doc_topic_best.clear();
    topic_word_best.clear();
    pt.clear();
    
    bool verbose=true;
    
    mapii hard_memberships;
    if (partition_file_.size()==0) {
        
        DI links1;
        DI links2;
        DD weights;
        // running null model (without parallelization going on)
        null_model(links1, links2, weights, 0, 0, 1, print_sig_words, out_dir);
        
        // collecting infomap initial partition
        if(links1.size()==0) {
            cerr<<"empty graph"<<endl;
            cerr<<"Please try again with a bigger p-value"<<endl;
            exit(-1);
        }
        
        get_infomap_partition_from_edge_list(Nruns, irand(100000000), \
                                             links1, links2, weights, \
                                             hard_memberships, verbose, out_dir);
        links1.clear();
        links2.clear();
        weights.clear();
        write_partition(hard_memberships, out_dir);
    } else {
        cout<<"reading partition from file: "<<partition_file_<<endl;
        int_matrix ten;
        get_partition_from_file(partition_file_, ten);
        RANGE_loop(i, ten) RANGE_loop(j, ten[i]) hard_memberships[ten[i][j]]=i;
    }
    
    // max-likelihood filter
    double eff_ntopics=optimal_filtering(hard_memberships, 
                                         pt, doc_topic_best,
                                         topic_word_best,
                                         verbose, step, out_dir);
    
    return eff_ntopics;
    
}





