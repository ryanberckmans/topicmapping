

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





int sample_topic(DD & probs, DI & topics, double & sum) {
    
    // check this
    //prints(topics);
    //prints(probs);
    int nn=lower_bound(probs.begin(), probs.end(), ran4()*sum)-probs.begin();
    //cout<<"nn:: "<<nn<<" "<<topics.size()<<endl;
    return topics[nn];
    
}

void word_corpus::gibbs_sampling(deque<mapii> & doc_assignments) {
    
    
    // you probably need arrays
    
    
    double alpha_hyper=0.01;
    double beta_hyper=0.01;
    
    deque<mapii> n_d_topic;
    map<int, mapii> n_w_topic;
    mapii n_topic;
    
    if(doc_assignments.size()!=docs_.size()) {  cerr<<"doc_assignments size does not match"<<endl; exit(-1);  }
    
    
    // initialization
    RANGE_loop(doc_number, doc_assignments) {
        
        mapii d_topic;
        IT_loop(mapii, itm, doc_assignments[doc_number]) {            
            
            //cout<<"-------------- "<<doc_number<<endl;
            //prints(doc_assignments[i]);
            int_histogram(itm->second, d_topic);
            if(n_w_topic.count(itm->first)==0) {
                mapii new_mapii;
                n_w_topic[itm->first]=new_mapii;
            }
            int_histogram(itm->second, n_w_topic[itm->first]);
            int_histogram(itm->second, n_topic);
        }
        n_d_topic.push_back(d_topic);
    }
    
    // sampling
    for(int iter=0; iter<1000; iter++) RANGE_loop(doc_number, doc_assignments) {
        
        //cout<<"-------------- "<<doc_number<<" "<<n_d_topic[doc_number].size()<<endl;
        
        IT_loop(mapii, itm, doc_assignments[doc_number]) {            
            
            int topic= itm->second;
            int word= itm->first;
            // make sure they are positive
            n_d_topic[doc_number].at(topic)-=1;
            n_w_topic[word].at(topic)-=1;
            n_topic.at(topic)-=1;
            
            if(n_d_topic[doc_number][topic]==0) n_d_topic[doc_number].erase(topic);
            if(n_w_topic[word][topic]==0) n_w_topic[word].erase(topic);
            if(n_topic[topic]==0) n_topic.erase(topic);
            
            // looping over doc_topics
            DI topics;
            DD probs;
            double sum=0.;
            IT_loop(mapii, itm2, n_d_topic[doc_number]) {
                double pr= (itm2->second + alpha_hyper) * (n_w_topic[word][itm2->first]+ beta_hyper);
                probs.push_back(sum+pr);
                topics.push_back(itm2->first);
                sum+=pr;
            }
            // check  norm
            topic=sample_topic(probs, topics, sum);
            itm->second=topic;
            int_histogram(topic, n_d_topic[doc_number]);
            int_histogram(topic, n_w_topic[word]);
            int_histogram(topic, n_topic);            
        }
    }
    
    
    // topics can be empty... check this
    
    deque<mapid> n_d_topic_d;
    map<int, mapid> n_topic_w_d;
    
    
    
    RANGE_loop(doc_number, n_d_topic) {
        mapid new_mapid;
        IT_loop(mapii, itm, n_d_topic[doc_number]) new_mapid[itm->first]=itm->second;
        cout<<doc_number<<" ---- "<<n_d_topic[doc_number].size()<<endl;
        normalize_mapid(new_mapid);
        n_d_topic_d.push_back(new_mapid);
    }
    
    
    for (map<int, mapii>::iterator word_itm= n_w_topic.begin(); 
         word_itm!=n_w_topic.end(); word_itm++) {
        
        IT_loop(mapii, itm2, word_itm->second) {
            
            int word= word_itm->first;
            int topic= itm2->first;
            int count= itm2->second;
            if(n_topic_w_d.count(topic)==0) {
                mapid new_mapid;
                n_topic_w_d[topic]=new_mapid;
            }
            n_topic_w_d[topic][word]=count;
        }
    }
    
    for(map<int, mapid>::iterator topic_itm= n_topic_w_d.begin(); 
        topic_itm!=n_topic_w_d.end(); topic_itm++) {
        
        normalize_mapid(topic_itm->second);
        
    }
    
    
    write_beta_and_theta_files(n_d_topic_d, n_topic_w_d, "thetas_gs.txt", "betas_gs.txt");
    
}




string word_corpus::get_topic_title(const mapid & topic_distr) {
    
    // returning the most probable 20 words in this topic
    
    deque<pair<double, int> > pr_word;
    IT_const_loop(mapid, itm2, topic_distr) {
        pr_word.push_back(make_pair(-itm2->second, itm2->first));
    }
    sort(pr_word.begin(), pr_word.end());
    
    string title="";
    for(int i=0; i<min(20,int(pr_word.size())); i++){
        title+=word_strings_[pr_word[i].second]+" ";
    }
    return title;
    
}




double get_z_score(int N, double p, int real) {
    
    /* this is the z-score of a binomial distribution
     we could actually compute the right cumularive.
     but the approximation should hold pretty well */
    
    
    double expected_usage =  p * N;
    double std= sqrt(max(0., p * N * (1.-p) ));
    double z_score=0;
    if(std<1e-10) {
        z_score= 1e300;
    } else {
        z_score = (double(real) - expected_usage)/std;
    }
    return z_score;
    
}


void invert_doc_topic_newid(deque<mapii> & doc_topic_newid, map<int, mapii> & topic_old_docs) {
    
    topic_old_docs.clear();
    RANGE_loop(i, doc_topic_newid) {
        IT_loop(mapii, itm, doc_topic_newid[i]) {
            
            int topic=itm->first;
            int new_id=itm->second;
            
            if(topic_old_docs.count(topic)==0) {
                mapii mapii_;
                topic_old_docs[topic]=mapii_;
            }
            // in this topic, here is how the new id maps to the original doc
            topic_old_docs[topic][new_id]=i;
        }
    }
}





// it would be better to pass this ====================================================

/*mapid & var_gamma = gammas_ldav_map_[doc_number];
 mapid digamma_gam;
 IT_loop(mapid, topic_pr, var_gamma) digamma_gam[topic_pr->first]=digamma(topic_pr->second);
 */

/*
 bool verbose=false;
 // =========== independent =================== of var_gamma values
 double sum_alphas=0.;
 double sum_lg_alphas=0.;
 
 for (int k = 0; k < num_topics_ldav_; k++) {
 sum_alphas += alphas_ldav_[k];
 sum_lg_alphas += lgamma(alphas_ldav_[k]);
 }
 
 double var_gamma_sum=sum_alphas + docs_[doc_number].num_words_;
 double digsum2 = digamma(var_gamma_sum);
 
 if(verbose) {
 cout<<"sum_alphas local "<<sum_alphas<<" =="<<endl;
 cout<<"digsum... "<<digsum2<<" == "<<digsum<<endl;
 cout<<"alphas_terms:: "<<lgamma(sum_alphas) - sum_lg_alphas<<" == "<<endl; 
 }
 
 double lik_const=-lgamma(var_gamma_sum);
 double likelihood = lgamma(sum_alphas) 
 - sum_lg_alphas 
 - lgamma(var_gamma_sum);
 
 // this sum if over values which are not found in var_gamma
 // var_gamma[k]= alphas_ldav_[k]
 
 // It's probably good to store digamma(alphas_ldav_[k]) !!!!!!!!!!!!!!
 // digamma_gam[k] = digamma(alphas_ldav_[k])
 
 for (int k = 0; k < num_topics_ldav_; k++) if(var_gamma.count(k)==0) {
 
 likelihood +=
 (alphas_ldav_[k] - 1)*(digamma(alphas_ldav_[k]) - digsum) 
 + lgamma(alphas_ldav_[k])
 - (alphas_ldav_[k] - 1)*(digamma(alphas_ldav_[k]) - digsum);
 lik_const+=(alphas_ldav_[k] - 1)*(digamma(alphas_ldav_[k]) - digsum) 
 + lgamma(alphas_ldav_[k])
 - (alphas_ldav_[k] - 1)*(digamma(alphas_ldav_[k]) - digsum);
 
 }
 
 if(verbose) {
 cout<<"compute_likelihood_sparse_constants:: "<<lik_const<<endl;
 cout<<"check compute_likelihood_sparse_constants:: ";
 double digggsum;
 cout<<compute_likelihood_sparse_constants(doc_number, var_gamma, digggsum, sum_alphas)<<endl;
 cout<<"digggsum "<<digggsum<<endl;
 }
 
 assert_floats(lik_const, likelihood_const, "error with const");
 assert_floats(likelihood_alpha, lgamma(sum_alphas) - sum_lg_alphas, "error with lik_alpha");
 */

// =========== independent =================== of var_gamma values

