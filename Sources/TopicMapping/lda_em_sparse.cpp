
# define LIK_precision 1e-5
# define MAX_ITER 1000
# define DIGAMMA_precision -14.     // exp(-14) ~ 1e-6
# define VARGAMMA_precision 1e-6
// use the following values if you do not want 
// to use sparse data structures
//# define DIGAMMA_precision -1000000000.
//# define VARGAMMA_precision 0.
# define SMALL_LOG -100


double word_corpus::compute_likelihood_sparse_constants(int doc_number, \
                                                        mapid & var_gamma, \
                                                        double & digsum, \
                                                        const double & sum_alphas) {
    
    
    double likelihood=0.;
    // =========== likelihood terms which are independent of var_gamma values
    double var_gamma_sum= sum_alphas + docs_[doc_number].num_words_;
    digsum = digamma(var_gamma_sum);
    likelihood= -lgamma(var_gamma_sum);
    
    // this sum if over values which are not found in var_gamma
    // var_gamma[k]= alphas_ldav_[k]
    // It's probably good to store digamma(alphas_ldav_[k])
    // (micro-optimization?)
    // digamma_gam[k] = digamma(alphas_ldav_[k])
    
    for (int k = 0; k < num_topics_ldav_; k++) if(var_gamma.count(k)==0) {
        
        likelihood +=   \
                        (alphas_ldav_[k] - 1)*(digamma(alphas_ldav_[k]) - digsum)  \
                        + lgamma(alphas_ldav_[k]) \
                        - (alphas_ldav_[k] - 1)*(digamma(alphas_ldav_[k]) - digsum);
        
    }
    
    return likelihood;
    
    
}

double word_corpus::lda_inference_sparse(int doc_number, const double & sum_alphas,\
                                         const double & likelihood_alpha,\
                                         bool print_assignments, ostream & asgout) {
    
    
    /*
        this function is the equivalent of lda_inference
        but here we only update var_gamma for topic 
        which are actually used by the document
        for the other topics var_gamma is going to be simply the prior
        and this is not updated
        similarly, phis_ldav_map_[wn] is computed only for 
        topics used by the document
    */
    
    
        
    
    bool verbose=false;
    
    general_assert(betas_ldav_map_.size()>0, "empty betas");
    
    mapid & var_gamma = gammas_ldav_map_[doc_number];
    mapid digamma_gam;
    
    if(verbose) cout<<"doc topic space: "<<var_gamma.size()<<endl;
    
    // compute posterior dirichlet
    IT_loop(mapid, topic_pr, var_gamma) {
        // uniform
        topic_pr->second = (alphas_ldav_[topic_pr->first] + double(docs_[doc_number].num_words_)/num_topics_ldav_);
        digamma_gam[topic_pr->first]=digamma(topic_pr->second);
    }
    
    if(verbose) {
        cout<<"sparse var_gamma::"<<endl;
        prints(var_gamma);
    }
    
    IT_loop(deqii, wn_occ, docs_[doc_number].wn_occs_) {            
        // clearing the map phis_ldav_map_[wn] only
        // !!! phis_ldav_map_ MUST always be used
        // !!! as phis_ldav_map_wn with wn over this doc
        phis_ldav_map_.at(wn_occ->first).clear();
        deqid & phis_ldav_map_wn= phis_ldav_map_[wn_occ->first];
        
        // consider only the common topic between var_gamma and betas
        //double topics_for_this_word = double(betas_ldav_map_.at(wn_occ->first).size());
        IT_loop(mapid, topic_pr, var_gamma) if(betas_ldav_map_.at(wn_occ->first).count(topic_pr->first)>0) {
            phis_ldav_map_wn.push_back(make_pair(topic_pr->first, 1.0/num_topics_ldav_));
            
        }
    }
    
    
    double digsum=0.;
    double likelihood_const= compute_likelihood_sparse_constants(doc_number, var_gamma, digsum, sum_alphas);

    
    double converged=1;
    double likelihood=0;
    double likelihood_old=0;
    int var_iter=0;
    
    while(converged>LIK_precision and var_iter<MAX_ITER) {        
        
        ++var_iter;
        // looping over words in doc
        // updating phis_ldav_
        
        IT_loop(deqii, wn_occ, docs_[doc_number].wn_occs_) {
        
            //======================================================
            // all this refers to this particular word (wn_occ->first)
            mapid oldphi;
            double phisum=-1;
            deqid & phis_ldav_map_wn= phis_ldav_map_[wn_occ->first];
            
            // loop over topics for this word
            bool first_time_here=true;
            IT_loop(deqid, topic_pr, phis_ldav_map_wn) {
                
                oldphi[topic_pr->first] = topic_pr->second;
                topic_pr->second =  digamma_gam[topic_pr->first] + 
                                    betas_ldav_map_.at(wn_occ->first).at(topic_pr->first);
                if (first_time_here==false)
                    phisum = log_sum(phisum, topic_pr->second);
                else {
                    phisum = topic_pr->second; // note, phi is in log space
                    first_time_here=false;
                }
            }
            
            // the first time we update on this word,
            // I need to remove the prior from all topics
            if(var_iter==1) {
                IT_loop(mapid, topic_pr, var_gamma) {
                    topic_pr->second -= double(wn_occ->second) /num_topics_ldav_;
                    digamma_gam[topic_pr->first]=digamma(topic_pr->second);
                }
                IT_loop(deqid, topic_pr, phis_ldav_map_wn) {
                    topic_pr->second = exp(topic_pr->second - phisum);
                    var_gamma[topic_pr->first] += wn_occ->second * (topic_pr->second);
                    digamma_gam[topic_pr->first]=digamma(var_gamma[topic_pr->first]);
                }
                
            } else {
                // 
                IT_loop(deqid, topic_pr, phis_ldav_map_wn) {
                    topic_pr->second = exp(topic_pr->second - phisum);
                    var_gamma[topic_pr->first] += wn_occ->second * (topic_pr->second - oldphi[topic_pr->first]);
                    digamma_gam[topic_pr->first]=digamma(var_gamma[topic_pr->first]);
                }
            }
            
            if(verbose) {
                cout<<"wn:: "<<wn_occ->first<<endl;
                //prints(phis_ldav_map_wn);
            }
            
            //======================================================
        }
        
        likelihood = compute_likelihood_sparse(doc_number, var_gamma, digamma_gam, digsum,
                                               likelihood_alpha, likelihood_const);
        
        if(likelihood!=likelihood) {
            cerr<<"error in likelihood:: "<<likelihood<<" doc_number:: "<<doc_number<<endl;
            cerr<<"This should not happen. Please contact me: arg.lanci@gmail.com"<<endl;
            exit(-1);
        }
        converged = (likelihood_old - likelihood) / likelihood_old;
        likelihood_old = likelihood;
    }
    
    // updates for M steps
    // updating stats for betas_ldav_
    IT_loop(deqii, wn_occ, docs_[doc_number].wn_occs_) {
    
        deqid & phis_ldav_map_wn= phis_ldav_map_[wn_occ->first];
        mapid & class_word_ldav_map_wn = class_word_ldav_map_.at(wn_occ->first);
        
        IT_loop(deqid, topic_pr, phis_ldav_map_wn) {
            const int & k = topic_pr->first;
            int_histogram(k, class_word_ldav_map_wn,  wn_occ->second * topic_pr->second );
            int_histogram(k, class_total_ldav_map_, wn_occ->second * topic_pr->second);
        }
    }
    
    if(verbose) {
        cout<<"alphas_ldav_"<<endl;
        prints(alphas_ldav_);
    }
    
    /*
     if exp(digamma_gam[topic_num])<< 1 
     and var_gamma is very close to the prior,
     it means this doc is barely using topic_num
     (no words have been assigned to that topic almost at all)
     therefore, we can remove it from var_gammas,
     meaning we set var_gamma[topic_num] == alphas_ldav_[topic_num]
     in other words, we are making a further approximation:
     we are going to explore only the topic 
     space where this doc is not using topic_num.
     this constraint make the likelihood optimization a little bit worse
     but the algorithm is WAY faster.
    */
    
    if(verbose)
        cout<<"======== erasing stuff ============="<<endl;
    
    IT_loop(mapid, itm, digamma_gam) {
        if  (\
             (itm->second < DIGAMMA_precision) \
             and \
             (var_gamma.at(itm->first) - alphas_ldav_[itm->first] < VARGAMMA_precision ) \
             ) {
            // erasing this topic from the possible ones for this doc
            if(verbose) cout<<"erasing this::: "<<exp(itm->second) <<var_gamma.at(itm->first) - alphas_ldav_[itm->first]<<endl;
            var_gamma.erase(itm->first);
        }
        else {
            if(verbose)  cout<<"not erasing this::: "<<exp(itm->second) <<" "<<var_gamma.at(itm->first) - alphas_ldav_[itm->first]<<endl;
        }
    }
    
    
    if(verbose) cout<<"=========doc topic space: "<<var_gamma.size()<<endl;
        
    if(print_assignments) {
        
        // writing word assignments for this document
        IT_loop(deqii, wn_occ, docs_[doc_number].wn_occs_) {
            deqid & phis_ldav_map_wn= phis_ldav_map_[wn_occ->first];            
            // loop over topics for this word
            int best_topic=-1;  // topics are all positive
            double best_prob=0.;
            IT_loop(deqid, topic_pr, phis_ldav_map_wn) {
                if(best_topic==-1 or topic_pr->second>best_prob) {
                    best_topic=topic_pr->first;
                    best_prob=topic_pr->second;
                }
            }
            // wn word_str occ topic
            asgout<<wn_occ->first<<" "<<word_strings_.at(wn_occ->first)<<" ";
            asgout<<wn_occ->second<<" "<<best_topic<<" ";
        }
        asgout<<endl;
    //
    }

    return likelihood;
}



double word_corpus::compute_likelihood_sparse(int doc_number, \
                                              mapid & var_gamma, mapid & digamma_gam, 
                                              const double & digsum, 
                                              const double & likelihood_alpha,
                                              const double & likelihood_const) {
    
    
    double likelihood=0.;
    IT_loop(mapid, topic_pr, var_gamma) {
        
        int k = topic_pr->first;
        // the part with the digsum could be removed from here
        // probably a micro-optimization, though  
        likelihood += (alphas_ldav_[k] - 1) * (digamma_gam[k] - digsum) \
                        + lgamma(topic_pr->second) \
                        - (topic_pr->second - 1) * (digamma_gam[k] - digsum);
    }    
    
    
    IT_loop(deqii, wn_occ, docs_[doc_number].wn_occs_) {
        const int & n = wn_occ->first;
        deqid & phis_ldav_map_wn= phis_ldav_map_[n];
        IT_loop(deqid, topic_pr, phis_ldav_map_wn) {
            const int & k = topic_pr->first;
            if (topic_pr->second>0) {
                likelihood += wn_occ->second * (topic_pr->second*((digamma_gam[k] - digsum) - log(topic_pr->second)
                                                               +  get_from_mapid(betas_ldav_map_.at(n), k, SMALL_LOG) )    );
            }
        }
    }
    

    return likelihood + likelihood_alpha + likelihood_const;

}


double word_corpus::compute_likelihood_alpha_terms(double & sum_alphas) {
    
    // part of likelihood bound which does not depend on doc
    sum_alphas=0.;
    double sum_lg_alphas=0.;
    
    for (int k = 0; k < num_topics_ldav_; k++) {
        sum_alphas += alphas_ldav_[k];
        sum_lg_alphas += lgamma(alphas_ldav_[k]);
    }
        
    return (lgamma(sum_alphas) - sum_lg_alphas);
    
}


double word_corpus::E_step(bool verbose, string out_dir, int iter) {
    
    
    // performing E step and returning likelihood
    // if verbose, write likelihood of each doc to a file
    
    
    string iter_s;
    if (iter<0) {
        iter_s="final";
    } else{
        char iter_c[100];
        sprintf(iter_c, "%d", iter);
        iter_s=string(iter_c);
    }

    
    
    ofstream likout;
    ofstream asgout;
    if(verbose) {
        likout.open((out_dir+"/lda_log_likelihood_per_doc_"+iter_s+".txt").c_str());
        asgout.open((out_dir+"/lda_word_assignments_"+iter_s+".txt").c_str());
        cout<<"final E step"<<endl;
    }
    double sum_alphas=0.;
    double likelihood_alpha = compute_likelihood_alpha_terms(sum_alphas);
    
    set_class_words_to_zeros_map();
    
    double likelihood_all=0.;
    RANGE_loop(doc_number, docs_) {
        if (doc_number%10000==0 and doc_number>0)
            cout<<"running lda E-step for doc "<<doc_number<<endl;
        double likelihood_doc = lda_inference_sparse(doc_number, sum_alphas, \
                                                     likelihood_alpha,\
                                                     verbose, asgout);
        // writing likelihood of this doc to file
        if(verbose)
            likout<<likelihood_doc<<endl;
        likelihood_all+=likelihood_doc;
        
    }
    asgout.close();
    likout.close();
    return likelihood_all;
    
    
}



double word_corpus::run_em_sparse(bool skip_alpha_opt, bool infer_flag,\
                                  int print_lag, string out_dir) {
    
    cout<<"running EM"<<endl;    
    
    ofstream likvalue((out_dir+"/lda_log_likelihood.txt").c_str());    
    double likelihood_old=-1e300;
    // these are the variational parameters
    // which are the main output of the program
    
    int iter=0;
    while(infer_flag==false and iter<MAX_ITER) {
        
        ++iter;
        // E step
        cout<<"E step "<<iter<<endl;
        double likelihood_all = E_step(iter%print_lag==1, out_dir, iter);

        // M step
        // optimizing betas
        RANGE_loop(wn, word_occurrences_) {
        
            mapid & class_word_ldav_map_wn = class_word_ldav_map_.at(wn);
            mapid & betas_ldav_wn = betas_ldav_map_.at(wn);
            
            IT_loop(mapid, topic_pr, class_word_ldav_map_wn) {
                if(topic_pr->second > 0) {
                    betas_ldav_wn[topic_pr->first] = log(topic_pr->second) \
                                                    - log(class_total_ldav_map_.at(topic_pr->first));
                }
            }
        }
        // optimizing alphas
        if (skip_alpha_opt==false) {
            optimize_alpha_sparse();
        }
        
        cout<<"log likelihood "<<likelihood_all<<endl;
        likvalue<<iter<<" "<<likelihood_all<<endl;
        if( fabs( ( likelihood_all - likelihood_old ) / likelihood_old ) < LIK_precision or \
           likelihood_old>likelihood_all )
            break;
        likelihood_old=likelihood_all;
        
        if(iter%print_lag==1) {
            print_lda_results(out_dir, iter);
        }
    }
    
    
    //final E step
    double likelihood_all = E_step(true, out_dir, -1);

    cout<<"log likelihood "<<likelihood_all<<endl;
    likvalue<<iter+1<<" "<<likelihood_all<<endl;
    likvalue.close();
    // -1 because it's final
    print_lda_results(out_dir, -1);
    
    // printing lda_class_words (just at the very end)
    map<int, mapid> topic_word;
    from_class_to_topic_word(class_word_ldav_map_, topic_word, false);
    print_topic_sparse_format_short(topic_word, out_dir+"/lda_class_words.txt");
    

    return 0.;
}





void word_corpus::print_lda_results(string out_dir, int iter) {
    
    string iter_s;
    if (iter<0) {
        iter_s="final";
    } else{
        char iter_c[100];
        sprintf(iter_c, "%d", iter);
        iter_s=string(iter_c);
    }
    
    
    cout<<"printing LDA results"<<endl;
    deque<DD> gammas_ldav;
    compute_non_sparse_gammas(gammas_ldav); 
    ofstream pout_final((out_dir+"/lda_gammas_"+iter_s+".txt").c_str());
    printm(gammas_ldav, pout_final);
    pout_final.close();
    
    map<int, mapid> topic_word;
    from_beta_to_topic_word(betas_ldav_map_, topic_word);
    
    DD ptopic;
    int num_words_shown=100;
    // p(topic)
    get_ptopic_distr(ptopic, gammas_ldav);
    
    cout<<"topics printed: "<<topic_word.size()<<" words: "<<betas_ldav_map_.size()<<endl;
    print_topic_sparse_format_complete(topic_word, out_dir+"/lda_betas_sparse_"+iter_s+".txt",\
                                       out_dir+"/lda_summary_"+iter_s+".txt", \
                                       word_strings_, \
                                       ptopic, num_words_shown);
    
    ofstream alpha_out((out_dir+"/lda_alphas_"+iter_s+".txt").c_str());
    prints(alphas_ldav_, alpha_out);
    alpha_out.close();
    

}
