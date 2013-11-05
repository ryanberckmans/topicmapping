
# define LIK_precision 1e-5
// exp(-14) ~ 1e-6
# define DIGAMMA_precision -14.
# define VARGAMMA_precision 1e-6



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
    // It's probably good to store digamma(alphas_ldav_[k]) !!!!!!!!!!!!!!
    // digamma_gam[k] = digamma(alphas_ldav_[k])
    
    for (int k = 0; k < num_topics_ldav_; k++) if(var_gamma.count(k)==0) {
        
        likelihood +=   \
                        (alphas_ldav_[k] - 1)*(digamma(alphas_ldav_[k]) - digsum)  \
                        + lgamma(alphas_ldav_[k]) \
                        - (alphas_ldav_[k] - 1)*(digamma(alphas_ldav_[k]) - digsum);
        
    }
    
    return likelihood;
    
    
}

double word_corpus::lda_inference_sparse(int doc_number, const double & sum_alphas, const double & likelihood_alpha) {
    
    // phis_ldav_map_wn
    // should probably be a deque<pair<int, double> >
    // should be faster
    
    /*
        this function is the equivalent of lda_inference
        but here we only update var_gamma for topic 
        which are actually used by the document
        for the other topics var_gamma is going to be simply the prior
        and this is not updated
        similarly, phis_ldav_map_[wn] is computed only for 
        topics used by the document
    */
    
    // 
    
    /*
    RANGE_loop(wn, betas_ldav_map_) {
        cout<<">>> wn "<<wn<<endl;
        prints(betas_ldav_map_[wn]);        
    }*/
    
    
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
        // phi[word][topic]
    }
    
    if(verbose) {
        cout<<"sparse var_gamma::"<<endl;
        prints(var_gamma);
    }
    
    IT_loop(deqii, wn_occ, docs_[doc_number].wn_occs_) {            
        // clearing the map phis_ldav_map_[wn]
        // phis_ldav_map_ must always be used 
        // as phis_ldav_map_wn with loop over this doc
        phis_ldav_map_.at(wn_occ->first).clear();
        deqid & phis_ldav_map_wn= phis_ldav_map_[wn_occ->first];
        
        // consider only the common topic between var_gamma and betas
        //double topics_for_this_word = double(betas_ldav_map_.at(wn_occ->first).size());
        IT_loop(mapid, topic_pr, var_gamma) if(betas_ldav_map_.at(wn_occ->first).count(topic_pr->first)>0) {
            //phis_ldav_map_wn[topic_pr->first] = 1.0/num_topics_ldav_;
            phis_ldav_map_wn.push_back(make_pair(topic_pr->first, 1.0/num_topics_ldav_));
            
        }
        // this considers everything
        //IT_loop(mapid, topic_pr, var_gamma) {
        //    phis_ldav_map_wn[topic_pr->first] = 1.0/num_topics_ldav_;
        //}
    }
    
    
    double digsum=0.;
    double likelihood_const= compute_likelihood_sparse_constants(doc_number, var_gamma, digsum, sum_alphas);

    
    double converged=1;
    double likelihood=0;
    double likelihood_old=0;
    int var_iter=0;
    
    while(converged>1e-5 and var_iter<1000) {        
        
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
            
            // the first time, we update on this word I need to remove the prior
            // from all topics
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
            exit(-1);
        }
        converged = (likelihood_old - likelihood) / likelihood_old;
        likelihood_old = likelihood;
        /*
         cout<<"iter "<<var_iter<<endl;
         prints(var_gamma);
         cout<<"likelihood:: "<<likelihood<<" "<<converged<<endl;
         */
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
    
    //cout<<"likelihood::: "<<likelihood<<endl;
    //return likelihood;
    //exit(-1);
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
        // you should also use topic_pr->!!!!!!!!!!!!!!!
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
                                                               +  get_from_mapid(betas_ldav_map_.at(n), k, -100) )    );
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


double word_corpus::run_em_sparse() {
    
    cout<<"running EM"<<endl;    
    
    // tmp    
    double likelihood_old=-1e300;
    // this are the varaiational parameters which are the main output of the program
    deque<DD> gammas_ldav;
    
    while(true) {
        
        // E steps
        double sum_alphas=0.;
        double likelihood_alpha = compute_likelihood_alpha_terms(sum_alphas);
        //cout<<"likelihood_alpha==== "<<likelihood_alpha<<" sum_alphas "<<sum_alphas<<endl;
        //likelihood_alpha=0.;
        set_class_words_to_zeros_map();
        //gammas_ldav_.clear();
        double likelihood_all=0.;
        //ofstream infout("inf.txt");
        RANGE_loop(doc_number, docs_) {
            if (doc_number%100==0)
                cout<<"running lda-inference for doc "<<doc_number<<endl;
            double likelihood_doc = lda_inference_sparse(doc_number, sum_alphas, likelihood_alpha);
            likelihood_all+=likelihood_doc;
            //infout<<likelihood_doc<<endl;
            //exit(-1);
        }
        //infout.close();
        //exit(-1);
        
        // M step
        // optimizing betas
        
        
        RANGE_loop(wn, word_occurrences_) {
        
            mapid & class_word_ldav_map_wn = class_word_ldav_map_.at(wn);
            mapid & betas_ldav_wn = betas_ldav_map_.at(wn);
            
            IT_loop(mapid, topic_pr, class_word_ldav_map_wn) {
                if(topic_pr->second > 0) {
                    //cout<<"topic "<<topic_pr->first<<" "<<topic_pr->second<<endl;
                    betas_ldav_wn[topic_pr->first] = log(topic_pr->second) \
                                                    - log(class_total_ldav_map_.at(topic_pr->first));
                } else {
                    // this should never happen
                    //betas_ldav_wn[topic_pr->first] = -100;
                    cerr<<"Small value in class_word_ldav_map_wn "<<topic_pr->second<<endl;
                    cerr<<"This is likely a bug! Please contact me: arg.lanci@gmail.com Thanks!"<<endl;
                    exit(-1);
                }
            }
        }
        
        // check betas is  normalized!!!!!!!!!!
        // remove thissssssssssssss
        mapid topic_norm;
        RANGE_loop(wn, word_occurrences_) {
            mapid & betas_ldav_wn = betas_ldav_map_.at(wn);
            IT_loop(mapid, topic_pr, betas_ldav_wn) {
                int_histogram(topic_pr->first, topic_norm, exp(topic_pr->second));
            }        
        }
        cout<<"________________"<<endl;
        prints(topic_norm);
        optimize_alpha_sparse(gammas_ldav);

        cout<<"log likelihood "<<likelihood_all<<endl;
        // this is arbitrary
        if( fabs( ( likelihood_all - likelihood_old ) / likelihood_old ) < LIK_precision )
            break;
        likelihood_old=likelihood_all;
    }
    
    // this was set from optimize_alpha_sparse
    ofstream pout_final("lda_gammas.txt");
    printm(gammas_ldav, pout_final);
    
    return 0.;
}






