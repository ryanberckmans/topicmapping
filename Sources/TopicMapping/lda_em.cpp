
double word_corpus::lda_inference(int doc_number) {
    
    
    //cout<<"alphas_ldav_:: "<<num_topics_ldav_<<endl;
    //prints(alphas_ldav_);
    
    // assert here
    int num_topics= betas_ldav_[0].size();    
    //cout<<"num_topics:: "<<num_topics<<endl;
    
    DD var_gamma;
    var_gamma.assign(num_topics, 0.);
    // psi(var_gamma)
    double digamma_gam[num_topics];
    double oldphi[num_topics];
    
    
    // compute posterior dirichlet
    for(int k=0; k<num_topics; k++) {
        // uniform
        var_gamma[k]=(alphas_ldav_[k] + double(docs_[doc_number].num_words_)/num_topics);
        digamma_gam[k]=digamma(var_gamma[k]);
        // phi[word][topic]
        IT_loop(mapii, itm, docs_[doc_number].wn_occurences_) {            
            phis_ldav_[itm->first][k] = 1.0/num_topics;
        }
    }
    
    double converged=1;
    double likelihood=0;
    double likelihood_old=0;
    int var_iter=0;
    
    while(converged>1e-5 and var_iter<1000) {
        
        var_iter++;
        
        // looping over words in doc
        IT_loop(mapii, itm, docs_[doc_number].wn_occurences_) {
            
            // all this refers to this particular word
            int n= itm->first;
            double phisum;
            
            for (int k = 0; k < num_topics; k++) {
                
                oldphi[k]=phis_ldav_[n][k];
                phis_ldav_[n][k] = digamma_gam[k] + betas_ldav_[n][k];
                if (k > 0)
                    phisum = log_sum(phisum, phis_ldav_[n][k]);
                else
                    phisum = phis_ldav_[n][k]; // note, phi is in log space
            }
            
            for (int k = 0; k < num_topics; k++) {
                
                phis_ldav_[n][k] = exp(phis_ldav_[n][k] - phisum);
                var_gamma[k] += itm->second*(phis_ldav_[n][k] - oldphi[k]);
                digamma_gam[k] = digamma(var_gamma[k]);
            }            
        }
        
        likelihood = compute_likelihood(doc_number, var_gamma);
        
        if(likelihood!=likelihood) {
            cerr<<"error in likelihood:: "<<likelihood<<endl;
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
    gammas_ldav_.push_back(var_gamma);
    
    // updating stats for betas_ldav_
    IT_loop(mapii, itm, docs_[doc_number].wn_occurences_) for (int k = 0; k < num_topics; k++) {
        class_word_ldav_[k][itm->first] += itm->second * phis_ldav_[itm->first][k];
        class_total_ldav_[k] += itm->second * phis_ldav_[itm->first][k];
    }
    
    
    return(likelihood);
}


/*
 * compute likelihood bound
 *
 */


double word_corpus::compute_likelihood(int doc_number, DD & var_gamma) {
    
    
    
    //cout<<"var gamma"<<endl;
    //prints(var_gamma);
    // assert number of topics is consistent (remove this later)
    assert_ints(var_gamma.size(), alphas_ldav_.size());
    assert_ints(var_gamma.size(), betas_ldav_[0].size());
    
    int num_topics= betas_ldav_[0].size();
    
    double dig[num_topics];
    double var_gamma_sum=0.;
    double sum_l_gamma=0.;
    double sum_alphas=0.;
    for (int k = 0; k < num_topics; k++) {
        dig[k] = digamma(var_gamma[k]);
        var_gamma_sum += var_gamma[k];
        sum_l_gamma += lgamma(alphas_ldav_[k]);
        sum_alphas += alphas_ldav_[k];
    }
    double digsum = digamma(var_gamma_sum);
    
    /*
     double likelihood = lgamma(alpha * num_topics) 
     - num_topics * lgamma(alpha) 
     - lgamma(var_gamma_sum);
     */
    
    double likelihood = lgamma(sum_alphas) 
    - sum_l_gamma 
    - lgamma(var_gamma_sum);
    
    
    for (int k = 0; k < num_topics; k++) {
        likelihood +=
        (alphas_ldav_[k] - 1)*(dig[k] - digsum) 
        + lgamma(var_gamma[k])
        - (var_gamma[k] - 1)*(dig[k] - digsum);
        
        
        IT_loop(mapii, itm, docs_[doc_number].wn_occurences_) {
            
            int n=itm->first;
            if (phis_ldav_[n][k] > 0) {
                
                likelihood += itm->second * (phis_ldav_[n][k]*((dig[k] - digsum) - log(phis_ldav_[n][k])
                                                               +  betas_ldav_[n][k] )    );
            }
        }
    }
    return(likelihood);
}





double word_corpus::run_em() {
    
    cout<<"running EM"<<endl;    
    
    // tmp    
    ofstream infout("inf.txt");
    double likelihood_old=-1e300;
    
    while(true) {
        
        // E steps
        set_class_words_to_zeros();
        gammas_ldav_.clear();
        double likelihood_all=0.;
        RANGE_loop(doc_number, docs_) {
            if (doc_number%100==0)
                cout<<"running lda-inference for doc "<<doc_number<<endl;
            double likelihood_doc = lda_inference(doc_number);
            likelihood_all+=likelihood_doc;
            infout<<likelihood_doc<<endl;
        }
        infout.close();
        
        // M step
        // optimizing betas
        for(int k=0; k < num_topics_ldav_; k++) {
            RANGE_loop(wn, word_occurrences_) {
                if(class_word_ldav_[k][wn] > 0)
                    betas_ldav_[wn][k] = log(class_word_ldav_[k][wn]) - log(class_total_ldav_[k]);
                else
                    betas_ldav_[wn][k] = -100;
            }
        }    
        optimize_alpha();
        
        cout<<"log likelihood "<<likelihood_all<<endl;
        
        // this is arbitrary
        if( fabs( ( likelihood_all - likelihood_old ) / likelihood_old ) < 1e-4 )
            break;
        likelihood_old=likelihood_all;
    }
    
    
    ofstream pout_final("lda_gammas.txt");
    printm(gammas_ldav_, pout_final);
    
    return 0.;
}


