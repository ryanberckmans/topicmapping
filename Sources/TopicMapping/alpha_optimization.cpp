





void word_corpus::optimize_alpha(deque<DD> & gammas_ldav) {
    
    
    /*
        This function could be optimized further.
        For instance, we could use the previous alphas_ldav_ as initial 
        conditions, instead of getting them for scratch every time
        Also, we could store all the sufficient statistics
        rather than getting them from gammas_ldav
        However, the bottleneck of the algorithm is not here.
    */
    
    
    double sum_alphas=0.;
    RANGE_loop(k, alphas_ldav_) sum_alphas+=alphas_ldav_[k];
    
    // ps[i][k] is the probability of topic k in sample i
    deque<DD> ps;
    
    DD sum_gamma_over_topics;
    sum_gamma_over_topics.assign(gammas_ldav.size(), 0.);
    
    RANGE_loop(i, gammas_ldav) {
        ps.push_back(gammas_ldav[i]);
        double sum_=normalize_one(ps[i]);
        // this assert is true but there might be round-off problems
        //assert_floats((docs_[i].num_words_+sum_alphas)/sum_, 1., "error in gamma terms");
        sum_gamma_over_topics[i]=digamma(sum_);
    }
    
    if(gammas_ldav.size()==0) {
        cerr<<"gammas_ldav are empty in optimize_alpha"<<endl; exit(-1);
    }

    // this is only important the first time 
    // I should move this to initialization
    
    
    //-------------------------- BEGIN INITIALIZATION --------------------
    // average_ps[k] is the average over i of p[i][k]
    DD average_ps;
    average_ps.assign(num_topics_ldav_, 0.);
    DD average_square_ps;
    average_square_ps.assign(num_topics_ldav_, 0.);
    RANGE_loop(i, ps) {
        RANGE_loop(k, ps[i]) {
            average_ps[k]+=ps[i][k];
            average_square_ps[k]+=(ps[i][k])*(ps[i][k]);            
        }
    }
    RANGE_loop(k, average_ps) average_ps[k]/= gammas_ldav.size();
    RANGE_loop(k, average_ps) average_square_ps[k]/= gammas_ldav.size();
    //-------------------------- INITIALIZATION END --------------------


    // 1/gammas_ldav.size() x sum_i [ digamma(gammas_ldav[i][k]) - sum_gamma_over_topics[i]] 
    DD gamma_terms;
    gamma_terms.assign(num_topics_ldav_, 0.);
    
    RANGE_loop(i, ps) {
        RANGE_loop(k, ps[i]) {
            // gammas_ldav should not never be so small because of the prior
            if(gammas_ldav[i][k]>1e-10)
                gamma_terms[k]+= digamma(gammas_ldav[i][k]) - sum_gamma_over_topics[i];
            else {
                //gamma_terms[k]+= -1e10 -sum_gamma_over_topics[i];
                cerr<<"Very small gamma value:: "<<gammas_ldav[i][k]<<endl;
                cerr<<"This is likely a bug. Please contact me: arg.lanci@gmail.com Thanks!"<<endl;
                exit(-1);
            }
        }
    }
    
    
    RANGE_loop(k, gamma_terms) gamma_terms[k]/= gammas_ldav.size();
    
    
    /*cout<<"average_ps"<<endl;
    prints(average_ps);
    cout<<"average_square_ps"<<endl;
    prints(average_square_ps);
    cout<<"gamma_terms"<<endl;
    prints(gamma_terms);*/
    //cout<<"sum_alphas:::: "<<sum_alphas<<endl;
    
    
    double alpha_sum=   (average_ps[0]-average_square_ps[0])/
                        (  average_square_ps[0] - average_ps[0]*average_ps[0] );
    
    
    
    DD alphas;
    RANGE_loop(k, average_ps) alphas.push_back(average_ps[k]*alpha_sum);
    sum_alphas=0.;
    RANGE_loop(k, alphas) sum_alphas+=alphas[k];

    //cout<<"alpha_sum:: "<<alpha_sum<<" "<<sum_alphas<<endl;

    
    int iter=0;
    while(true) {
        double new_sum=0.;
        //cout<<"gamma terms"<<endl;
        //prints(gamma_terms);
        //cout<<"-->> "<<digamma(sum_alphas)<<endl;
        RANGE_loop(k, alphas) {
            alphas[k]= invert_digamma(digamma(sum_alphas) + gamma_terms[k] );
            new_sum+=alphas[k];
        }
        //cout<<"new_sum:: "<<new_sum<<endl;
        if(new_sum!=new_sum) {
            RANGE_loop(k, alphas) alphas[k]=1.;
            cerr<<"warning! alphas are too big!"<<endl;
            new_sum=alphas.size();
            break;
        }
        
        if(fabs(sum_alphas-new_sum)<1e-7 or iter>100) {
            break;
        }
        sum_alphas=new_sum;
        ++iter;
        
    }

    cout<<"new sum_alphas: "<<sum_alphas<<endl;
    //prints(alphas);
    alphas_ldav_=alphas;
    
    
}

void word_corpus::compute_non_sparse_gammas(deque<DD> & gammas_ldav) {
    
    // copying gammas_ldav_map_ in gammas_ldav
    // topics which are not gammas_ldav_map_ get just the prior (alphas_ldav_)
    
    gammas_ldav.clear();    
    
    DD void_dd_numtops;
    void_dd_numtops.assign(num_topics_ldav_, 0.);
    
    // initializing gammas
    RANGE_loop(doc_number, docs_) {
        gammas_ldav.push_back(void_dd_numtops);
    }
    
    RANGE_loop(doc_number, docs_) {
        IT_loop(mapid, itm, gammas_ldav_map_.at(doc_number)) {
            gammas_ldav[doc_number][itm->first]=itm->second;
        }
        RANGE_loop(k, gammas_ldav[doc_number]) {
            gammas_ldav[doc_number][k]=max(gammas_ldav[doc_number][k], alphas_ldav_[k]);
        }
    }

}


void word_corpus::optimize_alpha_sparse() {
    
    deque<DD> gammas_ldav;
    compute_non_sparse_gammas(gammas_ldav);
    optimize_alpha(gammas_ldav);
    
}



