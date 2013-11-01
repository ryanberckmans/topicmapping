


double invert_digamma(double y) {
    
    // returns x suxh that digamma(x)=y
        
    double xold=0.;
    if(y>-2.22){
        xold=exp(y)+0.5;
    } else {
        xold= -1./(y -digamma(1));
    }
    
    double xnew;
    for(int iter=0; iter<5; iter++) {
        xnew= xold - (digamma(xold)-y)/(trigamma(xold));
        if(fabs(xnew-xold)<1e-6) {
            return xnew;
        }
    }
    return xnew;
}



void word_corpus::optimize_alpha() {
    
    
    
    // ps[i][k] is the probability of topic k in sample i
    deque<DD> ps;
    
    DD sum_gamma_over_topics;
    // sum_gamma_over_topics[i] = digamma(sum_k) 
    sum_gamma_over_topics.assign(gammas_ldav_.size(), 0.);
    // assert that sum_gamma_over_topics is simply
    // the digamma(number of words + sum of priors)
    
    RANGE_loop(i, gammas_ldav_) {
        ps.push_back(gammas_ldav_[i]);
        double sum_=normalize_one(ps[i]);
        sum_gamma_over_topics[i]=digamma(sum_);
    }
    
    if(gammas_ldav_.size()==0) {
        cerr<<"gammas_ldav_ are empty in optimize_alpha"<<endl; exit(-1);
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
    
    RANGE_loop(k, average_ps) average_ps[k]/= gammas_ldav_.size();
    RANGE_loop(k, average_ps) average_square_ps[k]/= gammas_ldav_.size();
    //-------------------------- INITIALIZATION END --------------------



    // 1/gammas_ldav_.size() x sum_i [ digamma(gammas_ldav_[i][k]) - sum_gamma_over_topics[i]] 
    DD gamma_terms;
    gamma_terms.assign(num_topics_ldav_, 0.);

    
    RANGE_loop(i, ps) {
        RANGE_loop(k, ps[i]) {
            // gammas_ldav_ should not never be so small because of the prior
            if(gammas_ldav_[i][k]>1e-10)
                gamma_terms[k]+= digamma(gammas_ldav_[i][k]) - sum_gamma_over_topics[i];
            else {
                gamma_terms[k]+= -1e10 -sum_gamma_over_topics[i];
                cerr<<"very small gamma value:: "<<gammas_ldav_[i][k]<<endl;
            }

        }
    }
    
    
    RANGE_loop(k, average_ps) gamma_terms[k]/= gammas_ldav_.size();
    
    cout<<"average_ps"<<endl;
    prints(average_ps);
    cout<<"average_square_ps"<<endl;
    prints(average_square_ps);
    cout<<"gamma_terms"<<endl;
    prints(gamma_terms);

    

    
    double alpha_sum=   (average_ps[0]-average_square_ps[0])/
                        (  average_square_ps[0] - average_ps[0]*average_ps[0] );
    
    DD alphas;
    RANGE_loop(k, average_ps) alphas.push_back(average_ps[k]*alpha_sum);
    
    cout<<"alphas init"<<endl;
    prints(alphas);
    
    
    double sum_alphas=0.;
    RANGE_loop(k, alphas) sum_alphas+=alphas[k];
    
    
    for(int iter=0; iter<20; iter++) {
        
        double new_sum=0.;
        RANGE_loop(k, alphas) {
            //alphas[k]= invert_digamma(digamma(sum_alphas) + average_log_ps[k] );
            alphas[k]= invert_digamma(digamma(sum_alphas) + gamma_terms[k] );
            new_sum+=alphas[k];
        }
        cout<<"diff:: "<<sum_alphas-new_sum<<" "<<new_sum<<" "<<sum_alphas<<" "<<digamma(sum_alphas)<<endl;
        sum_alphas=new_sum;

    }

    // insert error check
    cout<<"alphas final"<<endl;
    prints(alphas);
    alphas_ldav_=alphas;
    
    
}



