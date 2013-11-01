
// adapting code from lda-c ( http://www.cs.princeton.edu/~blei/topicmodeling.html )

#include "lda_util.cpp"
#include "lda_model.cpp"
#include "alpha_optimization.cpp"

double word_corpus::lda_inference(int doc_number) {
    
    // remove infout2
    
    // assert here
    int num_topics= betas_ldav[0].size();    
    //cout<<"num_topics:: "<<num_topics<<endl;
    
    DD var_gamma;
    var_gamma.assign(num_topics, 0.);
    // psi(var_gamma)
    double digamma_gam[num_topics];
    double oldphi[num_topics];


    // compute posterior dirichlet
    for(int k=0; k<num_topics; k++) {
        // uniform
        var_gamma[k]=(alphas_ldav[k] + double(docs[doc_number].num_words)/num_topics);
        digamma_gam[k]=digamma(var_gamma[k]);
        // phi[word][topic]
        IT_loop(mapii, itm, docs[doc_number].wn_occurences) {            
            phis_ldav[itm->first][k] = 1.0/num_topics;
        }
    }
    
    double converged=1;
    double likelihood=0;
    double likelihood_old=0;
    int var_iter=0;
    
    while(converged>1e-5 and var_iter<1000) {
    
        var_iter++;
        
        // looping over words in doc
        IT_loop(mapii, itm, docs[doc_number].wn_occurences) {
        
            // all this refers to this particular word
            int n= itm->first;
            double phisum;

            for (int k = 0; k < num_topics; k++) {

                oldphi[k]=phis_ldav[n][k];
                phis_ldav[n][k] = digamma_gam[k] + betas_ldav[n][k];
                if (k > 0)
                    phisum = log_sum(phisum, phis_ldav[n][k]);
                else
                    phisum = phis_ldav[n][k]; // note, phi is in log space
            }
            
            for (int k = 0; k < num_topics; k++) {
            
                phis_ldav[n][k] = exp(phis_ldav[n][k] - phisum);
                var_gamma[k] += itm->second*(phis_ldav[n][k] - oldphi[k]);
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

        
        /*cout<<"iter "<<var_iter<<endl;
        prints(var_gamma);
        cout<<"likelihood:: "<<likelihood<<" "<<converged<<endl;*/        
    }

    // updates for M steps
    gammas_ldav.push_back(var_gamma);

    // updating stats for betas_ldav
    IT_loop(mapii, itm, docs[doc_number].wn_occurences) for (int k = 0; k < num_topics; k++) {
        class_word_ldav[k][itm->first] += itm->second * phis_ldav[itm->first][k];
        class_total_ldav[k] += itm->second * phis_ldav[itm->first][k];
    }

    
    return(likelihood);
}


/*
 * compute likelihood bound
 *
 */

    
double word_corpus::compute_likelihood(int doc_number, DD & var_gamma) {
    
    // assert number of topics is consistent (remove this later)
    assert_ints(var_gamma.size(), alphas_ldav.size());
    assert_ints(var_gamma.size(), betas_ldav[0].size());
    
    int num_topics= betas_ldav[0].size();
    
    double dig[num_topics];
    double var_gamma_sum=0.;
    double sum_l_gamma=0.;
    double sum_alphas=0.;
    for (int k = 0; k < num_topics; k++) {
        dig[k] = digamma(var_gamma[k]);
        var_gamma_sum += var_gamma[k];
        sum_l_gamma += lgamma(alphas_ldav[k]);
        sum_alphas += alphas_ldav[k];
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
                        (alphas_ldav[k] - 1)*(dig[k] - digsum) 
                        + lgamma(var_gamma[k])
                        - (var_gamma[k] - 1)*(dig[k] - digsum);
        
        
        IT_loop(mapii, itm, docs[doc_number].wn_occurences) {
        
            int n=itm->first;
            if (phis_ldav[n][k] > 0) {
        
                likelihood += itm->second * (phis_ldav[n][k]*((dig[k] - digsum) - log(phis_ldav[n][k])
                            +  betas_ldav[n][k] )    );
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
        gammas_ldav.clear();
        double likelihood_all=0.;
        RANGE_loop(doc_number, docs) {
            if (doc_number%100==0)
                cout<<"running lda-inference for doc "<<doc_number<<endl;
            double likelihood_doc = lda_inference(doc_number);
            likelihood_all+=likelihood_doc;
            infout<<likelihood_doc<<endl;
        }
        infout.close();
        
        // M step
        // optimizing betas
        for(int k=0; k < num_topics_ldav; k++) {
            RANGE_loop(wn, word_occurrences) {
                if(class_word_ldav[k][wn] > 0)
                    betas_ldav[wn][k] = log(class_word_ldav[k][wn]) - log(class_total_ldav[k]);
                else
                    betas_ldav[wn][k] = -100;
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
    printm(gammas_ldav, pout_final);

    return 0.;
}




double word_corpus::lda_model(deque<mapid> & doc_topic, 
                               map<int, mapid> & topic_word) {

    // doc_topic[doc] is p(t|doc)
    // topic_word[topic][word] is p(w|t)
    // getting all data structures ready
    initialize_lda_data(doc_topic, topic_word);
    
    // loop until convergence
    run_em();
    
    return 0.;

}


