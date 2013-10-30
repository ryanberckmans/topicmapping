
// adapting code from lda-c ( http://www.cs.princeton.edu/~blei/topicmodeling.html )

#include "lda_util.cpp"

double word_corpus::lda_inference(int doc_number, double alpha, ostream & infout2) {

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
        var_gamma[k]=(alpha + double(docs[doc_number].num_words)/num_topics);
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

        likelihood = compute_likelihood(doc_number, var_gamma, alpha);
        
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
    
    for (int k = 0; k < num_topics; k++)
        infout2<<var_gamma[k]<<" ";
    infout2<<endl;
    
    return(likelihood);
}



/*
 * compute likelihood bound
 *
 */
    
double word_corpus::compute_likelihood(int doc_number, DD & var_gamma, double alpha) {
    
    
    int num_topics= betas_ldav[0].size();
    
    double dig[num_topics];
    double var_gamma_sum=0.;

    for (int k = 0; k < num_topics; k++) {
        dig[k] = digamma(var_gamma[k]);
        var_gamma_sum += var_gamma[k];
    }
    double digsum = digamma(var_gamma_sum);
    
    double likelihood = lgamma(alpha * num_topics) 
                        - num_topics * lgamma(alpha) 
                        - lgamma(var_gamma_sum);
    
    for (int k = 0; k < num_topics; k++) {
        likelihood +=
                        (alpha - 1)*(dig[k] - digsum) 
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


void assert_consecutive(DI & a) {
    // this function is checking that a is range(a.size())
    RANGE_loop(i, a) assert_ints(i, a[i], "error in assert_consecutive");
}



double word_corpus::lda_model(deque<mapid> & doc_topic, 
                               map<int, mapid> & topic_word, 
                               deque<mapii> & doc_assignments) {
    
    
    
    // I should make sure that topics are in order (meaning sorted)
    // and I should rename the words so that they are also sorted
    
    // doc_topic[doc] is p(t|doc)
    // topic_word[topic][word] is p(w|t)
    // doc_assignments[doc][wn] is the topic to which the word has been assigned

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
    
    
    //lda data structures
    phis_ldav.clear();
    betas_ldav.clear();
    alphas_ldav.clear();

    
    DI all_words;
    IT_loop(mapii, itm, wn_occurences_global) all_words.push_back(itm->first);
    assert_consecutive(all_words);
    DI all_topics;
    for (map<int, mapid>::iterator topic_itm= topic_word.begin(); 
         topic_itm!=topic_word.end(); topic_itm++) {
        all_topics.push_back(topic_itm->first);
    }
    assert_consecutive(all_topics);
    
    int num_topics=topic_word.size();
    DD void_dd;
    void_dd.assign(num_topics, 0.);
    RANGE_loop(i, all_words) {
        betas_ldav.push_back(void_dd);
        // this could be made more efficient
        phis_ldav.push_back(void_dd);
    }
    alphas_ldav=void_dd;
    
    
    // copying topic_word in betas_ldav
    for (map<int, mapid>::iterator topic_itm= topic_word.begin(); 
         topic_itm!=topic_word.end(); topic_itm++) {
        
        int word_wn=0;
        IT_loop(mapid, itm2, topic_itm->second) { 
            assert_ints(itm2->first, word_wn);
            betas_ldav[word_wn][topic_itm->first]=log(itm2->second);
            ++word_wn;
        }
    }
    
    
    
    cout<<"running lda model"<<endl;

    
    double alpha = 0.0065538836;
    
    ofstream infout("inf.txt");
    ofstream infout2("inf2.txt");
    mapid fake;
    //cout<<lda_inference(732, topic_word, fake, alpha, infout2)<<endl;

    //exit(-1);
    
    RANGE_loop(i, docs) {
        if (i%100==0)
            cout<<"running lda-inference for doc "<<i<<endl;
        mapid fake;
        infout<<lda_inference(i, alpha, infout2)<<endl;
        //exit(-1);
    }
    
    infout.close();
    infout2.close();
        
    return 0.;

}
