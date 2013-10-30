
// adapting code from lda-c ( http://www.cs.princeton.edu/~blei/topicmodeling.html )

#include "lda_util.cpp"

double word_corpus::lda_inference(int doc_number, map<int, mapid> & topic_word,
                                  mapid & doc_topic, double alpha, ostream & infout2) {

    // TODO:
    // add smoothing on topic_word (1e-2 should be OK)
    // and make topic sparse.
    // doc_number does not deal with all topics.


    // inference is performed using only topics which are in doc_topic
    // this is a futher approximation, but enables sparser data structures
    // (each document is only using a limited set of topics)

    int num_topics= topic_word.size();
    
    //cout<<"num_topics:: "<<num_topics<<endl;
    
    // p(topic|doc) inferred -not normalized-
    DD var_gamma;
    var_gamma.assign(num_topics, 0.);
    // psi(var_gamma)
    double digamma_gam[num_topics];
    double oldphi[num_topics];


    // each word of doc in a (topic-probability)
    // DD[topic] is the probability
    map<int, DD> phi;
    IT_loop(mapii, itm, docs[doc_number].wn_occurences) {            
        DD void_dd;
        void_dd.assign(num_topics, 0.);
        phi[itm->first]=void_dd;
    }
    // compute posterior dirichlet
    for(int k=0; k<num_topics; k++) {
        // uniform
        var_gamma[k]=(alpha + double(docs[doc_number].num_words)/num_topics);
        //var_gamma[k]=alpha + get_from_mapid(doc_topic, k) * docs[doc_number].num_words;
        digamma_gam[k]=digamma(var_gamma[k]);
        // phi[word][topic]
        IT_loop(mapii, itm, docs[doc_number].wn_occurences) {            
            phi[itm->first][k] = 1.0/num_topics;
            //phi[itm->first][k] = digamma_gam[k] + log(max(get_from_mapid(topic_word[k], itm->first), 1e-100));
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
                // this needs to be fixed
                // 1e-100 is a very small number, just because this word cannot be found in topic k
                // ADD SMOOTING HERE

                oldphi[k]=phi[n][k];
                phi[n][k] = digamma_gam[k] + log(max(get_from_mapid(topic_word[k], n), 1e-100));
                if (k > 0)
                    phisum = log_sum(phisum, phi[n][k]);
                else
                    phisum = phi[n][k]; // note, phi is in log space
            }
            
            for (int k = 0; k < num_topics; k++) {
            
                phi[n][k] = exp(phi[n][k] - phisum);
                var_gamma[k] += itm->second*(phi[n][k] - oldphi[k]);
                digamma_gam[k] = digamma(var_gamma[k]);
            }            
        }

        likelihood = compute_likelihood(doc_number, topic_word, phi, var_gamma, alpha);
        //cout<<"likelihood "<<likelihood<<endl;
        //prints(var_gamma);
        
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
    
double word_corpus::compute_likelihood(int doc_number, map<int, mapid> & topic_word,
                                       map<int, DD> & phi, DD & var_gamma, double alpha) {
    
    
    int num_topics= topic_word.size();
    
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
            if (phi[n][k] > 0) {
        
                likelihood += itm->second * (phi[n][k]*((dig[k] - digsum) - log(phi[n][k])
                            +  log(max(get_from_mapid(topic_word[k], n), 1e-100)) )    );
                }
            }
    }
    return(likelihood);
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
    
    cout<<"running lda model"<<endl;
    //RANGE_loop(i, doc_topic) prints(doc_topic[i]);
    
    /*
    cout<<"topics:: "<<topic_word.size()<<endl;
    for (map<int, mapid>::iterator topic_itm= topic_word.begin(); 
         topic_itm!=topic_word.end(); topic_itm++) {
        cout<<"--------- "<<topic_itm->first<<endl;
        //prints(topic_itm->second);
    }*/
    
    
    //gibbs_sampling(doc_assignments);
    //exit(-1);
    
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
        infout<<lda_inference(i, topic_word, fake, alpha, infout2)<<endl;
        //exit(-1);
    }
    
    infout.close();
    infout2.close();
        
    return 0.;

}
