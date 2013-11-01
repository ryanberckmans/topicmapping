
// adapting code from lda-c ( http://www.cs.princeton.edu/~blei/topicmodeling.html )

#include "lda_util.cpp"
#include "lda_model.cpp"
#include "alpha_optimization.cpp"
#include "lda_em.cpp"



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


