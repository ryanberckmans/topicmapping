
// adapting code from lda-c ( http://www.cs.princeton.edu/~blei/topicmodeling.html )

#include "lda_util.cpp"
#include "lda_model.cpp"
#include "alpha_optimization.cpp"
//#include "lda_em.cpp"
#include "lda_em_sparse.cpp"



double word_corpus::lda_model(map<int, mapid> & topic_word, double alpha_init,\
                              bool skip_alpha_opt, bool infer_flag) {

    // topic_word[topic][word] is p(w|t)
    // getting all data structures ready
    initialize_lda_data(topic_word, alpha_init);
    
    // loop until convergence
    run_em_sparse(skip_alpha_opt, infer_flag);
    
    return 0.;

}


