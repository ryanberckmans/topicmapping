

#include "standard_package/standard_include.cpp"
#include "comparison_utilities.cpp"



double effective_number(DD & a) {

    double h=0;
    RANGE_loop(i, a) if(a[i]>0) {
        h-=a[i]*log(a[i]);
    }
    return exp(h);
}

pair<double, double> randomized_bestmatch(deque<DD> & pd_given_t1, \
                                          deque<DD> & pd_given_t2, DD & pt1, \
                                          int reals=20) {
    
    double raw_bm=bestmatch(pd_given_t1, pd_given_t2, pt1);
    double bm_random=0;
    for(int real=0; real<reals; real++) {
        deque<DD> shuffled_matrix;
        shuffle_colums(pd_given_t2, shuffled_matrix);
        bm_random+=bestmatch(pd_given_t1, shuffled_matrix, pt1);
    }
    bm_random/=reals;
    cout<<"raw "<<raw_bm<<" random: "<<bm_random<<endl;
    pair<double, double> rawbm_randbm(raw_bm, bm_random);
    return rawbm_randbm;
}


double compute_bm(string gammafile1, string gammafile2, string corpus_txt) {
    
    deque<DD> gamma_matrix1;
    deque<DD> gamma_matrix2;
    DD pd;
    
    
    if (get_gamma_matrix(gammafile1, gamma_matrix1)==-1 or \
        get_gamma_matrix(gammafile2, gamma_matrix2)==-1)
        return -1;
    
    get_pd(corpus_txt,  pd);
    if (pd.size()!=gamma_matrix1.size() or pd.size()!=gamma_matrix2.size() ) {
        cerr<<"document numbers do not match"<<endl;
        return -1;
    }
    
    deque<DD> pd_given_t1;
    deque<DD> pd_given_t2;
    DD pt1;
    DD pt2;
    
    bayes_theorem(pd, gamma_matrix1, pt1, pd_given_t1);
    bayes_theorem(pd, gamma_matrix2, pt2, pd_given_t2);
    
    cout<<"eff1: "<<effective_number(pt1)<<endl;
    cout<<"eff2: "<<effective_number(pt2)<<endl;
    
    pair<double, double> m1vsm2= randomized_bestmatch(pd_given_t1, pd_given_t2, pt1);
    pair<double, double> m2vsm1= randomized_bestmatch(pd_given_t2, pd_given_t1, pt2);
    double raw_average = 0.5 * (m1vsm2.first + m2vsm1.first);
    double rand_average = 0.5 * (m1vsm2.second + m2vsm1.second);
    
    if (fabs(1-rand_average)<1e-10)
        return 0;
    return (raw_average-rand_average)/(1.-rand_average);
    
}


int main(int argc, char * argv[]) {
    
    if(argc<4) {
        cout<<argv[0]<<" [theta_file1] [theta_file2] [corpus.txt]"<<endl;
        return -1;
    }
    
    srand_file();
    double bm=compute_bm(string(argv[1]), string(argv[2]),string(argv[3]));
    cout<<"bm similarity: "<<bm<<endl;
    
    return 0;
        
}
