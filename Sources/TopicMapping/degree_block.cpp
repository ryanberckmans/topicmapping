



double get_log_fact(DD & log_factorials, int requested_k) {
    
    if (requested_k < int(log_factorials.size()))
        return log_factorials[requested_k];
    
    int k=log_factorials.size()-1;
    double log_k_fact=log_factorials[k];
        
    for(int i = k+1; i<=requested_k; i++) {
        log_k_fact+=log(i);
        log_factorials.push_back(log_k_fact);
    }
    
    /*cout<<"*** factorials ***"<<endl;
    RANGE_loop(i, log_factorials){
        cout<<i<<" "<<exp(log_factorials[i])<<endl;
    }*/
    return log_factorials[requested_k];
    
}


double get_right_cum(int & k, double lambda, DD & log_factorials) {


    double log_k_fact = get_log_fact(log_factorials, k);
    double q= exp(-lambda + k*log(lambda) - log_k_fact);
    if (q<1e-300) {
        cout<<"value too small. this should never happen. BUG! please contact me: arg.lanci@gmail.com"<<endl;
        exit(-1);
    }
    double vk=q;
    
    while(true) {    
        double vk_plus1 = vk *  lambda / (k+1);
        q+=vk_plus1;
        if (vk_plus1/q<1e-6)
            break;
        ++k;
        vk=vk_plus1;
    }
    return q;
}





int poissian_quantile(double p, double lambda, DD & log_factorials) {
	
	// find the minimum k such that 1 - exp(-lambda) sum_[0,k] lambda**i/i! < p	
    if(p>=1.)
        return 0;
    if(lambda<1e-50)
        return 1;
    
    int starting_k= int(lambda);
    double q=1.1;
    while(true) {
        q = get_right_cum(starting_k, lambda, log_factorials);
        if (q<p)
            break;
        starting_k+=max(1,int(sqrt(lambda)));
    }
    
    double log_k_fact = get_log_fact(log_factorials, starting_k);
    double vstarting_k= exp(-lambda + starting_k*log(lambda) - log_k_fact);
    
    while(true) {
        
        double vk_minus1 = vstarting_k *  starting_k / lambda;
        q+=vk_minus1;
        --starting_k;
        if (q>p)
            break;
        vstarting_k=vk_minus1;
        
    }
	return starting_k;
}




int qgeneral(int k1, int k2, double p, map<pair<int, int> , int> & significant_p_percent, const double & pfact_poisson, DD & log_factorials){
	
	pair<int, int> dds(k1,k2);
    
    
	if(significant_p_percent.find(dds)==significant_p_percent.end()) {
        significant_p_percent[dds]=poissian_quantile(p, pfact_poisson*k1*k2, log_factorials);
	}
	return significant_p_percent[dds];
	
}



class degree_block {
	
public:
	
	degree_block(){ log_factorials.push_back(0); log_factorials.push_back(0); };
	~degree_block(){};
	
	int qfive(int k1, int k2,double pvalue);
	void random_similarities_poissonian_set_data(DI & doc_sizes);
	
	
private:
	
	map<pair<int, int> , int> significant_five_percent;
	double pfact_poisson;
    DD log_factorials;
    

	
};


void degree_block::random_similarities_poissonian_set_data(DI & doc_sizes) {
	
	/* 
	 this function set a few things necessary for return_quantiles function.
	 Things set:
		
	1. doc_sizes_poisson
	2. pfact_poisson
	3. dN_poisson
	
	 */
	
	pfact_poisson=0.;
	significant_five_percent.clear();
	DI doc_sizes_poisson;
	
	int N=0;
	RANGE_loop(i, doc_sizes) {
		N+=doc_sizes[i];
		doc_sizes_poisson.push_back(doc_sizes[i]);
	}
	double dN_poisson=N;
	RANGE_loop(i, doc_sizes_poisson) {
		pfact_poisson+= doc_sizes_poisson[i]*doc_sizes_poisson[i]/ (dN_poisson*dN_poisson);
	}
    cout<<"dN_poisson "<<dN_poisson<<" pfact_poisson "<<pfact_poisson<<endl;
}


int degree_block::qfive(int k1, int k2, double pvalue) {
	
	return qgeneral(k1,k2,pvalue,significant_five_percent,pfact_poisson, log_factorials);
}











