
#define log_table_pr 1e-5

class log_fact_table {

public:
    log_fact_table() {};
    ~log_fact_table(){};
    
    void _set_(int);
    double cum_binomial_right(int x, int N, double prob);
    double cum_binomial_left(int x, int N, double prob);

private:
    
    inline double binom(int x, int N, double p) { return exp(  log_choose(N, x) + x*log(p) + (N-x) * log(1-p)    ); };
    inline double log_choose(int tm, int degree_node) {  return lnf.at(tm) - lnf.at(tm - degree_node) - lnf.at(degree_node);  };

    DD lnf;
};


void log_fact_table::_set_(int size) {

	//cout<<"allocating "<<size<<" factorials..."<<endl;
	lnf.clear();
		
	double f=0;
	lnf.push_back(0);
	for(int i=1; i<=size; i++) {		
		f+=log(i);
		lnf.push_back(f);
	}
	
	//cout<<"done"<<endl;
}




double log_fact_table::cum_binomial_right(int x, int N, double prob) {
    
	// this is bigger  or equal p(x >= kin_node)   *** EQUAL ***
	
	//cout<<"x "<<x<<" N "<<N <<"  prob "<<prob<<endl;
	
	
	if(x<=0) 
		return 1;
	
	if(x>N)
		return 0;
	
	
	if(prob-1> - 1e-11)
		return 1;
	
	if(x<N*prob)
		return 1-cum_binomial_left(x, N, prob);
	
	
	double pzero= binom(x, N, prob);
	
	
	if(pzero<=1e-40)
		return 0;
	
	double z_zero= 1.;
	double sum= z_zero;
	
	
	while(true) {
        
		z_zero *=  prob * double(N-x) / ((x+1)*(1-prob));
		x++;
		//cout<<"zzero sum "<<z_zero<<" "<<sum<<" "<<endl;
		
		if(z_zero< log_table_pr * sum)
			break;
		
		sum+=z_zero;	
	}
	
	
	return pzero * sum;
    
    
    
    
}





double log_fact_table::cum_binomial_left(int x, int N, double prob) {
    
	// this is less strictly p(x < kin_node)   *** NOT EQUAL ***
	
	
	if(x<=0)
		return 0;
	
	if(x>N)
		return 1;
	
	
	if(prob<1e-11)
		return 1;
	
	if(x>N*prob)
		return 1-cum_binomial_right(x, N, prob);
	
	--x;
	double pzero= binom(x, N, prob);
	
	
	if(pzero<=1e-40)
		return 0;
	
	
	double z_zero= 1.;
	double sum= z_zero;
	
	while(true) {
        
		
		
		--x;
		z_zero *=  (1-prob) * double(x+1) / ((N-x) *prob);
		
		//cout<<"zzero sum "<<z_zero<<" "<<sum<<" "<<(ga + x)<<endl;
		
		if(z_zero< log_table_pr * sum)
			break;
		
		sum+=z_zero;	
	}
	
	
	return pzero * sum;
}


