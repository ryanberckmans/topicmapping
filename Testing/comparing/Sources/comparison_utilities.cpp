

int shuffle_and_set(DI & indices, int dim) {		

	indices.clear();
    multimap <double, int> a;
	for (int i=0; i<dim; i++)
		a.insert(make_pair(ran4(), i));	
	for (multimap<double, int>::iterator it=a.begin(); it!=a.end(); it++)
		indices.push_back(it->second);
	
	return 0;
    
}


void shuffle_colums(const deque<DD> & m1, deque<DD> & m2) {

    DI indices;
    shuffle_and_set(indices, m1[0].size());
    RANGE_loop(i, m1) {
        DD d;
        RANGE_loop(j, indices) {
            d.push_back(m1[i][indices[j]]);
        }
        m2.push_back(d);
    }
    
}


double find_best_match(DD & d, deque<DD> & matrix) {

    double bestm=0.;
    
    RANGE_loop(i, matrix) {

        double similarity=1.-diff_norm_one(d, matrix[i])*0.5;        
        if (similarity>bestm) {
            bestm=similarity;
        }
    }

    return bestm;
}


double bestmatch(deque<DD> & pd_given_t1, deque<DD> & pd_given_t2, DD & pt1) {
    
    double global_bestm=0.;
    
    RANGE_loop(i, pd_given_t1) {
        
        double bestm = find_best_match(pd_given_t1[i], pd_given_t2);
        global_bestm+=pt1[i] * bestm;
    }
    
    return global_bestm;
}



void set_matrix_to_zero(int rows, int columns, deque<DD> & matrix) {
    
    matrix.clear();
    DD zeros;
    zeros.assign(columns, 0.);
    for(int z=0; z<rows; z++) {
        matrix.push_back(zeros);
    }
}



int bayes_theorem(const DD & pa, const deque<DD> & pb_given_a, \
                  DD & pb, deque<DD> & pa_given_b) {
    
    // this function should be self-explanatory
    
    if (pa.size()==0 or pb_given_a.size()==0) {
        return -1;
    }
    
    pb.clear();
    pa_given_b.clear();
    
    int na= pa.size();
    int nb= pb_given_a[0].size();
    
    pb.assign(nb, 0.);
    set_matrix_to_zero(nb, na, pa_given_b);
    
    RANGE_loop(a, pb_given_a) RANGE_loop(b, pb_given_a[a]) {
        pa_given_b[b][a]=pa[a]*pb_given_a[a][b];
    }
    
    RANGE_loop(b, pa_given_b) {
        pb[b]=norm_one(pa_given_b[b]);
        normalize_one(pa_given_b[b]);
    }
    
    
    
    return 0;
    
}




int get_gamma_matrix(string filename, deque<DD> & matrix) {
    
    ifstream gin(filename.c_str());
    
    string s;
    matrix.clear();
    
    while(getline(gin, s)) {
        DD d;
        cast_string_to_doubles(s, d);
        if (norm_one(d)<1e-10) {
            RANGE_loop(i, d) d[i]=1.;
        }
        normalize_one(d);
        matrix.push_back(d);
    }
    
    for(int i=0; i<int(matrix.size())-1; i++) {
        if (matrix[i].size()!=matrix[i+1].size()) {
            cerr<<"something is wrong in "<<filename<<endl;
            return -1;
        }
    }
    
    return 0;
    
}




void get_pd(string corpus_file_txt, DD & pd) {

    pd.clear();
    
    ifstream gin(corpus_file_txt.c_str());    
    string s;
    while(getline(gin, s)) {
        
        stringstream ss(s);
        deque<string> tokens;
        string buf;
        while (ss >> buf)
            tokens.push_back(buf);
        
        if (tokens.size()>0)
            pd.push_back(tokens.size());
    }
    normalize_one(pd);

}

