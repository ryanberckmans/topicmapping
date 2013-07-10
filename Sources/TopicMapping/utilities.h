


void dotpr_similarity_of_connected_words(int_matrix & corpus_matrix, map<pair<int, int> , int> & cooc) {
	
    cooc.clear();
	RANGE_loop(i, corpus_matrix) {
    
        mapii hist;
        DI ws;
        RANGE_loop(j, corpus_matrix[i]) int_histogram(corpus_matrix[i][j], hist);        
        IT_loop(mapii, itm, hist) ws.push_back(itm->first);
        
        RANGE_loop(j, ws) for(UI k=0; k<j; k++) {
            
            int w1=min(ws[j], ws[k]);
            int w2=max(ws[j], ws[k]);
            
            if(w1!=w2) {
                pair<int, int> pp(w1, w2);
                if (cooc.find(pp)==cooc.end())
                    cooc[pp]=0;
                cooc[pp]+=hist[w1]*hist[w2];
            }
        }
    }
}



int argmax(mapid & topic_distr) {
    
    int argm=-1;
    double maxe=0.;
    IT_loop(mapid, itm, topic_distr) {
        if(itm->second>maxe) {
            argm=itm->first;
            maxe=itm->second;
        }
    }
    return argm;
    
}


double normalize_mapid(mapid & a) {

    double norm=0.;
    IT_loop(mapid, itm, a) norm+=itm->second;
    if(norm<1e-10) {
        a.clear();
    }
    else {
        IT_loop(mapid, itm, a) itm->second/=norm;
    }
    return norm;
}


void set_matrix_to_zero(int rows, int columns, deque<DD> & matrix) {
    
    matrix.clear();
    DD zeros;
    zeros.assign(columns, 0.);
    for(int z=0; z<rows; z++) {
        matrix.push_back(zeros);
    }
}

void print_matrix_to_file_or_screen(deque<DD> & matrix, string filename) {
    
    if (filename.size()==0)
        printm(matrix);
    else {
        ofstream pout(filename.c_str());
        printm(matrix, pout);
    }
    
}


inline double get_from_mapid(const mapid & m, int key) {
    
    /* this function returns m.at(key) if it exists
     and 0 if it does not */
    
    mapid::const_iterator pos = m.find(key);
    if (pos == m.end()) {
        return 0.;
    } else {
        return pos->second;
    }
    
}






int max_element_mapid(mapid & m) {
	
	int mx=0;
	int mc=-1;
	
	IT_loop(mapid, itm, m ) if(itm->second>mx) {
		mx=itm->second;
		mc=itm->first;
	}
	
	return mc;
    
}



void make_topic_names_consecutive(deque<mapid> & doc_topic, \
                                  map<int, mapid> & topic_word, \
                                  deque<mapii> & doc_assignments, mapid & pt) {
    
    // for each topic, I use consecutive names
    mapii old_to_new_labels;
    
    for(map<int, mapid>::iterator itm = topic_word.begin(); itm!=topic_word.end(); itm++) {
        int tp=old_to_new_labels.size();
        old_to_new_labels[itm->first]=tp;
    }

    // fixing doc_topic
    RANGE_loop(i, doc_topic) {
        mapid new_topic_distr;
        IT_loop(mapid, itm, doc_topic[i]) if(itm->second>1e-10) {
            new_topic_distr[old_to_new_labels.at(itm->first)]=itm->second;
        }
        doc_topic[i]=new_topic_distr;
    }
    
    // fixing topic_word
    map<int, mapid> new_topic_word;
    for(map<int, mapid>::iterator itm = topic_word.begin(); itm!=topic_word.end(); itm++) {
        new_topic_word.insert(make_pair(old_to_new_labels.at(itm->first), itm->second));
    }
    topic_word= new_topic_word;
    
    // fixing doc_assignments
    RANGE_loop(i, doc_assignments) {        
        IT_loop(mapii, itm, doc_assignments[i]) {
            itm->second=old_to_new_labels.at(itm->second);
        }
    }

    // fixing pt
    mapid new_pt;
    for(mapid::iterator itm = pt.begin(); itm!=pt.end(); itm++) {
        new_pt.insert(make_pair(old_to_new_labels.at(itm->first), itm->second));
    }
    pt= new_pt;
    

    
    
}


