

void check_folder(string folder_name) {
    
    cout<<"making directory:: "<<folder_name<<endl;
    string make_folder="nohup mkdir "+folder_name+" > mkdir.tmp";
    int sy_err=system(make_folder.c_str());
    cout<<"mkdir "<<folder_name<<" returned:: "<<sy_err<<endl;
    
    ifstream gin("mkdir.tmp");
    deque<string> lines;
    string line;
    while(getline(gin, line)) {
        lines.push_back(line);
    }
    
    if(lines.size()>0) {
        cerr<<"folder \""<<folder_name<<"\" already exists!"<<endl;
        cerr<<"please remove it or run this program with a different out-directory"<<endl;
        exit(-1); 
    }
    
    sy_err=system("rm mkdir.tmp");
    
    
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
        pout.close();
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



void assert_consecutive(DI & a) {
    // this function is checking that a is range(a.size())
    RANGE_loop(i, a) assert_ints(i, a[i], "error in assert_consecutive");
}



void update_map_wordnum(mapsi & word_number, mapis & number_word, string s) {
	
	if(word_number.find(s)==word_number.end()) {
		int num=word_number.size();
		word_number[s]=num;
		number_word[num]=s;
	}
}



double compute_eff_num_topics(const mapid & pt) {
    double h=0.;
    IT_const_loop(mapid, itm, pt) h+=itm->second*log2(itm->second);
    return pow(2,-1*h);
}



/*
 
 for in-out, we only use topic_word format
 first number if the topic
 everything with a # is ignored
 
 How it looks:
 #topic No. 3: star sky blue
 3 1 0.74647 45 0.73636 78 10.536
 
 NB1. There can be multiple lines associated with the same topic
 NB2. there is no need for normalization
 
 topic_word[wn][topic]
 beta_kind[topic][wn] and it is in  lig space
 
 */


void from_beta_to_topic_word(deque<mapid> & beta_kind, map<int, mapid> & topic_word) {
    
    // take exp and inserting in topic_word
    
    topic_word.clear();
    RANGE_loop(wn, beta_kind) {
        
        mapid & betas_ldav_wn = beta_kind.at(wn);
        IT_loop(mapid, topic_pr, betas_ldav_wn) { 
            int topic=topic_pr->first;
            // inserting topic
            if(topic_word.count(topic)==0) {
                mapid void_mapid;
                topic_word[topic]=void_mapid;
            }
            // inserting pr for this wn
            double pr=exp(topic_pr->second);
            if(pr>0)
                topic_word[topic][wn]=pr;
        }
    }
}


void from_class_to_topic_word(deque<mapid> & class_kind, map<int, mapid> & topic_word, bool normalize) {
    
    // same function as above but without exp
    
    topic_word.clear();
    RANGE_loop(wn, class_kind) {        
        mapid & class_wn = class_kind.at(wn);
        IT_loop(mapid, topic_pr, class_wn) { 
            int topic=topic_pr->first;
            // inserting topic
            if(topic_word.count(topic)==0) {
                mapid void_mapid;
                topic_word[topic]=void_mapid;
            }
            // inserting pr for this wn
            double pr=topic_pr->second;
            if(pr>0)
                topic_word[topic][wn]=pr;
        }
    }
    
    // normalize it
    if(normalize) for (map<int, mapid>::iterator topic_itm= topic_word.begin(); 
         topic_itm!=topic_word.end(); topic_itm++) {
        normalize_mapid(topic_itm->second);
    }
}


void from_topic_word_to_beta(map<int, mapid> & topic_word, deque<mapid> & beta_kind, double sparse_limit) {
    
    // converting topic_word in beta_kind
    // the function assumes that beta_kind
    // was already initialized
    // and words in topic_word are present in beta_kind    
    
    RANGE_loop(wn, beta_kind) beta_kind[wn].clear();
    
    for (map<int, mapid>::iterator topic_itm= topic_word.begin(); 
         topic_itm!=topic_word.end(); topic_itm++) {
        IT_loop(mapid, itm2, topic_itm->second) { 
            int wn=itm2->first;
            // words which are below this threshold should never 
            // be relevant in this topic
            if(itm2->second>sparse_limit)
                beta_kind.at(wn)[topic_itm->first] = log(itm2->second);
        }
    }
}



void print_topic_sparse_format_complete(map<int, mapid> & topic_word, string outfile,\
                                        string outfile_short, const deque<string> & word_strings,
                                        DD & ptopic, int num_words_shown) {
    
    // writes print_topic_sparse_format
    
    if(outfile_short.size()>0 and word_strings.size()>0 and ptopic.size()>0)
        cout<<"printing results. summary can be found in: "<<outfile_short<<endl;
    ofstream pout1;
    ofstream pout2;
    
    pout1.open(outfile.c_str());
    // the file is opened only if outfile_short and word_strings are provided
    if(outfile_short.size()>0 and word_strings.size()>0 and ptopic.size()>0) {
        pout2.open(outfile_short.c_str());
    }
    for (map<int, mapid>::iterator topic_itm= topic_word.begin(); 
         topic_itm!=topic_word.end(); topic_itm++) {
        // printing topic distribution over words
        mapid & topic_distr = topic_itm->second;
        deque<pair<double, int> > pr_word;
        IT_loop(mapid, itm2, topic_distr) {
            pr_word.push_back(make_pair(-itm2->second, itm2->first));
        }
        sort(pr_word.begin(), pr_word.end());
        if(pout2.is_open()) {
            pout2<<"topic: "<<topic_itm->first<<" #words: "<<pr_word.size()<<" pt: "<<ptopic[topic_itm->first]<<endl;
        }
        pout1<<topic_itm->first<<" ";
        RANGE_loop(i, pr_word) {
            pout1<<pr_word[i].second<<" "<<-pr_word[i].first<<" ";
            if (int(i)<num_words_shown and pout2.is_open()) {
                pout2<<word_strings.at(pr_word[i].second)<<" ";
            } 
        }
        pout1<<endl;
        if(pout2.is_open()) pout2<<endl;
    }
    pout1.close();
    pout2.close();
    
}



void print_topic_sparse_format_short(map<int, mapid> & topic_word, string outfile) {
    
    string outfile_short;
    deque<string> word_strings;
    DD ptopic;
    print_topic_sparse_format_complete(topic_word, outfile, outfile_short, word_strings, ptopic, 0);
}



void draw_a_random_model(map<int, mapid> & topic_word, int K, int num_terms) {
    topic_word.clear();
    for(int k=0; k<K; k++) {
        mapid new_mapid;
        for(int wn=0; wn<num_terms; wn++) {
            new_mapid[wn]=1./num_terms+ran4();
        }
        normalize_mapid(new_mapid);
        topic_word[k]=new_mapid;
    }
}    


void read_topic_model_from_file(string infile, map<int, mapid> & topic_word) {
    
    
    // !!!this function is a bit silent!!!
    // !!!add some format checks!!!
    
    // getting topic_word from file similar to "topic_words.txt"
    // aggregates lines which start with the same topic
    // and normalizes the mapid
    
    // topic_word[topic][wn]
    topic_word.clear();
    
    ifstream gin(infile.c_str());
    string gins;    
    
    while(getline(gin, gins)) if(gins.size()>0 and gins[0]!='#') {
        deque<string> vs;
        separate_strings(gins, vs);
        
        int wn=-1;
        int topic=-1;
        double pr=0.;
        RANGE_loop(i, vs) {
            if(i==0) {
                topic=atoi(vs[i].c_str());
            } else if(i%2==1) {
                wn=atoi(vs[i].c_str());
            } else {
                pr=atof(vs[i].c_str());
                // inserting
                if(topic_word.count(topic)==0) {
                    mapid new_mapid;
                    topic_word[topic]=new_mapid;
                }
                int_histogram(wn, topic_word[topic], pr);
            }
        }
    }
    gin.close();
    // asserting topics are consecutive
    DI all_topics;    
    for(map<int, mapid>::iterator itm = topic_word.begin(); itm!=topic_word.end(); itm++) {
        all_topics.push_back(itm->first);
    }
    assert_consecutive(all_topics);

    
    // normalize it
    for (map<int, mapid>::iterator topic_itm= topic_word.begin(); 
         topic_itm!=topic_word.end(); topic_itm++) {
        normalize_mapid(topic_itm->second);
    }
    
}



void get_ptopic_distr(DD & ptopic, const deque<DD> & gammas) {
    
    // ptopic(t) ~ sum_doc gammas[doc][t]
    
    int size=-1;
    ptopic.clear();
    RANGE_loop(doc, gammas) {
        const DD & vs = gammas[doc];
        
        general_assert(size==-1 or int(vs.size())==size, \
                       "ERROR: gammas should all have same length");
        if(size==-1) {
            RANGE_loop(t, vs) ptopic.push_back(vs[t]);
        } else {
            RANGE_loop(t, vs) ptopic[t]+=vs[t];
        }
        if(size==-1)
            size=vs.size();
    }
    normalize_one(ptopic);
}

