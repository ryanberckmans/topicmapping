

#include "../standard_package/standard_include.cpp"
#include "utilities.cpp"




void compute_word_strings(string infile, DI & word_occurrences, deque<string> & word_strings) {
	
    word_strings.clear();
    word_occurrences.clear();
    map<int, pair<string, int> > wn_word_occ_map;
    string word_str;
    int wn;
    int occ;
    ifstream gin_wn(infile.c_str());
    while(gin_wn>>word_str) {
        gin_wn>>wn; gin_wn>>occ;
        wn_word_occ_map[wn]=make_pair(word_str, occ);
    }
    int counter=0;
    mapsi word_wn_all;
    for(map<int, pair<string, int> >::iterator itm= wn_word_occ_map.begin(); \
        itm!=wn_word_occ_map.end(); itm++) {
        assert_ints(itm->first, counter, "error in word_wn_count.txt. word_ids are not consecutive!");
        // itm->second.first is a string
        word_strings.push_back(itm->second.first);
        // itm->second.second is the occs
        word_occurrences.push_back(itm->second.second);
        ++counter;
        word_wn_all[itm->second.first]=itm->first;
    }
    gin_wn.close();

}


void compute_p_topic(string infile, DD & ptopic) {
    // ptopic(t) ~ sum_doc gammas[doc][t]
    
    ifstream gin(infile.c_str());
    general_assert(gin.is_open(), "gamma_file not found!");
    ptopic.clear();
    int size=-1;
    string gins;

    while(getline(gin, gins)) {
        DD vs;
        cast_string_to_doubles(gins, vs);
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
    gin.close();
    normalize_one(ptopic);
}




void print_topic_per_words(deque<mapid> & beta_kind, deque<string> & word_strings, \
                           DI & word_occurrences, string outfile, DD & ptopic, double threshold) {

    ofstream pout(outfile.c_str());
    deqii occ_words;
    RANGE_loop(wn, word_occurrences) {
        occ_words.push_back(make_pair(-word_occurrences[wn], wn));
    }
    sort(occ_words.begin(), occ_words.end());
    pout<<"#word topic:p(topic|word)"<<endl;
    pout<<"#threshold for printing topics:: "<<threshold<<endl;
    
    RANGE_loop(occ_wn, occ_words) {
        // wn
        int wn=occ_words[occ_wn].second;
        mapid & betas_ldav_wn = beta_kind.at(wn);
        deque<pair<double, int> > pr_topic_given_wn;
        double norm=0.;
        IT_loop(mapid, topic_pr, betas_ldav_wn) {
            // inserting pr for this wn
            double ptopic_givenw = exp(topic_pr->second)*ptopic[topic_pr->first];
            pr_topic_given_wn.push_back(make_pair(-ptopic_givenw, topic_pr->first));
            norm+=ptopic_givenw;
        }
        sort(pr_topic_given_wn.begin(), pr_topic_given_wn.end());
        pout<<word_strings[wn]<<" ";
        RANGE_loop(i, pr_topic_given_wn) if(-pr_topic_given_wn[i].first/norm>threshold) {
            pout<<pr_topic_given_wn[i].second<<":"<<-pr_topic_given_wn[i].first/norm<<" ";
        }
        pout<<endl;
    }


}



int main(int argc, char * argv[]) {
    
    if(argc<6) {
        cout<<argv[0]<<" [betas_sparse] [word_wn_count] [gammas] [num_words_shown] [topic_threshold]"<<endl;
        cout<<"provide general information about the model"<<endl;
        return -1;
    }
    
    
    // ==== getting word_strings =============
    DI word_occurrences;           // word_occurrences[wn] = occurences
    deque<string> word_strings;          // word_strings[wn] = "word_in_text"

    string word_wn_infile(argv[2]);    
    compute_word_strings(word_wn_infile, word_occurrences, word_strings);

    cout<<"#words: "<<word_strings.size()<<endl;
    deque<mapid> beta_kind;
    RANGE_loop(wn, word_strings) {
        mapid void_mapid;
        beta_kind.push_back(void_mapid);
    }
    
    // ==== getting topic_word =============
    cout<<"reading model"<<endl;
    map<int, mapid> topic_word;
    string betas_sparse_infile(argv[1]);    
    read_topic_model_from_file(betas_sparse_infile, topic_word);
    
    // ============= betas =========
    cout<<"computing betas"<<endl;
    from_topic_word_to_beta(topic_word, beta_kind, 0.);

    
    //============ ptopic[topic] is p(topic)
    cout<<"computing p(topic)"<<endl;
    DD ptopic;
    string gammas_infile(argv[3]);
    compute_p_topic(gammas_infile, ptopic);
    
    // ==== printing =============
    cout<<"recording"<<endl;
    string num_words_shown(argv[4]);
    string threshold_shown(argv[5]);
    cout<<"threshold_shown "<<threshold_shown<<endl;
    print_topic_sparse_format_complete(topic_word, "betas_sparse.txt",\
                                       "summary.txt", word_strings, ptopic, \
                                       cast_int(cast_string_to_double(num_words_shown)));

    print_topic_per_words(beta_kind, word_strings, word_occurrences,\
                          "word_summary.txt", ptopic, cast_string_to_double(threshold_shown));
    
    
    return 0;
}
