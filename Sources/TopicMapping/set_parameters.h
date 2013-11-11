

class parameters {
	
public:
	
	parameters(){};
	~parameters(){};
	
	/* if required==true, the program will not run if the parameter is not set by the user */
	void set_double(string flag, double default_p, bool required, string statement);	
	void set_int(string flag, int default_p, bool required, string statement);	
	void set_string(string flag, string default_p, bool required, string statement);	
	void set_bool(string flag, bool default_p, bool required, string statement);	
	void printing(ostream & pout);
	void set_from_argv(int argc, char * argv[]);
	
	/* maps the flag to the parameters */
	map<string, double> double_ps;	
	map<string, int> int_ps;
	map<string, string> string_ps;
	map<string, bool> bool_ps;
	/* maps the flag to the parameters */

	
private:
	
	void set_the_rest(string flag, bool required, string statement);
	void required_flags_are_set(char * argv[]);
	void error_statement(char * argv[]);
	
	deque<pair<string, string> > flag_statement;
	set<string> required_flags;
	
};


void parameters::set_the_rest(string flag, bool required, string statement) {
	
	flag_statement.push_back(make_pair(flag,statement));
	if(required)
		required_flags.insert(flag);
	
	

}

void parameters::set_double(string flag, double default_p, bool required, string statement) {

	double_ps[flag]=default_p;
	set_the_rest(flag, required, statement);
}


void parameters::set_int(string flag, int default_p, bool required, string statement) {
	
	int_ps[flag]=default_p;
	set_the_rest(flag, required, statement);
}



void parameters::set_string(string flag, string default_p, bool required, string statement) {
	
	string_ps[flag]=default_p;
	set_the_rest(flag, required, statement);
}



void parameters::set_bool(string flag, bool default_p, bool required, string statement) {
	
	bool_ps[flag]=default_p;
	set_the_rest(flag, required, statement);
}



void parameters::error_statement(char * argv[]) {

	// here there could be a front statement
    cerr<<"This program is to find topics in a set of documents using network ";
    cerr<<"clustering (Infomap) and the LDA model (variational optimization)."<<endl<<endl;
    
	for(deque<pair<string, string> >::iterator itm=flag_statement.begin(); itm!=flag_statement.end(); itm++) {
		cerr<<itm->first<<" "<<itm->second<<endl;
	}
    
    cerr<<endl<<"Basic example: "<<endl;
    cerr<<argv[0]<<" -f quantum-and-granular-large-stemmed"<<endl;
    cerr<<"Please look at ReadMe.pdf for more info."<<endl;
	exit(-1);
}


void parameters::set_from_argv(int argc, char * argv[]) {
	
	int i=1;
	
	if(i==argc)
		return error_statement(argv);
	
	while(i<argc) {
		
		string temp(argv[i]);
		
		
		if(required_flags.find(temp)!=required_flags.end())
			required_flags.erase(temp);

		/***************** find what temp is *****************/
		if(double_ps.find(temp)!=double_ps.end()) {
			
			++i;
			if(i==argc)
				error_statement(argv);
			string temp2(argv[i]);
			double_ps[temp]=cast_string_to_double(temp2);
			
		} else if(int_ps.find(temp)!=int_ps.end()) {
			
			++i;
			if(i==argc)
				error_statement(argv);
			string temp2(argv[i]);
			int_ps[temp]=cast_int(cast_string_to_double(temp2));
			
		} else if(string_ps.find(temp)!=string_ps.end()) {
			
			++i;
			if(i==argc)
				error_statement(argv);
			string temp2(argv[i]);
			string_ps[temp]=temp2;
			
		} else if(bool_ps.find(temp)!=bool_ps.end()) {
			
			bool_ps[temp]= !bool_ps[temp];
		} else {
			cerr<<temp<<" is an unknown flag!\n\n"<<endl;
			error_statement(argv);
		}
		
		/***************** find what temp is *****************/
		++i;		
	}

	required_flags_are_set(argv);
	
}


void parameters::required_flags_are_set(char * argv[]) {
	
	if(required_flags.size()!=0) {		
		cerr<<"flag "<<*required_flags.begin()<<" needs to be set\n\n\n"<<endl;
		error_statement(argv);
	}

}


void parameters::printing(ostream & pout) {
	
	for(map<string, int>::iterator itm=int_ps.begin(); itm!=int_ps.end(); itm++) {
		pout<<itm->first<<" "<<itm->second<<endl;
	}

	for(map<string, double>::iterator itm=double_ps.begin(); itm!=double_ps.end(); itm++) {
		pout<<itm->first<<" "<<itm->second<<endl;
	}
	
	for(map<string, string>::iterator itm=string_ps.begin(); itm!=string_ps.end(); itm++) {
		pout<<itm->first<<" "<<itm->second<<endl;
	}
	
	for(map<string, bool>::iterator itm=bool_ps.begin(); itm!=bool_ps.end(); itm++) {
		if(itm->second)
			pout<<itm->first<<" true"<<endl;
		else
			pout<<itm->first<<" false"<<endl;
			
	}
}



void set_parameters_for_docmap(parameters & P, int argc, char * argv[]) {


    P.set_string("-f", "nofile", true, "[string]: name of the corpus file (plain txt file)");
	P.set_double("-p", 0.05, false, "[float]: the p-value, any number bigger than 0 and smaller than 1. Default is 0.05. Bigger the p-value, fewer the topics.");
    P.set_int("-r", 10, false, "[int]: number of runs for Infomap. Default is 10");
    P.set_int("-t", 0., false, "[int]: minimum number of documents per topic. Default is 0, but 10 is recommended for big corpuses.");    
    P.set_double("-minf", 0., false, "[double] minimum value for the likelihood filter. Default is 0.");
    P.set_double("-maxf", 0.51, false, "[double] maximum value for the likelihood filter.");
    P.set_double("-alpha", 0.01, false, "[double] initial value of alpha.");
    P.set_double("-step", 0.01, false, "[double] step in PLSA filtering.");
	P.set_int("-seed", -1, false, "[int]: seed for random number generator. default is read from file time_seed.dat.");
	P.set_string("-part", "", false, "[string]: a file like \"infomap.part\" saved from a previous run.");
    P.set_bool("-skip_opt_al", false, false, " : use this option to skip alpha optimization");
    P.set_string("-model", "", false, "[string]: a file like \"topic_words.txt\" saved from a previous run");
    P.set_bool("-write_net", false, false, " : writes a file called \"sig_words.edges\" in the format \"wn wn weight\"");
    P.set_bool("-infer", false, false, " : performs one single E step for inferring gammas -No alpha optimization is involved. As -alpha is just a scalar, this is symmetric LDA-");
    P.set_string("-parall", "", false, "[string= \"i:j:n\"] : this is for building sig_words.edges with multiple (n^2) jobs. i,j must be >=0 and <n.");
    P.set_int("-lag", 20, false, "[int] : lda model is printed every [-lag] EM steps.");
    
	P.set_from_argv(argc, argv);
	P.printing(cout);

	if(P.int_ps.at("-seed")==-1) {
		cout<<"setting random number seed from file"<<endl;
		srand_file();
	} else
		srand5(P.int_ps["-seed"]);
    
    general_assert(P.double_ps.at("-alpha")>0, "-alpha error:: alpha should be positive");
    general_assert(P.double_ps.at("-step")>0, "-step error:: step should be positive");
    general_assert(P.string_ps.at("-parall").size()==0 or P.string_ps.at("-part").size()==0, \
                   "error: you cannot set -parall and -part together");
    
    
    
    if (P.bool_ps.at("-infer")) {
        // if program is in infer mode, you need to specify -alpha
        // or better, -alpha_file 
        if(P.string_ps.at("-alpha_file").size()==0 \
           or P.string_ps.at("-model").size()==0 \
           or P.string_ps.at("-word_wn").size()==0) {
            
            cout<<"ERROR: you are running -infer without providing alphas from file. Current alpha is: ";
            cout<<P.double_ps.at("-alpha")<<endl;
            cout<<"OR you are not providing -model OR -word_wn"<<endl;
            
        }
    }
    


}

