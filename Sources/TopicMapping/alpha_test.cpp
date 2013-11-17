
# include "../standard_package/standard_include.cpp"
# include "lda_util.cpp"
# include "alpha_optimization.cpp"


int main(int argc, char * argv[]) {



    ifstream gin("sample.txt");
    
    
    deque<DD> prs;
    string gins;
    while(getline(gin, gins)) {
        DD ginv;
        cast_string_to_doubles(gins, ginv);
        prs.push_back(ginv);
    }
    gin.close();
    //printm(prs);
    
    cout<<"size:: "<<prs.size()<<endl;
    
    double inv_=invert_digamma(12.64675656);
    cout<<inv_<<" "<<digamma(inv_)<<endl;
    optimize_alpha(prs);


    
    return 0;
}