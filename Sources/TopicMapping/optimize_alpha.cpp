
# include "../standard_package/standard_include.cpp"

#include "lda_util.cpp"
#include "alpha_optimization.cpp"


int main(int argc, char * argv[]) {

    
    if(argc<2) {
        
        cout<<argv[0]<<" [gamma_file]"<<endl;
        cout<<"computes alphas and writes to file:: ";
        cout<<"alphas.txt"<<endl;
        return -1;
    }
    
    string gins;
    ifstream gin(argv[1]);
    general_assert(gin.is_open(), "gamma_file not found!");
    
    deque<DD> gammas;
    int size=-1;
    while(getline(gin, gins)) {
        DD vs;
        cast_string_to_doubles(gins, vs);
        general_assert(size==-1 or int(vs.size())==size, \
                       "ERROR: gammas should all have same length");
        size=vs.size();
        gammas.push_back(vs);
    }
    gin.close();
    
    DD alphas;
    optimize_alpha(gammas, alphas);
    ofstream pout("alphas.txt");
    prints(alphas, pout);
    pout.close();
    
    

    return 0;

}