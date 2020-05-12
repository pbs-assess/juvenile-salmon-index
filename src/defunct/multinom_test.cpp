# include <TMB.hpp >
template < class Type >
Type objective_function <Type >:: operator
() ()
{

DATA_MATRIX (X);
DATA_INTEGER(nobs);


PARAMETER_VECTOR(p);

Type nll = 0;

vector<Type> pest;

pest = exp(p)/sum(exp(p));

for(int i=0; i<nobs; ++i){

	vector<Type> x=X.row(i);

    nll += -dmultinom(x,pest,true);

}


ADREPORT(pest);
return nll;
}

