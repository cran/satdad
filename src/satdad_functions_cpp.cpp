// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <boost/math/special_functions/binomial.hpp>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;


template < typename OutputIterator >
struct SolutionWriter
{
  OutputIterator outputIt;
  template < typename SolBeginIterator, typename SolEndIterator >
  void operator() ( SolBeginIterator first, SolEndIterator const& end )
  {
    for ( ; first != end; ++first )
      *outputIt++ = *first;
  }
};

struct SolutionPrinter
{
  template < typename SolBeginIterator, typename SolEndIterator >
  void operator() ( SolBeginIterator first, SolEndIterator const& end ) const
  {
    std::cout << "[ ";
    for ( ; first != end; ++first )
      std::cout << *first << " ";
    std::cout << "]" << std::endl;
  }
};

template < typename SetBeginIterator, typename SetEndIterator, typename Value, typename Observer >
void genSubset( SetBeginIterator setBegin, SetEndIterator setEnd, std::vector<Value> & currentSol, std::size_t length, Observer & observer )
{
  if ( length == 0 )
  {
    observer( currentSol.cbegin(), currentSol.cend() );
    return;
  }
  while ( setBegin != setEnd )
  {
    currentSol.push_back( *setBegin );
    genSubset( ++setBegin, setEnd, currentSol, length-1, observer );
    currentSol.pop_back();
  }
}


IntegerMatrix matsubset( IntegerVector const& A, int length ){
  std::vector<int> currentSol;
  currentSol.reserve( length );
  IntegerMatrix M( length, boost::math::binomial_coefficient<double>( A.length(), length ) );
  SolutionWriter<IntegerMatrix::iterator> observer={ M.begin() };
  genSubset( A.begin(), A.end(), currentSol, length, observer );
  return M;
}

// [[Rcpp::export(.subsets_cpp)]]
List subsets_cpp(unsigned int const& d){
  List result(pow(2,d)-1);
  unsigned int ListIndice=0;
  for (unsigned int j=0; j<d; ++j){
    IntegerMatrix M=matsubset(seq(1,d),j+1);
    unsigned int n=M.nrow();
    unsigned int N=M.ncol();
    for(unsigned int i=0; i<N; ++i){
      IntegerVector V=rep(0,n);
      for (unsigned int k=0; k<n; ++k){
        V(k)=M(k,i);
      }
      result[ListIndice]=V;
      ++ListIndice;
    }
  }
  return result;
}

// [[Rcpp::export(.sort_sub_cpp)]]
std::vector<std::vector<int>> sort_sub_cpp(std::vector<std::vector<int>> sub){
  std::sort(sub.begin(), sub.end(), [](const std::vector<int> & a, const std::vector<int> & b){ return a.size() < b.size(); });
  return sub;
}


std::vector<double> runif_cpp(int n) {
  std::vector<double> values(n);
	std::generate(values.begin(), values.end(), [&](){ return R::runif(0.0,1.0); });
	return values;
}


std::vector<std::vector<double>> matunif(unsigned int const&  n, unsigned int const&  m){
  std::vector<std::vector<double>> mat(n, std::vector<double>(m));
  for (unsigned int i=0; i<n; ++i){
    mat[i]=runif_cpp(m);
  }
  return mat;
}


// [[Rcpp::export(.generate_asy_sub_cpp)]]
std::vector<std::vector<double>> generate_asy_sub_cpp(unsigned int  const& d, std::vector<std::vector<int>> const& sub){
	std::vector<std::vector<double>> result;
	std::vector<std::vector<double>> Unif(d, std::vector<double>(sub.size()));
	for (unsigned int j=0; j<sub.size(); ++j){
		for(unsigned int k=0; k<sub[j].size(); ++k){
			Unif[sub[j][k]-1][j]=as<double>(Rcpp::runif(1));
		}
	}
	std::vector<std::vector<double>> Tnif(sub.size(),std::vector<double>(d));
	for (unsigned int i=0; i<d; ++i){
		double S=std::accumulate(Unif[i].begin(), Unif[i].end(), 0.0);
		std::transform(Unif[i].begin(), Unif[i].end(), Unif[i].begin(),
                 std::bind(std::multiplies<double>(), std::placeholders::_1,1./S ));
		for (int j=0;j<sub.size(); j++){
			Tnif[j][i] = Unif[i][j];
		}
	}
	for(unsigned int j=0; j<sub.size(); ++j){
		std::vector<double> V=Tnif[j];
		V.erase(std::remove(V.begin(), V.end(), 0.0),V.end());
		result.push_back(V);
	}
	return result;
}

// [[Rcpp::export(.find_missing_indices_cpp)]]
std::vector<int> find_missing_indices_cpp(int d, std::vector<std::vector<int>> sub)
{
	std::vector<int> res(d);
	for(int k=0; k<d;++k){res[k]=k+1;}
	for (int i = 0; i < sub.size(); ++i)
	{
		for (int j = 0; j < sub[i].size(); ++j)
		{
			res[sub[i][j]-1]=0;
		}
	}
	res.erase(remove(res.begin(), res.end(), 0),res.end());
	if(res.size()==0){res={0};}
return res;
}



double qexpo(double u) {
	return -log(1-u);
}

double sinu(double v) {
	return sin(v);
}


double sim_d_stable(double alpha){
	double U=M_PI*R::runif(0.0,1.0);
	return pow(sinu(U*(1.-alpha))/qexpo(R::runif(0.0,1.0)),1./alpha-1.)*sinu(U*alpha)/pow(sin(U),1./alpha);
}


double invqexpo(double u, double alpha) {
	double v=-log(1-u);
	return 1./pow(v,alpha);
}


// [[Rcpp::export(.algo21_numMat_cpp)]]
NumericMatrix algo21_numMat_cpp(int n,  int d, double alpha){
	NumericMatrix X(n,d);
	for(unsigned int k=0; k<n ; ++k){
		NumericVector W=runif(d);
		std::transform(W.begin(), W.end(),W.begin(),[=](double w){ return invqexpo(w, alpha); });
		double S=sim_d_stable(alpha);
		std::transform(W.begin(), W.end(), W.begin(),
                 std::bind(std::multiplies<double>(), std::placeholders::_1,pow(S,alpha)));
		X.row(k)=W;
	}
	return X;
}



NumericVector max_vectNum_evd( NumericVector x, NumericVector y){
	int n=x.size();
	NumericVector z(n);
	for(unsigned int k=0; k<n ; ++k){
		z[k]=std::max(x[k],y[k]);
	}
	return z;
}



// [[Rcpp::export(.algo22_numMat_cpp)]]
NumericMatrix algo22_numMat_cpp(int n, int d, std::vector<std::vector<int>> sub,  std::vector<double> dep, std::vector<std::vector<double>> asy){
	NumericMatrix X(n,d);
	for(unsigned int k=0; k<sub.size() ; ++k){
		NumericMatrix aux=algo21_numMat_cpp(n,sub[k].size(),dep[k]);
		for(unsigned int j=0; j<sub[k].size() ; ++j){
					std::transform(aux.column(j).begin(), aux.column(j).end(), aux.column(j).begin(),
			              std::bind(std::multiplies<double>(), std::placeholders::_1,asy[k][j] ));
			int indice=sub[k][j]-1;
			NumericVector bux=max_vectNum_evd(X.column(indice),aux.column(j));
			X.column(indice)=bux;
		}
	}
	return X;
}


// [[Rcpp::export(.ellmevlogv_cpp)]]
double ellmevlogv_cpp(NumericVector v,  int d, std::vector<std::vector<int>> sub, std::vector<double> dep, std::vector<std::vector<double>> asy){
	double res = 0;
	for(unsigned int k=0; k<sub.size() ; ++k){
		double aux =0 ;
		for(unsigned int j=0; j<sub[k].size() ; ++j){
			aux += pow(asy[k][j]*v[sub[k][j]-1],1/dep[k]);
		}
		res += pow(aux,dep[k]);
	}
	return res;
}

double norm_subsetb_cpp(NumericVector v, int d, std::vector<int> b, double alpha, std::vector<double> asyb){
	double bux_norm=0;
	int db=b.size();
	for(unsigned int j=0; j<db ; ++j){
		bux_norm += pow(asyb.at(j)*v[b.at(j)-1],1/alpha);
	}
	return pow(bux_norm,alpha);
}

double diff_subsets_weighted_norm_cpp(NumericVector v,  int d,std::vector<int> b ,std::vector<int> sub, double dep, std::vector<double> asy){
	double res=1;
	double norm=norm_subsetb_cpp(v,d,sub,dep,asy);
	for(unsigned int i=0; i<b.size(); ++i){
		auto  pos=std::find(sub.begin(), sub.end(), b[i]);
		int index = std::distance(sub.begin(), pos);
		if( pos != sub.end()) {
			/* sub contains b(i) */
			if(i==0){res *= (dep-i)/dep*pow(asy.at(index),1/dep)*pow(v[sub[index]-1],1/dep-1)*pow(norm,(dep-1)/dep);
			}else{res *= (dep-i)/dep*pow(asy.at(index),1/dep)*pow(v[sub[index]-1],1/dep-1)*pow(norm,-1/dep);}
		} else {
			/* sub does not contain b(i) */
			res *= 0;
			i=b.size();
		}
	}
	return res;
}



double diff_ellmevlogv_cpp(NumericVector v,  int d, std::vector<int> b, std::vector<std::vector<int>> sub, std::vector<double> dep, std::vector<std::vector<double>> asy){
	double res = 0;
	for(unsigned int k=0; k<sub.size() ; ++k){
		res += diff_subsets_weighted_norm_cpp(v,d,b,sub[k],dep[k],asy[k]);
	}
	return res;
}

NumericVector invers(NumericVector v){
	NumericVector res(v.size());
	for(unsigned int k=0; k<v.size(); ++k){
		res[k]=1/v[k];
	}
	return res;
}

NumericMatrix invers_m(NumericMatrix v){
	NumericMatrix res(v.nrow(),v.ncol());
	for(unsigned int k=0; k<v.nrow(); ++k){
		res(k,_)=invers(v(k,_));
	}
	return res;
}

NumericVector inversmar(NumericVector v, NumericVector mar){
	// vérifier domaine où c'est positif
	NumericVector res(v.size());
	for(unsigned int k=0; k<v.size(); ++k){
		res[k]=pow((1+mar[2]*(v[k]-mar[0])/mar[1]),-1/mar[2]);
	}
	return res;
}

NumericMatrix inversmar_m(NumericMatrix v, NumericVector mar){
	// vérifier domaine où c'est positif
	NumericMatrix res(v.nrow(),v.ncol());
	for(unsigned int k=0; k<v.nrow(); ++k){
		res(k,_)=inversmar(v(k,_),mar);
	}
	return res;
}



NumericVector inversmarm(NumericVector v, NumericMatrix mar){
	// vérifier domaine où c'est positif
	NumericVector res(v.size());
	for(unsigned int k=0; k<v.size(); ++k){
		res[k]=pow((1+mar(k,2)*(v[k]-mar(k,0))/mar(k,1)),-1/mar(k,2));
	}
	return res;
}

NumericMatrix inversmarm_m(NumericMatrix v, NumericMatrix mar){
	// vérifier domaine où c'est positif
	NumericMatrix res(v.nrow(),v.ncol());
	for(unsigned int k=0; k<v.nrow(); ++k){
		res(k,_)=inversmarm(v(k,_),mar);
	}
	return res;
}

double pis(NumericVector v,std::vector<int> b){
	double res=1;
	for(unsigned int k=0; k<b.size(); ++k){
		res *= 1/pow(v[b.at(k)-1],2.0);
	}
	return res;
}

double pismar(NumericVector v,std::vector<int> b,NumericVector mar){
	double res=1;
	if(mar[2]!=0){
		for(unsigned int k=0; k<b.size(); ++k){
		res *= pow(1+mar[2]*(v[b.at(k)-1]-mar[0])/mar[1],-1/mar[2]-1)/mar[1];
		}
	}else{
		for(unsigned int k=0; k<b.size(); ++k){
			res *= exp(-(v[b.at(k)-1]-mar[0])/mar[1])/mar[1];
		}
	}
	return res;
}

double pismarm(NumericVector v,std::vector<int> b,NumericMatrix mar){
	double res=1;
	for(unsigned int k=0; k<b.size(); ++k){
		if(mar(k,2)!=0){
			res *= pow(1+mar(k,2)*(v[b.at(k)-1]-mar(k,0))/mar(k,1),-1/mar(k,2)-1)/mar(k,1);
		}else{
			res *= exp(-(v[b.at(k)-1]-mar(k,0))/mar(k,1))/mar(k,1);
		}
	}
	return res;
}


IntegerMatrix setparts_cpp(int d){
  Rcpp::Function setparts = Rcpp::Environment("package:partitions")["setparts"];
  return setparts(d);
}


std::vector<std::vector<int>> transf_cpp(IntegerVector v) {
  int n = v.size();
  std::vector<std::vector<int>> result(max(v));
  for(int i = 0; i < n; i++) {
    int idx = v[i] - 1;
    std::vector<int>& vec = result[idx];
    vec.push_back(i+1);
  }
  return result;
}


List ListPart_cpp(int d){
  IntegerMatrix m = setparts_cpp(d);
  List result(m.ncol());
  for(int i = 0; i<m.ncol(); i++){
    std::vector<std::vector<int>> aux=transf_cpp(m.column(i));
    result[i]=aux;
  }
  return result;
}



// [[Rcpp::export(.ratio_diff_pmevlogv_cpp)]]
double ratio_diff_pmevlogv_cpp(NumericVector v,  int d,  std::vector<std::vector<int>> sub, std::vector<double> dep, std::vector<std::vector<double>> asy){
	List P=ListPart_cpp(d);
	NumericVector vinv=invers(v);
	double aux=0;
	for(unsigned int k=0; k<P.size(); ++k){
		List p=P[k];
		double auxp=1;
		for(unsigned int i=0; i<p.size(); ++i){
			std::vector<int> b=p[i];
			auxp *= diff_ellmevlogv_cpp(vinv,d,b,sub,dep,asy) * pow(-1,b.size()) * pis(v,b);
		}
		aux += pow(-1,p.size())*auxp;
	}
	return aux;
}


std::vector<double> sumweight(std::vector<double> aux,std::vector<double> auxp,double sign){
	std::vector<double> res(aux.size());
	for(unsigned int i=0; i<aux.size(); ++i){
		res[i]=aux[i]+sign*auxp[i];
	}
	return res;
}

// [[Rcpp::export(.ratio_diff_pmevlogm_cpp)]]
std::vector<double> ratio_diff_pmevlogm_cpp(NumericMatrix v,  int d,  std::vector<std::vector<int>> sub, std::vector<double> dep, std::vector<std::vector<double>> asy){
	List P=ListPart_cpp(d);
	NumericMatrix vinv=invers_m(v);
	std::vector<double> aux(v.nrow());
	for(unsigned int k=0; k<P.size(); ++k){
		List p=P[k];
		std::vector<double> auxp(v.nrow());
		std::fill(auxp.begin(), auxp.end(), 1);
		for(unsigned int i=0; i<p.size(); ++i){
			std::vector<int> b=p[i];
			for(unsigned int j=0; j<v.nrow(); ++j){
				auxp.at(j) *= diff_ellmevlogv_cpp(vinv(j,_),d,b,sub,dep,asy) * pow(-1,b.size()) * pis(v(j,_),b);
			}
		}
		aux = sumweight(aux,auxp,pow(-1,p.size()));
	}
	return aux;
}


// [[Rcpp::export(.ratio_diff_pmevlogv_marv_cpp)]]
double ratio_diff_pmevlogv_marv_cpp(NumericVector v,  int d,  std::vector<std::vector<int>> sub, std::vector<double> dep, std::vector<std::vector<double>> asy, NumericVector mar){
	List P=ListPart_cpp(d);
	NumericVector vinv=inversmar(v,mar);
	double aux=0;
	for(unsigned int k=0; k<P.size(); ++k){
		List p=P[k];
		double auxp=1;
		for(unsigned int i=0; i<p.size(); ++i){
			vector<int> b=p[i];
			auxp *= diff_ellmevlogv_cpp(vinv,d,b,sub,dep,asy) * pow(-1,b.size()) * pismar(v,b,mar);
		}
		aux += pow(-1,p.size())*auxp;
	}
	return aux;
}

// [[Rcpp::export(.ratio_diff_pmevlogm_marv_cpp)]]
std::vector<double> ratio_diff_pmevlogm_marv_cpp(NumericMatrix v,  int d,  std::vector<std::vector<int>> sub, std::vector<double> dep, std::vector<std::vector<double>> asy, NumericVector mar){
	List P=ListPart_cpp(d);
	NumericMatrix vinv=inversmar_m(v,mar);
	std::vector<double> aux(v.nrow());
	for(unsigned int k=0; k<P.size(); ++k){
		List p=P[k];
		std::vector<double> auxp(v.nrow());
		std::fill(auxp.begin(), auxp.end(), 1);
		for(unsigned int i=0; i<p.size(); ++i){
			std::vector<int> b=p[i];
			for(unsigned int j=0; j<v.nrow(); ++j){
				auxp.at(j) *= diff_ellmevlogv_cpp(vinv(j,_),d,b,sub,dep,asy) * pow(-1,b.size()) * pismar(v(j,_),b,mar);
			}
		}
		aux = sumweight(aux,auxp,pow(-1,p.size()));
	}
	return aux;
}



// [[Rcpp::export(.ratio_diff_pmevlogv_marm_cpp)]]
double ratio_diff_pmevlogv_marm_cpp(NumericVector v,  int d,  std::vector<std::vector<int>> sub, std::vector<double> dep, std::vector<std::vector<double>> asy, NumericMatrix mar){
	List P=ListPart_cpp(d);
	NumericVector vinv=inversmarm(v,mar);
	double aux=0;
	for(unsigned int k=0; k<P.size(); ++k){
		List p=P[k];
		double auxp=1;
		for(unsigned int i=0; i<p.size(); ++i){
			std::vector<int> b=p[i];
			auxp *= diff_ellmevlogv_cpp(vinv,d,b,sub,dep,asy) * pow(-1,b.size()) * pismarm(v,b,mar);
		}
		aux += pow(-1,p.size())*auxp;
	}
	return aux;
}

// [[Rcpp::export(.ratio_diff_pmevlogm_marm_cpp)]]
std::vector<double> ratio_diff_pmevlogm_marm_cpp(NumericMatrix v,  int d,  std::vector<std::vector<int>> sub, std::vector<double> dep, std::vector<std::vector<double>> asy, NumericMatrix mar){
	List P=ListPart_cpp(d);
	NumericMatrix vinv=inversmarm_m(v,mar);
	vector<double> aux(v.nrow());
	for(unsigned int k=0; k<P.size(); ++k){
		List p=P[k];
		vector<double> auxp(v.nrow());
		std::fill(auxp.begin(), auxp.end(), 1);
		for(unsigned int i=0; i<p.size(); ++i){
			std::vector<int> b=p[i];
			for(unsigned int j=0; j<v.nrow(); ++j){
				auxp.at(j) *= diff_ellmevlogv_cpp(vinv(j,_),d,b,sub,dep,asy) * pow(-1,b.size()) * pismarm(v(j,_),b,mar);
			}
		}
		aux = sumweight(aux,auxp,pow(-1,p.size()));
	}
	return aux;
}


// [[Rcpp::export(.ellmevlogm_cpp)]]
NumericVector ellmevlogm_cpp(NumericMatrix v,  int d, std::vector<std::vector<int>> sub, std::vector<double> dep, std::vector<std::vector<double>> asy){
	int m=v.nrow();
	NumericVector res(m);
	for(unsigned int i=0; i<m ; ++i){
		res(i)=ellmevlogv_cpp(v.row(i),d,sub,dep,asy);
	}
	return res;
}



// [[Rcpp::export(.ecdsmevlog_cpp)]]
NumericVector ecdsmevlog_cpp ( int d, std::vector<std::vector<int>> sub, std::vector<double> dep, std::vector<std::vector<double>> asy, List ind)
{
	NumericVector res(ind.size());
	for (int i=0; i<ind.size(); ++i)
	{
		NumericVector v(d);
		for(int j=0; j<as<IntegerVector>(ind[i]).size(); ++j){
			v(as<IntegerVector>(ind[i])[j]-1)=1;
		}
		res(i)=ellmevlogv_cpp(v,d,sub,dep,asy);
	}
	return res;
}



NumericVector subs(NumericVector U, int i,double u){
	int d=U.size();
	NumericVector Ui(d);
	for(unsigned int j=0; j<d; ++j){
		if(j==(i-1)){Ui[j]=u;
		}else{Ui[j]=U[j];}
	}
	return Ui;
}

NumericVector substit(NumericVector U, NumericVector V, IntegerVector pos){
	NumericVector W(V.size());
	std::copy(V.begin(), V.end(), W.begin() ) ;
	int l=pos.size();
	for(unsigned int j=0; j<l; ++j){
		W[pos[j]-1]=U[pos[j]-1];
	}
	return W;
}

double int_Dijs_cpp(NumericVector U, double v, double w, int i, int j, int d, std::vector<std::vector<int>> sub, std::vector<double> dep, std::vector<std::vector<double>> asy){
	NumericVector Ui=subs(U,i,v);
	NumericVector Uj=subs(U,j,w);
	NumericVector Uij=subs(Ui,j,w);
	double l=ellmevlogv_cpp(U,d,sub,dep,asy);
	double li=ellmevlogv_cpp(Ui,d,sub,dep,asy);
	double lj=ellmevlogv_cpp(Uj,d,sub,dep,asy);
	double lij=ellmevlogv_cpp(Uij,d,sub,dep,asy);
	double aux=l-li-lj+lij;
	return aux*aux;
}



double Dijs_cpp(int N,int i, int j, int d, std::vector<std::vector<int>> sub, std::vector<double> dep, std::vector<std::vector<double>> asy){
	double sum=0;
	for(int k=0;k<N;++k){
		NumericVector U=runif(d);
		double v=R::runif(0.0,1.0);
		double w=R::runif(0.0,1.0);
		sum += int_Dijs_cpp(U,v,w,i,j,d,sub,dep,asy);
	}
	return sum/(4*N);
}


NumericVector tsicdsmevlog_pairs_cpp(int N, int d,std::vector<std::vector<int>> sub,std::vector<double> dep, std::vector<std::vector<double>> asy){
	IntegerMatrix M=matsubset(seq(1,d),2);
	int m=M.cols();
	NumericVector res=NumericVector(m);
	for(int k=0; k<m; ++k){
		res[k]=Dijs_cpp(N,M(0,k),M(1,k),d,sub,dep,asy);
	}
	return res;
}

double tsicdsmevlog_cpp (int N, int d,std::vector<std::vector<int>> sub,std::vector<double> dep,std::vector<std::vector<double>> asy,IntegerVector U){
	int u=U.size();
	double integrant=0;
	for(int i=0;i<N;++i){
		NumericVector X=runif(d);
		NumericVector Z=runif(d);
		double aux= pow(-1.0,u-0)*ellmevlogv_cpp(Z,d,sub,dep,asy);
		for(int k=1; k<u+1;++k){
		IntegerMatrix V=matsubset(U,k);
		double sign=pow(-1.0,u-k);
		for(int j=0; j<V.cols();++j){
			IntegerVector v=V(_,j);
			NumericVector x=substit(X,Z,v);
			aux += sign*ellmevlogv_cpp(x,d,sub,dep,asy);
			}
		}
		integrant +=  pow(aux,2.0);
	}
	integrant /= (pow(2,u)*N);
	return integrant;
}


// [[Rcpp::export(.tsicdsmevlog_list_cpp)]]
NumericVector tsicdsmevlog_list_cpp(int N, int d,std::vector<std::vector<int>> sub,std::vector<double> dep,std::vector<std::vector<double>> asy,List L){
	int l=L.size();
	NumericVector res(l);
	for(int i=0;i<l;++i){
		res(i)=tsicdsmevlog_cpp(N,d,sub,dep,asy,as<IntegerVector>(L[i]));
	}
	return res;
}

NumericMatrix runifmat(int n, int m) {
	NumericVector v = runif(n * m);
	v.attr("dim") = Dimension(n, m);
	return as<NumericMatrix>(v);
}

double meanpow2(NumericVector aux, int N){
	double res=0;
	for(int i=0;i<N;++i){
		res += pow(aux(i),2.0)/N;
	}
	return res;
}

	NumericVector tsicdsmevlog_paires_cpp (int N, int d,std::vector<std::vector<int>> sub,std::vector<double> dep,std::vector<std::vector<double>> asy){
		NumericVector integrants(d*(d-1)/2);
		int compteur=0;
		NumericMatrix const X=runifmat(N,d);
		NumericVector RES=ellmevlogm_cpp(X,d,sub,dep,asy);
		NumericMatrix const Z=runifmat(N,d);
			for(int i=0; i<d-1;++i){
				NumericMatrix X0(Rcpp::clone(X));
				X0(_,i)=Z(_,i);
				NumericVector RES0=ellmevlogm_cpp(X0,d,sub,dep,asy);
				for(int j=i+1; j<d;++j){
					NumericMatrix X01(Rcpp::clone(X));
					X01(_,i)=Z(_,i);
					NumericMatrix X1(Rcpp::clone(X));
					X1(_,j)=Z(_,j);
					X01(_,j)=Z(_,j);
					NumericVector aux=RES-RES0-ellmevlogm_cpp(X1,d,sub,dep,asy)+ellmevlogm_cpp(X01,d,sub,dep,asy);
					integrants(compteur)=meanpow2(aux,N)/4.0;
					compteur += 1;
				}
			}
			return integrants;
		}

NumericVector tsicdsmevlog_singletons_cpp (int N, int d,std::vector<std::vector<int>> sub,std::vector<double> dep,std::vector<std::vector<double>> asy){
	NumericVector integrants(d);
	NumericMatrix const X=runifmat(N,d);
	NumericVector RES=ellmevlogm_cpp(X,d,sub,dep,asy);
	NumericMatrix const Z=runifmat(N,d);
	for(int i=0; i<d;++i){
		NumericMatrix X0(Rcpp::clone(X));
		X0(_,i)=Z(_,i);
		NumericVector RES0=ellmevlogm_cpp(X0,d,sub,dep,asy);
		NumericVector aux=RES-RES0;
		integrants(i)=meanpow2(aux,N)/2.0;
	}
	return integrants;
}


// [[Rcpp::export(.tsicdsmevlog_empty_cpp)]]
double tsicdsmevlog_empty_cpp(int N, int d,std::vector<std::vector<int>> sub,std::vector<double> dep,std::vector<std::vector<double>> asy){
  NumericMatrix const X=runifmat(N,d);
  NumericVector RES0=ellmevlogm_cpp(X,d,sub,dep,asy);
  double m=mean(RES0);
  NumericMatrix const Z=runifmat(N,d);
  NumericVector RES=ellmevlogm_cpp(Z,d,sub,dep,asy);
  double va=meanpow2(RES,N)-pow(m,2);
  return va;
}
/////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Empirical ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


NumericVector EllEmp( double k, arma::mat const & x, arma::mat const & sample )
{
  const std::size_t n = sample.n_rows;
  const std::size_t d = sample.n_cols;
  const std::size_t nx = x.n_rows;
  NumericVector res(nx);
  for ( std::size_t jx = 0; jx < nx; ++jx ) {
    arma::uvec vec(n, arma::fill::zeros);
    for ( std::size_t j = 0; j < d; ++j ) {
    int tmpjx = std::floor(n + 0.5 - k*x(jx,j));
    if(tmpjx < 0) {
      vec.fill(1);
      break;
    } else if(tmpjx < n) {
      const arma::uvec idx = arma::sort_index( sample.col(j) );
      vec(idx.subvec(tmpjx, n-1)).fill(1);
    }
    }
    double aux = arma::accu(vec) / k;
    res(jx) =  std::min( std::max( aux, (x.row(jx)).max()), arma::accu(x.row(jx)) );
  }
  return res;
}

// [[Rcpp::export(.ellEmp_cpp)]]
NumericMatrix ellEmp_cpp( arma::vec const &  kvec, arma::mat const & x, arma::mat const & sample )
{
  const std::size_t nk = kvec.size();
  NumericMatrix res(nk,x.n_rows);
  for (unsigned int ik = 0; ik < nk; ++ik ) {
    res(ik,_)=EllEmp(kvec(ik),x,sample );
  }
  return res;
}





// [[Rcpp::export(.rankbar_cpp)]]
arma::mat rankbar_cpp(arma::mat const & x, unsigned int k){
  const std::size_t n = x.n_rows;
  const std::size_t d = x.n_cols;
  arma::mat U = arma::zeros(n, d);
  for(unsigned int j =0; j<d ; j++){
    const arma::uvec ord = arma::sort_index(x.col(j))+1;
    const arma::uvec rank = arma::sort_index(ord)+1;
    for(unsigned int i =0; i<n ; i++){
      double res=rank[i];
      double aux1=(n+1.-res)/k;
      double aux2=std::min(aux1,1.);
      U(i,j)=aux2;
    }
  }
  return U;
}


arma::mat ratiop1surp2(arma::mat const & R,int const & i){
  const std::size_t n = R.n_rows;
  arma::mat U(n, n);
  for(unsigned int s=0; s<n ; s++){
    for(unsigned int sprime=0; sprime<n ; sprime++){
      U(s,sprime)=(std::min(R(s,i),R(sprime,i))-R(s,i)*R(sprime,i))/std::min(R(s,i),R(sprime,i));
    }
  }
  return U;
}


arma::mat ratiop1surp2sobol(arma::mat const & R,int const & i){
  const std::size_t n = R.n_rows;
  arma::mat U(n, n);
  for(unsigned int s=0; s<n ; s++){
    for(unsigned int sprime=0; sprime<n ; sprime++){
      U(s,sprime)=(std::min(R(s,i),R(sprime,i))-R(s,i)*R(sprime,i))/(R(s,i)*R(sprime,i));
    }
  }
  return U;
}


// [[Rcpp::export(.ratioall_cpp)]]
List ratioall_cpp(arma::mat const & R){
  const std::size_t d = R.n_cols;
  List Uall(d);
  for(unsigned int i=0; i<d ; i++){
    Uall[i]=ratiop1surp2(R,i);
  }
  return Uall;
}

// [[Rcpp::export(.ratioallsobol_cpp)]]
List ratioallsobol_cpp(arma::mat const & R){
  const std::size_t d = R.n_cols;
  List Uall(d);
  for(unsigned int i=0; i<d ; i++){
    Uall[i]=ratiop1surp2sobol(R,i);
  }
  return Uall;
}

// [[Rcpp::export(.prodmin_cpp)]]
arma::mat prodmin_cpp(arma::mat const & R){
  const std::size_t n = R.n_rows;
  const std::size_t d = R.n_cols;
  arma::mat U=arma::ones(n, n);
  for(unsigned int s=0; s<n ; s++){
    for(unsigned int sprime=0; sprime<n ; sprime++){
      for(unsigned int i=0; i<d ; i++){
        U(s,sprime) *= std::min(R(s,i),R(sprime,i));
      }
    }
  }
  return U;
}


// [[Rcpp::export(.prodprod_cpp)]]
arma::mat prodprod_cpp(arma::mat const & R){
  const std::size_t n = R.n_rows;
  const std::size_t d = R.n_cols;
  arma::mat U=arma::ones(n, n);
  for(unsigned int s=0; s<n ; s++){
    for(unsigned int sprime=0; sprime<n ; sprime++){
      for(unsigned int i=0; i<d ; i++){
        U(s,sprime) *= R(s,i)*R(sprime,i);
      }
    }
  }
  return U;
}
