#include <TMB.hpp>


typedef TMBad::ad_aug ad;
struct model {
  matrix<ad> Xrow;
  vector<ad> dfcost_row; // Data
  ad N, rep, st1;// alpha0, alpha1,  logsigma_a;
  vector<ad> sigma;
  ad operator()(vector<ad> a1) {
    //type individual likelihood here
    ad logsigma_a = sigma(0);
    int sigma_size =sigma.size();
    int alpha_size = sigma_size - 1;
    vector<ad> alpha = sigma.segment(1,alpha_size);
    //ad alpha0 = sigma(1);
    //ad alpha1 = sigma(2);
    int st2 = CppAD::Integer(st1);
    //int w1 = CppAD::Integer(w);
    int rep2 = CppAD::Integer(rep);
    int repind = 0;
    
    matrix<ad> Xrow1(rep2,alpha_size);
    for(int j=st2; j<st2+rep2; j++){
      for(int k=0; k<alpha_size; k++){
        Xrow1(repind,k) = Xrow(j,k); //
      }
      repind = repind + 1;
    }
    
    vector<ad> eta_row =  Xrow1 * alpha + a1(0);
    
    ad nll = -dnorm(a1(0),ad(0),exp(logsigma_a), true); 
    vector<ad> pi_i = exp(eta_row)/(1+exp(eta_row));
    
    //ad l_a = 0; 
    for(int h=0; h<rep; h++){
      if(dfcost_row(h) > 0){
        nll -= log(pi_i(h));
        
      }
      if(dfcost_row(h) == 0){
        nll -= log(1-pi_i(h)); 
        
      }
    }
    return nll;
  }
  template<class Type>
  Type eval_nldens(vector<Type> &start) {
    //vector<Type> start(1); //start.setZero();
    newton::newton_config cfg;
    //cfg.trace=true;
    Type res = newton::Laplace(*this, start, cfg);
    return res;
  }
};

// for gradient fn


 struct model_evaluator {
 model obj;
 ad operator()(vector<ad> sigma) {
 vector<ad> a1(1);
 a1.setZero();
 obj.sigma = sigma;
 return obj.eval_nldens<ad>(a1);
 }
 };
 


// Joint likelihood  
template<class Type>
Type objective_function<Type>::operator() ()
{ 
  
  /* Data */
  DATA_MATRIX(Xrow);
  DATA_INTEGER(N);
  DATA_VECTOR(replicates); 
  DATA_MATRIX(dfcost);
  DATA_VECTOR(weights);
  
  DATA_VECTOR(a1); //randeff initial guess
  DATA_VECTOR(start);
  /* Parameters */   
  //PARAMETER(logsigma_a);
  //PARAMETER(alpha0);
  //PARAMETER(alpha1);
  PARAMETER_VECTOR(sigma);
  
  /* Joint likelihood */
  
  
  
  
  Type nll = 0;
  
  //R.E initial guess
  //vector<Type> a1(N);
  //a1.setZero();
  vector<Type> a2(N);
  
  //Type logsigma_a = sigma(0);
  //Type alpha0 = sigma(1);
  //Type alpha1 = sigma(2);
  //Type alpha2 = sigma(3);
  //Type alpha3 = sigma[4];
  
  int sigma_size = sigma.size();
  matrix<Type> g(N,sigma_size);
  
  //the joint negative log-likelihood
  for(int i=0; i<N; i++){
    
    //vector<Type> X_row = Xrow.row(i);
    
    //vector<Type> eta_row = eta.row(i);
    vector<Type> dfcost_row = dfcost.row(i);
    Type rep = replicates(i);
    Type st1 = start(i);
    Type w = weights(i);
    

    
 
    model obj = {Xrow, dfcost_row, N, rep, st1, sigma}; // put data vectors here
    //nll += obj.eval_nldens<Type>();
    nll += w * obj.eval_nldens<Type>(a1);
    a2(i) = a1(0);
    model_evaluator F = {obj};
    g.row(i) = autodiff::gradient(F, sigma);
    
  }
  //model_evaluator F = {obj};
  //vector<Type> g = autodiff::gradient(F, sigma);
  //matrix<Type> h = autodiff::hessian(F, sigma);
  ADREPORT(a2);
  ADREPORT(g);
  return nll;
}
