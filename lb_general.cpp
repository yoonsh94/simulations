#include <TMB.hpp>

// Likelihood for the amount of +ve medical cost
typedef TMBad::ad_aug ad;
struct model {
  matrix<ad> Xrow;
  vector<ad> dfcost_row; // Data
  ad a_1, N, rep, costindexsum1, costsum1, st1;
  vector<ad> sigma;
  ad operator()(vector<ad> b1) {
    //type individual likelihood here
    //int sigma_size = sigma.size();

    ad logsigma_b = sigma(0);
    ad delta = sigma(1);
    ad lambda1 = sigma(2);
    
    int sigma_size = sigma.size();
    int beta_size = sigma_size - 3; //3
    //vector<ad> beta(beta_size);
    vector<ad> beta = sigma.segment(3,beta_size); //modular betas
      
    int st2 = CppAD::Integer(st1);
    int rep2 = CppAD::Integer(rep);
    int repind = 0;
    
    matrix<ad> Xrow1(rep2,beta_size);
    for(int j=st2; j<st2+rep2; j++){
      for(int k=0; k<beta_size; k++){
        Xrow1(repind,k) = Xrow(j,k); 
      }
      repind = repind + 1;
    }

    //vector<ad> mu_row = 1/ (Xrow1 * beta + lambda1*a_1 + b1(0)); // returns mu with length rep
    vector<ad> mu_row = exp(Xrow1 * beta + lambda1*a_1 + b1(0));
    
    
    // ad sigma1 = delta;
    ad sigma_exp = exp(delta); 
    
    ad nll = -dnorm(b1(0),ad(0),exp(logsigma_b), true); 
    
    //Type lik_b = 0;
    
    if (costindexsum1 == 0) {
      nll -= 0;
    }

    if (costindexsum1 > 0) {
      
      for(int i=0; i<rep2; i++){
        if (dfcost_row(i) > 0){
          
          nll -= dgamma(dfcost_row(i),sigma_exp, mu_row(i)/sigma_exp,TRUE);
          // nll -=
          
          //nll -= 0;
        }
        if (dfcost_row(i) == 0){
          nll -= 0;
        }
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

 
struct model_evaluator {
model obj;
ad operator()(vector<ad> sigma) {
     vector<ad> b1(1);
     b1.setZero();
     obj.sigma = sigma;
     return obj.eval_nldens<ad>(b1);
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
  
  DATA_VECTOR(a1);
  DATA_VECTOR(b1); //randeff initial guess


  DATA_VECTOR(costindexsum);
 DATA_VECTOR(costsum);
  DATA_VECTOR(start);

  //DATA_INTEGER(maxrep);
  //DATA_VECTOR(weight);
  /* Parameters */
  

  //PARAMETER_VECTOR(b1); 

  //PARAMETER(logsigma_b);
  //PARAMETER(delta);
  //PARAMETER(beta0);
  //PARAMETER(beta1);
  //PARAMETER(lambda1);
  //PARAMETER(delta);
 
  PARAMETER_VECTOR(sigma);  
  /* Joint likelihood */
  
  Type nll = 0;
  
  vector<Type> b2(N); //randeff
  
  int sigma_size = sigma.size();
  matrix<Type> g(N,sigma_size); //grad

  //the joint negative log-likelihood
  for(int i=0; i<N; i++){
    //vector<Type> X_row = Xrow.row(i);
    vector<Type> dfcost_row = dfcost.row(i);
    Type rep = replicates(i);
    //vector<Type> mu_row = mu.row(i);
    Type costindexsum1 = costindexsum(i);
    Type costsum1 = costsum(i);
    Type a_1 = a1(i);
    Type st1 = start(i);
    Type w = weights(i);
    //Type weight1 = weight(i);
    model obj = {Xrow, dfcost_row, a_1, N, rep,costindexsum1,costsum1, st1,sigma};

    nll += w * obj.eval_nldens<Type>(b1);
    
    b2(i) = b1(0);
    model_evaluator F = {obj};
    g.row(i) = autodiff::gradient(F, sigma);
  } 
  ADREPORT(b2);
  ADREPORT(g);
  return nll;
}




