#include <TMB.hpp>

template <class Type>
Type cumhaz(vector<Type> Qpartition, Type survt_1, Type h1, Type h2, Type h3, Type h4, Type h5, Type h6, Type h7, Type h8, Type h9, Type h10, vector<Type> cumhaz_int,
            Type index){
  
  vector<Type> htilde(10);
  htilde(0)=h1;
  htilde(1)=h2;
  htilde(2)=h3;
  htilde(3)=h4;
  htilde(4)=h5;
  htilde(5)=h6;
  htilde(6)=h7;
  htilde(7)=h8;
  htilde(8)=h9;
  htilde(9)=h10;
  vector<Type> H = cumhaz_int * htilde;
  Type retval = 0;
  if (index == 0){
    retval = H(0);}
  if (index == 1){
    retval = H(0) + H(1);}
  if (index == 2){
    retval = H(0) + H(1) + H(2);}
  
  if (index == 3){
    retval = H(0) + H(1) + H(2) + H(3);}
  
  if (index == 4){
    retval = H(0) + H(1) + H(2) + H(3) + H(4);}
  
  if (index == 5){
    retval = H(0) + H(1) + H(2) + H(3) + H(4) + H(5);}
  
  if (index == 6){
    retval = H(0) + H(1) + H(2) + H(3) + H(4) + H(5) + H(6);}
  
  if (index == 7){
    retval = H(0) + H(1) + H(2) + H(3) + H(4) + H(5) + H(6) + H(7);}
  
  if (index == 8){
    retval = H(0) + H(1) + H(2) + H(3) + H(4) + H(5) + H(6) + H(7) + H(8);}
  
  if (index == 9){
    retval = H(0) + H(1) + H(2) + H(3) + H(4) + H(5) + H(6) + H(7) + H(8) + H(9);}
  
  
  return retval;}



// Likelihood for the survival submodel
template <class Type>
Type llc(Type h1, Type h2, Type h3, Type h4, Type h5, Type h6, Type h7, Type h8, Type h9, Type h10, Type lambda2, Type lambda3, vector<Type> gamma, vector<Type> Z, Type a_1, Type b_1, Type survt_1, Type di_1, vector<Type> Qpartition, vector<Type> cumhaz_int,
         Type index){
  
  

  Type zeta1 = (Z * gamma).sum() + lambda2 * a_1 + lambda3 * b_1; // z* gamma is a multiplied vector
  Type lik_c = 0;
  Type basehaz = 0;
  
  if(survt_1 <= Qpartition(1)){
    basehaz = h1;}
  if(survt_1 > Qpartition(1) && survt_1 <= Qpartition(2)){
    basehaz = h2;}
  if(survt_1 > Qpartition(2) && survt_1 <= Qpartition(3)){
    basehaz = h3;}
  
  if(survt_1 > Qpartition(3) && survt_1 <= Qpartition(4)){
    basehaz = h4;}
  if(survt_1 > Qpartition(4) && survt_1 <= Qpartition(5)){
    basehaz = h5;}
  
  if(survt_1 > Qpartition(5) && survt_1 <= Qpartition(6)){
    basehaz = h6;}
  if(survt_1 > Qpartition(6) && survt_1 <= Qpartition(7)){
    basehaz = h7;}
  
  if(survt_1 > Qpartition(7) && survt_1 <= Qpartition(8)){
    basehaz = h8;}
  
  if(survt_1 > Qpartition(8) && survt_1 <= Qpartition(9)){
    basehaz = h9;}
  
  if(survt_1 > Qpartition(9) && survt_1 <= Qpartition(10)){
    basehaz = h10;}
  
  lik_c = di_1 * (log(basehaz) + zeta1) - cumhaz(Qpartition, survt_1, h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,cumhaz_int,index) * exp(zeta1);
  return lik_c;  } 

template<class Type>
struct func {
  Type h1;
  Type h2;
  Type h3;
  Type h4;
  Type h5;
  Type h6;
  Type h7;
  Type h8;
  Type h9;
  Type h10;
  Type lambda2;
  Type lambda3;
  vector<Type> gamma;
  vector<Type> Z;
  Type a_1;
  Type b_1;
  Type index;
  Type survt_1; 
  Type di_1;
  Type w;
  vector<Type> Qpartition;
  vector<Type> cumhaz_int;
  
  func(Type h1_,
       Type h2_,
       Type h3_,
       Type h4_,
       Type h5_,
       Type h6_,
       Type h7_,
       Type h8_,
       Type h9_,
       Type h10_,
       Type lambda2_,
       Type lambda3_,
       vector<Type> gamma_,
       vector<Type> Z_,
       Type a_1_,
       Type b_1_,
       Type index_,
       Type survt_1_,
       Type di_1_,
       Type w_,
       vector<Type> Qpartition_,
       vector<Type> cumhaz_int_) :
    h1(h1_), h2(h2_), h3(h3_), h4(h4_), h5(h5_),
    h6(h6_), h7(h7_), h8(h8_), h9(h9_), h10(h10_), lambda2(lambda2_), lambda3(lambda3_), gamma(gamma_), Z(Z_), a_1(a_1_), b_1(b_1_),index(index_), survt_1(survt_1_),
    di_1(di_1_), w(w_), Qpartition(Qpartition_), cumhaz_int(cumhaz_int_) {}

  
  template <class T>
  T operator()(vector<T> theta){  // Evaluate function
    //T zeta1 = T(this->zeta1); // constant
    T h1 = theta(0); // variable
    T h2 = theta(1); // variable
    T h3 = theta(2); // variable
    T h4 = theta(3); // variable
    T h5 = theta(4); // variable
    T h6 = theta(5); // variable 
    T h7 = theta(6); // variable
    T h8 = theta(7); // variable
    T h9 = theta(8); // variable
    T h10 = theta(9); // variable
    
    T lambda2 = theta(10);
    T lambda3 = theta(11);
    vector<T> gamma = theta.segment(12,3);
    vector<T> Z = this->Z.template cast<T>();
    T a_1 = T(this->a_1);
    T b_1 = T(this->b_1);
    T index = T(this->index);
    T survt_1 = T(this->survt_1);
    T di_1 = T(this->di_1);
    T w = T(this-> w);
    vector<T> Qpartition = this->Qpartition.template cast<T>(); // constant
    vector<T> cumhaz_int = this->cumhaz_int.template cast<T>(); // constant
  
    return llc(h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,lambda2,lambda3,gamma,Z,a_1,b_1,survt_1,di_1,Qpartition, cumhaz_int,index) ;
  }
};


// Joint likelihood  
template<class Type>
Type objective_function<Type>::operator() ()
{ 
  
  /* Data */
  DATA_VECTOR(a1);
  DATA_VECTOR(b1);
  DATA_MATRIX(Z0);
  DATA_INTEGER(N);
  DATA_VECTOR(Qpartition);
  DATA_VECTOR(survt1);
  DATA_VECTOR(di1);
  DATA_VECTOR(cumhaz_int);
  DATA_VECTOR(index_vec);
  DATA_VECTOR(weights);

  /* Parameters */
  
  PARAMETER_VECTOR(theta);
  

  /* Joint likelihood */
  
  
  
  
  Type jnll = Type(0.0);
  
  Type h1 = theta(0);
  Type h2 = theta(1);
  Type h3 = theta(2);
  Type h4 = theta(3);
  Type h5 = theta(4);
  Type h6 = theta(5);
  Type h7 = theta(6);
  Type h8 = theta(7);
  Type h9 = theta(8);
  Type h10 = theta(9);
  Type lambda2 = theta(10);
  Type lambda3 = theta(11);
  int theta_size = theta.size();
  int gamma_size = theta_size - 12;
  vector<Type> gamma = theta.segment(12,gamma_size);
  
  //vector<Type> zeta = Z0 * gamma + lambda2*a1 + lambda3*b1;
   
  
  
  matrix<Type> g(N,theta_size);
  //the joint negative log-likelihood
  for(int i=0; i<N; i++){
    
    
    Type index = index_vec(i);
    Type survt_1 = survt1(i); 
    Type di_1 = di1(i);
    Type a_1 = a1(i);
    Type b_1 = b1(i);
    //Type zeta1 = zeta(i);
    vector<Type> Z = Z0.row(i);
    Type w = weights(i);
   
    
    //func grad;
    func<Type> f(h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, lambda2, lambda3,gamma,Z,a_1,b_1,index,survt_1,di_1,w,Qpartition,cumhaz_int);
    // put theta in the first arg?
    jnll -= w * llc(h1,h2,h3,h4,h5,h6,h7,h8,h9,h10, lambda2, lambda3, gamma, Z, a_1, b_1, survt_1,di_1,Qpartition, cumhaz_int, index);
  
    g.row(i) = autodiff::gradient(f, theta);
  } 
  ADREPORT(g);
  return jnll;
}

