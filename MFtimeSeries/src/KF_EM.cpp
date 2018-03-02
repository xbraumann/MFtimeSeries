// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
double glogdet (arma::mat X) {
  
  vec eigval;
  mat eigvec;
  
  eig_sym(eigval, eigvec, X);
  
  int i;
  int k=0;
  double gdet=0.0;
  for (i=0; i<X.n_rows; i++ ) {
    
    if (eigval(i)>1e-8) {
      k+=1;
      gdet += log(eigval(i));
    }
  }
  //Rprintf("\n %d \n", k);
  return(gdet);
}



// [[Rcpp::export]]
void showMatrix(arma::mat X, const char* name) {
    Rcout << name << std::endl << X << std::endl;
}

//' The Kalman filter equations for state space systems
//' 
//' @param Y The data matrix
//' @param A The transiation matrix
//' @param B The matrix corresponding to the transition error, i.e. \eqn{\Sigma=BB'}
//' @param C The observation matrix
//' @param H The observation error covariance matrix
//' @param S The covariance between observation and transition error
//' @param a1 The mean of state x1
//' @param P1 The covariance of state x1
//' @param N The sampling frequency of the slow component
//' @param nf The number of fast components
//' 
//' @return A list with the usual Kalman filter variables.
//' 
// [[Rcpp::export]]
Rcpp::List kfilter (const arma::mat Y,
                    const arma::mat A,
                    const arma::mat B0,
                    const arma::mat C,
                    const arma::mat H,
                    const arma::mat S,
                    const arma::mat a1,
                    const arma::mat P1,
                    const int N,
                    const int nf) {
  
   int j;
   int n = Y.n_cols;
   int m = A.n_rows;
   int T = Y.n_rows;
   int q = B0.n_cols;
   mat at = zeros(T+1, m);             
   cube Pt = zeros(m, m, T+1); 
   cube Kt = zeros(m, n, T);
   double ll = 0;         
               
   //  Rprintf("\n %d \n", q);

   mat B = zeros(m, q);
   B(span(), span(0, q-1)) = B0;
      
   mat C2 (C.n_rows, C.n_cols, fill::zeros);
   mat H2 (H.n_rows, H.n_cols, fill::zeros);
   mat S2 (S.n_rows, S.n_cols, fill::zeros);
   
   if (nf<n) {
    C2(span(0,nf-1), span()) = C( span(0,nf-1), span());
    H2(span(0,nf-1), span(0,nf-1)) = H(span(0,nf-1), span(0,nf-1));
    H2(span(nf, H.n_rows-1), span(nf, H.n_rows-1) ) = eye (n-nf,n-nf);
    S2(span(),span(0,nf-1)) = S(span(),span(0,nf-1));
   }
   
   at(0, span()) = a1;
   Pt(span(), span(), span(0)) = P1;   
   
   mat tmpPt = zeros(m, m);
   mat tmpmxm = zeros(m,m);
   mat tmpnxn = zeros(n,n);
   mat tmpAPCBS = zeros(m, n);
   mat tmpFt = zeros(n, n);
   mat tmpKt = zeros(m, n);
          
   cube Ft = zeros(n, n, T);
   mat vt = zeros(T, n);
   
   mat tmpll = zeros(1,1);
   mat tmpvt = zeros(1, n);
       
 
   for (j=0; j<T; j++) {

     // Rprintf("\n %d", j);

     tmpPt = Pt(span(), span(), span(j));
               
     // NA values 
     if ((j+1) % N) {
       
       vt(j, span()) = Y(j, span()) - at(j, span()) * C2.t();
       tmpnxn = C2*tmpPt*C2.t() + H2;
       tmpFt = (tmpnxn + tmpnxn.t())/2;
       Ft(span(), span(), span(j)) = tmpFt;
       
       tmpAPCBS = A*tmpPt*C2.t() + B*S2;
   
       // showMatrix(tmpFt, "Ft");

       tmpKt = tmpAPCBS*pinv(tmpFt);
       Kt(span(), span(), span(j)) = tmpKt;
         
     } else {
             
       vt(j, span()) = Y(j, span()) - at(j, span()) * C.t();
       tmpnxn = C*tmpPt*C.t() + H;
       tmpFt = (tmpnxn + tmpnxn.t())/2;  
       Ft(span(), span(), span(j)) = tmpFt;
       
       tmpAPCBS = A*tmpPt*C.t() + B*S;
       tmpKt = tmpAPCBS*pinv(tmpFt);
       Kt(span(), span(), span(j)) = tmpKt;
   
     }  
     
     at(j+1, span()) = at(j, span())*A.t() + vt(j, span())*tmpKt.t();
     // NUMERICAL PROBLEMS!! AUSLOESCHUNG --> offenbar hat Sigma = (Sigma + t(Sigma))/2 geholfen
     
     Pt(span(), span(), span(j+1)) = A*tmpPt*A.t() + B*B.t() - tmpKt*tmpFt*tmpKt.t();
   
     tmpvt = vt(j, span());
     tmpll = tmpvt*pinv(tmpFt)*tmpvt.t();
     //Rprintf("\n %d \n", q);
     
     // vec eigval = eigs_sym (tmpFt, q);
     
     ll = ll + glogdet(tmpFt) + tmpll(0,0);
   }              
   
   // Rprintf("\n After filter equations \n");
   
   /* The smoothing equations */
   mat atT = zeros(T+1, m);             
   cube PtT = zeros(m, m, T+1); 
   
   cube Pttm1 = zeros(m, m, T+1);
   cube Nt = zeros(m, m, T);
   mat rt = zeros(T, m);
   mat tmpLt = zeros(m, m);
   mat tmprt = zeros(1, m);
   mat tmpNt = zeros(m,m);
   mat tmpPtp1 = zeros(m,m);
   
   // r_T = N_T = 0 by definition;
   tmprt = rt(T-1,span());
   tmpNt = Nt(span(), span(), span(T-1));
   
   // 
   atT(T,span(0,n-1)) = at(T,span(0,n-1));
   PtT(span(), span(), span(T)) = Pt(span(), span(), span(T));
   // and PtT is zero;
   for (j=T-1; j>=0; j--) {
     
   // Rprintf("\n %d \n", j);
    
    tmpPt = Pt(span(), span(), span(j));
    tmpFt = Ft(span(), span(), span(j));
    tmpFt = pinv(tmpFt);
    tmpKt = Kt(span(), span(), span(j));
     
    if ((j+1) % N) {
     
      tmpLt = A - tmpKt*C2;
      // note that tmpFt is the inverse here!
      tmprt = vt(j,span())*tmpFt*C2 + tmprt*tmpLt;
      
      tmpNt = C2.t()*tmpFt*C2 + tmpLt.t()*tmpNt*tmpLt;
     
     } else {
       
      tmpLt = A - tmpKt*C;
      // note that tmpFt is the inverse here!
      tmprt = vt(j,span())*tmpFt*C + tmprt*tmpLt;
      
      tmpNt = C.t()*tmpFt*C + tmpLt.t()*tmpNt*tmpLt;
     }
     
     atT(j,span()) = at(j, span()) + tmprt*tmpPt.t();
     PtT(span(), span(), span(j)) = tmpPt - tmpPt*tmpNt*tmpPt;
     
     if (j<T-1) {
       tmpPtp1 = Pt(span(), span(), span(j+1));
       Pttm1(span(), span(), span(j+1)) = (eye(m,m)-tmpPtp1*tmpNt.t())*tmpLt*tmpPt;
     }
     
     
   }

 //   showMatrix(Pt(span(),span(),span(2)));
//    showMatrix(Pt(span(),span(),span(3)));

   return Rcpp::List::create(Rcpp::Named("at")=at,
                             Rcpp::Named("Pt")=Pt,
                             Rcpp::Named("vt")=vt,
                             Rcpp::Named("Ft")=Ft,
                             Rcpp::Named("Kt")=Kt,
                             Rcpp::Named("atT")=atT,
                             Rcpp::Named("PtT")=PtT,
                             Rcpp::Named("Pttm1")=Pttm1,
                             Rcpp::Named("rt")=rt,
                             Rcpp::Named("Nt")=Nt,
                             Rcpp::Named("loglik")=ll,
                             Rcpp::Named("Cred")=C2,
                             Rcpp::Named("Hred")=H2,
                             Rcpp::Named("Sred")=S2);
                      
}



*/

//' The Expectation step in the EM Algorithm
//' 
//' Compute the sufficient statistics in the EM Algorithm for state space models. Ignore the starting values.
//' 
//' @rdname EMalg
//' 
// [[Rcpp::export]]
Rcpp::List Estep0002 (const arma::mat Y,
                  const arma::mat A,
                  const arma::mat B,
                  const arma::mat C,
                  const arma::mat H,
                  const arma::mat S,
                  const arma::mat a1,
                  const arma::mat P1,
                  const int N,
                  const int nf) {
  
  Rcpp::List reskf;
  reskf = kfilter (Y, A, B, C, H, S, a1, P1, N, nf);
    
  int T;
  int n;
  int m;
  int j;
  int p;
  T = Y.n_rows;
  n = Y.n_cols;
  m = A.n_rows;
  p = floor(C.n_cols/C.n_rows);
  
  SEXP at0 = reskf[5];
  NumericVector at1(at0);
  mat at(at1);
  at.set_size(T+1,m);
  
  SEXP Pt0 = reskf[6];
  NumericVector Pt1(Pt0);
  cube Pt(Pt1.begin(), m, m, T+1);
  
  SEXP Pcs0 = reskf[7];
  NumericVector Pcs1(Pcs0);
  cube Pcs(Pcs1.begin(), m, m, T+1);

  SEXP ll0 = reskf[10];
  NumericVector ll(ll0);
  
  mat Sx = zeros(m, m);
  mat Syx = zeros(n, m);
  mat Sy = zeros(n, n);
  
  mat tmp = zeros(m, m);
  mat tmp2 = zeros(m,m);
  
  // tmp = Pt(span(), span(), span(1));
  // S00 = at(1, span()).t()*at(1, span()) + tmp;
  
  //Rprintf("\n %d \n", p);
  
  for (j=p; j<T; j++) {
  
//  Rprintf("\n %d \n", j);
    tmp = Pt(span(), span(), span(j));
    tmp2 = Pt(span(0, n-1), span(0, n-1), span(j+1));
       
    Sy = Sy + at(j+1, span(0,n-1)).t()*at(j+1, span(0,n-1)) + tmp2;
    Sx = Sx + at(j, span()).t()*at(j, span()) + tmp;
    
    tmp = Pcs(span(0, n-1), span(), span(j+1));
    
    Syx = Syx + at(j+1, span(0,n-1)).t()*at(j, span()) + tmp;
  
  }

/*
  tmp = Pt(span(), span(), span(T));
  S11 = S11 + at(T, span()).t()*at(T, span()) + tmp;

  tmp = Pcs(span(), span(), span(T-1)); 
  S10 = S10 + at(T, span()).t()*at(T-1, span()) + tmp;
*/
  mat V1 = zeros(m, m);
 // V1 = Pt(span(), span(), span(0));
  mat mu1 = zeros(1, m);
//  mu1 = at(0, span());
  
  return Rcpp::List::create(Rcpp::Named("at")=at,
                            Rcpp::Named("Pt")=Pt,
                            Rcpp::Named("Pttm1")=Pcs,
                            Rcpp::Named("Sx")=Sx,
                            Rcpp::Named("Syx")=Syx,
                            Rcpp::Named("Sy")=Sy,
			    Rcpp::Named("loglik")=ll,
			    Rcpp::Named("mu1")=mu1,
			    Rcpp::Named("V1")=V1);
  
}





//' The Expectation step in the EM Algorithm
//' 
//' Compute the sufficient statistics in the EM Algorithm for state space models. Including starting values x1 with mean mu1 and covariance V1. See Shumway and Stoffer 1982.
//' 
//' @rdname EMalg
//'
// [[Rcpp::export]]
Rcpp::List Estep0003 (const arma::mat Y,
                      const arma::mat A,
                      const arma::mat B,
                      const arma::mat C,
                      const arma::mat H,
                      const arma::mat S,
                      const arma::mat a1,
                      const arma::mat P1,
                      const int N,
                      const int nf) {
  
  Rcpp::List reskf;
  reskf = kfilter (Y, A, B, C, H, S, a1, P1, N, nf);
    
  int T;
  int n;
  int m;
  int j;
  int p;
  T = Y.n_rows;
  n = Y.n_cols;
  m = A.n_rows;
  p = floor(C.n_cols/C.n_rows);
  
  SEXP at0 = reskf[5];
  NumericVector at1(at0);
  mat at(at1);
  at.set_size(T+1,m);
  
  SEXP Pt0 = reskf[6];
  NumericVector Pt1(Pt0);
  cube Pt(Pt1.begin(), m, m, T+1);
  
  SEXP Pcs0 = reskf[7];
  NumericVector Pcs1(Pcs0);
  cube Pcs(Pcs1.begin(), m, m, T+1);

  SEXP ll0 = reskf[10];
  NumericVector ll(ll0);
  
  mat Sx = zeros(m, m);
  mat Syx = zeros(n, m);
  mat Sy = zeros(n, n);
  
  mat tmp = zeros(m,m);
  mat tmp2 = zeros(m,m);
  
  //tmp = Pt(span(), span(), span(0));
  //Sx = at(0, span()).t()*at(0, span()) + tmp;
  
//  Rprintf("\n\n New STEP \n");
  
  for (j=0; j<T; j++) {
  
    if (j<=4) {
  //     Rprintf("\n %d \n", j);
    }
    tmp = Pt(span(), span(), span(j));
    Sx = Sx + at(j, span()).t()*at(j, span()) + tmp;
   
    if (j<=4) {
 //    showMatrix(tmp);
    }
   
    tmp2 = Pt(span(0, n-1), span(0, n-1), span(j+1));
    Sy = Sy + at(j+1, span(0,n-1)).t()*at(j+1, span(0,n-1)) + tmp2;
       
    tmp = Pcs(span(0, n-1), span(), span(j+1));
    if (j<=4) {
 //   Rprintf("\n Pcs \n");
 //    showMatrix(tmp);
    }
    Syx = Syx + at(j+1, span(0,n-1)).t()*at(j, span()) + tmp;
  }

/*
  tmp = Pt(span(), span(), span(T));
  S11 = S11 + at(T, span()).t()*at(T, span()) + tmp;

  tmp = Pcs(span(), span(), span(T-1)); 
  S10 = S10 + at(T, span()).t()*at(T-1, span()) + tmp;
*/
  mat mu1 = zeros(1, m);
  mu1 = at(0, span());
  mat V1 = mu1.t()*mu1;
  tmp = Pt(span(), span(), span(0));
  V1 = V1 + tmp;
  
  return Rcpp::List::create(Rcpp::Named("at")=at,
                            Rcpp::Named("Pt")=Pt,
                            Rcpp::Named("Pttm1")=Pcs,
                            Rcpp::Named("Sx")=Sx,
                            Rcpp::Named("Syx")=Syx,
                            Rcpp::Named("Sy")=Sy,
			    Rcpp::Named("loglik")=ll,
			    Rcpp::Named("mu1")=mu1,
			    Rcpp::Named("V1")=V1);
  
}


