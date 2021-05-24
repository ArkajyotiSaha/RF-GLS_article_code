//------------------------------------------------
//  updated: 3/19/2019
//  Author: Arkajyoti Saha
//  Email: arkajyotisaha93@gmal.com
//------------------------------------------------
#include <cstdlib>
#include <time.h>
#include <Rmath.h>
#include <iostream>     // std::cout
#include <algorithm>    // std::min
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

// // --------------------------------//--------------------------------
// #include <omp.h>
// //[[Rcpp::plugins(openmp)]]
// // --------------------------------//--------------------------------

using namespace arma;
using namespace Rcpp;
using namespace std;
using namespace stats;
using namespace R;

#define swapInt(a, b) ((a ^= b), (b ^= a), (a ^= b))
#define NODE_TERMINAL -1
#define NODE_TOSPLIT  -2
#define NODE_INTERIOR -3
#define a 1.0e-10
#define a_jump 1
#define inf_temp 1/1e-22
//--------------------------------
//--------------------------------
// [[Rcpp::export]]
double dist2(double a1, double a2, double b1, double b2){
  return(sqrt(pow(a1 - b1,2)+pow(a2 - b2,2)));
}

// [[Rcpp::export]]
double spCor(double D, double phi, double nu, int covModel){
  
  //0 and default exponential
  //1 spherical
  //2 matern
  //3 gaussian
  //4 Identity
  
  if(covModel == 1){//spherical
    
    if(D > 0 && D <= 1.0/phi){
      return 1.0 - 1.5*phi*D + 0.5*pow(phi*D,3);
    }else if(D >= 1.0/phi){
      return 0.0;
    }else{
      return 1.0;
    }
  }else if(covModel == 2){//matern
    
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*pi/2*(besselI(d*phi,-nu)-besselI(d*phi, nu))/sin(nu*pi), or
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi, nu=nu)
    
    if(D*phi > 0.0){
      return pow(D * phi, nu)/(pow(2, nu-1)*gammafn(nu))*bessel_k(D * phi, nu, 1.0);//thread safe bessel
    }else{
      return 1.0;
    }
  }else if(covModel == 3){//gaussian
    
    return exp(-1.0*(pow(phi*D,2)));
    
  }else if(covModel == 4){//gaussian
    
    if(D == 0){
      return 1;
    }else{
      return 0;
    }
    
  }else{
    //exponential
    
    return exp(-phi*D);
  }
}


// [[Rcpp::export]]
mat SigmaCholinv(double sigmasq, double tausq, double phi, double nu, int covModel, int n, mat coords){
  mat covar(n,n), covarchol(n,n);
  double D;
  
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      D = dist2(coords(i,0),coords(i,1),coords(j,0),coords(j,1));
      covar(i,j) = sigmasq * spCor(D, phi, nu, covModel);
      if(i == j){
        covar(i,j) = covar(i,j) + tausq;
      }
    }
  }
  covarchol = (inv(trans(chol(covar))));
  return(covarchol);
}


// [[Rcpp::export]]
mat PSigmaCholinv(uvec &P_index, mat Sigmacholinv){
  return(Sigmacholinv.rows(P_index));
}

// [[Rcpp::export]]
double costnew(uvec &Z_index, mat PSigmacholInv, vec y, int n, int rc, uvec P_index){
  mat Z(n,rc);
  mat res_temp;
  double res_temp_min;
  mat res;
  double new_test = 100000;
  for(int i = 0; i < n ; i++){
    for(int j = 0; j < rc; j++){
      Z(i,j) = 0;
    }
    Z(i,Z_index(i)) = 1;
  }
  vec beta(rc);
  vec residual(rc);
  mat temp(n,n);
  temp = trans(PSigmacholInv) * PSigmacholInv;
  res_temp_min = 100000;
  for(int lpower = 0; lpower < a_jump; lpower++){
    beta = pinv(trans(Z)*temp*Z, a/(10^lpower))*trans(Z)*temp*y;
    residual = y - Z * beta;
    res_temp = residual.t() * temp * residual; 
    if(res_temp_min > res_temp(0,0)){
      res_temp_min = res_temp(0,0);
    }
  }
  res = sum(Z.rows(P_index));
  for(int test = 0; test < rc; test++)
  {
    if(new_test > res(0,test) ){
      new_test = res(0,test);
    }
  }
  if(new_test > 0){
    return(res_temp_min);
  }else{
    return(inf_temp);
  }
}

// [[Rcpp::export]]
double rightbeta(uvec &Z_index, mat PSigmacholInv, vec y, int n, int rc){
  mat Z(n,rc);
  mat res_temp;
  double res_temp_min;
  for(int i = 0; i < n ; i++){
    for(int j = 0; j < rc; j++){
      Z(i,j) = 0;
    }
    Z(i,Z_index(i)) = 1;
  }
  vec beta(rc);
  vec residual(rc);
  mat temp(n,n);
  temp = trans(PSigmacholInv) * PSigmacholInv;
  res_temp_min = 10000000;
  int lpower_min = 0;
  for(int lpower = 0; lpower < a_jump; lpower++){
    beta = pinv(trans(Z)*temp*Z, a/(10^lpower))*trans(Z)*temp*y;
    residual = y - Z * beta;
    res_temp = residual.t() * temp * residual; 
    if(res_temp_min > res_temp(0,0)){
      res_temp_min = res_temp(0,0);
      lpower_min = lpower;
    }
  }
  beta = pinv(trans(Z)*temp*Z, a/(10^lpower_min))*trans(Z)*temp*y;
  return(beta[rc-1]);
}

// [[Rcpp::export]]
double leftbeta(uvec &Z_index, mat PSigmacholInv, vec y, int n, int rc, int nnumber){
  mat Z(n,rc);
  mat res_temp;
  double res_temp_min;
  for(int i = 0; i < n ; i++){
    for(int j = 0; j < rc; j++){
      Z(i,j) = 0;
    }
    Z(i,Z_index(i)) = 1;
  }
  vec beta(rc);
  vec residual(rc);
  mat temp(n,n);
  temp = trans(PSigmacholInv) * PSigmacholInv;
  res_temp_min = 10000000;
  int lpower_min = 0;
  for(int lpower = 0; lpower < a_jump; lpower++){
    beta = pinv(trans(Z)*temp*Z, a/(10^lpower))*trans(Z)*temp*y;
    residual = y - Z * beta;
    res_temp = residual.t() * temp * residual; 
    if(res_temp_min > res_temp(0,0)){
      res_temp_min = res_temp(0,0);
      lpower_min = lpower;
    }
  }
  beta = pinv(trans(Z)*temp*Z, a/(10^lpower_min))*trans(Z)*temp*y;
  return(beta[nnumber]);
}

// [[Rcpp::export]]
List findBestSplit(mat &x, uvec &jdex, vec &y, int mdim, int nsample,
                   int ndstart, int ndend, int &msplit, double &decsplit,
                   double &ubest, int &ndendl, int &jstat, int mtry,
                   int nodecnt, uvec &Z_index, int rc, mat PSigmaCholinvmat, uvec P_index) {
  int last, nl, nr;
  int i, j, kv, tieVal, tieVar;
  vec xt(nsample), ut(nsample), yl(nsample); 
  
  double ubestt;
  double crit, critmax, critvar;
  vec v(ndend - ndstart + 1), sortv(ndend - ndstart + 1);
  uvec mind(mdim), ncase2(nsample), ncase(nsample), ncase3(nsample);
  uvec ndind(ndend - ndstart + 1), ncase2_short(ndend - ndstart + 1), ncase_short(ndend - ndstart + 1);
  uvec Z_index_local(nsample);
  
  
  /* START BIG LOOP */
  msplit = -1;
  decsplit = 0.0;
  //critmax = arma::datum::inf;
  critmax = 100000;
  ubestt = 0.0;
  
  for (i=0; i < mdim; ++i) mind[i] = i;
  last = mdim - 1;
  
  
  for(int ilp = 0; ilp < mtry; ++ilp){
    for(int il = 0; il < nsample; il++){
      Z_index_local[il] = Z_index[il];
    }
    
    //critvar = arma::datum::inf;
    critvar = 100000;
    j = (int) (R::unif_rand() * (last+1));
    kv = mind[j];
    swapInt(mind[j], mind[last]);
    last--;
    
    /* numeric variable */
    for (j = ndstart; j <= ndend; ++j) {
      xt(j) = x(jdex[j] - 1, kv);
      yl(j) = y(jdex[j] - 1);
    }
    /* copy the x data in this node. */
    for (j = ndstart; j <= ndend; ++j) v(j  - ndstart) = xt(j);
    for (j = 1; j <= nsample; ++j) ncase[j - 1] = j;
    for (j = 1; j <= nsample; ++j) ncase2[j - 1] = jdex[j-1];
    for(int index = 0; index < (ndend - ndstart + 1); index++){
      ncase2_short[index]  = ncase2[ndstart+index];
      ncase_short[index]  = ncase[ndstart+index];
    }
    ndind = sort_index(v);
    sortv = sort(v);
    for(int index = 0; index < (ndend - ndstart + 1); index++){
      ncase2[ndstart+index] = ncase2_short[ndind[index]];
      ncase[ndstart+index] = ncase_short[ndind[index]];
    }
    //crit = arma::datum::inf;
    crit = 100000;
    tieVal = 1;
    tieVar = 1;
    for(j = ndstart; j <= ndend - 1; ++j) {
      Z_index_local[ncase2[j] - 1] = rc;//new mode is right child
      if (sortv[j - ndstart] < sortv[j+1 - ndstart]) {
        crit = costnew(Z_index_local, PSigmaCholinvmat, y, nsample, rc+1, P_index);
        if (crit == critvar) {
          tieVal++;
          if(R::unif_rand() < 1.0/tieVal){
            ubestt =  (sortv[j - ndstart] + sortv[j+1 - ndstart]) / 2.0;
            critvar = crit;
          }
        }
        if (crit < critvar) {
          ubestt = (sortv[j - ndstart] + sortv[j+1 - ndstart]) / 2.0;
          critvar = crit;
          tieVal = 1;
        }
      }
    }
    if(critvar < 100000){
      if (critvar == critmax) {
        tieVar++;
        if (R::unif_rand() < 1.0 / tieVar) {
          ubest = ubestt;
          msplit = kv + 1;
          ncase3 = ncase2;
          for (j = ndstart; j <= ndend; ++j) {
            ut[j] = xt[j];
          }
          critmax = critvar;
        }
      }
      if (critvar < critmax) {
        ubest = ubestt;
        msplit = kv + 1;
        critmax = critvar;
        ncase3 = ncase2;
        for (j = ndstart; j <= ndend; ++j) {
          ut[j] = xt[j];
        }
        tieVar = 1;
      } 
    }
  }
  
  
  decsplit = critmax;
  
  /* If best split can not be found, set to terminal node and return. */
  if (msplit != -1) {
    nl = ndstart;
    for (j = ndstart; j <= ndend; ++j) {
      if (ut[j] <= ubest) {
        nl++;
        Z_index[jdex(j) - 1] = rc;
      }
    }
    ndendl = imax2(nl - 1, ndstart);
    nr = ndendl + 1;
    for (j = ndstart; j <= ndend; ++j) {
      if (ut[j] > ubest) {
        if (nr >= nsample) break;
        nr++;
      }
    }
    if (ndendl >= ndend) ndendl = ndend - 1;
    jdex = ncase3;
    
  }
  else jstat = 1;
  
  
  mat Z(nsample,rc+1);
  mat res;
  double new_test = 100000;
  for(int i = 0; i < nsample ; i++){
    for(int j = 0; j < (rc+1); j++){
      Z(i,j) = 0;
    }
    Z(i,Z_index(i)) = 1;
  }
  res = sum(Z.rows(P_index));
  for(int test = 0; test < rc; test++)
  {
    if(new_test > res(0,test) ){
      new_test = res(0,test);
    }
  }
  //Rprintf("the value of col_min : %f and value %f and msplit is %i jstat is %i critvar is %f critmax is %f\n", new_test, decsplit, msplit, jstat, critvar, critmax);
  List ret;
  ret ["jdex"] = jdex;
  ret ["ubest"] = ubest;
  ret ["Z_index"] = Z_index;
  ret ["ndendl"] = ndendl;
  ret ["msplit"] = msplit;
  ret ["critmax"] = critmax;
  ret ["jstat"] = jstat;
  return ret;
}



// [[Rcpp::export]]
uvec regTree(mat &x, vec &y, int mdim, int nsample, NumericVector &lDaughter,
             NumericVector &rDaughter,
             vec &upper, vec &avnode, NumericVector &nodestatus, int nrnodes,
             int &treeSize, int nthsize, int mtry, NumericVector &mbest, NumericVector &varUsed, mat &PSigmacholInv, uvec P_index) {
  int i, k, ncur;
  NumericVector nodestart(nrnodes), nodepop(nrnodes);
  uvec jdex(nsample);
  int ndstart, ndend, ndendl, nodecnt, jstat, msplit = 0;
  double decsplit, ubest;
  uvec Z_index(nsample);
  for (i = 1; i <= nsample; ++i){
    jdex[i-1] = i;
    Z_index[i-1] = 0;
  }
  ncur = 0;
  nodestart[0] = 0;
  nodepop[0] = nsample;
  nodestatus[0] = NODE_TOSPLIT;
  
  int rc = 1;
  avnode[0] = rightbeta(Z_index, PSigmacholInv, y, nsample, rc);
  
  for (k = 0; k < nrnodes - 2; ++k){
    if (k > ncur || ncur >= nrnodes - 2) break;
    if (nodestatus[k] != NODE_TOSPLIT) continue;
    /* initialize for next call to findbestsplit */
    ndstart = nodestart[k];
    ndend = ndstart + nodepop[k] - 1;
    nodecnt = nodepop[k];
    jstat = 0;
    decsplit = 0.0;
    ndstart = nodestart[k];
    ndend = ndstart + nodepop[k] - 1;
    nodecnt = nodepop[k];
    jstat = 0;
    decsplit = 0.0;
    
    List result = findBestSplit(x, jdex, y, mdim, nsample,
                                ndstart, ndend, msplit, decsplit,
                                ubest, ndendl, jstat, mtry,
                                nodecnt, Z_index, rc, PSigmacholInv, P_index);
    if (jstat == 1) {
      /* Node is terminal: Mark it as such and move on to the next. */
      nodestatus[k] = NODE_TERMINAL;
    }
    if (jstat != 1){
      mbest[k] = msplit;
      varUsed[msplit - 1] = 1;
      upper[k] = ubest;
      nodestatus[k] = NODE_INTERIOR;
      //Rprintf("the value of nodestatus[%u] : %i \n", k, nodestatus[k]);
      /* leftnode no.= ncur+1, rightnode no. = ncur+2. */
      nodepop[ncur + 1] = ndendl - ndstart + 1;
      nodepop[ncur + 2] = ndend - ndendl;
      nodestart[ncur + 1] = ndstart;
      nodestart[ncur + 2] = ndendl + 1;
      
      
      avnode[ncur + 1] = rightbeta(Z_index, PSigmacholInv, y, nsample, rc+1);
      nodestatus[ncur + 1] = NODE_TOSPLIT;
      if (nodepop[ncur + 1] <= nthsize) {
        nodestatus[ncur + 1] = NODE_TERMINAL;
      }
      
      avnode[ncur + 2] = leftbeta(Z_index, PSigmacholInv, y, nsample, rc+1, Z_index[jdex[ndend] - 1]);
      nodestatus[ncur + 2] = NODE_TOSPLIT;
      if (nodepop[ncur + 2] <= nthsize) {
        nodestatus[ncur + 2] = NODE_TERMINAL;
      }
      
      lDaughter[k] = ncur + 1 + 1;
      rDaughter[k] = ncur + 2 + 1;
      /* Augment the tree by two nodes. */
      ncur += 2;
      rc = rc+1;
      /* Found the best split. */
    }
  }
  treeSize = nrnodes;
  for (k = nrnodes - 1; k >= 0; --k) {
    if (nodestatus[k] == 0) (treeSize)--;
    if (nodestatus[k] == NODE_TOSPLIT) {
      nodestatus[k] = NODE_TERMINAL;
    }
  }
  return Z_index;
}


// [[Rcpp::export]]
vec predictRegTree(mat &x, int nsample, int mdim, NumericVector &lDaughter,
                   NumericVector &rDaughter, NumericVector &nodestatus, vec &upper, vec &avnode, 
                   NumericVector &mbest){
  int i, k, m;
  vec ypred(nsample);
  for(int kj = 0; kj < nsample; kj++){
    ypred[kj] = 0;
  }
  for (i = 0; i < nsample; ++i) {
    k = 0;
    while (nodestatus[k] != NODE_TERMINAL) { /* go down the tree */
      m = mbest[k] - 1;
      
      if(x(i,m) <= upper[k]){
        k = lDaughter[k] - 1;
      }else{
        k = rDaughter[k] - 1;
      }
    }
    /* terminal node: assign prediction and move on to next */
    ypred[i] = avnode[k];
  }
  return(ypred);
}




// [[Rcpp::export]]
List regRFtestnew(mat &x, vec &y, int mdim, int nsample, int nthsize, int nrnodes, int ntree, int mtry, 
                  double sigmasq, double tausq, double phi, double nu, int covModel, mat &coords, mat &xtest, int ntest){
  mat PSigmaCholinvmat(nsample, nsample), SigmaCholinvmat(nsample, nsample);
  uvec P_index(nsample);
  umat P_index_record(nsample, ntree);
  unsigned int k;\
  int l;
  double xrand;
  NumericVector lDaughter(nrnodes), rDaughter(nrnodes);
  NumericVector varUsed(mdim), nodestatus(nrnodes), mbest(nrnodes);
  vec upper(nrnodes), avnode(nrnodes);
  SigmaCholinvmat = SigmaCholinv(sigmasq, tausq, phi, nu, covModel, nsample, coords);
  uvec test;
  mat prediction(nsample, ntree), predictiontest(ntest, ntree);
  IntegerVector treesize(ntree);
  
  for(int j =0; j < ntree; ++j){
    for (int n = 0; n < nsample; ++n) {
      xrand = R::unif_rand();
      k = xrand * nsample;
      P_index(n) = k;
      P_index_record(n,j) = k;
    }
    for(l = 0; l < nrnodes; l++){
      lDaughter[l]  = 0;
      rDaughter[l] = 0;
      upper[l]  = 0;
      avnode[l] = 0;
      nodestatus[l] = 0;
      mbest[l] = 0;
    }
    PSigmaCholinvmat= PSigmaCholinv(P_index, SigmaCholinvmat);
    test = regTree(x, y, mdim, nsample, lDaughter,
                   rDaughter,
                   upper, avnode, nodestatus, nrnodes,
                   treesize[j], nthsize, mtry, mbest, varUsed, PSigmaCholinvmat, P_index);
    prediction.col(j) = predictRegTree(x, nsample, mdim, lDaughter, rDaughter, nodestatus, upper, avnode, mbest);
    predictiontest.col(j) = predictRegTree(xtest, ntest, mdim, lDaughter, rDaughter, nodestatus, upper, avnode, mbest);
  }
  arma::mat betahat ;
  mat regX(nsample,2);
  
  for(int hj = 0; hj < nsample; hj++){
    regX(hj,0) = 1;
  }
  regX.col(1) = trans(sum(prediction.t())/ntree);
  betahat = (regX.t() * regX ).i() * regX.t() * y ;
  List ret;
  ret ["insample_prediction_uncor"] = prediction;
  ret ["outsample_rediction_uncor"] = predictiontest;
  ret ["3"] = P_index_record;
  ret ["4"] = ntree;
  ret ["6"] = regX;
  //ret ["7"] = trans(sum(prediction.t())/ntree) * betahat[1] + betahat[0];
  //ret ["8"] = trans(sum(predictiontest.t())/ntree) * betahat[1] + betahat[0];
  ret ["7"] = trans(sum(prediction.t())/ntree);
  ret ["8"] = trans(sum(predictiontest.t())/ntree);
  
  
  return ret;
}
