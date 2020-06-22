# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;


arma::mat mrdivide(arma::mat A, arma::mat B)
{
  return (solve(B.t(), A.t(), arma::solve_opts::fast )).t();
}

// [[Rcpp::export()]]
arma::mat makePSM (arma::umat myMat) {
  
  arma::uword nIndividuals = myMat.n_cols;
  arma::uword nIterations  = myMat.n_rows;
  std::cout << nIndividuals << "   " << nIterations << "\n";
  
  arma::umat myPSM(nIndividuals, nIndividuals);;
  myPSM.zeros();

  for( int k = 0; k < nIterations; k++ ){
    //if(k % 10 == 0)
    //{
    std::cout << "Iteration " << k << " of " << nIterations << "\n";
    //}
    //std::cout << k << "\n";
    
    arma::urowvec currentClustering = myMat.row(k);
    for( int i = 0; i < nIndividuals; i++ ){
      // if(i % 1000 == 0)
      // {
      //   std::cout << "Individual " << i << " of " << nIndividuals << "\n";
      // }
      
      for( int j = (i+1); j < nIndividuals; j++ ){
        if(currentClustering(i) == currentClustering(j))
        {
          myPSM(i,j)++;
        }
      }
    }
  }
  
  myPSM = myPSM + myPSM.t();
  //ret.diag = nIterations;
  //arma::mat ret = myPSM/nIndividuals;
  //arma::mat ret(nIndividuals, nIndividuals);
  //ret.zeros();
  
  arma::mat ret = arma::conv_to<arma::mat>::from(myPSM);;
  ret = ret/nIterations;
  ret.diag().ones();
  return(ret) ;
}