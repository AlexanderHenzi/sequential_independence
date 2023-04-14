// headers ---------------------------------------------------------------------
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// pd_ds headers for ordered sets
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;

// useful definitions and functions --------------------------------------------
// ordered set
typedef tree<double, null_type, less<double>, rb_tree_tag,
             tree_order_statistics_node_update>
  ordered_set;

// Rcpp implementation of findInterval search (from http://adv-r.had.co.nz/Rcpp.html)
IntegerVector find_interval(NumericVector x, NumericVector breaks) {
  IntegerVector out(x.size());
  NumericVector::iterator it, pos;
  IntegerVector::iterator out_it;
  for(it = x.begin(), out_it = out.begin(); it != x.end(); ++it, ++out_it) {
    pos = std::upper_bound(breaks.begin(), breaks.end(), *it);
    *out_it = std::distance(breaks.begin(), pos);
  }
  return out - 1;
}

// Sinkhorn's algorithm for normalizing rows/columns of matrix
NumericMatrix sinkhorn(NumericMatrix x, int maxit, double eps) {
  int count = 0;
  int d = x.nrow();
  double eps1 = 1.0 / eps;
  double mc, Mc, mr, Mr;
  while (count < maxit) {
    count += 1;
    NumericVector cols(d);
    NumericVector rows(d);
    mc = 1;
    Mc = 1;
    mr = 1;
    Mr = 1;
    for (int i = 0; i < d; i++) {
      for (int j = 0; j < d; j++) {
        cols[i] += x(j, i);
      }
      if (cols[i] < mc) mc = cols[i];
      if (cols[i] > Mc) Mc = cols[i];
    }
    for (int i = 0; i < d; i++) {
      for (int j = 0; j < d; j++) {
        x(i, j) = x(i, j) / cols[j];
      }
    }
    for (int i = 0; i < d; i++) {
      for (int j = 0; j < d; j++) {
        rows[i] += x(i, j);
      }
      if (rows[i] < mr) mr = rows[i];
      if (rows[i] > Mr) Mr = rows[i];
    }
    for (int i = 0; i < d; i++) {
      for (int j = 0; j < d; j++) {
        x(i, j) = x(i, j) / rows[i];
      }
    }
    if (mc > eps1 & Mc < eps & mr > eps1 & Mr < eps) break;
    
  }
  return x / sum(x);
}

// naive sequential test by binning --------------------------------------------
//[[Rcpp::export]]
List sequential_histogram(
    NumericVector X,
    NumericVector Y,
    IntegerVector d,
    double n0
  ) {
  int N = X.length();
  int K = d.length();

  // containers for sequential ranks
  NumericVector FX(N); // unnormalized sequential ranks for X and Y
  NumericVector GY(N);
  ordered_set OX; // ordered sets, to compute the sequential ranks
  ordered_set OY;

  // compute sequential ranks
  for (int i = 0; i < N; i++) {
    OX.insert(X[i]);
    OY.insert(Y[i]);
    FX[i] = (OX.order_of_key(X[i]) + 1.0) / (1.0 + i);
    GY[i] = (OY.order_of_key(Y[i]) + 1.0) / (1.0 + i);
  }

  // positions in grid
  NumericVector rpos(N); // positions of ranks in the grid
  NumericVector spos(N);
  int r, s; // for specific positions (rpos, spos) in the loop to reduce indexing
  double r1, s1; // grid positions corresponding to r,s
  bool r_same_cell, s_same_cell; // to determine if derandomization is necessary
  double width; // width of grid cell (1/d[j])
  double rand_size, rand_size2; // width of randomization area (1/n for nth observation) and its square
  double d2; // d[j]^2
  double d2n0; // d[j]^2*n0
  double tmp; // to store weights in derandomization
  List out(2*K); // to store output

  for (int j = 0; j < K; j++) { // loop over all grid sizes
    // containers for histogram and martingale
    NumericMatrix bkldn(d[j], d[j]);
    NumericVector grid(d[j]);
    NumericVector Mn(N);

    // grid, grid size, and related
    d2 = pow(d[j], 2);
    d2n0 = d2 * n0;
    width = 1.0/d[j];
    for (int k = 0; k < d[j]; k++) {
      grid[k] = 1.0 * k/d[j];
    }
    
    // compute usual ranks for first d observations...
    for (int k = 0; k < d[j]; k++) {
      FX[k] = 1.0;
      GY[k] = 1.0;
      for (int l = 0; l < d[j]; l++) {
        if (X[l] < X[k]) {
          FX[k] = FX[k] + 1.0;
        }
        if (Y[l] < Y[k]) {
          GY[k] = GY[k] + 1.0;
        }
        // ... and include n0 artificial observations in each cell
        bkldn(k,l) = n0;
      }
      FX[k] = FX[k] / d[j];
      GY[k] = GY[k] / d[j];
    }
    rpos = find_interval(FX, grid);
    spos = find_interval(GY, grid);
    
    // initialize histogram and martingale (first d[j] observations)
    for (int k = 0; k < d[j]; k++) {
      // martingale equal to 1, no testing
      Mn[k] = 1.0;
      
      // update cells (derandomized strategy)
      r = rpos[k];
      s = spos[k];
      r1 = grid[r];
      s1 = grid[s];
      r_same_cell = (FX[k] - width >= r1);
      s_same_cell = (GY[k] - width >= s1);
      
      if (r_same_cell & s_same_cell) {
        // if whole area for randomization is in single cell, no derandomization is necessary
        bkldn(r,s) += 1.0;
      } else if ((!r_same_cell) & s_same_cell) {
        // otherwise compute the areas of rectangles in each relevant grid cell
        bkldn(r,s) += (FX[k] - r1) * d[j];
        bkldn(r-1,s) += (r1 - FX[k] + width) * d[j];
      } else if (r_same_cell & (!s_same_cell)) {
        bkldn(r,s) += (GY[k] - s1) * d[j];
        bkldn(r,s-1) += (s1 - GY[k] + width) * d[j];
      } else {
        bkldn(r,s) += (FX[k] - r1) * (GY[k] - s1) * d2;
        bkldn(r-1,s) += (r1 - FX[k] + width) * (GY[k] - s1) * d2;
        bkldn(r,s-1) += (FX[k] - r1) * (s1 - GY[k] + width) * d2;
        bkldn(r-1,s-1) += (r1 - FX[k] + width) * (s1 - GY[k] + width) * d2;
      }
    }

    // sequential update for observations d[j]+1,...,N
    for (int i = d[j]; i < N; i++) {
      rand_size = 1.0/(1.0 + i);
      r = rpos[i];
      s = spos[i];
      r1 = grid[r];
      s1 = grid[s];
      r_same_cell = (FX[i] - rand_size >= r1);
      s_same_cell = (GY[i] - rand_size >= s1);
      if (r_same_cell & s_same_cell) {
        Mn[i] = bkldn(r,s);
        bkldn(r,s) += 1.0;
      } else if (!r_same_cell & s_same_cell) {
        tmp = (FX[i] - r1) / rand_size;
        Mn[i] += tmp*bkldn(r,s);
        bkldn(r,s) += tmp;
        tmp = (r1 - FX[i] + rand_size) / rand_size;
        Mn[i] += tmp*bkldn(r-1,s);
        bkldn(r-1,s) += tmp;
      } else if (r_same_cell & !s_same_cell) {
        tmp = (GY[i] - s1) / rand_size;
        Mn[i] += tmp*bkldn(r,s);
        bkldn(r,s) += tmp;
        tmp = (s1 - GY[i] + rand_size) / rand_size;
        Mn[i] += tmp*bkldn(r,s-1);
        bkldn(r,s-1) += tmp;
      } else {
        rand_size2 = pow(rand_size, 2);
        tmp = (FX[i] - r1) * (GY[i] - s1) / rand_size2;
        Mn[i] += tmp*bkldn(r,s);
        bkldn(r,s) += tmp;
        tmp = (r1 - FX[i] + rand_size) * (GY[i] - s1) / rand_size2;
        Mn[i] += tmp*bkldn(r-1,s);
        bkldn(r-1,s) += tmp;
        tmp = (FX[i] - r1) * (s1 - GY[i] + rand_size) / rand_size2;
        Mn[i] += tmp*bkldn(r,s-1);
        bkldn(r,s-1) += tmp;
        tmp = (r1 - FX[i] + rand_size) * (s1 - GY[i] + rand_size) / rand_size2;
        Mn[i] += tmp*bkldn(r-1,s-1);
        bkldn(r-1,s-1) += tmp;
      }
      Mn[i] = Mn[i] * d2 / (i + d2n0);
    }
    out[j] = Mn;
    out[j + K] = bkldn/(d2n0 + N);
  }

  return(out);
}

// histogram with Sinkhorn normalization for fors and colums -------------------
//[[Rcpp::export]]
List sequential_histogram_sinkhorn(
    NumericVector X,
    NumericVector Y,
    IntegerVector d,
    double n0,
    double eps,
    int maxit
) {
  int N = X.length();
  int K = d.length();
  
  // containers for sequential ranks
  NumericVector FX(N); // unnormalized sequential ranks for X and Y
  NumericVector GY(N);
  ordered_set OX; // ordered sets, to compute the sequential ranks
  ordered_set OY;
  
  // compute sequential ranks
  for (int i = 0; i < N; i++) {
    OX.insert(X[i]);
    OY.insert(Y[i]);
    FX[i] = (OX.order_of_key(X[i]) + 1.0) / (1.0 + i);
    GY[i] = (OY.order_of_key(Y[i]) + 1.0) / (1.0 + i);
  }
  
  // positions in grid
  NumericVector rpos(N); // positions of ranks in the grid
  NumericVector spos(N);
  int r, s; // for specific positions (rpos, spos) in the loop to reduce indexing
  double r1, s1; // grid positions corresponding to r,s
  bool r_same_cell, s_same_cell; // to determine if derandomization is necessary
  double width; // width of grid cell (1/d[j])
  double rand_size, rand_size2; // width of randomization area (1/n for nth observation) and its square
  double d2; // d[j]^2
  double d2n0; // d[j]^2*n0
  double tmp; // to store weights in derandomization
  List out(2*K); // to store output
  
  for (int j = 0; j < K; j++) { // loop over all grid sizes
    // containers for histogram and martingale
    NumericMatrix bkldn(d[j], d[j]);
    NumericMatrix probs(d[j], d[j]);
    NumericVector grid(d[j]);
    NumericVector Mn(N);
    
    // grid, grid size, and related
    d2 = pow(d[j], 2);
    d2n0 = d2 * n0;
    width = 1.0/d[j];
    for (int k = 0; k < d[j]; k++) {
      grid[k] = 1.0 * k/d[j];
    }
    
    // compute usual ranks for first d observations...
    for (int k = 0; k < d[j]; k++) {
      FX[k] = 1.0;
      GY[k] = 1.0;
      for (int l = 0; l < d[j]; l++) {
        if (X[l] < X[k]) {
          FX[k] = FX[k] + 1.0;
        }
        if (Y[l] < Y[k]) {
          GY[k] = GY[k] + 1.0;
        }
        // ... and include n0 artificial observations in each cell
        bkldn(k,l) = n0;
      }
      FX[k] = FX[k] / d[j];
      GY[k] = GY[k] / d[j];
    }
    rpos = find_interval(FX, grid);
    spos = find_interval(GY, grid);
    
    // initialize histogram and martingale (first d[j] observations)
    for (int k = 0; k < d[j]; k++) {
      // martingale equal to 1, no testing
      Mn[k] = 1.0;
      
      // update cells (derandomized strategy)
      r = rpos[k];
      s = spos[k];
      r1 = grid[r];
      s1 = grid[s];
      r_same_cell = (FX[k] - width >= r1);
      s_same_cell = (GY[k] - width >= s1);
      
      if (r_same_cell & s_same_cell) {
        // if whole area for randomization is in single cell, no derandomization is necessary
        bkldn(r,s) += 1.0;
      } else if ((!r_same_cell) & s_same_cell) {
        // otherwise compute the areas of rectangles in each relevant grid cell
        bkldn(r,s) += (FX[k] - r1) * d[j];
        bkldn(r-1,s) += (r1 - FX[k] + width) * d[j];
      } else if (r_same_cell & (!s_same_cell)) {
        bkldn(r,s) += (GY[k] - s1) * d[j];
        bkldn(r,s-1) += (s1 - GY[k] + width) * d[j];
      } else {
        bkldn(r,s) += (FX[k] - r1) * (GY[k] - s1) * d2;
        bkldn(r-1,s) += (r1 - FX[k] + width) * (GY[k] - s1) * d2;
        bkldn(r,s-1) += (FX[k] - r1) * (s1 - GY[k] + width) * d2;
        bkldn(r-1,s-1) += (r1 - FX[k] + width) * (s1 - GY[k] + width) * d2;
      }
    }
    probs = clone(bkldn);
    probs = sinkhorn(probs, maxit, eps);
    
    // sequential update for observations d[j]+1,...,N
    for (int i = d[j]; i < N; i++) {
      rand_size = 1.0/(1.0 + i);
      r = rpos[i];
      s = spos[i];
      r1 = grid[r];
      s1 = grid[s];
      r_same_cell = (FX[i] - rand_size >= r1);
      s_same_cell = (GY[i] - rand_size >= s1);
      if (r_same_cell & s_same_cell) {
        Mn[i] = probs(r,s);
        bkldn(r,s) += 1.0;
      } else if (!r_same_cell & s_same_cell) {
        tmp = (FX[i] - r1) / rand_size;
        Mn[i] += tmp*probs(r,s);
        bkldn(r,s) += tmp;
        tmp = (r1 - FX[i] + rand_size) / rand_size;
        Mn[i] += tmp*probs(r-1,s);
        bkldn(r-1,s) += tmp;
      } else if (r_same_cell & !s_same_cell) {
        tmp = (GY[i] - s1) / rand_size;
        Mn[i] += tmp*probs(r,s);
        bkldn(r,s) += tmp;
        tmp = (s1 - GY[i] + rand_size) / rand_size;
        Mn[i] += tmp*probs(r,s-1);
        bkldn(r,s-1) += tmp;
      } else {
        rand_size2 = pow(rand_size, 2);
        tmp = (FX[i] - r1) * (GY[i] - s1) / rand_size2;
        Mn[i] += tmp*probs(r,s);
        bkldn(r,s) += tmp;
        tmp = (r1 - FX[i] + rand_size) * (GY[i] - s1) / rand_size2;
        Mn[i] += tmp*probs(r-1,s);
        bkldn(r-1,s) += tmp;
        tmp = (FX[i] - r1) * (s1 - GY[i] + rand_size) / rand_size2;
        Mn[i] += tmp*probs(r,s-1);
        bkldn(r,s-1) += tmp;
        tmp = (r1 - FX[i] + rand_size) * (s1 - GY[i] + rand_size) / rand_size2;
        Mn[i] += tmp*probs(r-1,s-1);
        bkldn(r-1,s-1) += tmp;
      }
      probs = clone(bkldn);
      probs = sinkhorn(probs, maxit, eps);
      Mn[i] = Mn[i] * d2;
    }
    out[j] = Mn;
    out[j + K] = probs;
  }
  
  return(out);
}

// histogram with Sinkhorn normalization for fors and colums, no deramdomization
//[[Rcpp::export]]
List sequential_histogram_sinkhorn_rand(
    NumericVector X,
    NumericVector Y,
    IntegerVector d,
    double n0,
    double eps,
    int maxit
) {
  int N = X.length();
  int K = d.length();
  
  // containers for sequential ranks
  NumericVector FX(N); // unnormalized sequential ranks for X and Y
  NumericVector GY(N);
  ordered_set OX; // ordered sets, to compute the sequential ranks
  ordered_set OY;
  
  // compute sequential ranks
  for (int i = 0; i < N; i++) {
    OX.insert(X[i]);
    OY.insert(Y[i]);
    FX[i] = (OX.order_of_key(X[i]) + 1.0) / (1.0 + i) - R::runif(0, 1.0 / (1.0 + i)); 
    GY[i] = (OY.order_of_key(Y[i]) + 1.0) / (1.0 + i) - R::runif(0, 1.0 / (1.0 + i));
  }
  
  // positions in grid
  NumericVector rpos(N); // positions of ranks in the grid
  NumericVector spos(N);
  int r, s; // for specific positions (rpos, spos) in the loop to reduce indexing
  double r1, s1; // grid positions corresponding to r,s
  bool r_same_cell, s_same_cell; // to determine if derandomization is necessary
  double width; // width of grid cell (1/d[j])
  double rand_size, rand_size2; // width of randomization area (1/n for nth observation) and its square
  double d2; // d[j]^2
  double d2n0; // d[j]^2*n0
  double tmp; // to store weights in derandomization
  List out(2*K); // to store output
  
  for (int j = 0; j < K; j++) { // loop over all grid sizes
    // containers for histogram and martingale
    NumericMatrix bkldn(d[j], d[j]);
    NumericMatrix probs(d[j], d[j]);
    NumericVector grid(d[j]);
    NumericVector Mn(N);
    
    for (int k = 0; k < d[j]; k++) {
      for (int l = 0; l < d[j]; l++) {
        bkldn(k, l) = n0;
      }
    }
    
    // grid, grid size, and related
    d2 = pow(d[j], 2);
    d2n0 = d2 * n0;
    width = 1.0/d[j];
    for (int k = 0; k < d[j]; k++) {
      grid[k] = 1.0 * k/d[j];
    }
    rpos = find_interval(FX, grid);
    spos = find_interval(GY, grid);
    
    for (int k = 0; k < d[j]; k++) {
      Mn[k] = 1.0;
      r = rpos[k];
      s = spos[k];
      bkldn(r,s) += 1.0;
    }
    probs = clone(bkldn);
    probs = sinkhorn(probs, maxit, eps);
    
    for (int i = d[j]; i < N; i++) {
      r = rpos[i];
      s = spos[i];
      Mn[i] = probs(r,s) * d2;
      bkldn(r,s) += 1.0;
      probs = clone(bkldn);
      probs = sinkhorn(probs, maxit, eps);
    }
    out[j] = Mn;
    out[j + K] = probs;
  }
  
  return(out);
}

// sequential KS-test by Shekhar & Ramdas --------------------------------------
//[[Rcpp::export]]
List sequential_ks(
  NumericVector X,
  NumericVector Y,
  NumericVector xgrid,
  NumericVector ygrid
) {
  int N = X.length();
  int Nx = xgrid.length();
  int Ny = ygrid.length();
  
  if (floor(N/2.0) != N/2.0) {
    // only even number of observations can be used
    N = N - 1;
    X.erase(N);
    Y.erase(N);
  }
  int N05 = N/2;
  
  // containers for lambda, z, a, and test function g (evaluated on grid)
  NumericVector lambda_plus(N05+2);
  NumericVector lambda_minus(N05+2);
  NumericVector a_plus(N05+1);
  a_plus[0] = 1.0;
  NumericVector a_minus(N05+1);
  a_minus[0] = 1.0;
  NumericVector z_plus(N05+1);
  NumericVector z_minus(N05+1);
  NumericVector Mn_plus(N05);
  NumericVector Mn_minus(N05);
  NumericMatrix err(Nx,Ny);
  IntegerVector pos_x = find_interval(X,xgrid);
  IntegerVector pos_y = find_interval(Y,ygrid);
  IntegerVector xp;
  IntegerVector yp;
  NumericVector xx;
  NumericVector yy;
  NumericVector g_plus = {xgrid[Nx-1], ygrid[Ny-1]};
  NumericVector g_minus = {xgrid[Nx-1], ygrid[Ny-1]};
  double add_to_err = 0;
  bool sign_equal = false;
  double max_plus = 0;
  double max_minus = 0;
  double v_plus;
  double v_minus;
  
  for (int i = 0; i < N05; i++) {
    // evaluate g_plus, g_minus at new observations and get Mn
    sign_equal = ((X[2*i] > X[2*i+1]) == (Y[2*i] > Y[2*i+1]));
    if (X[2*i] < X[2*i+1]) {
      xp = {pos_x[2*i], pos_x[2*i+1]};
      xx = {X[2*i], X[2*i+1]};
    } else {
      xp = {pos_x[2*i+1], pos_x[2*i]};
      xx = {X[2*i+1], X[2*i]};
    }
    if (Y[2*i] < Y[2*i+1]) {
      yp = {pos_y[2*i], pos_y[2*i+1]};
      yy = {Y[2*i], Y[2*i+1]};
    } else {
      yp = {pos_y[2*i+1], pos_y[2*i]};
      yy = {Y[2*i+1], Y[2*i]};
    }
    if (xx[0] <= g_plus[0] & xx[1] > g_plus[0] & 
        yy[0] <= g_plus[0] & yy[1] > g_plus[1]) {
      if (sign_equal) {
        v_plus = 1.0;
      } else {
        v_plus = -1.0;
      }
    } else {
      v_plus = 0.0;
    }
    if (xx[0] <= g_minus[0] & xx[1] > g_minus[0] &
        yy[0] <= g_minus[0] & yy[1] > g_minus[1]) {
      if (sign_equal) {
        v_minus = -1.0;
      } else {
        v_minus = 1.0;
      }
    } else {
      v_minus = 0.0;
    }
    Mn_plus[i] = 1.0-lambda_plus[i+1]*v_plus;
    Mn_minus[i] = 1.0-lambda_minus[i+1]*v_minus;
    
    // compute next lambdas
    z_plus[i+1] = v_plus/(1.0-v_plus*lambda_plus[i+1]);
    a_plus[i+1] = a_plus[i] + pow(z_plus[i+1], 2);
    lambda_plus[i+2] = lambda_plus[i+1] - 
      2.0/(2.0-log(3.0))*z_plus[i+1]/a_plus[i+1];
    if (lambda_plus[i+2] > 0.5) {
      lambda_plus[i+2] = 0.5;
    } else if (lambda_plus[i+2] < -0.5) {
      lambda_plus[i+2] = -0.5;
    }
    z_minus[i+1] = v_minus/(1.0-v_minus*lambda_minus[i+1]);
    a_minus[i+1] = a_minus[i] + pow(z_minus[i+1], 2);
    lambda_minus[i+2] = lambda_minus[i+1] - 
      2.0/(2.0-log(3.0))*z_minus[i+1]/a_minus[i+1];
    if (lambda_minus[i+2] > 0.5) {
      lambda_minus[i+2] = 0.5;
    } else if (lambda_minus[i+2] < -0.5) {
      lambda_minus[i+2] = -0.5;
    }
    
    // update err, g_plus, g_minus
    if ((xp[0] < xp[1] or xx[0] == xgrid[xp[0]]) & 
        (yp[0] < yp[1] or yy[0] == ygrid[yp[0]])) {
      if (sign_equal) {
        add_to_err = 1.0;
      } else {
        add_to_err = -1.0;
      }
      for (int k = xp[0]; k < xp[1]; k++) {
        for (int l = yp[0]; l < yp[1]; l++) {
          err(k,l) += add_to_err;
          if (err(k,l) > max_plus) {
            max_plus = err(k,l);
            g_plus = {xgrid[k],ygrid[l]};
          }
          if ((-err(k,l)) > max_minus) {
            max_minus = (-err(k,l));
            g_minus = {xgrid[k],ygrid[l]};
          }
        }
      }
    }
  }
  List out = List::create(Mn_plus, Mn_minus, err);
  return out;
}

// sequential binning approach with uniform marginals --------------------------

// functions for gradient descent
double loglikelihood(
    NumericVector x,
    NumericMatrix basis_mat,
    NumericVector q,
    double c0) {
  double out = 0.0;
  int d2 = x.length();
  int m = q.length();
  double tmp;
  for (int i = 0; i < d2; i++) {
    tmp = c0;
    for (int j = 0; j < m; j++) {
      tmp += basis_mat(i, j) * q[j];
    }
    if (tmp < 0) return R_PosInf;
    out -= x[i] * log(tmp);
  }
  return out;
}

NumericVector gradient_loglikelihood(
    NumericVector x,
    NumericMatrix basis_mat,
    NumericVector q,
    double c0) {
  int d2 = x.length();
  int m = q.length();
  NumericVector out(m);
  NumericVector bq(d2);
  for (int i = 0; i < d2; i++) {
    for (int j = 0; j < m; j++) {
      bq[i] += basis_mat(i, j) * q[j];
    }
  }
  for (int j = 0; j < m; j++) {
    for (int i = 0; i < d2; i++) {
      out[j] -= x[i] * basis_mat(i, j) / (c0 + bq[i]);
    }
  }
  return out;
}

NumericVector gradient_descent(
    NumericVector q0,
    NumericVector x,
    NumericMatrix B,
    double c0,
    double eps, // relative tolerance on log-likelihood for stopping
    double delta, // absolute tolerance for improvement in stepsize correction 
    int maxit, // maximum number of iterations
    int d // dimension
) {
  double current_ll;
  double stepsize;
  double current_improvement = 4.0 * eps;
  double next_improvement = 0.0;
  int iterations = 0;
  int m = q0.length();
  NumericVector gradient(m);
  NumericVector qnew(m);
  qnew = clone(q0);
  
  current_ll = loglikelihood(x, B, q0, c0);
  while ((current_improvement > eps) & (iterations < maxit)) {
    iterations++;
    stepsize = 1.0;
    gradient = gradient_loglikelihood(x, B, qnew, c0);
    current_improvement = 
      current_ll - loglikelihood(x, B, qnew - stepsize * gradient, c0);
    next_improvement = 
      current_ll - loglikelihood(x, B, qnew - stepsize * gradient / 2.0, c0);
    while ((next_improvement < -delta) or 
             (next_improvement > delta + current_improvement)) {
      stepsize = stepsize / 2.0;
      current_improvement = next_improvement;
      next_improvement = 
        current_ll - loglikelihood(x, B, qnew - stepsize * gradient / 2.0, c0);
    }
    qnew = qnew - stepsize * gradient;
    current_ll = loglikelihood(x, B, qnew, c0);
  }
  return qnew;
}

// sequential histogram with uniformity constraints on marginals
//[[Rcpp::export]]
List sequential_histogram_constrained(
    NumericVector X,
    NumericVector Y,
    List basis_mat,
    NumericVector c0,
    int dmax,
    double n0,
    double eps,
    double delta,
    int maxit
) {
  int N = X.length();
  int K = dmax;

  NumericVector FX(N);
  NumericVector GY(N);
  ordered_set OX;
  ordered_set OY;
  
  for (int i = 0; i < N; i++) {
    OX.insert(X[i]);
    OY.insert(Y[i]);
    FX[i] = (OX.order_of_key(X[i]) + 1.0) / (1.0 + i);
    GY[i] = (OY.order_of_key(Y[i]) + 1.0) / (1.0 + i);
  }
  
  NumericVector rpos(N);
  NumericVector spos(N);
  int r, s;
  double r1, s1;
  bool r_same_cell, s_same_cell;
  double width;
  double rand_size, rand_size2;
  double d2;
  double d2n0;
  double tmp;
  int d;
  int m;
  List out(2*K);

  for (int j = 0; j < dmax; j++) {
    d = pow(2, j + 1);
    NumericMatrix B = basis_mat[j];
    NumericMatrix bkldn(d, d);
    NumericMatrix probs(d, d);
    NumericVector grid(d);
    NumericVector Mn(N);
    
    d2 = pow(d, 2);
    m = d2 - 2*d + 1;
    d2n0 = d2 * n0;
    width = 1.0/d;
    for (int k = 0; k < d; k++) {
      grid[k] = 1.0 * k/d;
    }
    
    for (int k = 0; k < d; k++) {
      FX[k] = 1.0;
      GY[k] = 1.0;
      for (int l = 0; l < d; l++) {
        if (X[l] < X[k]) {
          FX[k] = FX[k] + 1.0;
        }
        if (Y[l] < Y[k]) {
          GY[k] = GY[k] + 1.0;
        }
        bkldn(k,l) = n0;
      }
      FX[k] = FX[k] / d;
      GY[k] = GY[k] / d;
    }
    rpos = find_interval(FX, grid);
    spos = find_interval(GY, grid);
    
    for (int k = 0; k < d; k++) {
      Mn[k] = 1.0;
      
      r = rpos[k];
      s = spos[k];
      r1 = grid[r];
      s1 = grid[s];
      r_same_cell = (FX[k] - width >= r1);
      s_same_cell = (GY[k] - width >= s1);
      
      if (r_same_cell & s_same_cell) {
        bkldn(r,s) += 1.0;
      } else if ((!r_same_cell) & s_same_cell) {
        bkldn(r,s) += (FX[k] - r1) * d;
        bkldn(r-1,s) += (r1 - FX[k] + width) * d;
      } else if (r_same_cell & (!s_same_cell)) {
        bkldn(r,s) += (GY[k] - s1) * d;
        bkldn(r,s-1) += (s1 - GY[k] + width) * d;
      } else {
        bkldn(r,s) += (FX[k] - r1) * (GY[k] - s1) * d2;
        bkldn(r-1,s) += (r1 - FX[k] + width) * (GY[k] - s1) * d2;
        bkldn(r,s-1) += (FX[k] - r1) * (s1 - GY[k] + width) * d2;
        bkldn(r-1,s-1) += (r1 - FX[k] + width) * (s1 - GY[k] + width) * d2;
      }
    }
    NumericVector q0;
    q0 = rep(0, m);
    q0 = gradient_descent(q0, bkldn, B, c0[j], eps, delta, maxit, d);
    for (int k = 0; k < d2; k++) {
      probs[k] = c0[j];
      for (int l = 0; l < m; l++) {
        probs[k] += B(k, l) * q0[l];
      }
    }

    for (int i = d; i < N; i++) {
      rand_size = 1.0/(1.0 + i);
      r = rpos[i];
      s = spos[i];
      r1 = grid[r];
      s1 = grid[s];
      r_same_cell = (FX[i] - rand_size >= r1);
      s_same_cell = (GY[i] - rand_size >= s1);
      if (r_same_cell & s_same_cell) {
        Mn[i] = probs(r,s);
        bkldn(r,s) += 1.0;
      } else if (!r_same_cell & s_same_cell) {
        tmp = (FX[i] - r1) / rand_size;
        Mn[i] += tmp*probs(r,s);
        bkldn(r,s) += tmp;
        tmp = (r1 - FX[i] + rand_size) / rand_size;
        Mn[i] += tmp*probs(r-1,s);
        bkldn(r-1,s) += tmp;
      } else if (r_same_cell & !s_same_cell) {
        tmp = (GY[i] - s1) / rand_size;
        Mn[i] += tmp*probs(r,s);
        bkldn(r,s) += tmp;
        tmp = (s1 - GY[i] + rand_size) / rand_size;
        Mn[i] += tmp*probs(r,s-1);
        bkldn(r,s-1) += tmp;
      } else {
        rand_size2 = pow(rand_size, 2);
        tmp = (FX[i] - r1) * (GY[i] - s1) / rand_size2;
        Mn[i] += tmp*probs(r,s);
        bkldn(r,s) += tmp;
        tmp = (r1 - FX[i] + rand_size) * (GY[i] - s1) / rand_size2;
        Mn[i] += tmp*probs(r-1,s);
        bkldn(r-1,s) += tmp;
        tmp = (FX[i] - r1) * (s1 - GY[i] + rand_size) / rand_size2;
        Mn[i] += tmp*probs(r,s-1);
        bkldn(r,s-1) += tmp;
        tmp = (r1 - FX[i] + rand_size) * (s1 - GY[i] + rand_size) / rand_size2;
        Mn[i] += tmp*probs(r-1,s-1);
        bkldn(r-1,s-1) += tmp;
      }
      Mn[i] = Mn[i] * d2;
      q0 = gradient_descent(q0, bkldn, B, c0[j], eps, delta, maxit, d);
      for (int k = 0; k < d2; k++) {
        probs[k] = c0[j];
        for (int l = 0; l < m; l++) {
          probs[k] += B(k, l) * q0[l];
        }
      }
    }
    out[j] = Mn;
    out[j + K] = clone(probs);
  }
  return out;
}

// sequential version of the BET by Kai Zhang ----------------------------------
//[[Rcpp::export]]
NumericMatrix sequential_bet(
    NumericVector X,
    NumericVector Y,
    double n0,
    List zero_one,
    IntegerVector d_groups
) {
  int N = X.length();
  int K = zero_one.length();
  int dmax = d_groups.length() - 2;
  
  NumericVector FX(N);
  NumericVector GY(N);
  ordered_set OX;
  ordered_set OY;
  
  for (int i = 0; i < N; i++) {
    OX.insert(X[i]);
    OY.insert(Y[i]);
    FX[i] = (OX.order_of_key(X[i]) + 1.0) / (1.0 + i);
    GY[i] = (OY.order_of_key(Y[i]) + 1.0) / (1.0 + i);
  }
  
  NumericVector rpos(N);
  NumericVector spos(N);
  int r, s, ind1, ind2, ind3, ind4; // positions in grid
  double r1, s1; // grid positions corresponding to r,s
  bool r_same_cell, s_same_cell;
  double width;
  double rand_size, rand_size2;
  double d, d2; 
  double tmp1, tmp2, tmp3, tmp4;
  NumericMatrix out(N, K);
  NumericMatrix probs(K, 2);
  
  for (int u = 0; u <= dmax; u++) { // loop over all orders of interaction terms
    d = pow(2.0, u + 1);
    NumericVector grid(d);
    d2 = pow(d, 2);
    width = 1.0/d;
    for (int k = 0; k < d; k++) {
      grid[k] = 1.0 * k/d;
    }
    for (int k = 0; k < d; k++) {
      // compute usual ranks for first d observations
      FX[k] = 1.0;
      GY[k] = 1.0;
      for (int l = 0; l < d; l++) {
        if (X[l] < X[k]) {
          FX[k] = FX[k] + 1.0;
        }
        if (Y[l] < Y[k]) {
          GY[k] = GY[k] + 1.0;
        }
      }
      FX[k] = FX[k] / d;
      GY[k] = GY[k] / d;
    }
    rpos = find_interval(FX, grid);
    spos = find_interval(GY, grid);
    
    for (int j = d_groups[u]; j < d_groups[u + 1]; j++) {
      // run sequential test for all interaction terms
      NumericMatrix zo = zero_one[j];
      probs(j, 0) = n0;
      probs(j, 1) = n0;

      for (int k = 0; k < d; k++) {
        // initialize martingale
        out(k, j) = 1.0;
        r = rpos[k];
        s = spos[k];
        r1 = grid[r];
        s1 = grid[s];
        r_same_cell = (FX[k] - width >= r1);
        s_same_cell = (GY[k] - width >= s1);

        if (r_same_cell & s_same_cell) {
          ind1 = zo(r, s);
          probs(j, ind1) += 1.0;
        } else if ((!r_same_cell) & s_same_cell) {
          ind1 = zo(r, s);
          probs(j, ind1) += (FX[k] - r1) * d;
          ind1 = zo(r - 1, s);
          probs(j, ind1) += (r1 - FX[k] + width) * d;
        } else if (r_same_cell & (!s_same_cell)) {
          ind1 = zo(r, s);
          probs(j, ind1) += (GY[k] - s1) * d;
          ind1 = zo(r, s - 1);
          probs(j, ind1) += (s1 - GY[k] + width) * d;
        } else {
          ind1 = zo(r, s);
          probs(j, ind1) += (FX[k] - r1) * (GY[k] - s1) * d2;
          ind1 = zo(r - 1, s);
          probs(j, ind1) += (r1 - FX[k] + width) * (GY[k] - s1) * d2;
          ind1 = zo(r, s - 1);
          probs(j, ind1) += (FX[k] - r1) * (s1 - GY[k] + width) * d2;
          ind1 = zo(r - 1, s - 1);
          probs(j, ind1) += (r1 - FX[k] + width) * (s1 - GY[k] + width) * d2;
        }
      }

      for (int i = d; i < N; i++) {
        rand_size = 1.0/(1.0 + i);
        r = rpos[i];
        s = spos[i];
        r1 = grid[r];
        s1 = grid[s];
        r_same_cell = (FX[i] - rand_size >= r1);
        s_same_cell = (GY[i] - rand_size >= s1);
        if (r_same_cell & s_same_cell) {
          ind1 = zo(r, s);
          out(i, j) = probs(j, ind1);
          probs(j, ind1) += 1.0;
        } else if (!r_same_cell & s_same_cell) {
          tmp1 = (FX[i] - r1) / rand_size;
          ind1 = zo(r, s);
          out(i, j) += tmp1*probs(j, ind1);
          tmp2 = (r1 - FX[i] + rand_size) / rand_size;
          ind2 = zo(r - 1, s);
          out(i, j) += tmp2*probs(j, ind2);
          probs(j, ind1) += tmp1;
          probs(j, ind2) += tmp2;
        } else if (r_same_cell & !s_same_cell) {
          tmp1 = (GY[i] - s1) / rand_size;
          ind1 = zo(r, s);
          out(i, j) += tmp1*probs(j, ind1);
          tmp2 = (s1 - GY[i] + rand_size) / rand_size;
          ind2 = zo(r, s - 1);
          out(i, j) += tmp2*probs(j, ind2);
          probs(j, ind1) += tmp1;
          probs(j, ind2) += tmp2;
        } else {
          rand_size2 = pow(rand_size, 2);
          tmp1 = (FX[i] - r1) * (GY[i] - s1) / rand_size2;
          ind1 = zo(r, s);
          out(i, j) += tmp1*probs(j, ind1);
          tmp2 = (r1 - FX[i] + rand_size) * (GY[i] - s1) / rand_size2;
          ind2 = zo(r - 1, s);
          out(i, j) += tmp2*probs(j, ind2);
          tmp3 = (FX[i] - r1) * (s1 - GY[i] + rand_size) / rand_size2;
          ind3 = zo(r, s - 1);
          out(i, j) += tmp3*probs(j, ind3);
          tmp4 = (r1 - FX[i] + rand_size) * (s1 - GY[i] + rand_size) / rand_size2;
          ind4 = zo(r - 1, s - 1);
          out(i, j) += tmp4*probs(j, ind4);
          probs(j, ind1) += tmp1;
          probs(j, ind2) += tmp2;
          probs(j, ind3) += tmp3;
          probs(j, ind4) += tmp4;
        }
        out(i, j) = out(i, j) * 2.0 / (2.0 * n0 + i);
      }
    }
  }
  
  return out;
}

// get first time of rejection with test martingale ----------------------------
//[[Rcpp::export]]
int get_first_rejection(NumericVector martingale, double threshold) {
  int ind = 0;
  martingale.push_back(R_PosInf);
  while (martingale[ind] < threshold) {
    ind++;
  }
  if (ind == martingale.length() - 1) {
    ind = 0;
  } else {
    ind++;
  }
  return ind;
}

