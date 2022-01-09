#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;


// [[Rcpp::export]]
arma::colvec cc_mpd(arma::umat cdm, arma::mat dist, bool ab = false) {
    arma::colvec mpd(cdm.n_rows);
    arma::umat w;
    arma::uvec ind;
    if (ab) {
        for (int ii = 0; ii < cdm.n_rows; ii++) {
            w = cdm.row(ii).t() * cdm.row(ii);
            mpd(ii) = accu(w % dist)/accu(w);
        }
    }
    else {
        for (int ii = 0; ii < cdm.n_rows; ii++) {
            ind = find(cdm.row(ii));
            if (ind.size() > 1) {
                mpd(ii) = mean(nonzeros(trimatl(dist.submat(ind, ind))));
            }
            else
                mpd(ii) = NA_REAL;
        }
    }
    return mpd;
}


// [[Rcpp::export]]
double cc_mpd_v(arma::uvec cdm, arma::mat dist, bool ab = false) {
    double mpd;
    arma::umat w;
    arma::uvec ind;
    if (ab) {
        w = cdm * cdm.t();
        mpd = accu(w % dist)/accu(w);
    }
    else {
        ind = find(cdm);
        if (ind.size() > 1) {
            mpd = mean(nonzeros(trimatl(dist.submat(ind, ind))));
        }
        else
            mpd = NA_REAL;
    }
    return mpd;
}


// [[Rcpp::export]]
arma::colvec cc_rand_mpd(arma::umat cdm, arma::mat dist, bool ab = false) {
    cdm = shuffle(cdm, 1);
    return cc_mpd(cdm, dist, ab);
}


// [[Rcpp::export]]
arma::umat make_cdm(arma::umat id, int size = 8, int S = 500) {
    int n = std::floor(id.n_rows/size);
    arma::umat cdm(n*n, S, arma::fill::zeros);
    //arma::umat temp;
    int counter = 0;

    for (int ii = 0; ii < n; ii++) {
        for (int jj = 0; jj < n; jj++) {
            //temp = id.submat(ii*size, jj*size, (ii+1)*size-1, (jj+1)*size-1);

            for (int uu = 0; uu < size; uu++) {
                for (int ee = 0; ee < size; ee++) {
                    cdm(counter, id(ii*size + uu, jj*size + ee)-1)++;
                }
            }
            counter++;
        }
    }

    return cdm;
}


// [[Rcpp::export]]
arma::ucube make_cda(arma::mat id, int size = 8, int S = 500) {
    int n = std::floor(id.n_rows/size);
    arma::ucube cda(n, n, S, arma::fill::zeros);

    for (int ii = 0; ii < n; ii++) {
        for (int jj = 0; jj < n; jj++) {
            for (int uu = 0; uu < size; uu++) {
                for (int ee = 0; ee < size; ee++) {
                    cda(ii, jj, id(ii*size + uu, jj*size + ee)-1)++;
                }
            }
        }
    }
    return cda;
}

// [[Rcpp::export]]
arma::mat scaling_grain(arma::umat id, arma::mat dis, bool ab = true, int iter = 100) {

    int rmax = std::floor((id.n_rows/8 - 1)/2);
    int S = dis.n_rows;
    int ii, jj, r, itloc, uu;
    int n = rmax + 1;

    arma::umat cdm;
    arma::vec mpd;
    arma::mat null;
    arma::mat res(n, 2, arma::fill::zeros);
    int counter = 0;

    for (r = 0; r < rmax+1; r++) {
        //Rcout << r << " " << std::endl;

        cdm = make_cdm(id, (2*r + 1) * 8, S);

        if (cdm.n_rows > 50) itloc = 10; else itloc = iter;

        null.set_size(cdm.n_rows, itloc);

        for (ii = 0; ii < itloc; ii++) {
            //Rcout << ii << " ";
            null.col(ii) = cc_rand_mpd(cdm, dis, ab);
        }
        //Rcout << std::endl;

        mpd = cc_mpd(cdm, dis, ab);

        res(r, 0) = r;
        res(r, 1) = arma::mean( (mpd - arma::mean(null, 1))/arma::stddev(null, 0, 1) );
    }

    return res;
}

// [[Rcpp::export]]
arma::mat scaling_grain2(arma::umat id, arma::mat dis, bool ab = true, int iter = 100) {

    int rmax = std::floor((id.n_rows/8 - 1)/2);
    int S = dis.n_rows;
    int ii, jj, r, itloc, uu;
    int n = rmax + 1;

    arma::urowvec pool;
    arma::uvec index;
    arma::umat cdm;
    arma::vec mpd;
    arma::mat null, pool_dist;
    arma::mat res(n, 2, arma::fill::zeros);
    int counter = 0;

    for (r = 0; r < rmax+1; r++) {
        //Rcout << r << " " << std::endl;

        cdm = make_cdm(id, (2*r + 1) * 8, S);
        pool = sum(cdm);
        index = sum(cdm, 1);
        pool_dist = dis( find(pool), find(pool) );
        cdm = cdm( find(index), find(pool) );
        if (cdm.n_rows > 50) itloc = 10; else itloc = iter;

        null.set_size(cdm.n_rows, itloc);

        for (ii = 0; ii < itloc; ii++) {
            //Rcout << ii << " ";
            null.col(ii) = cc_rand_mpd(cdm, pool_dist, ab);
        }
        //Rcout << std::endl;

        mpd = cc_mpd(cdm, pool_dist, ab);

        res(r, 0) = r;
        res(r, 1) = arma::mean( (mpd - arma::mean(null, 1))/arma::stddev(null, 0, 1) );
    }

    return res;
}


// [[Rcpp::export]]
arma::mat scaling_pool(arma::ucube cda, arma::mat dis, bool ab = true, int iter = 100) {
    int K = cda.n_rows;

    int rmax = std::floor((K-1)/2);
    int ii, jj, r, nloc, uu;
    int n = 0;
    for (ii = 1; ii <= rmax; ii++) n += (K - 2*ii)*(K - 2*ii);

    arma::ucube temp;
    arma::uvec pool, target;
    arma::vec null(iter);
    arma::mat res(n, 5, arma::fill::zeros), pool_dist;
    int counter = 0;

    for (r = 1; r < rmax+1; r++) {
        nloc = K - 2*r;
        for (ii = 1; ii < nloc + 1; ii++) {
            for (jj = 1; jj < nloc + 1; jj++) {
                res(counter, 0) = r;
                res(counter, 1) = r+ii;
                res(counter, 2) = r+jj;

                temp = sum(sum(cda.tube(ii-1, jj-1, ii + 2*r - 1, jj + 2*r - 1)), 1);
                pool = temp(arma::span(0), arma::span(0), arma::span::all);
                target = cda(arma::span(r + ii - 1), arma::span(r + jj - 1), arma::span::all);
                target = target( find(pool) );
                pool_dist = dis( find(pool), find(pool));

                res(counter, 3) = cc_mpd_v(target, pool_dist, ab);

                for (uu = 0; uu < iter; uu++) {
                    null(uu) = cc_mpd_v(shuffle(target), pool_dist, ab);
                }
                res(counter, 4) = (res(counter, 3) - mean(null)) / stddev(null);


                counter++;
            }
        }
    }
    return res;
}

// [[Rcpp::export]]
arma::mat scaling_coherent(arma::ucube cda, arma::mat dis, bool ab = true, int iter = 100) {
    int K = cda.n_rows;

    int rmax = std::floor((K-3)/6);
    int ii, jj, r, nloc, uu;
    int n = 0;
    for (ii = 0; ii <= rmax; ii++) n += (K - 2*ii - 2*(2*ii + 1))*(K - 2*ii - 2*(2*ii + 1));

    arma::ucube temp;
    arma::uvec pool, target;
    arma::vec null(iter);
    arma::mat res(n, 5, arma::fill::zeros), pool_dist;
    int counter = 0;

    for (r = 0; r < rmax+1; r++) {
        nloc = K - 2*r - 2*(2*r + 1);

        for (ii = 1; ii < nloc + 1; ii++) {
            for (jj = 1; jj < nloc + 1; jj++) {
                res(counter, 0) = r;
                res(counter, 1) = r + ii + 2*r + 1;
                res(counter, 2) = r + jj + 2*r + 1;

                temp = sum(sum(cda.tube(ii-1, jj-1, ii + 6*r + 1, jj + 6*r + 1)), 1);
                pool = temp(arma::span(0), arma::span(0), arma::span::all);

                temp = sum(sum(cda.tube(2*r + ii, 2*r + jj, 4*r + ii, 4*r + jj)), 1);
                target = temp(arma::span(0), arma::span(0), arma::span::all);
                target = target( find(pool) );
                pool_dist = dis( find(pool), find(pool));

                res(counter, 3) = cc_mpd_v(target, pool_dist, ab);

                for (uu = 0; uu < iter; uu++) {
                    null(uu) = cc_mpd_v(shuffle(target), pool_dist, ab);
                }
                res(counter, 4) = (res(counter, 3) - mean(null)) / stddev(null);

                counter++;
            }
        }
    }
    return res;
}

