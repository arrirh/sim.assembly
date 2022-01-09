#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;



// [[Rcpp::export]]
double shannon(NumericVector n) {
    NumericVector p = n/sum(n);
    double out = 0;
    for (int ii = 0; ii < p.size(); ii++) {
        if (p[ii] > 0) out -= p[ii]*log(p[ii]);
    }
    return out;
}


// [[Rcpp::export]]
List createCom(NumericMatrix env,
               NumericVector niche,
               double nu = .001,
               int range = 3,
               double betaEnv = 1,
               double betaComp = 1,
               double betaAbund = 1,
               double sigma = 5) {

    int L = env.ncol(), richness = niche.size();
    NumericMatrix id(L,L);
    NumericVector n(richness);

    NumericVector nn(max(id), 0);

    for (int ii = 0; ii < L; ++ii) {
        for (int jj = 0; jj < L; ++jj) {
            id(ii,jj) = std::ceil(R::runif(0.00000001, richness));
            n[id(ii,jj)-1]++;
        }
    }

    NumericVector B, S, H, N;

    List out = List::create(_["env"] = env,
                            _["id"] = id,
                            _["n"] = wrap(n),
                            _["niche"] = niche,
                            _["L"] = L,
                            _["nu"] = nu,
                            _["range"] = range,
                            _["betaEnv"] = betaEnv,
                            _["betaComp"] = betaComp,
                            _["betaAbund"] = betaAbund,
                            _["sigma"] = sigma,
                            _["S"] = S,
                            _["H"] = H,
                            _["N"] = N,
                            _["t"] = 0);

    return out;
}

// [[Rcpp::export]]
List cc_simCom(List com, int itNum, int report) {

    int curt = com["t"];
    NumericVector S(itNum + curt), H(itNum + curt), N(itNum + curt);
    NumericVector oldS = com["S"], oldH = com["H"], oldN = com["N"];

    std::copy( oldS.begin(), oldS.end(), S.begin() ) ;
    std::copy( oldH.begin(), oldH.end(), H.begin() ) ;
    std::copy( oldN.begin(), oldN.end(), N.begin() ) ;

    NumericMatrix env = com["env"], id = com["id"];
    NumericVector n = com["n"], niche = com["niche"];
    double nu = com["nu"], betaEnv = com["betaEnv"], betaComp = com["betaComp"], betaAbund = com["betaAbund"], sigma = com["sigma"];
    int L = com["L"], range = com["range"];

    int richness = n.size();

    double pEnvMax = R::dnorm(0, 0, sigma, 0);
    //Rcout << "Max env fit: " << pEnvMax << std::endl;
    NumericMatrix overlap(richness, richness);
    for (int ii = 0; ii < richness; ii++)
    {
        for (int jj = 0; jj < richness; jj++)
        {
            overlap(ii, jj) = 2 * R::pnorm(-abs(niche[ii] - niche[jj])/2, 0, sigma, true, false);
        }
    }

    int betaType = 0;
    if (betaEnv > 0 & betaComp > 0) betaType = 1;
    else if (betaEnv > 0 & betaComp == 0) betaType = 2;
    else if (betaComp > 0 & betaEnv == 0) betaType = 3;
    else if (betaComp == 0 & betaEnv == 0) betaType = 4;

    bool flag = false;

    int counter = 0, siz = L*L, x, y, nx, ny;
    int fromX, fromY, toX, toY, K;
    std::map<int, int> counts;
    std::vector<double> pEnv, pComp, pAbund, pTot;

    double r, pd_new, bag;
    int ii, jj, uu, ee;

    for (ii = 0; ii < itNum; ii++) {
        //Rcout << "iteration: " << ii << std::endl;
        //Rcout << "initial id matrix: " << std::endl << id << std::endl;
        while (counter < siz) {
            // choose a position on a lattice
            x = (int)std::floor(R::runif(0,L));
            y = (int)std::floor(R::runif(0,L));

            r = R::runif(0,1);

            if (r > nu) {
                // choose new individual locally

                // compute table of local abundances

                fromX = x - range;
                fromY = y - range;
                toX = x + range;
                toY = y + range;
                if (fromX < 0) fromX = 0;
                if (fromY < 0) fromY = 0;
                if (toX > L - 1) toX = L-1;
                if (toY > L - 1) toY = L-1;
                K = (toX-fromX+1) * (toY-fromY+1);
                //Rcout << "Counting the range from " << fromX << " to " << toX << " and from " << fromY << " to " << toY << ":" << std::endl;
                counts.clear();
                for (uu = fromX; uu <= toX; uu++)
                {
                    for (ee = fromY; ee <= toY; ee++)
                    {
                        counts[id(uu, ee)]++;
                    }
                }

                // compute filters for species in table
                pEnv.clear();
                pComp.clear();
                pAbund.clear();
                pTot.clear();
                for(auto& it : counts) {
                    //Rcout << it.first << " " << it.second << " " << niche[it.first-1] << " ";
                    if (betaEnv > 0) {
                        pEnv.push_back(R::dnorm(env(x, y), niche[it.first-1], sigma, 0)/pEnvMax);
                    } else pEnv.push_back(1);
                    //Rcout << pEnv.back() << " ";
                    if (betaComp > 0) {
                        bag = 0;
                        for (auto& it2 : counts) {
                            bag += overlap(it.first-1, it2.first-1) * it2.second;
                        }
                        pComp.push_back(1 - bag/K);
                    } else pComp.push_back(1);
                    //Rcout << pComp.back() << " ";
                    pAbund.push_back((double)it.second/K);
                    //Rcout << pAbund.back() << " ";
                    //pTot.push_back(pEnv.back() * pComp.back() * pAbund.back());
                    switch (betaType) {
                    case 1:
                        pTot.push_back(std::exp(betaEnv * std::log(pEnv.back()) + betaComp * std::log(pComp.back()) + betaAbund * std::log(pAbund.back())));
                        break;
                    case 2:
                        pTot.push_back(std::exp(betaEnv * std::log(pEnv.back()) + betaAbund * std::log(pAbund.back())));
                        break;
                    case 3:
                        pTot.push_back(std::exp(betaComp * std::log(pComp.back()) + betaAbund * std::log(pAbund.back())));
                        break;
                    case 4:
                        pTot.push_back(pAbund.back());
                        break;
                    default:
                        Rcout << "Error: handling betaType" << std::endl;
                    }
                }

                // Final weights
                bag = 0;
                for (auto val : pTot) bag += val;
                //bag = std::accumulate(pTot.begin(), pTot.end(), 0);
                //Rcout << "Total pTot:" << bag << std::endl;
                for (int uu = 0; uu < pTot.size(); uu++) pTot[uu] = pTot[uu]/bag;
                //Rcout << "Final pTot:" << std::accumulate(pTot.begin(), pTot.end(), 0.0) << std::endl;

                // Choose parent and update the lattice
                r = R::runif(0,1);
                auto it = counts.begin();
                uu = 0;
                bag = pTot[uu];
                while (bag < r) {
                    uu++;
                    bag += pTot[uu];
                    it++;
                }
                //Rcout << "Parent choice:" << it->first << " / " << uu << std::endl;
                n[id(x,y)-1]--;
                id(x, y) = it->first;
                n[id(x,y)-1]++;

            }
            else {
                // choose new individual from pool
                n[id(x,y)-1]--;
                id(x, y) = std::ceil(R::runif(0.00000001, richness));
                n[id(x,y)-1]++;
            }

            counter++;

        }


        //B[ii] <- sum(lat)
        S[ii + curt] = sum(n > 0);
        N[ii + curt] = sum(n);
        H[ii + curt] = shannon(n);

        if ((ii+1) % report == 0) Rcout << ii+1 << " " << N[ii + curt] << " " << S[ii + curt] << " " << H[ii + curt] << std::endl;

        counter = 0;
    }

    List out = List::create(_["env"] = env,
                            _["id"] = id,
                            _["n"] = wrap(n),
                            _["niche"] = niche,
                            _["L"] = L,
                            _["nu"] = nu,
                            _["range"] = range,
                            _["betaEnv"] = betaEnv,
                            _["betaComp"] = betaComp,
                            _["betaAbund"] = betaAbund,
                            _["sigma"] = sigma,
                            _["S"] = S,
                            _["H"] = H,
                            _["N"] = N,
                            _["t"] = curt + itNum);
    return out;
}



// [[Rcpp::export]]
List cSAR(NumericMatrix id, NumericVector q) {
    int L = id.ncol();
    //Rcout << L << std::endl;
    int leg = (int)round(log(L)/log(2));
    //Rcout << leg << std::endl;
    int sizer = (leg - 4)*3 + 1;
    //Rcout << sizer << std::endl;
    int sp = (int)max(id);


    NumericVector siz(sizer+7), siz1 = NumericVector::create(2,3,4,6,8,10,13);
    std::copy(siz1.begin(), siz1.end(), siz.begin());
    double a = 0;
    for (int ii = 0; ii < sizer; ii++) {
        siz[ii+7] = std::ceil(pow(2.0, 4 + a));
        a += 0.333333333333333;
    }

    sizer = siz.size();

    NumericVector p(sp);

    NumericMatrix moment(sizer,q.size());

    int cur, sqNum, sqAb, qsize = q.size();

    double momBag = 0;

    for (int aa = 0; aa < sizer; aa++) {
        cur = (int)floor((double)L/siz[aa]);
        sqNum = 0;
        for (int ii = 0; ii < cur; ii++) {
            for (int jj = 0; jj < cur; jj++) {
                p.fill(0);
                sqAb = 0;
                for (int iii = ii*siz[aa]; iii < (ii+1)*siz[aa]; iii++) {
                    for (int jjj = jj*siz[aa]; jjj < (jj+1)*siz[aa]; jjj++) {
                        if (id(iii,jjj) > 0) {
                            p[id(iii,jjj)-1]++;
                            sqAb++;
                        }
                    }
                }
                if (sqAb > 0) {
                    //p = p/sqAb;
                    for (int uu = 0; uu < qsize; uu++) {
                        momBag = 0;
                        for (int oo = 0; oo < sp; oo++) {
                            if (p[oo] > 0)
                                momBag += pow(p[oo]/sqAb, q[uu]);
                        }
                        moment(aa,uu) += momBag;
                    }
                    sqNum++;
                }
            }
        }
        for (int uu = 0; uu < qsize; uu++)
            moment(aa,uu) /= sqNum;
    }

    List out = List::create(_["mom"] = moment,
                            _["A"] = pow(siz, 2),
                            _["q"] = q);


    return out;
}

/*

 // [[Rcpp::export]]
 List cc_simCom_opt1(List com, int itNum, int report) {

 int curt = com["t"];
 NumericVector S(itNum + curt), H(itNum + curt), N(itNum + curt);
 NumericVector oldS = com["S"], oldH = com["H"], oldN = com["N"];

 std::copy( oldS.begin(), oldS.end(), S.begin() ) ;
 std::copy( oldH.begin(), oldH.end(), H.begin() ) ;
 std::copy( oldN.begin(), oldN.end(), N.begin() ) ;

 NumericMatrix env = com["env"], id = com["id"];
 NumericVector n = com["n"], niche = com["niche"];
 double nu = com["nu"], betaEnv = com["betaEnv"], betaComp = com["betaComp"], betaAbund = com["betaAbund"];
 int L = com["L"], range = com["range"];

 int richness = n.size();

 double pEnvMax = R::dnorm(0, 0, 5, 0);
 //Rcout << "Max env fit: " << pEnvMax << std::endl;
 NumericMatrix overlap(richness, richness);
 for (int ii = 0; ii < richness; ii++)
 {
 for (int jj = 0; jj < richness; jj++)
 {
 overlap(ii, jj) = 2 * R::pnorm(-abs(niche[ii] - niche[jj])/2, 0, 5, true, false);
 }
 }


 bool flag = false;

 int counter = 0, siz = L*L, x, y, nx, ny;
 int fromX, fromY, toX, toY, K;
 std::map<int, int> counts;
 std::vector<double> pEnv, pComp, pAbund, pTot;

 double r, pd_new, bag;
 int ii, jj, uu, ee;

 if (betaEnv > 0 & betaComp > 0) {
 Rcout << "betaEnv > 0 & betaComp > 0" << std::endl;
 for (ii = 0; ii < itNum; ii++) {
 //Rcout << "iteration: " << ii << std::endl;
 //Rcout << "initial id matrix: " << std::endl << id << std::endl;
 while (counter < siz) {
 // choose a position on a lattice
 x = (int)std::floor(R::runif(0,L));
 y = (int)std::floor(R::runif(0,L));

 r = R::runif(0,1);

 if (r > nu) {
 // choose new individual locally

 // compute table of local abundances

 fromX = x - range;
 fromY = y - range;
 toX = x + range;
 toY = y + range;
 if (fromX < 0) fromX = 0;
 if (fromY < 0) fromY = 0;
 if (toX > L - 1) toX = L-1;
 if (toY > L - 1) toY = L-1;
 K = (toX-fromX+1) * (toY-fromY+1);
 counts.clear();
 for (uu = fromX; uu <= toX; uu++)
 {
 for (ee = fromY; ee <= toY; ee++)
 {
 counts[id(uu, ee)]++;
 }
 }

 // compute filters for species in table
 pEnv.clear();
 pComp.clear();
 pAbund.clear();
 pTot.clear();
 for(auto& it : counts) {
 pEnv.push_back(R::dnorm(env(x, y), niche[it.first-1], 5, 0)/pEnvMax);
 bag = 0;
 for (auto& it2 : counts) {
 bag += overlap(it.first-1, it2.first-1) * it2.second;
 }
 pComp.push_back(1 - bag/K);
 pAbund.push_back((double)it.second/K);
 pTot.push_back(std::exp(betaEnv * std::log(pEnv.back()) + betaComp * std::log(pComp.back()) + betaAbund * std::log(pAbund.back())));
 }

 // Final weights
 bag = 0;
 for (auto val : pTot) bag += val;
 for (int uu = 0; uu < pTot.size(); uu++) pTot[uu] = pTot[uu]/bag;

 // Choose parent and update the lattice
 r = R::runif(0,1);
 auto it = counts.begin();
 uu = 0;
 bag = pTot[uu];
 while (bag < r) {
 uu++;
 bag += pTot[uu];
 it++;
 }
 n[id(x,y)-1]--;
 id(x, y) = it->first;
 n[id(x,y)-1]++;

 }
 else {
 // choose new individual from pool
 n[id(x,y)-1]--;
 id(x, y) = std::ceil(R::runif(0.00000001, richness));
 n[id(x,y)-1]++;
 }

 counter++;

 }

 S[ii + curt] = sum(n > 0);
 N[ii + curt] = sum(n);
 H[ii + curt] = shannon(n);

 if ((ii+1) % report == 0) Rcout << ii+1 << " " << N[ii + curt] << " " << S[ii + curt] << " " << H[ii + curt] << std::endl;

 counter = 0;
 }

 } else if (betaComp == 0 & betaEnv > 0) {
 Rcout << "betaComp == 0 & betaEnv > 0" << std::endl;
 for (ii = 0; ii < itNum; ii++) {
 //Rcout << "iteration: " << ii << std::endl;
 //Rcout << "initial id matrix: " << std::endl << id << std::endl;
 while (counter < siz) {
 // choose a position on a lattice
 x = (int)std::floor(R::runif(0,L));
 y = (int)std::floor(R::runif(0,L));

 r = R::runif(0,1);

 if (r > nu) {
 // choose new individual locally

 // compute table of local abundances

 fromX = x - range;
 fromY = y - range;
 toX = x + range;
 toY = y + range;
 if (fromX < 0) fromX = 0;
 if (fromY < 0) fromY = 0;
 if (toX > L - 1) toX = L-1;
 if (toY > L - 1) toY = L-1;
 K = (toX-fromX+1) * (toY-fromY+1);
 counts.clear();
 for (uu = fromX; uu <= toX; uu++)
 {
 for (ee = fromY; ee <= toY; ee++)
 {
 counts[id(uu, ee)]++;
 }
 }

 // compute filters for species in table
 pEnv.clear();
 //pComp.clear();
 pAbund.clear();
 pTot.clear();
 for(auto& it : counts) {
 pEnv.push_back(R::dnorm(env(x, y), niche[it.first-1], 5, 0)/pEnvMax);
 pAbund.push_back((double)it.second/K);
 pTot.push_back(std::exp(betaEnv * std::log(pEnv.back()) + betaAbund * std::log(pAbund.back())));
 }

 // Final weights
 bag = 0;
 for (auto val : pTot) bag += val;
 for (int uu = 0; uu < pTot.size(); uu++) pTot[uu] = pTot[uu]/bag;

 // Choose parent and update the lattice
 r = R::runif(0,1);
 auto it = counts.begin();
 uu = 0;
 bag = pTot[uu];
 while (bag < r) {
 uu++;
 bag += pTot[uu];
 it++;
 }
 n[id(x,y)-1]--;
 id(x, y) = it->first;
 n[id(x,y)-1]++;

 }
 else {
 // choose new individual from pool
 n[id(x,y)-1]--;
 id(x, y) = std::ceil(R::runif(0.00000001, richness));
 n[id(x,y)-1]++;
 }

 counter++;

 }

 S[ii + curt] = sum(n > 0);
 N[ii + curt] = sum(n);
 H[ii + curt] = shannon(n);

 if ((ii+1) % report == 0) Rcout << ii+1 << " " << N[ii + curt] << " " << S[ii + curt] << " " << H[ii + curt] << std::endl;

 counter = 0;
 }



 } else if (betaEnv == 0 & betaComp > 0) {
 Rcout << "betaEnv == 0 & betaComp > 0" << std::endl;
 for (ii = 0; ii < itNum; ii++) {
 //Rcout << "iteration: " << ii << std::endl;
 //Rcout << "initial id matrix: " << std::endl << id << std::endl;
 while (counter < siz) {
 // choose a position on a lattice
 x = (int)std::floor(R::runif(0,L));
 y = (int)std::floor(R::runif(0,L));

 r = R::runif(0,1);

 if (r > nu) {
 // choose new individual locally

 // compute table of local abundances

 fromX = x - range;
 fromY = y - range;
 toX = x + range;
 toY = y + range;
 if (fromX < 0) fromX = 0;
 if (fromY < 0) fromY = 0;
 if (toX > L - 1) toX = L-1;
 if (toY > L - 1) toY = L-1;
 K = (toX-fromX+1) * (toY-fromY+1);
 counts.clear();
 for (uu = fromX; uu <= toX; uu++)
 {
 for (ee = fromY; ee <= toY; ee++)
 {
 counts[id(uu, ee)]++;
 }
 }

 // compute filters for species in table
 //pEnv.clear();
 pComp.clear();
 pAbund.clear();
 pTot.clear();
 for(auto& it : counts) {
 bag = 0;
 for (auto& it2 : counts) {
 bag += overlap(it.first-1, it2.first-1) * it2.second;
 }
 pComp.push_back(1 - bag/K);
 pAbund.push_back((double)it.second/K);
 pTot.push_back(std::exp(betaComp * std::log(pComp.back()) + betaAbund * std::log(pAbund.back())));
 }

 // Final weights
 bag = 0;
 for (auto val : pTot) bag += val;
 for (int uu = 0; uu < pTot.size(); uu++) pTot[uu] = pTot[uu]/bag;

 // Choose parent and update the lattice
 r = R::runif(0,1);
 auto it = counts.begin();
 uu = 0;
 bag = pTot[uu];
 while (bag < r) {
 uu++;
 bag += pTot[uu];
 it++;
 }
 n[id(x,y)-1]--;
 id(x, y) = it->first;
 n[id(x,y)-1]++;

 }
 else {
 // choose new individual from pool
 n[id(x,y)-1]--;
 id(x, y) = std::ceil(R::runif(0.00000001, richness));
 n[id(x,y)-1]++;
 }

 counter++;

 }

 S[ii + curt] = sum(n > 0);
 N[ii + curt] = sum(n);
 H[ii + curt] = shannon(n);

 if ((ii+1) % report == 0) Rcout << ii+1 << " " << N[ii + curt] << " " << S[ii + curt] << " " << H[ii + curt] << std::endl;

 counter = 0;
 }
 } else if (betaEnv == 0 & betaComp == 0) {
 Rcout << "betaEnv == 0 & betaComp == 0" << std::endl;
 for (ii = 0; ii < itNum; ii++) {
 //Rcout << "iteration: " << ii << std::endl;
 //Rcout << "initial id matrix: " << std::endl << id << std::endl;
 while (counter < siz) {
 // choose a position on a lattice
 x = (int)std::floor(R::runif(0,L));
 y = (int)std::floor(R::runif(0,L));

 r = R::runif(0,1);

 if (r > nu) {
 // choose new individual locally

 // compute table of local abundances

 fromX = x - range;
 fromY = y - range;
 toX = x + range;
 toY = y + range;
 if (fromX < 0) fromX = 0;
 if (fromY < 0) fromY = 0;
 if (toX > L - 1) toX = L-1;
 if (toY > L - 1) toY = L-1;
 K = (toX-fromX+1) * (toY-fromY+1);
 counts.clear();
 for (uu = fromX; uu <= toX; uu++)
 {
 for (ee = fromY; ee <= toY; ee++)
 {
 counts[id(uu, ee)]++;
 }
 }

 // compute filters for species in table
 //pEnv.clear();
 //pComp.clear();
 //pAbund.clear();
 pTot.clear();
 for(auto& it : counts) {
 pTot.push_back((double)it.second/K);
 }

 // Final weights

 // Choose parent and update the lattice
 r = R::runif(0,1);
 auto it = counts.begin();
 uu = 0;
 bag = pTot[uu];
 while (bag < r) {
 uu++;
 bag += pTot[uu];
 it++;
 }
 n[id(x,y)-1]--;
 id(x, y) = it->first;
 n[id(x,y)-1]++;

 }
 else {
 // choose new individual from pool
 n[id(x,y)-1]--;
 id(x, y) = std::ceil(R::runif(0.00000001, richness));
 n[id(x,y)-1]++;
 }

 counter++;

 }

 S[ii + curt] = sum(n > 0);
 N[ii + curt] = sum(n);
 H[ii + curt] = shannon(n);

 if ((ii+1) % report == 0) Rcout << ii+1 << " " << N[ii + curt] << " " << S[ii + curt] << " " << H[ii + curt] << std::endl;

 counter = 0;
 }

 }





 List out = List::create(_["env"] = env,
 _["id"] = id,
 _["n"] = wrap(n),
 _["niche"] = niche,
 _["L"] = L,
 _["nu"] = nu,
 _["range"] = range,
 _["betaEnv"] = betaEnv,
 _["betaComp"] = betaComp,
 _["betaAbund"] = betaAbund,
 _["S"] = S,
 _["H"] = H,
 _["N"] = N,
 _["t"] = curt + itNum);
 return out;
 }


 // [[Rcpp::export]]
 List cc_simCom(List com, int itNum, int report) {

 int curt = com["t"];
 NumericVector S(itNum + curt), H(itNum + curt), N(itNum + curt);
 NumericVector oldS = com["S"], oldH = com["H"], oldN = com["N"];

 std::copy( oldS.begin(), oldS.end(), S.begin() ) ;
 std::copy( oldH.begin(), oldH.end(), H.begin() ) ;
 std::copy( oldN.begin(), oldN.end(), N.begin() ) ;

 NumericMatrix env = com["env"], id = com["id"];
 NumericVector n = com["n"], niche = com["niche"];
 double nu = com["nu"], betaEnv = com["betaEnv"], betaComp = com["betaComp"], betaAbund = com["betaAbund"];
 int L = com["L"], range = com["range"];

 int richness = n.size();

 double pEnvMax = R::dnorm(0, 0, 5, 0);
 //Rcout << "Max env fit: " << pEnvMax << std::endl;
 NumericMatrix overlap(richness, richness);
 for (int ii = 0; ii < richness; ii++)
 {
 for (int jj = 0; jj < richness; jj++)
 {
 overlap(ii, jj) = 2 * R::pnorm(-abs(niche[ii] - niche[jj])/2, 0, 5, true, false);
 }
 }


 bool flag = false;

 int counter = 0, siz = L*L, x, y, nx, ny;
 int fromX, fromY, toX, toY, K;
 std::map<int, int> counts;
 std::vector<double> pEnv, pComp, pAbund, pTot;

 double r, pd_new, bag;
 int ii, jj, uu, ee;

 for (ii = 0; ii < itNum; ii++) {
 //Rcout << "iteration: " << ii << std::endl;
 //Rcout << "initial id matrix: " << std::endl << id << std::endl;
 while (counter < siz) {
 // choose a position on a lattice
 x = (int)std::floor(R::runif(0,L));
 y = (int)std::floor(R::runif(0,L));

 r = R::runif(0,1);

 if (r > nu) {
 // choose new individual locally

 // compute table of local abundances

 fromX = x - range;
 fromY = y - range;
 toX = x + range;
 toY = y + range;
 if (fromX < 0) fromX = 0;
 if (fromY < 0) fromY = 0;
 if (toX > L - 1) toX = L-1;
 if (toY > L - 1) toY = L-1;
 K = (toX-fromX+1) * (toY-fromY+1);
 //Rcout << "Counting the range from " << fromX << " to " << toX << " and from " << fromY << " to " << toY << ":" << std::endl;
 counts.clear();
 for (uu = fromX; uu <= toX; uu++)
 {
 for (ee = fromY; ee <= toY; ee++)
 {
 counts[id(uu, ee)]++;
 }
 }

 // compute filters for species in table
 pEnv.clear();
 pComp.clear();
 pAbund.clear();
 pTot.clear();
 for(auto& it : counts) {
 //Rcout << it.first << " " << it.second << " " << niche[it.first-1] << " ";
 pEnv.push_back(R::dnorm(env(x, y), niche[it.first-1], 5, 0)/pEnvMax);
 //Rcout << pEnv.back() << " ";
 bag = 0;
 for (auto& it2 : counts) {
 bag += overlap(it.first-1, it2.first-1) * it2.second;
 }
 pComp.push_back(1 - bag/K);
 //Rcout << pComp.back() << " ";
 pAbund.push_back((double)it.second/K);
 //Rcout << pAbund.back() << " ";
 //pTot.push_back(pEnv.back() * pComp.back() * pAbund.back());
 pTot.push_back(std::exp(betaEnv * std::log(pEnv.back()) + betaComp * std::log(pComp.back()) + betaAbund * std::log(pAbund.back())));
 //Rcout << pTot.back() << std::endl;
 }

 // Final weights
 bag = 0;
 for (auto val : pTot) bag += val;
 //bag = std::accumulate(pTot.begin(), pTot.end(), 0);
 //Rcout << "Total pTot:" << bag << std::endl;
 for (int uu = 0; uu < pTot.size(); uu++) pTot[uu] = pTot[uu]/bag;
 //Rcout << "Final pTot:" << std::accumulate(pTot.begin(), pTot.end(), 0.0) << std::endl;

 // Choose parent and update the lattice
 r = R::runif(0,1);
 auto it = counts.begin();
 uu = 0;
 bag = pTot[uu];
 while (bag < r) {
 uu++;
 bag += pTot[uu];
 it++;
 }
 //Rcout << "Parent choice:" << it->first << " / " << uu << std::endl;
 n[id(x,y)-1]--;
 id(x, y) = it->first;
 n[id(x,y)-1]++;

 }
 else {
 // choose new individual from pool
 n[id(x,y)-1]--;
 id(x, y) = std::ceil(R::runif(0.00000001, richness));
 n[id(x,y)-1]++;
 }

 counter++;
 }
 //Rcout << "final id matrix: " << std::endl << id << std::endl;
 //Rcout << "final n vector: " << std::endl << n << std::endl;
 //Rcout << "final emptyID vector: " << std::endl << idList << std::endl;
 //check(id, n);


 //B[ii] <- sum(lat)
 S[ii + curt] = sum(n > 0);
 N[ii + curt] = sum(n);
 H[ii + curt] = shannon(n);

 if ((ii+1) % report == 0) Rcout << ii+1 << " " << N[ii + curt] << " " << S[ii + curt] << " " << H[ii + curt] << std::endl;

 counter = 0;
 }


 //     com["t"] = itNum + as<int>com["t"];
 //     //com$B <- c(com$B,B)
 //     com["S"] = com["S"].insert(com["S"].end(), S.begin(), s.end());
 //     com["H"] = com["H"].insert(com["H"].end(), H.begin(), H.end());
 //     com["N"] = com["N"].insert(com["N"].end(), N.begin(), N.end());
 //     com["id"] = id;
 //     com["n"] = n;
 //     com["pd"] = pd;


 List out = List::create(_["env"] = env,
                         _["id"] = id,
                         _["n"] = wrap(n),
                         _["niche"] = niche,
                         _["L"] = L,
                         _["nu"] = nu,
                         _["range"] = range,
                         _["betaEnv"] = betaEnv,
                         _["betaComp"] = betaComp,
                         _["betaAbund"] = betaAbund,
                         _["S"] = S,
                         _["H"] = H,
                         _["N"] = N,
                         _["t"] = curt + itNum);
                         return out;
                         }



// [[Rcpp::export]]
List c_createCom(int L = 64, int initS = 100, double nu = .01, int dist = 1) {

    NumericMatrix id(L,L);
    //NumericMatrix lat(L,L);

    NumericVector n(initS);

    for (int ii = 0; ii < L; ++ii) {
        for (int jj = 0; jj < L; ++jj) {
            if (R::runif(0,1) < 0.5)
            {
                //lat(ii,jj) = 0;
                id(ii,jj) = std::ceil(R::runif(0.00000001,initS));
                n[id(ii,jj)-1]++;
            }
            else
            {
                id(ii,jj) = 0;
            }
        }
    }

    NumericVector pd = runif(initS, 0, 0.5);
    NumericVector B, S, H, N;

    List out = List::create(_["id"] = id,
                            //_["lat"] = lat,
                            _["n"] = wrap(n),
                            _["pd"] = pd,
                            _["L"] = L,
                            _["nu"] = nu,
                            _["dist"] = dist,
                            _["B"] = B,
                            _["S"] = S,
                            _["H"] = H,
                            _["N"] = N,
                            _["t"] = 0);

    return out;
}




// [[Rcpp::export]]
double wmean(NumericVector n, NumericVector p) {
    double out = sum(n*p)/sum(n);
    return out;
}

// [[Rcpp::export]]
void check(NumericMatrix id, NumericVector n) {
    NumericVector nn(max(id), 0);
    int L = id.ncol();

    for (int ii = 0; ii < L; ++ii) {
        for (int jj = 0; jj < L; ++jj) {
            n[id(ii,jj)-1]++;
        }
    }

    if (shannon(n) != shannon(nn)) throw std::range_error("boom");
}


// [[Rcpp::export]]
List cc_simCom01(List com, int itNum, int report) {

    int curt = com["t"];
    NumericVector S(itNum + curt), H(itNum + curt), N(itNum + curt);
    NumericVector oldS = com["S"], oldH = com["H"], oldN = com["N"];

    std::copy( oldS.begin(), oldS.end(), S.begin() ) ;
    std::copy( oldH.begin(), oldH.end(), H.begin() ) ;
    std::copy( oldN.begin(), oldN.end(), N.begin() ) ;

    NumericMatrix id = com["id"];
    NumericVector n = com["n"], pd = com["pd"];
    double nu = com["nu"];
    int L = com["L"], dist = com["dist"];

    NumericVector idList;
    for (int ii = 0; ii < n.size(); ii++) {
        if (n[ii] == 0) idList.push_back(ii);
    }

    //Rcout << "list of zeros: " << idList << std::endl;

    bool flag = false;

    int counter = 0, siz = L*L, x, y, nx, ny;

    double r, pd_new;


    //IntegerMatrix events(itNum,6);

    for (int ii = 0; ii < itNum; ii++) {
        //Rcout << "iteration: " << ii << std::endl;
        //Rcout << "initial id matrix: " << std::endl << id << std::endl;
        while (counter < siz) {
            // choose a position on a lattice
            x = (int)std::floor(R::runif(0,L));
            y = (int)std::floor(R::runif(0,L));
            // check if it is occupied

            while (id(x,y) == 0) {
                x = (int)std::floor(R::runif(0,L));
                y = (int)std::floor(R::runif(0,L));

            }

            r = R::runif(0,1);
            if (r < pd[id(x,y)-1]) {
                //Rcout << "death at: " << x << " " << y << std::endl;
                // death
                n[id(x,y)-1]--;
                if (n[id(x,y)-1] == 0)
                    idList.push_back(id(x,y)-1);
                id(x,y) = 0;

                counter++;
                //events(ii,0)++;
            }
            else if (r < pd[id(x,y)-1] + 0.5) {
                //Rcout << "reproduction from: " << x << " " << y << std::endl;
                //events(ii,1)++;
                if (R::runif(0,1) < nu)
                {
                    //Rcout << "speciation" << std::endl;
                    // speciation
                    pd_new = R::runif(0, 0.5);
                    flag = true;
                }
                else
                {
                    pd_new = pd[id(x,y)-1];
                }

                // choose site of competition
                nx = x + (int)round(R::rnorm(0,dist));
                ny = y + (int)round(R::rnorm(0,dist));
                while (nx == x && ny == y) {
                    nx = x + (int)round(R::rnorm(0,dist));
                    ny = y + (int)round(R::rnorm(0,dist));

                }
                nx = nx % L;
                ny = ny % L;
                if (nx < 0) nx += L;
                if (ny < 0) ny += L;

                //Rcout << "competition at: " << nx << " " << ny << std::endl;

                if (id(nx,ny) == 0) {
                    //events(ii,2)++;
                    //Rcout << "site is empty" << std::endl;
                    if (flag)
                    {
                        // new species establishes
                        if (idList.size() == 0)
                        {
                            // new id
                            n.push_back(1);
                            pd.push_back(pd_new);
                            id(nx,ny) = n.size();
                        }
                        else
                        {
                            // empty id
                            n[idList[idList.size()-1]] = 1;
                            pd[idList[idList.size()-1]] = pd_new;
                            id(nx,ny) = idList[idList.size()-1]+1;
                            idList.erase(idList.size()-1);
                        }
                        flag = false;
                    }
                    else
                    {

                        // newborn establishes
                        id(nx,ny) = id(x,y);
                        n[id(x,y)-1]++;
                    }
                }
                else if (pd_new > pd[id(nx,ny)-1])
                {
                    //events(ii,3)++;
                    //Rcout << "successfull replace" << std::endl;
                    if (flag)
                    {
                        // new species establishes
                        if (idList.size() == 0)
                        {
                            // new id
                            n.push_back(1);
                            pd.push_back(pd_new);

                            n[id(nx,ny)-1]--;
                            if (n[id(nx,ny)-1] == 0) idList.push_back(id(nx,ny)-1);
                            id(nx,ny) = n.size();
                        }
                        else
                        {
                            // empty id
                            n[id(nx,ny)-1]--;
                            if (n[id(nx,ny)-1] == 0) idList.push_back(id(nx,ny)-1);

                            n[idList[idList.size()-1]] = 1;
                            pd[idList[idList.size()-1]] = pd_new;

                            id(nx,ny) = idList[idList.size()-1]+1;
                            idList.erase(idList.size()-1);
                        }
                        flag = false;
                    }
                    else
                    {
                        // newborn establishes
                        n[id(nx,ny)-1]--;
                        if (n[id(nx,ny)-1] == 0) idList.push_back(id(nx,ny)-1);

                        id(nx,ny) = id(x,y);
                        n[id(x,y)-1]++;
                    }
                }
                else
                {
                    flag = false;
                    //events(ii,4)++;
                    //Rcout << "Failed replace" << std::endl;
                }

                counter++;
            } else
            {
                counter++;
                //events(ii,5)++;
                //Rcout << "Nothing happened" << std::endl;
            }
        }
        //Rcout << "final id matrix: " << std::endl << id << std::endl;
        //Rcout << "final n vector: " << std::endl << n << std::endl;
        //Rcout << "final emptyID vector: " << std::endl << idList << std::endl;
        //check(id, n);

        //B[ii] <- sum(lat)
        S[ii + curt] = sum(n > 0);
        N[ii + curt] = sum(n);
        H[ii + curt] = shannon(n);

        if (ii % report == 0) Rcout << ii << " " << N[ii + curt] << " " << S[ii + curt] << " " << H[ii + curt] << std::endl;

        counter = 0;
    }

    //     com["t"] = itNum + as<int>com["t"];
    //     //com$B <- c(com$B,B)
    //     com["S"] = com["S"].insert(com["S"].end(), S.begin(), s.end());
    //     com["H"] = com["H"].insert(com["H"].end(), H.begin(), H.end());
    //     com["N"] = com["N"].insert(com["N"].end(), N.begin(), N.end());
    //     com["id"] = id;
    //     com["n"] = n;
    //     com["pd"] = pd;


    List out = List::create(_["id"] = id,
                            _["n"] = wrap(n),
                            _["pd"] = pd,
                            _["L"] = L,
                            _["nu"] = nu,
                            _["dist"] = dist,
                            _["S"] = S,
                            _["H"] = H,
                            _["N"] = N,
                            _["t"] = curt + itNum);


    return out;
}



// [[Rcpp::export]]
double genpf(double F) {
    return 1/R::runif(1/F, 1) + R::runif(0, 1);
}


*/
