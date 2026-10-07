#include <random>
#include <chrono>

#include <Rcpp.h>

using namespace Rcpp;

 // Use R random number generator to control seeds
 // Otherwise use the following:
 // static std::random_device rd;
 // initialize Mersennes' twister using rd to generate the seed
 // static std::mt19937 gen{rd()};
 // std::uniform_real_distribution<double> dist(0, 1);

 // calls an arbitrary R function `f` on a numeric vector `x` and returns its result.
 NumericVector callFunction(NumericVector x, Function f) {
   NumericVector res = f(x);
   return res;
 }

 // Draws a single random number from Uniform(0,1). Separated to control from R.
 double drawUniform01() {
   double randomnumber = R::runif(0,1);
   return randomnumber;
 }

 // Draws a single value from `choiceset`, with sampling probabilities given
 // by `probavec`. Separated to control from R.
 double drawSample(IntegerVector choiceset, NumericVector probavec) {
   double generatedrandom = Rcpp::sample(choiceset, 1, false, probavec)[0];
   return generatedrandom;
 }

 //' indexPatterns
 //'
 //' @param U IntegerVector with factor values corresponding to the patterns of Partially Identifying Variables (PIVs) observed among records in the concerned source
 //' @param nvals integer for the total number of possible patterns (among all sources)
 //'
 //' @return List: for each pattern in key, a vector of indices of the records having the pattern in the concerned source
 //' @export
 // [[Rcpp::export]]
 List indexPatterns(IntegerVector U, int nvals)
 {
   // Create one empty "bucket" (integer vector) for each pattern
   List out(nvals);
   for (int i = 0; i < nvals; i++)
   {
     IntegerVector tmp;
     out[i] = tmp;
   }
   // Which value does observation i have
   // Assign it to the appropriate "bucket" in the list
   for (int i = 0; i < U.length(); i++)
   {
     IntegerVector tmpX = out[U[i]-1];
     tmpX.push_back(i+1);
     out[U[i]-1] = tmpX;
   }
   return out;
 }

 //' pairPatterns
 //'
 //' @param A List with for each pattern in value, count of the records having the pattern in the concerned source
 //' @param B List with for each pattern in value, count of the records having the pattern in the concerned source
 //' @param nvals integer for the total number of possible patterns (among all sources)
 //'
 //' @return IntegerMatrix: nrow=nbr of potential links, ncol = 2, indicates the indices (of records from A, records from B) of records with matching patterns to be considered in the likelihood as potential links. It represents the cartesian product of both lists passed as parameters to represent all possible linked pairs of records.
 //' @export
 // [[Rcpp::export]]
 IntegerMatrix pairPatterns(List A, List B, int nvals)
 {
   // Figure out how many candidate pairs there will be in total
   int ntotal = 0;
   for (int k = 0; k < nvals; k++)
   {
     // Check a bucket to see whether and how many combos are possible
     IntegerVector tmpA = A[k];
     IntegerVector tmpB = B[k];
     ntotal += tmpA.length() * tmpB.length();
   }
   IntegerMatrix tmpC(ntotal,2);
   int counter = 0;
   // For each pattern (index in list A or B), give all possible pairs of indices carrying that pattern
   for (int k = 0; k < nvals; k++)
   {
     IntegerVector tmpA = A[k];
     IntegerVector tmpB = B[k];
     for (int i = 0; i < tmpA.length(); i++)
     {
       for (int j = 0; j < tmpB.length(); j++)
       {
         tmpC(counter,0) = tmpA[i];
         tmpC(counter,1) = tmpB[j];
         counter += 1;
       }
     }
   }
   return tmpC;
 }

 //' pasteIntoPattern
 //'
 //' @param A IntegerMatrix with values to form patterns
 //'
 //' @return CharacterVector: "pasteIntoPattern(A)" is a faster version in Rcpp of "do.call(paste, c(as.data.frame(A), list(sep="_")))"
 //' @export
 // [[Rcpp::export]]
 CharacterVector pasteIntoPattern(IntegerMatrix A)
 {
   int i = 0, j = 0;
   int sz = A.nrow();
   CharacterVector res(sz);
   for (std::ostringstream oss; i < sz; i++, oss.str(""))
   {
     oss << A(i,0);
     for (j = 1; j < A.ncol(); j++)
     {
       oss << "_" << A(i,j);
     }
     res[i] = oss.str();
   }
   return res;
 }

 std::map<int, std::set<int>> _DeltaMap;

 //' initDeltaMap
 //'
 //' @return void: Initialise the cpp map _DeltaMap representing the sparse linkage matrix Delta.
 //' @export
 // [[Rcpp::export]]
 void initDeltaMap()
 {
   _DeltaMap = {};
 }

 //' Deltafind()
 //'
 //' @return IntegerMatrix: find the indices of the elements in the cpp map _DeltaMap representing the sparse linkage matrix Delta.
 //' @export
 // [[Rcpp::export]]
 IntegerMatrix Deltafind()
 {
   // Count total number of links
   int length = 0;
   for(const auto & rowColSet: _DeltaMap)
   {
     length += rowColSet.second.size();
   }
   IntegerMatrix found(length, 2);
   int index = 0;
   // Flatten the map
   for(const auto & rowColSet: _DeltaMap)
   {

     for (int c : rowColSet.second)
     {
       found(index,0)   = rowColSet.first;
       found(index++,1) = c;
     }
   }
   if(index != length)
   {
     Rcpp::stop("Something went wrong creating Δ.");
   }
   return found;
 }

 //' sampleD
 //'
 //' Performs one Gibbs-sampling sweep over all candidate record pairs in `S`,
 //' proposing for each pair to either:
 //'  - add a link, by comparing the log-likelihood of the "linked" vs "not linked" state, or
 //'  - remove an existing link (if the pair is currently linked)
 //' Each proposal is accepted stochastically.
 //'
 //' @param S IntegerMatrix where each row correspond to the indices (from source A and source B) of records for which the true values matches (representing the potential links)
 //' @param LLA NumericVector gives the likelihood contribution of each non linked record from A
 //' @param LLB NumericVector gives the likelihood contribution of each non linked record from B
 //' @param LLL NumericVector gives the likelihood contribution of each potential linked records (from select)
 //' @param gamma NumericVector repeats the value of the parameter gamma (proportion of linked records) number of potential linked records (nrow of S) times
 //' @param loglik double for the value of the current complete log likelihood of the model
 //' @param nlinkrec integer for the current number of linked records
 //' @param sumRowD A LogicalVector vector indicating, for each row of the linkage matrix, i.e. for each record in the smallest file A, whether the record has a link in B or not.
 //' @param sumColD A LogicalVector vector indicating, for each column of the linkage matrix, i.e. for each record in the largest file B, whether the record has a link in A or not.
 //'
 //' @return
 //' List:
 //' - new set of links
 //' - new sumRowD
 //' - new sumColD
 //' - new value of the complete log likelihood
 //' - new number fo linked records
 //' @export
 // [[Rcpp::export]]
 List sampleD(const IntegerMatrix & S,
              const NumericVector & LLA,
              const NumericVector & LLB,
              const NumericVector & LLL,
              const NumericVector & gamma,
              double loglik,
              int nlinkrec,
              LogicalVector & sumRowD,
              LogicalVector & sumColD)
 {
   // Iterate over every candidate pair
   for (int q = 0; q < S.nrow(); q++)
   {
     int i = S(q,0)-1;
     int j = S(q,1)-1;
     // If non-linked -> possibly becomes linked
     if((sumRowD(i)==false) && (sumColD(j)==false))
     {
       // Log-likelihood if this pair were linked:
       // remove the two "unlinked" contributions, add the "linked" contribution,
       // and adjust for the bipartite-matching
       double loglikNew = loglik
       // Comparisons
       - LLB(j) - LLA(i)
       + LLL(q)
       // Bipartite matching
       - log(1-gamma(i)) + log(gamma(i))
       - log(LLB.length() - nlinkrec);
       double sumlogdensity = log(1 + exp(loglik-loglikNew)) + loglikNew;
       double pswitch = exp(loglikNew - sumlogdensity);
       // Random number smaller than pswitch -> generates binomial value
       // Accept the new link with probability pswitch
       bool link = drawUniform01()  < pswitch;
       if(link)
       {
         loglik = loglikNew;
         _DeltaMap[i].insert(j);
         sumRowD(i) = true;
         sumColD(j) = true;
         nlinkrec = nlinkrec + 1;
       }
     }else if(_DeltaMap.count(i) && _DeltaMap.at(i).count(j))
     {
       // If linked -> possibly becomes non-linked (analog)
       double loglikNew = loglik
       // Comparisons
       + LLB(j) + LLA(i)
       - LLL(q)
       // Bipartite matching
       + log(1-gamma(i)) - log(gamma(i))
       + log(LLB.length() - nlinkrec+1);
       double sumlogdensity = log(1 + exp(loglik-loglikNew)) + loglikNew;
       double pswitch = exp(loglikNew - sumlogdensity);
       // Accept breaking the link with probability pswitch
       bool nolink = drawUniform01()  < pswitch;
       if(nolink)
       {
         loglik = loglikNew;
         if(_DeltaMap.count(i) && _DeltaMap[i].count(j))
           _DeltaMap[i].erase(j);
         sumRowD(i) = false;
         sumColD(j) = false;
         nlinkrec = nlinkrec - 1;
       }
     }
   }
   // Convert the updated sparse linkage map into the (idxA, idxB) matrix form
   IntegerMatrix links=Deltafind();
   // Return to R
   List ret;
   ret["links"] = links;
   ret["sumRowD"] = sumRowD;
   ret["sumColD"] = sumColD;
   ret["loglik"] = loglik;
   ret["nlinkrec"] = nlinkrec;
   return ret;
 }

 //' sampleNL
 //'
 //' Samples the latent "true" value of a single PIV for non-linked records

 //' @param G IntegerVector of registered values for a certain Partially Identifying Variable (PIV) for non linked records
 //' @param eta NumericVector parameter for the distribution of the PIV concerned
 //' @param phi NumericVector parameter for the registration errors for the PIV concerned
 //'
 //' @return IntegerVector: of latent true values underlying G
 //' @export
 // [[Rcpp::export]]
 IntegerVector sampleNL(IntegerVector G, NumericVector eta, NumericVector phi)
 {
   IntegerVector H(G.length());
   // Number of possible values
   int nval = eta.length();
   // Create the possible values to sample from
   IntegerVector choice_set = seq_len(nval);
   // Possible registration errors
   double pMissing = phi[1];
   double pTypo = (1-pMissing) * (1-phi[0]) / (eta.length()-1);
   double pAgree = (1-pMissing) * phi[0];
   // Iterate over all elements
   for(int i = 0; i < G.length(); i++)
   {
     // Create a vector indicating P(Registered=X|True)
     // First value is for the missings
     // What happens if missing:
     if(G(i) == 0)
     {
       // Registered value is missing
       // Mistakes equally likely
       // Sample from the PIV distribution
       H(i) = drawSample(choice_set, eta);
     }else
     {
       // Registered value is observed
       NumericVector help1(nval, pTypo);
       help1(G(i)-1) = pAgree;
       // Check the joint probability of having the registered and true value
       NumericVector prob = eta * help1;
       H(i) = drawSample(choice_set, prob);
     }
   }
   return H;
 }

 //' sampleL
 //'
 //' Samples the latent joint "true" values of a single PIV for linked records
 //'
 //' @param GA IntegerVector of registered values for a certain Partially Identifying Variable (PIV) for linked records from A
 //' @param GB IntegerVector of registered values for a certain PIV for linked records from B
 //' @param survivalpSameH NumericVector of probabilities that the concerned PIV values coincide between file A and file B
 //' @param choice_set IntegerMatrix of 2 columns (for A and for B) with possible joint true values underlying GA and GB
 //' @param choice_equal IntegerVector of booleans indicating whether the 2 true values (from A and B) in the choice set are equal
 //' @param nval integer for the number of unique values in the PIV concerned
 //' @param phikA NumericVector parameter for the registration errors in A for the PIV concerned
 //' @param phikB NumericVector parameter for the registration errors in B for the PIV concerned
 //' @param eta NumericVector parameter for the distribution of the PIV concerned
 //'
 //' @return IntegerVector: of indices from the joint latent true values choice set underlying GA and GB
 //' @export
 // [[Rcpp::export]]
 IntegerVector sampleL(IntegerVector GA, IntegerVector GB, NumericVector survivalpSameH,
                       IntegerMatrix choice_set, IntegerVector choice_equal,
                       int nval, NumericVector phikA, NumericVector phikB, NumericVector eta)
 {
   IntegerVector H(GA.length());
   int size_choice_set = choice_set.nrow();
   IntegerVector choice_index = seq_len(size_choice_set);
   // Possible actions for file A
   double pMissingA = phikA[1];
   double pTypoA = (1-pMissingA) * (1-phikA[0]) / (nval-1);
   double pAgreeA = (1-pMissingA) * phikA[0];
   // Possible actions for file B
   double pMissingB = phikB[1];
   double pTypoB = (1-pMissingB) * (1-phikB[0]) / (nval-1);
   double pAgreeB = (1-pMissingB) * phikB[0];
   // Iterate over all matches
   for(int i = 0; i < GA.length(); i++)
   {
     // Prob that both TRUE values are the same
     double pSameH = survivalpSameH[i];
     // Dummy Vectors
     // Define P(Hb|Ha) (taking account of dynamic PIV that change with time)
     NumericVector probH(size_choice_set, pSameH);
     // Define P(Ga|Ha) and P(Gb|Hb) (taking account of registration errors)
     NumericVector helpA(size_choice_set, pTypoA);
     NumericVector helpB(size_choice_set, pTypoB);
     for(int j = 0; j < size_choice_set; j++)
     {
       // Prob to observe HB=b|HA=a where a!=b
       if(choice_equal(j)==0)
       {
         probH(j) = (1-pSameH)/(nval-1); //
       }
       if(GA(i) == choice_set(j,0))
         helpA(j) = pAgreeA;
       if(GB(i) == choice_set(j,1))
         helpB(j) = pAgreeB;
     }
     // There are four options possible, to generate true values
     // Both missing
     if(GA(i)==0 && GB(i)==0)
     {
       NumericVector prob = eta * probH;
       H(i) = drawSample(choice_index, prob);
     }else if(GA(i)>0 && GB(i)==0)
     {
       // Joint probability to have the registered and true value
       NumericVector prob = eta * probH * helpA;
       H(i) = drawSample(choice_index, prob);
     }else if(GA(i)==0 && GB(i)>0)
     {
       // Joint probability to have the registered and true value
       NumericVector prob = eta * probH * helpB;
       H(i) = drawSample(choice_index, prob);
     }else if(GA(i)>0 && GB(i)>0)
     {
       // None missing
       // Create vectors indicating P(Registered=X|True)
       // Joint probability to have the registered and true value
       NumericVector prob = eta * probH * helpA * helpB;
       H(i) = drawSample(choice_index, prob);
     }
   }
   return H;
 }

 //' ExpandGrid
 //'
 //' Compute the cartesian product of vectors
 //'
 //' @param vector1 first IntegerVector of values to compute the cartesian product
 //' @param vector2 second IntegerVector of values to compute the cartesian product
 //'
 //' @return IntegerMatrix: of 2 columns with the cartesian product of vec1 and vec2
 //' @export
 // [[Rcpp::export]]
 IntegerMatrix ExpandGrid(IntegerVector vec1, IntegerVector vec2) {
   int n = vec1.size();
   int m = vec2.size();
   IntegerMatrix result(n * m, 2);
   for (int i = 0; i < n; i++) {
     for (int j = 0; j < m; j++) {
       result(i * m + j, 0) = vec1[i];
       result(i * m + j, 1) = vec2[j];
     }
   }
   return result;
 }

 //' generateSequence
 //'
 //' @param n integer superior to 1
 //'
 //' @return IntegerVector: with n values from 1 to n
 //' @export
 // [[Rcpp::export]]
 IntegerVector generateSequence(int n) {
   IntegerVector result(n);
   for (int i = 0; i < n; i++) {
     result[i] = i + 1;
   }
   return result;
 }

 //' sampleH
 //'
 //' @param nA IntegerVector of dimensions of registered values of the Partially Identifying Variables (PIVs) in A
 //' @param nB IntegerVector of dimensions of registered values of the PIVs in B
 //' @param links IntegerMatrix of 2 columns with the indices of the linked records
 //' @param survivalpSameH NumericMatrix with for each PIV the probability that true values coincide (if stable: filled with 1)
 //' @param pivs_stable LogicalVector indicating for each PIV whether it is stable of not (if not we expect survivalpSameH for that same element to not be filled with 1 but with lower values)
 //' @param pivsA List ith registered data from A
 //' @param pivsB List with registered data from B
 //' @param nvalues IntegerVector with number of unique values of each PIV
 //' @param nonlinkedA LogicalVector indicating for all records in A whether they are linked or not
 //' @param nonlinkedB LogicalVector indicating for all records in B whether they are linked or not
 //' @param eta List parameters of the PIVs distributions
 //' @param phi List parameters of the PIVs registration errors
 //'
 //' @return List:
 //' - truePIVsA, true values underlying data in A
 //' - truePIVsB, true values underlying data in B
 //' @export
 // [[Rcpp::export]]
 List sampleH(IntegerVector nA, IntegerVector nB, IntegerMatrix links, NumericMatrix survivalpSameH, LogicalVector pivs_stable, List pivsA, List pivsB, IntegerVector nvalues, LogicalVector nonlinkedA, LogicalVector nonlinkedB, List eta, List phi)
 {
   // Output matrices of true PIV values, same shape as the registered data
   IntegerMatrix truepivsA(nA[0],nA[1]);
   IntegerMatrix truepivsB(nB[0],nB[1]);
   int nphi = 2;
   for (int k = 0; k < nvalues.length(); k++) {
     IntegerVector truepivsA_k = truepivsA(_,k);
     IntegerVector truepivsB_k = truepivsB(_,k);
     NumericVector eta_k = eta[k];
     NumericVector phi_k = phi[k];
     // Unpack the registration-error parameters for this PIV:
     // phi_k = [agree_A, agree_B, missing_A, missing_B]
     NumericVector phi_k_A(nphi);
     phi_k_A[0] = phi_k[0];
     phi_k_A[1] = phi_k[2];
     NumericVector phi_k_B(nphi);
     phi_k_B[0] = phi_k[1];
     phi_k_B[1] = phi_k[3];
     IntegerVector pivsA_k = pivsA[k];
     IntegerVector pivsB_k = pivsB[k];
     IntegerVector pivsA_k_L = pivsA_k[links(_,0)];
     IntegerVector pivsB_k_L = pivsB_k[links(_,1)];
     IntegerVector pivsA_k_NL = pivsA_k[nonlinkedA];
     IntegerVector pivsB_k_NL = pivsB_k[nonlinkedB];
     // MArginally sample true values for non-linked pairs
     IntegerVector truepivsA_k_NL = sampleNL(pivsA_k_NL, eta_k, phi_k_A);
     IntegerVector truepivsB_k_NL = sampleNL(pivsB_k_NL, eta_k, phi_k_B);
     truepivsA_k[nonlinkedA] = truepivsA_k_NL;
     truepivsB_k[nonlinkedB] = truepivsB_k_NL;
     IntegerMatrix choice_set;
     IntegerVector choice_equal;
     NumericVector eta_choice;
     if (links.nrow()>0)
     {
       IntegerVector values = generateSequence(nvalues[k]);
       choice_set = ExpandGrid(values, values);
       choice_equal = choice_set(_,1) == choice_set(_,0);
       eta_choice = eta_k[choice_set(_,1) - 1];
       NumericVector survivalpSameH_k = survivalpSameH(_,k);
       // Jointly sample true values for linked pairs
       IntegerVector out = sampleL(pivsA_k_L, pivsB_k_L, survivalpSameH_k, choice_set, choice_equal, nvalues[k], phi_k_A, phi_k_B, eta_choice);
       IntegerVector choice_set_A = choice_set(_,0);
       IntegerVector choice_set_B = choice_set(_,1);
       truepivsA_k[links(_,0)] = choice_set_A[out - 1];
       truepivsB_k[links(_,1)] = choice_set_B[out - 1];
     }
     truepivsA(_,k) = truepivsA_k;
     truepivsB(_,k) = truepivsB_k;
   }
   List ret;
   ret["truepivsA"] = truepivsA;
   ret["truepivsB"] = truepivsB;
   return ret;
 }
