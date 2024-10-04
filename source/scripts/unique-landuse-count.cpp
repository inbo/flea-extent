#include <Rcpp.h>
#include <unordered_set>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector count_unique_landuse(NumericVector x, size_t ni, size_t nw) {
    IntegerVector out(ni);
    
    // loop over cells
    size_t start = 0;
    for (size_t i = 0; i < ni; i++) {
        size_t end = start + nw;
        
        // Use an unordered_set to automatically handle uniqueness
        std::unordered_set<int> unique_values;
        bool all_na = true;
        
        // loop over the values of a window
        for (size_t j = start; j < end; j++) {
            if (!R_IsNA(x[j])) {
                // insert() only adds the value if it's not already in the set
                unique_values.insert(static_cast<int>(x[j]));
                all_na = false;
            }
        }
        
        // If all values in the window are NA, set the output to NA
        // Otherwise, set it to the count of unique values
        if (all_na) {
            out[i] = NA_INTEGER;
        } else {
            out[i] = unique_values.size();
        }
        
        start = end;
    }
    
    return out;
}
