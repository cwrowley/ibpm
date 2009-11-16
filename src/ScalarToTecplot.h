#ifndef _SCALARTOTECPLOT_H_
#define _SCALARTOTECPLOT_H_

#include "Scalar.h"
#include <string>
#include <vector>
using std::string;

namespace ibpm {
   
class ScalarToTecplot {
    
public:
    bool write( const Scalar& var, string varName, string filename, string title );
    
    bool write( const vector<Scalar>& varVec, vector<string> varNameVec, string filename, string title );
    
};

} // namespace ibpm

#endif /* _SCALARTOTECPLOT_H_ */
