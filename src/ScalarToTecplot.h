#ifndef _SCALARTOTECPLOT_H_
#define _SCALARTOTECPLOT_H_

#include "Scalar.h"
#include <string>
#include <vector>
using std::string;

namespace ibpm {

bool ScalarToTecplot( const Scalar* var, string varName, string filename, string title );
    
bool ScalarToTecplot( vector<const Scalar*> varVec, vector<string> varNameVec, string filename, string title );
    
} // namespace ibpm

#endif /* _SCALARTOTECPLOT_H_ */
