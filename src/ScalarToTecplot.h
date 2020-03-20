#ifndef _SCALARTOTECPLOT_H_
#define _SCALARTOTECPLOT_H_

#include "Scalar.h"
#include <string>
#include <vector>
using std::string;
using std::vector;

namespace ibpm {

bool ScalarToTecplot( const Scalar* var, string varName, string filename, string title );
    
bool ScalarToTecplot( vector<const Scalar*> varVec, vector<string> varNameVec, string filename, string title );
    
bool ScalarToTecplot( const Scalar* var, string varName, string filename, string title, int lev );

bool ScalarToTecplot( vector<const Scalar*> varVec, vector<string> varNameVec, string filename, string title, int lev );    
    
} // namespace ibpm

#endif /* _SCALARTOTECPLOT_H_ */
