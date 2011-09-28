#ifndef _OUTPUTENERGY_H_
#define _OUTPUTENERGY_H_

#include "Output.h"
#include <stdio.h>
#include <string>
using std::string;

namespace ibpm {

class OutputEnergy : public Output {
public:
    /// \brief Constructor
    /// \param[in] filename to which force data will be written.
    OutputEnergy(string filename);

    /// \brief Open the file for writing.
    /// If a file with the same name is already present, it is overwritten.
    /// Returns true if successful.
    bool init();

    /// \brief Close the file.
    /// Returns true if successful
    bool cleanup();

    /// \brief Write data to the energy file.
    bool doOutput(const State& x);
    
    /// \brief Write data to the energy file, making use of the baseflow.
    bool doOutput(const BaseFlow& q, const State& x);
    
private:
    string _filename;
    FILE* _fp;
};

} // namespace ibpm

#endif /* _OUTPUTENERGY_H_ */
