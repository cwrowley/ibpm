#ifndef _OUTPUT_H_
#define _OUTPUT_H_

class State;

/*!
\file Output.h
\class Output

\brief Abstract base class for output routines.

Subclasses provide a callback routine for writing output in various forms.

\author Clancy Rowley
\author $LastChangedBy: $
\date 21 Aug 2008
\date $LastChangedDate: $
\version $Revision: $
*/
class Output {
public:
    /// \brief Provide initialization, if needed (e.g. opening a file).
    /// Returns true if successful
    virtual bool init() { return true; }

    /// \brief Clean up, if needed (e.g. close a file)
    /// Returns true if successful
    virtual bool cleanup() { return false; }

    /// \brief Callback for performing the actual output.
    /// Returns true if successful; otherwise false.
    /// Pure virtual: must be defined by subclasses.
    virtual bool doOutput(const State& x) = 0;
    
};

#endif /* _OUTPUT_H_ */
