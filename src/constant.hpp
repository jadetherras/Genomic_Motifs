#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <tclap/CmdLine.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <list>
#include <limits.h>
#include <algorithm>

/*! \file constant.hpp
 */
 
/*! 
    \brief the max of nucleotide in a fasta line return by lire_fasta
 */
int const _MAX_FASTA_(1000);

/*! 
    \brief Maximum of iteration for EM algorithm
 */
int const _MAX_EM_(100);

/*! \class ProjectError
    \brief Class that heritate from runtime_error,
    *allow to throw special error, in order to implement
    * error gestion
 */
class ProjectError : public std::runtime_error {
public:
    ProjectError(const char *c, int v=0) : std::runtime_error(c), code(v) {}
    ProjectError(const std::string &s, int v=0) : std::runtime_error(s), code(v) {}
    int value() const {return code;}
protected:
    const int code;
};

#define _SIMULERR_(_N, _id) class _N : public ProjectError { \
    public: _N(const char *c) : ProjectError(c,_id) {} \
            _N(const std::string &s) : ProjectError(s,_id) {} };


/*! 
    \brief specifit errors for the code
 */
_SIMULERR_(TCLAP_ERROR, 10)
_SIMULERR_(FASTA_ERROR, 20)
_SIMULERR_(BEDFILE_ERROR, 30)
_SIMULERR_(MATRIXFILE_ERROR, 40)
_SIMULERR_(BEDGRAPH_ERROR, 50)
_SIMULERR_(OUTPUT_ERROR, 60)
_SIMULERR_(ZERO_ERROR, 70)
_SIMULERR_(TOOSHORT_ERROR, 80)
