#include "production.hpp"

/*! \file Interface.hpp
 */

/*! \class Interface
    \brief allow to read different type of files, (matrix and bed) 
 * and to product this files using c++ data
 * using as the interface with the user*/
class Interface {
public :

/*!
 * @brief parametrer for TCLAP and redirecting to the useful production method
 * @param argc
 * @param argv
*/
Interface(int argc, char** argv);

/*!
 * @brief parametrer for TCLAP and redirecting to the useful production method
*/
~Interface() {}

/*! @brief Methods in order to produce the matrix in EM_algorithm. using methods of the Matrix class
 * and of the production class
 * @param length : length of the motif
 * @param bed : name of the bedfile in order to produce the first matrix
 * @param fasta : name of the fasta file 
 * @param FastaOutput : name of the fasta file with only sequence of interest (by default, interest_sequence.fasta)
 * @param out : output name (by default, std::cout)
 */
void EM_algorith(int length,const string& bed, const string& fasta, const string& FastaOutput, const string& output);


};

