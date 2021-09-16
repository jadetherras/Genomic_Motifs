#include "constant.hpp"

using namespace std;

/*! \file traitement.hpp
 */

/*! \class Traitement
    \brief Class for traitement of fasta (read, inversion, conversion) and for the ouput
    *allow to do all traitement important for the two major class
 */
 

class Traitement
{
	public :
/*!
    \brief Constructor initializes the class Traitement with the sequence
    \param fasta : name of the file containing the sequence of nucleotides
*/

	Traitement(const string& f);
/*!
   \brief  Default constructor, if haven't got a fasta
 */
	Traitement() {}

    /*!
     * \brief : Destructor : iberates the memory after the process
     */
	~Traitement();

/*!
   \brief Gives the reversed sequence to iterate on the negative strand
   @param seq : the sequence of nucleotide
   @return sequence of negative strand of the sequence in parameter in the order 5'->3'
 */
	
    string inversion(const string& seq);
/*!
 * \brief Returns the number of the column corresponding to the nucleotide in the Position Weight Matrix
 * @param c : Nucleotide to be converted
 * @return column of the nucleotide in the PWM
 */
	double conversion(char c);
	
	/*!
	 * \brief Read the fasta 1000 nucleotides at a time
	 * \brief Make change the chromosome if necessary
	 * @param line : string containing the nucleotides read a thousand at a time
	 * @param changement : boolean in reference to be changed by the methoed (parameter and output)
	 * @return the lecture is a success or not
	 */
	bool lire_fasta(string& line, bool& changement);
	
	/*!
     * @brief set the initial fasta file
     * @param f : name of the fasta file
     */
    void set_fasta(const string& f);
    
    /*!
     * @brief set the traitement class in EM mode :
     * Lire_fasta return one by one the sequence of interest,
     * and not 1000 by 1000 nucleotides of the same chromosome
     */
    void set_EM() {EM_algorith = true;}
    
    /*!
     * @brief : change fasta if we have a new fasta file,
     * close the old folder and open the new
     * @param f : name of the new fasta file
     */
    void new_loop(const string& f);
    
    /*!
     * @brief Change the outfile for a specific file, not default (terminal)
     * @param name_ : name of the outfile
     */
	void set_outfile (const std::string& name_);
	
	/*!
	* @brief return outname
	* @return outname
	*/
	string get_out() const {return out_name;}
	
	/*!
	* @brief return fasta_name
	* @return fasta name
	*/
	string get_fasta_name() const {return fasta_name;}
	
	/*!
	* @brief Output stream
	*/
	ofstream out;
	
	private :

    /*!
     * @brief imput stream corresponding to the fasta where we iterate
     */
	ifstream* fasta;
	/*!
	 @brief name of the file containing the fasta
	 */
	string fasta_name;
    /*!
     @brief Name of the chromosome under the treatment
     */
	string chr_name;

	/*!
	 @brief Do we deal with another chromosome during the treatment?
	 */
	bool changement_chr = false;
	
	/*!
     @brief Are we in the EM_algorithm ?
     */
	bool EM_algorith = false;

	/*!
	*  @brief name of output fil
	*/
	string out_name;
	
	
};

