#include "traitement.hpp"

/*! \file Matrix.hpp
 */
  
 /*! @struct Bedline
 *  @brief Stores one sequence's data from a Bed file 
 *  A Bed line contains : 
 *  @param chromosome : number of the chromosome 
 *  @param pos_ini : position of the begining of the sequence
 *  @param pos_end : position of the end of the sequence 
 *  @param sign : Orientation of the sequence : 
* -  if the sequence belongs to the forward strand : +
* -  if the sequence belongs to the reverse strand :  -  
 */
 
 
struct Bedline {

	//!Name of the chromosome.
	int chromosome; 
	//!Start position of the region of interest.
	size_t pos_ini; 
	//!Stop position of the region of interest.
	size_t pos_end; 
	//! Orientation of the region: + for forward and - for reverse strand.
	char sign; 
};

/*! \class Matrix
     \brief Main class for part 2 of the project. Proceeds to fill a matrix thanks
      to a fasta file, where we search specific sequences we obtain from a Bed file.
      * public vector of vector of double
  */

class Matrix: public std::vector<std::vector<double>>{
	
	public:
	
	/*!
     * \brief operator !=
     * @param m2 : the matrix to compare with the current matrix
     * @return true if matrix are different, false otherwise
     */
	bool operator!=(const Matrix& m2);
	
	/*!
     * \brief operator =
     * @param m : all the value of the current matrix will be egal to the value of
     * this matrix
     */
	void operator=(const Matrix& m);
	
    /*!
     * \brief : Destructor : iberates the memory after the process
     */
	~Matrix();

	
	/*!
	 * \brief Constructor which takes as arguments : 
	 * @param fasta :   The fasta we want to read the sequences from 
	 * @param sequence_length : The length of the motifs we want to analyse 
	 * @param initialisation : Boolean initialisation set at true to signify we're initializing our matrix 
	 * @param filled : Boolean filled set at true to signify we need to fill the first tab of 1000 chars 
	 */
	 
	Matrix(const string& fasta, double,bool initialisation = true,bool filled = true);
	
	
	
	/*!
	 * @brief Construct a matrix using a file and that takes as parameters :
	 *  @param name_mat : The filestream of the file to use to fill the matrix 
	 *  @param fasta : The fasta we want to read the sequences from
	 */
	Matrix(const string& name_mat,const string& fasta = "");
	
	/*!
	 * @brief set the outfile for the sequence of interest of the fasta
	 * @param out : name of the outfile
	 */
	void set_FastaOutput(const string& Out);
	
	/*!
	 * @brief method which opens the bedfile and proceeds to stock the smallest "next" one, as well as the following line.
	 *
	 * @param Bed : filestream of the bedfile to open and use.
	 */
	Bedline read_bed(const std::string&);
	
	/*!
	 * @brief Main method of the Matrix class that, with the help of other functions, produces the PWM matrix needed.
	 * if asked, can also add the analyzed sequences read from the fasta file, to an output file.
	 * @param bed : Bed file that tells us the sequences to analyze
	 * @param WantFasta : Boolean set to true if we want to have an output file with the sequences analyzed
	 * @param fastaOutput : File where we note down the sequences analyzed, if WantFasta is true
	 * 
	 */
	void Make_Matrix( const std::string &bed, const std::string& FastaOutput, bool WantFasta) ;
	
	/*!
	 * @brief The prime objective of the method is to find the next sequence in the fasta_file : this is done through numerous checks about first_line, second_line, and about the 1 000 chars line from the fasta.
	 * @param bed : Bed file that tells us the sequences to analyze
	 * @param WantFasta : Boolean set to true if we want to have an output file with the sequences analyzed
	 * @param fastaOutput : File where we note down the sequences analyzed, if WantFasta is true
	 * @param fasta_line : the fasta line to be traited
	 */
	
	void find_sequence(bool WantFasta, const std::string &Bed, const std::string &FastaOutput,const std::string& fasta_line);
	
	
	/*!
	 * @brief fills the PWM matrix for a given sequence which is taken as parameter
	   @param sequence : the sequence of interest
	 */
	  void fill_matrix(const std::string& sequence);
	  
	/*!
	 * @brief Print the matrix values in an output file
	 */
	 void print_matrix();
	
	
	/*!
	 * @brief Returns the number of motifs possible for each sequence to analyze
	 */	
	   double motifs_number();
	
	/*!
	  * @brief Returns true if the bedline a we take as the first parameter is located before the second parameter b on the genome
		@param a : to be compare 
		@param b : to be compare
		@return true if a is before b
	  */
	  bool compare_lower(const Bedline&a, const Bedline&b);
	
	 /*!
	  * @brief Returns true if the bedline a we take as the first parameter is located after the second parameter b on the genome
		@param a : to be compare 
		@param b : to be compare
		@return true if a if after b
	  */
	  bool compare_higher(const Bedline&a, const Bedline&b);
	  

	/*!
	 * @brief Adds +0.25 to each compartment of the matrix as well as dividing them by (( length of the sequences) +1 ), for the sum of probabilities of each row to be equal to 1.
	 */	
	  void last();
	  
	
	
    /*!
     * @brief Method needed for the 2nd extension : allows to directly add a motif ( adding +1 where we need to) into our matrix
     * @param motif : the motif to add to the matrix
     */
	void add_motif(const string& motif);
	
	/*!
     * @brief Empties the matrix ( fills it with Os) + resets nb_seq ( see below )
     */
	void empty_matrix();
	
	/*!
     * @brief Traitement element
     */
	Traitement* traitement; 
	
	private:
		
	 /*!
	  * @brief Returns true if the chromosome of the Bedfile's line we're reading is indeed the one we aim to analyze
	  */
	  bool good_chromosome();
	
	 /*!
	   * @brief Prints the sequence that has previously been used to fill the matrix
		@param out : the outfile to print the sequence
	   */
	  void print_seq(const std::string& Out);
	  
	 /*!
	  * @brief Simple method that is used to set two of our Booleans , KeepLine & SameTab. 
		@param keep
		@param same
	  */
	  void setStay(const bool & Keep, const bool &Same);
	  
	  /*!
	   * @brief Resets PositionOnChromosome and gets ready for a change of chromosome
		* @param line : the name of the new chromosome
	   */ 
	  
	  void reset(std::string & line);

	
	  /*!
	   * @brief Stocks the first / previous smallest sequence we analyzed, tthus working as a threshold for second_line to take the sequence we need to analyze now. Is therefore replaced by second_line each time read_bed is called.
	   */
	   Bedline first_line;
	   
	   /*!
	   * @brief Stocks the characteristics of the sequence we are going to analyze 
	   */
	   Bedline second_line;
	   
	   /*!
	   * @brief Stocks the characteristics of the next sequence from the bedfile.
	   */
	   Bedline next_line;
	   
	   /*!
	   * @brief  Number of sequences from the Bed file we've analyzed, plus one.
	   */
	   int nb_seq = 0;
	   
	   /*!
	   * @brief Current chromosome from which we want to analyze all the sequences given in the bedfile.
	   */
	   int current_chr;
	   
	    /*!
	   * @brief Length of the motifs we need to analyze, which is also the number of rows of our Matrix.
	   */
	   double length;
	   
	   /*!
	   * @brief  Sequence of nucleotides from the fasta we use to fill the final matrix
	   */
	   std::string temp_seq;
	   
	   
      /*!
	   * @brief Position, in the fasta file, on the chromosome we're on, while the file is being read.
	   */
	   size_t PositionOnChromosome=0;
	  
	  /*!
	   * @brief current position on each line of the fasta
	   */
	   size_t PositionOnTab = 0;
	    
	  /*!
	   * @brief Marks the initial position in the current line of 1 000 chars of the sequence we're currently analyzing 
	   */
	    unsigned int ini = 0;
	   
	  /*!
	   * @brief  Marks the final position in the current line of 1 000 chars of the sequence we're currently analyzing 
	   */
	   unsigned int end = 0;
	  
	  /*!
	   * @brief Checks if we keep the same line of 1 000 chars from the fasta that we used for the previous sequence
	   */
	   bool SameTab = false;
	  
	  /*!
	   * @brief  Checks if we've changed tab of 1000 chars: if true, allows to recalculate ini and end.
	   */
	   bool NextTab = false;
	   
	  /*!
	   * @brief Checks if we need to keep the same line from the Bedline for further calculations : for instance, if the sequence is split between two tabs. 
	   */
	   bool keepLine =false;
	   
	   /*!
	   * @brief For the method read_bed, checks if we've previously used it or if we need to initialize the first_line / second_line attributes.
	   */
	   bool Initialisation = true;
	   
	   /*!
	   * @brief returns true if the current tab is filled
	   */
	   bool filled;
	   
	   
	   /*!
	   * @brief For print_seq, checks whether the next chromosome to analyze is a new one or not
	   */
	    bool new_chr = false;
	    
		/*!
	   * @brief Name of the fasta_file where we want to print our analyzed sequences
	   */
		ofstream FastaOutput;
};
