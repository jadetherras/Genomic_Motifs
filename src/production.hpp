#include "Matrix.hpp"

using namespace std;

/*! \file production.hpp
 */

/*! @struct line_bedgraph
 *  @brief Store the datas of a line of the bedgraph
 *  A line contains :
 *  @param chromosome : number of the chromosome
 *  @param begin : position of the begining of the sequence corresponding to the data
 *  @param end : position of the end of the sequence corresponding to the data
 *  @param data : experimental value of the measurement of protein binding
 */

struct line_bedgraph
{
	//!number of the chromosome
    int chromosome;
    //!begin position
    size_t begin;
    //!end position
    size_t end;
    //!the score allowed to a nucleotide between begin and end
    double data;
};

/*! \class Production
    \brief Class for fonctionnality one, product of a bed file with a fasta and a matrix
    Extension : If the user demand it, this class can also display a supplementary output with quantitative values from a bedgraph
 */

class Production
{
    public :

    /*!
     * \brief Constructor of the class
     * @param name_mat : name of the file containing the matrix
     * @param name_seq : name of the file containing the fasta
     * @param T : treshold sequence beyond which the motif is displayed
     */
    Production(const std::string& name_mat,const std::string& name_seq, double T);
    
    /*!
     * \brief Constructor of the class in case of EM algorith
     * @param M : the matrix we use to compute the sequence with best score
     * @param a_remplir : the matrix to fill with the motif with the best score
     */
    Production(Matrix* M, Matrix* a_remplir);

    /*!
     * \brief : Destructor : liberates the memory after the process
     */
    ~Production();
	
    /*!
     * \brief Calls the function changement_chr() or traitement_ligne() in function of the situation
     * extension : if find_max mode (find the motif with best score), reset for all new sequence and do "last" to fill the matrix
     * 
     */
    void production_file();

    /*! \brief Changes the chromosome by changing the attribute number_chr and reset pos and position
     * @param line : name of the new chromosome
     */
    void changement_chr(string& line);

    /*! \brief Calls calcul_score and checks if motifs are beyond the treshold
     * @param line : we iterate with to check all the motif (fasta line)
    */
    void traitement_ligne(const string& line);

    /*!
     * \brief Writes the line of the bedgraph in the output file
     * @param sens : positive or negative strand?
     * @param sequence : sequence with the score beyond the treshold
     * @param s : score
     */
    void write(const std::string& sens, const std::string& sequence, double s);

    /*!
     * \brief Calculate the value of protein binding equal to the sum of the elements in the attribute table
     * @return Value of protein binding (if there is a Bedgraph) corresponding to a motif
     */
    double surface();

    /*!
     * @brief calcul the score of a motif
     * @param seq : the motif if interest
     * @return score of the motif
     */
    double calcul_score(const std::string& seq);

    /*!
     * \brief Gets the number of the chromosome
     * @return number of the chromosome
     */
    int get_name_chrm() const {return number_chr;}

    /*!
     * \brief Gets the attribute pos
     * @return Attribute pos
     */
    size_t get_position() const {return pos;}

    /*!
     * \brief Gets the size of the attribute position
     * @return size of attribute position
     */
    size_t get_position_size() const {return position.size();}

    /*!
    * \brief Gets the attribute act_line
    * @return attribute act_line
    */
    line_bedgraph get_act_line() const {return *act_line;}
   
    /*!
    *\brief Gets the attribute bedgragh
    *@return attribute bedgraph
    */
    bool get_bedgraph() const {return bedgraph;}

    /*!
    *\brief Gets the attribute bed_name
    *@return attribute bed_name
    */
    string get_bed_name() const {return bed_name;}

    /*!
    *\brief Gets the attribute get_table
    *@return attribute get_table
    */
    list<double> get_table() const {return table;}

    /*!
     * @brief Makes the bool bedgraph true
     * @param name_ : name of the file containing the bedgraph
     */
	void set_bedfile(const std::string& name_);

	/*!
	 * \brief Updates the attribute table
	 * @param p : position to reach in the bedgraph
	 */
	void edit_table(size_t p);

	/*!
	 * \brief Initialisation of the list table
	 */
	void init();

    /*!
     * \brief Read the line of the bedgraph to store the data in the address act_line
     */
	void read_line();
	
	/*!
     * \brief reset the position and the vector of motif
     *  */
	void reset();
	
	/*!
 * @brief PWM
 */
	Matrix* M;
	
	/*!
 * @brief PWM to fill
 */
	Matrix* a_remplir;
    
private :
/*!
* @brief datas of bedgraph
*/
list<double> table;

/*!
 * \brief adress contaning the line_bedgraph of the line of the bedgraph to be read
 */
line_bedgraph* act_line;

/*!
 * \brief is there a bedgraph
 */
bool bedgraph = false;

/*!
 * @brief input stream of bedgraph
 */
ifstream* file;

/*!
 * \brief Reading of the bedgraph over?
 */

bool fin_bedfile =false;

/*!
 * @brief Name of the bed graph
 */
string bed_name;

/*!
 * \brief Sequences stored during the iteration/treatment
 */
std::string line, inverse;

/*!
 * @brief number of the chromosome
 */
int number_chr;

/*!
 * \brief tool to detect the motifs
 */
vector<string> position;

/*!
 * \brief current position on the chromosome
 */
size_t pos = 0;

/*!
 * @brief treshold
 */
double T;

/*!
 * @brief find the sequence with best score
 */
bool find_max;

/*!
 * @brief score maximum
 */
double max_score = 0;

/*!
 * @brief motif with best score
 */
string best_seq = "";
};
