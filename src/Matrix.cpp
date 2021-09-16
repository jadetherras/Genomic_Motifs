#include "Matrix.hpp"
#include <sstream>
#include <fstream>


using namespace std;
	
void Matrix::operator=(const Matrix& m) {
	for ( size_t i(0) ;  i < length ; ++i ) {
       for ( size_t j(0) ; j < 4 ; ++j) { 
		   (*this)[i][j] = m[i][j]; 
		  }
    }
	}
	
bool Matrix::operator!=(const Matrix& m2) {
	for ( size_t i(0) ;  i < length; ++i ) {
       for ( size_t j(0) ; j < 4 ; ++j) { 
		   if ((*this)[i][j] != m2[i][j]) {
			   return true;
			   }
		  }
    }
    return false;
}
	

Matrix::Matrix(const string& fasta, double sequence_length,bool init,bool fill): length(sequence_length),Initialisation(init) , filled(fill)
{
	std::vector < double > temp;
	traitement = new Traitement(fasta);
	
	for ( size_t i(0) ;  i < length ; ++i ) {
	
       for ( size_t j(0) ; j < 4 ; ++j) {
		   
		   temp.push_back(0);
	   }
	   
	   this->push_back(temp);
	   temp.clear();
	
	}
	
	first_line.chromosome=0;
	first_line.pos_ini=0;
	first_line.pos_end=0;
	
	second_line.chromosome=INT_MAX;
	second_line.pos_ini=INT_MAX;
	second_line.pos_end=INT_MAX;
	
	next_line.chromosome = INT_MAX;
	next_line.pos_ini=INT_MAX;
	next_line.pos_end=INT_MAX;
	
	Initialisation = true;
	filled = true;
	SameTab = false;
	keepLine = false;

	
}

void Matrix::empty_matrix() {
	for(size_t i(0); i<length; ++i) {
		for(size_t j(0) ; j < 4 ; ++j) {
			(*this)[i][j] = 0.0;
	 }}
	 nb_seq = 0;
}
	
Matrix::~Matrix(){
	delete traitement;
	
	traitement=nullptr;
	

}

Matrix::Matrix(const string& name_mat,const string& fasta)
{
	if(fasta != "") {
	traitement = new Traitement(fasta);
	}
	
    std::ifstream file(name_mat);
    std::string value, line;
    std::vector< double> new_line(4);
    
    try {
		if (file.is_open())
    {
        while (std::getline(file, line))
        {
            std::stringstream ss(line);
            for (size_t n=0; n<4; n++) {
				std::getline(ss, value, ' ');
			std::stringstream ss(line);
            for(size_t n(0); n <4; ++ n) {
				std::getline(ss, value, '	');
                new_line[n] = stod(value);
                }
               }
                
            this->push_back(new_line);
			}
        file.close();
	} else throw(MATRIXFILE_ERROR("Could not open matrix file " + name_mat));
    } catch(std::ifstream::failure &e) {
        throw(FASTA_ERROR("Error with matrix file " + name_mat + ": " + e.what()));
    }

}


Bedline Matrix::read_bed (const std::string& Bed) {
	
	
	ifstream BedFile(Bed);
	
	string key, line;
	
	Bedline actual;
	
try {
	if(BedFile.is_open()){
			
		while(std::getline(BedFile, line)){
				
			if (!line.empty()){
				
				std::stringstream ss(line);
									 
				std::getline ( ss , key ,'\t');
					
				try{
				   
				    	
					key=key.substr(3); 
					
					}catch(out_of_range&){
						throw(TOOSHORT_ERROR("String is too short"));
					}
			
					
				actual.chromosome=stoi(key);
				
				   
				std::getline(ss, key, '\t');
				actual.pos_ini=stoi(key);
				
				 
				std::getline(ss, key, '\t');
				actual.pos_end=stoi(key);
			
				       
				 while ( (std::getline( ss, key ,'\t')) ){
	                         
				     if(key=="+" or key.empty()){  
		
						actual.sign= '+';
				
							
					}else{
						
						actual.sign = '-';		
						
					}
				}
				
				if ( !(std::getline(ss,key,'\n')) and actual.sign != '-'){
					
					actual.sign = '+' ;
				}
				
				
				
				if(Initialisation){
			
					if(compare_lower(actual, second_line) && compare_higher(actual,first_line)){  
							
						second_line = actual;
						Initialisation=false;
						
					}
			
						}else{
							
							
							if( compare_lower(actual, second_line) &&  compare_higher(actual, first_line) ){
								
								second_line = actual;
								
								}
						
							if ( compare_lower ( actual,next_line) && compare_higher (actual,second_line) ) {
								
								next_line = actual;
							}	
								
							if(second_line.chromosome==actual.chromosome   
							 and second_line.pos_ini==actual.pos_ini 
							 and second_line.pos_end==actual.pos_end){
								 	
								
								second_line = actual;
								
							}
							
							
						}
			        }
				}
			        
			        first_line = second_line;
			        second_line = next_line;
				
	                
	                next_line.chromosome=INT_MAX;
	                next_line.pos_ini=INT_MAX;
	                next_line.pos_end=INT_MAX;
	                
	                ++nb_seq;
	             
	               
	           

				}else throw(BEDFILE_ERROR("Could not open Bed file " + Bed));
				} catch(std::ifstream::failure &e) {
        throw(BEDFILE_ERROR("Error with Bed file " + Bed + ": " + e.what()));
    }			    
			
	BedFile.close();
	return first_line;
	
}




void Matrix::Make_Matrix(const string &bed, const string& FastaOutput, bool WantFasta) {
	
	std::string line;	
	
     bool changement(true); 

	
	while(traitement->lire_fasta(line,changement)){     
	        

		if ((changement) and !line.empty()) { 
				if(filled){
				

	  reset(line);

	 }else{
			
			find_sequence(WantFasta, bed, FastaOutput, line);

			reset(line);
				
			
		}
		
	}else{

		if (!SameTab) { 
			
		PositionOnTab += line.size();
		PositionOnChromosome += line.size();
			}
		}
		 if(PositionOnTab == line.size()){    
			
			filled = true;
			
			find_sequence(WantFasta, bed, FastaOutput, line);
			
			
			if (!SameTab) 
			{
			PositionOnTab = 0; 
			
			
					}
			

			}

		 }
		 try{

		 last();
		
				}catch(std::ifstream::failure &e) {
			throw(ZERO_ERROR(string("Division by zero: ") + e.what()));
	}

 }


void Matrix::find_sequence(bool WantFasta, const string &Bed, const string& FastaOutput,const string& fasta_line){
	
	do{
	
				if(!keepLine) { 
				
					first_line=read_bed(Bed);

				}

                 
				if(first_line.pos_ini > PositionOnChromosome && first_line.pos_end > PositionOnChromosome){  

				   filled=true;
				   setStay(true, false);


				    }else if(first_line.pos_ini < PositionOnChromosome && first_line.pos_end>PositionOnChromosome){ 
					
	
					 for(size_t i(first_line.pos_ini); i < PositionOnChromosome; ++i)
					 {
					
				        temp_seq.push_back(fasta_line[i-(PositionOnChromosome - fasta_line.size())] );    
				     }
				    
                     ini = 0;
					 end = (first_line.pos_end-PositionOnChromosome);
		
	                 NextTab =true;
					 filled=true;
					 setStay(true, false);
					
		
				}else{
					
				
				     if (!NextTab ){  
	                 				
					 ini = (first_line.pos_ini - (PositionOnChromosome - fasta_line.size()));
					 end = (first_line.pos_end - (PositionOnChromosome - fasta_line.size()));
				     }
					
					for (size_t i(ini); i< end ; ++i){
						

		temp_seq.push_back(fasta_line[i]);
            
					}
										
					if(!temp_seq.empty() and WantFasta) {
						
					  print_seq(FastaOutput);
				
				    }		
				    
					
					if(first_line.sign=='-'){
						
						temp_seq=traitement->inversion(temp_seq);
						
						}
					fill_matrix(temp_seq);
		
					
					keepLine = false;
					NextTab = false;
					temp_seq.clear();	
					
		           
					if (second_line.pos_ini < PositionOnChromosome && second_line.pos_end < PositionOnChromosome)  
					{                                                                                                
						setStay(true, true);
						read_bed(Bed);
					

																				
		            }else if( second_line.pos_ini < PositionOnChromosome && second_line.pos_end > PositionOnChromosome){  
						                                                                                                  
						keepLine = false;
						SameTab = true;
						}else{                                                                                             
						
						keepLine = false;
						SameTab = false;
					}
					
				}
					
		}while(SameTab);
}
	
	
void Matrix::fill_matrix(const std::string& sequence){
	
	for ( size_t i(0) ; i < motifs_number(); ++i )
	{
	
		for ( size_t j(0) ; j < length ; ++j )
		{
			try{
				
			(*this)[j][traitement->conversion(sequence[i+j])] += (1.0/motifs_number());
			
		}catch(std::ifstream::failure &e) {
			throw(ZERO_ERROR(string("Division by zero: ") + e.what()));
			}
		}
	} 	
	
}



void Matrix::last(){
	
	double n(nb_seq);
	for ( size_t i(0) ; i < length ; ++i)
	{
		for ( size_t j(0) ; j < 4 ; ++j)
		{
			(*this)[i][j] += 0.25;
			if(n==0){
				throw(ZERO_ERROR("Number of sequences is zero"));
				}else{
            (*this)[i][j] /= (n);
			} 			
	    }
	}
	
}

double Matrix::motifs_number(){
	double seq_nb(0);
	seq_nb += ( first_line.pos_end-first_line.pos_ini - length +1);
	if(seq_nb==0){
		throw(ZERO_ERROR("Number of sequences is zero"));
	}else{
	return seq_nb;
	}
}

void Matrix::print_matrix(){
	
	std::ostream *outstr = &std::cout;
	if(traitement->out.is_open()) {
		outstr = &traitement->out;
		}
try{
	  for(size_t i(0); i<length; ++i){
		  for(size_t j(0); j<4; ++j){
	  
			*outstr<< (*this)[i][j] << ' ';
		} *outstr <<endl;
	}
} catch(std::ofstream::failure &e) {
        throw(OUTPUT_ERROR("Error with outfile file " + traitement->get_out() + ": " + e.what()));
    }
	
}

	
bool Matrix::compare_lower(const Bedline&a, const Bedline&b){
	
	return (a.chromosome<b.chromosome or (a.chromosome==b.chromosome and a.pos_ini<b.pos_ini));
	
}

bool Matrix::compare_higher(const Bedline&a, const Bedline&b){
	
	return (a.chromosome>b.chromosome or (a.chromosome==b.chromosome and a.pos_ini>b.pos_ini));
	
}	

bool Matrix::good_chromosome(){
	
	return first_line.chromosome==current_chr;

}
	


void Matrix::set_FastaOutput(const string& Out) {
	try {
		if (!Out.empty()) {
		FastaOutput.open(Out, std::ios_base::out);
	}
	} catch(std::ofstream::failure &e) {
        throw(OUTPUT_ERROR("Could not open output file " + Out));
    }
 }

void Matrix::print_seq(const string& Out){
								
		std::ostream *outstr = &std::cout;
	if(FastaOutput.is_open()) {
		outstr = &FastaOutput;
		}
		try{				
			if(new_chr) {
			*outstr<< ">chr" << current_chr << endl;
			new_chr = false;
			}
				
			for(size_t i(0); i<temp_seq.size(); ++i){
			
				*outstr<<temp_seq[i];
			}
				
				*outstr<<endl;
						
				
				}catch(std::ifstream::failure &e) {
        throw(OUTPUT_ERROR("Error with outfile file " + Out + ": " + e.what()));
    }
				
	
}

void Matrix::setStay(const bool & Keep, const bool & Same){
	
	keepLine=Keep;
	SameTab=Same;
	
}

void Matrix::reset( string & line){
	
	current_chr = stoi(line.substr(3,line.size()-3 ));
	new_chr = true;
	       
	           PositionOnChromosome=0;
	          
}

void Matrix::add_motif(const string& motif) {
	
		for ( size_t j(0) ; j < length ; ++j )
		{
			try{
				
			(*this)[j][traitement->conversion(motif[j])] += 1.0;

		}catch(std::ifstream::failure &e) {
			throw(ZERO_ERROR(string("Division by zero: ") + e.what()));
			}
		}
		++nb_seq;
} 	


