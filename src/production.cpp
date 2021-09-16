#include "production.hpp"

Production::Production(const std::string& name_mat,const std::string& name_seq, double T)
:
T(T)
{
	 M = new Matrix(name_mat, name_seq);
}

Production::Production(Matrix* M, Matrix* a_remplir) 
:
M(M),
a_remplir(a_remplir)
{
	find_max = true;
}

Production::~Production() {
	if(bedgraph) {
	delete act_line;
	act_line = nullptr;
	delete file;
	file = nullptr;
	}
	
	if(find_max) {
	a_remplir = nullptr;
	M = nullptr;
	} else {
	delete M;
	M = nullptr;
	}
		
}

void Production::production_file()
{
		bool change_chr(false);
		string line;
		
		while(M->traitement->lire_fasta(line, change_chr)) {
			if (change_chr == true) {
				changement_chr(line);
				} else {
					if(find_max) {
						reset();
						 }
					traitement_ligne(line);
					}
				}
				if (!line.empty()) {
					if(find_max) {
						reset();
						 }
					traitement_ligne(line);
				}
				
		if(find_max) {
			a_remplir->last();
		}
		
}


void Production::changement_chr(string& line) {
    try {
		line = line.substr(3);
		number_chr = stoi(line);
    } catch (out_of_range &) {
			throw(TOOSHORT_ERROR("string is too short"));
    }
	reset();
	if(bedgraph) {
		init();
	}
}

void Production::reset() {
	pos = 0;			   
	position.clear();
	for(size_t i(0); i < M->size(); ++i) {
	position.push_back("");
	}
}

void Production::traitement_ligne(const string& line) {
	
    double first = true;
    
	for(size_t i(0); i < line.size(); ++i) {
						  if(bedgraph) {
							edit_table(pos+10);
						  }
						  
						  for(size_t j(0); j<std::min(M->size(),pos+1); ++j) {
						  position[j] += line[i];
						  if (position[j].size() == M->size()) {
							  
							 double S(calcul_score(position[j]));
							 
							 if(!find_max) {
								  if (S >= T) {
									  write("+",position[j], S);
									  }
									  
							  inverse = M->traitement->inversion(position[j]);
							  S =calcul_score(inverse);
								  if (S >= T) {
									  write("-",position[j], S);
									  }
							}else {
								 if(S > max_score or first) {
									max_score = S;
									best_seq = position[j];
									first = false;
								}
							}
							     position[j].clear();
							  }
							  }
							  ++pos;
	}
	if(find_max) {
	a_remplir->add_motif(best_seq);
	}
}

void Production::write(const std::string& sens, const std::string& sequence, double s) {
	std::ostream *outstr = &std::cout;
		if (M->traitement->out.is_open()) {
			outstr = &M->traitement->out; 
		}
		
		try {
	*outstr << "chr" << number_chr << " " << pos-M->size()+2 << " " << sens << " " << sequence << " " << s;
	if (bedgraph) {
		*outstr<<" "<< surface();
	}
	*outstr<<endl;
	} catch(std::ofstream::failure &e) {
        throw(OUTPUT_ERROR("Error with outfile file " + M->traitement->get_out() + ": " + e.what()));
    }
	
}
	
double Production::calcul_score(const std::string& seq) {
	double score(0);
	for(size_t k(0); k<seq.size(); ++k) {
		double N(M->traitement->conversion(seq[k]));
		if (N != -1) {
			score += 2 + log2((*M)[k][N]);
		} 
		}
		return score;
}

double Production::surface()
{
    double return_(0);
    for (auto I = table.begin(); I!= table.end();I++)
        return_ += *I;
    return return_;

}


void Production::init()
{
	table.clear();
	
	for(size_t i(0); i < 20+M->size(); ++i) {
		table.push_back(0.0);
		}
		
	if(act_line->chromosome == 0) {
	read_line();
	}
	for(size_t i(0); i < 10; ++i) {
		edit_table(i);
	}
	
}

void Production::read_line() {
    string key, line_;
    try {
    if (file->is_open()) {
        if (std::getline(*file, line_)) {
        if (!line_.empty()) {
            std::stringstream ss(line_);
            std::getline(ss, key, '\t');
            try {
                key = key.substr(3);
                act_line->chromosome = stoi(key);
            } catch (out_of_range &) {
				throw(TOOSHORT_ERROR("string is too short"));
            }
            std::getline(ss, key, '\t');
            act_line->begin = stoi(key);
            std::getline(ss, key, '\t');
            act_line->end = stoi(key);
            std::getline(ss, key, '\n');
            act_line->data = stod(key);
        } else {fin_bedfile = true;}
       } else {fin_bedfile = true;}
    } else throw(BEDGRAPH_ERROR("Could not open bedgraph file " + bed_name));
    } catch(std::ifstream::failure &e) {
        throw(BEDGRAPH_ERROR("Error with bedgraph file " + bed_name + ": " + e.what()));
    }
}

void Production::set_bedfile (const std::string& name_) {
	bedgraph = true;
	bed_name = name_;
	file = new ifstream(name_);
	act_line = new line_bedgraph{0,0,0,0};
}
	
void Production::edit_table(size_t p) {
	table.pop_front();
	double to_push(0);
	while((number_chr!= act_line->chromosome or p>=act_line->end) and !fin_bedfile) {
		read_line();
	}
	if(number_chr == act_line->chromosome) {
	if (p < act_line->begin) {
		to_push = 0;
	} else if(p < act_line->end and p >= act_line->begin) {
		to_push = act_line->data;
	} else if(!fin_bedfile) throw(BEDGRAPH_ERROR("try to find bedgraph value for a past nucleotide "));
	}
	table.push_back(to_push);
}
