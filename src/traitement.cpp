#include "traitement.hpp"
	
Traitement::Traitement(const string& f) 
{
	set_fasta(f);
}

Traitement::~Traitement(){
fasta->close();
out.close();
}

void Traitement::new_loop(const string& f) {
	fasta->close();
	set_fasta(f);
}

double Traitement::conversion(char c) {
	double n(0);
	switch (c) {
	case 'A':
	case 'a':
	n = 0;
	break;
	case 'C':
	case 'c':
	n = 1;
	break;
	case 'G':
	case 'g':
	n = 2;
	break;
	case 'T':
	case 't':
	n = 3;
	break;
	default:
	n = -1;
	break;
	}
	return n;
}

std::string Traitement::inversion(const string& seq) {
	std::string inverse;
	for(int u(seq.size()-1); u >= 0; --u) {
		switch (seq[u]) {
	case 'A':
	case 'a':
	inverse.push_back('T');
	break;
	case 'C':
	case 'c':
	inverse.push_back('G');
	break;
	case 'G':
	case 'g':
	inverse.push_back('C');
	break;
	case 'T':
	case 't':
	inverse.push_back('A');
	break;
	default :
	inverse.push_back('N');
	break;
	}
		}
		return inverse;
}


bool Traitement::lire_fasta(string& line, bool& changement) {
		line.clear();
		char n[_MAX_FASTA_];
		changement = false;
		try {
		if(fasta->is_open()) {
			if(changement_chr == true) {
				changement = true;
				changement_chr = false;
				line = chr_name;
				return true;
		} else {
	     while(!fasta->eof()) {
			if(fasta->fail() && size_t(fasta->gcount()) == _MAX_FASTA_ -1) 
			{fasta->clear();}			
			fasta->getline(n, _MAX_FASTA_ - line.size()+1, '\n');
			if(n[0] == '>') {
				if(!line.empty()) {
				chr_name = n;
				chr_name.erase(0,1);
					changement_chr = true;
					return true;
				} else {
					line = n;
					line.erase(0,1);
					changement = true;
					return true;
				}
				} else {
					if(EM_algorith) {
						line = n;
						return true;
					} else {
					line += n;
					if (line.size() == _MAX_FASTA_) {
						return true;
					}
				}
			}
		}
		}
			return false;
	}
	 else throw(FASTA_ERROR("Could not open fasta file " + fasta_name));
    } catch(std::ifstream::failure &e) {
        throw(FASTA_ERROR("Error with fasta file " + fasta_name + ": " + e.what()));
    }
}

void Traitement::set_fasta(const string& f) {
 fasta = new ifstream(f);
 fasta_name = f;
}

void Traitement::set_outfile (const std::string& name_) {
	try {
		if (!name_.empty()) {
		out_name =name_;
		out.open(name_, std::ios_base::out);
	}
	} catch(std::ofstream::failure &e) {
        throw(OUTPUT_ERROR("Could not open output file " + out_name));
    }
	}


