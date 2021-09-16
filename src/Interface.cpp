#include "Interface.hpp"

Interface::Interface(int argc, char** argv){
	try{
	TCLAP::CmdLine cmd("Motif genomic PWM");
	
	TCLAP::SwitchArg file("F", "o_file", "output list of motifs");    
	TCLAP::SwitchArg PWM("M", "o_matrix", "output position weight matrix");
	cmd.xorAdd(file,PWM);
	
	
	TCLAP::ValueArg<std::string> fasta("f", "name_fasta", "name of fasta file", true, "", "string");     
	cmd.add(fasta);
	TCLAP::ValueArg<std::string> output("o", "output", "output file", false, "", "string");     
	cmd.add(output);
	TCLAP::ValueArg<std::string> bed("b", "bed", "bed file", false, "", "string");     
	cmd.add(bed);
	
	TCLAP::SwitchArg EM_algo("E", "o_matrix_optimisee", "output position weight matrix with optimisation");
	cmd.add(EM_algo);
	TCLAP::ValueArg<std::string> matrix("m", "name_matrix", "name of matrix file", false, "", "string");     
	cmd.add(matrix);
	TCLAP::ValueArg< double > length("l", "length", "sequence size",  false, 0.0, "double"); 
	cmd.add(length);
	TCLAP::SwitchArg WantFasta("W", "WantFasta", "Have a fasta folder"); 
	cmd.add(WantFasta);
	TCLAP::ValueArg<std::string> FastaOutput("s", "FastaOutput", "Fasta output file", false, "interest_sequence.fasta", "string");     
	cmd.add(FastaOutput);
	TCLAP::ValueArg<std::string> BedGraph("B", "BedGraph", "name of the bedgraph file in input", false, "", "string");
	cmd.add(BedGraph);
	TCLAP::ValueArg< double > T("T", "threshold", "threshold of the score of motif",  false, 0.0, "double"); 
	cmd.add(T);
	
	cmd.parse(argc, argv);
	
	if (PWM.isSet() and length.isSet()) {
		if(EM_algo.isSet()) {
			EM_algorith(length.getValue(),bed.getValue(), fasta.getValue(),FastaOutput.getValue(), output.getValue());
		} else {
		Matrix matrix(fasta.getValue(), length.getValue());
		if (output.isSet()) {
			matrix.traitement->set_outfile(output.getValue());
		}
		matrix.set_FastaOutput(FastaOutput.getValue());
		matrix.Make_Matrix(bed.getValue(), FastaOutput.getValue(), WantFasta.getValue());
		matrix.print_matrix();
		}
		
	} else if (file.isSet() and T.isSet() and matrix.isSet()) {
		Production P(matrix.getValue(), fasta.getValue(), T.getValue());
		
		if (output.isSet()) {
			P.M->traitement->set_outfile(output.getValue());
		}
		if (BedGraph.isSet()) {
                if (BedGraph.getValue() != "")
					P.set_bedfile(BedGraph.getValue());
            }
			P.production_file();
    }
    }
	
	catch (TCLAP:: ArgException &e){
			std::cerr<<"error: "<<e.error()<<" for arg "<<e.argId() << std::endl;
	}
}

void Interface::EM_algorith(int length, const string& bed, const string& fasta, const string& FastaOutput, const string& output) {
	Matrix matrix(fasta, length);
	if (!output.empty()) {
			matrix.traitement->set_outfile(output);
		}
	matrix.set_FastaOutput(FastaOutput);
	Matrix remplir(fasta, length);
	
	matrix.Make_Matrix(bed, FastaOutput, true);
	matrix.traitement->set_fasta(FastaOutput);
	matrix.traitement->set_EM();
	
	Production P(&matrix,&remplir);
	for(size_t i(0); remplir !=matrix and i<_MAX_EM_; ++i) {
		if(i !=0 ) {
			matrix = remplir;
		}
	matrix.traitement->new_loop(FastaOutput);
	remplir.empty_matrix();
	P.production_file();
	}
	
	matrix.print_matrix();
}
