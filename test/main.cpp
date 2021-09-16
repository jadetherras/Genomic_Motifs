#include <gtest/gtest.h>
#include "../src/production.hpp"


Traitement traitement("../test/read_fasta.fasta");
Traitement traitement_EM("../test/promoters.fasta");
Production P("../test/DBP.mat","../test/promoters.fasta", 0);
Production P_extension("../test/DBP.mat","../test/promoters.fasta", 0);
std::string seq_test("NACGT");
std::string seq_score_connu("AGCTGAT");
std::string seq_chr("chr7");
std::string ligne("TCCAATCCACCCTCTGACTATTCCACATCCC");
line_bedgraph bedg {7, 40, 65, 0.955};

Matrix matrix("../test/promoters.fasta", 3); 
Matrix matrix2("../test/promoters.fasta", 3); 
Matrix M ("../test/DBP.mat","../test/promoters.fasta");
Matrix testMat_last ("../test/last.mat","../test/promoters.fasta");
Matrix testMat("../test/testMat.mat","../test/promoters.fasta");

struct AllBedlines{
		
		Bedline first;
		Bedline second;
		Bedline third;
		Bedline fourth;
};

TEST(Traitement,conversion) {
	double test(-1);
	for(size_t i(0); i<seq_test.size(); ++i){
		EXPECT_EQ(test,traitement.conversion(seq_test[i]));
		++test;
		}
}

TEST(Traitement,inversion) {
	std::string inverse(traitement.inversion(seq_test));
	EXPECT_TRUE(inverse == "ACGTN");
}

TEST(Traitement, lire_fasta) {
	bool changement(false);
	string lecture("");
	EXPECT_TRUE(traitement.lire_fasta(lecture, changement));
	EXPECT_TRUE(lecture == "chr7");
	EXPECT_TRUE(changement);
	EXPECT_FALSE(traitement.lire_fasta(lecture, changement));
	EXPECT_TRUE(lecture == "ATGCCCTC");
	EXPECT_FALSE(changement);
	EXPECT_FALSE(traitement.lire_fasta(lecture, changement));
	
	//en mode EM, on lit ligne par ligne
	changement = false;
	lecture = "";
	traitement_EM.set_EM();
	EXPECT_TRUE(traitement_EM.lire_fasta(lecture, changement));
	EXPECT_TRUE(lecture == "chr7");
	EXPECT_TRUE(changement);
	EXPECT_TRUE(traitement_EM.lire_fasta(lecture, changement));
	EXPECT_TRUE(lecture == "TCCAATCCACCCTCTGACTATTCCACATCCCATACCTCCTCCTCATCCCCTTGTCTCCATGAGGATTTCTCCACCCCACC");
}

TEST(Production,changement_chr) {
        string lecture;
        bool changement (false);
        P.M->traitement->lire_fasta(lecture, changement);
        EXPECT_TRUE(changement = true);
        P.changement_chr(lecture);
        EXPECT_EQ(0, P.get_position());
        EXPECT_EQ(7, P.get_name_chrm());
}

TEST(Production,calcul_score) {
	EXPECT_NEAR(0.387621, P.calcul_score(seq_score_connu), pow(10,-5));
}

void check_line(const line_bedgraph &ref,const Production*P){
    EXPECT_EQ(ref.chromosome, P->get_act_line().chromosome);
    EXPECT_EQ(ref.begin, P->get_act_line().begin);
    EXPECT_EQ(ref.end, P->get_act_line().end);
    EXPECT_NEAR(ref.data, P->get_act_line().data, pow(10,-5));
}


TEST(Production, set_bedfile){
    P_extension.set_bedfile("../test/test_bedgraph.bedgraph");
    EXPECT_TRUE(P_extension.get_bedgraph());
    EXPECT_EQ("../test/test_bedgraph.bedgraph", P_extension.get_bed_name());
    line_bedgraph ref{0, 0, 0, 0};
    check_line( ref, &P_extension);
}

TEST(Production, edit_table) {
    P_extension.init();
    bool changement;
    string line;
    P_extension.M->traitement->lire_fasta(line, changement);
    P_extension.changement_chr(line);
    for (size_t i(0); i < 65; ++i) {
		P_extension.edit_table(i+10);
		if(i == 0) {
    EXPECT_EQ(0, P_extension.surface());
	}
    if (i == 70) {
    EXPECT_EQ(39.46, P_extension.surface());
	}
	}
}
	
TEST(Production, read_line) {
	P_extension.set_bedfile("../test/test_bedgraph.bedgraph");
	P_extension.read_line();
    check_line(bedg, &P_extension);
}

void check_bedline( const Bedline &test1,const Bedline &test2)
{
    EXPECT_EQ(test1.chromosome, test2.chromosome);	//checking the header of the chromosome
    EXPECT_EQ(test1.pos_ini, test2.pos_ini);   //checking the start position
    EXPECT_EQ(test1.pos_end, test2.pos_end); 	//checking the stop position
    EXPECT_EQ(test1.sign, test2.sign);  //checking the sign
}

TEST(Matrix,read_bed) {
	
	//For a bedfile of 6 columns:
	AllBedlines All = { {1,120,125,'+'},
		                {1,525,621,'-'},
		                {2,109276,110001,'-'},
		                {3,1,12,'+'}
					  };
		
	Bedline test_bed(matrix2.read_bed("../test/bed6cols.bed"));
	check_bedline(test_bed,All.first);

	Bedline test_bed2(matrix2.read_bed("../test/bed6cols.bed"));
    check_bedline(test_bed2,All.second);
	
	Bedline test_bed3(matrix2.read_bed("../test/bed6cols.bed"));
    check_bedline(test_bed3,All.third);
	
	Bedline test_bed4(matrix2.read_bed("../test/bed6cols.bed"));
    check_bedline(test_bed4,All.fourth);
	
}

void check_matrix(const Matrix &mat1, const Matrix &mat2)
{
    for ( size_t i(0) ; i < 3 ; ++i){

        for ( size_t j(0) ; j < 4 ; ++j) {

            EXPECT_EQ( mat1[i][j],mat2[i][j]);
        }
    }
}

TEST(Matrix,fill_matrix) {
	matrix.read_bed("../test/testfill_matrix.bed");
	matrix.fill_matrix(seq_score_connu);

	check_matrix(matrix, testMat);
}

TEST(Matrix,last) {
	matrix.last();
	check_matrix(matrix, testMat_last);
}


TEST(Matrix,motifs_number) {
	matrix.read_bed("../test/sequence_number.bed");
	EXPECT_EQ(matrix.motifs_number(),6);
}


TEST(Matrix,compare_lower) {
	Bedline test1{1,12,15, '+'}; 
	Bedline test2{1,13,16, '+'}; 
	EXPECT_TRUE(matrix.compare_lower(test1,test2));
}

TEST(Matrix,compare_higher) {
	Bedline test1{2,12,15, '+'}; 
	Bedline test2{1,13,16, '+'}; 
	EXPECT_TRUE(matrix.compare_higher(test1,test2));
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
