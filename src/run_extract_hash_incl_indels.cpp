#include "run_extract_hash_incl_indels.h"

#include <getopt.h>
#include <iostream>


//** COPIED FROM ORIGINAL LINKER CODE ************************

/////////////// namespaces ///////////////////
namespace opt {
  static std::string chr_choice = "chr20";
  static std::string technology = "tenx";
  static std::string input_vcf_file = "./sample_data/BL1954_PCRFree.hets.recalibrated.vcf";
	static bool vcf_load = false;
  static std::string input_bam_file = "./sample_data/BL1954_10xG_chr20.bam";
  static bool multiple_bams = false;
  static std::string second_input_bam_file = "./sample_data/BL1954_10xG_chr20.bam";
  static std::string id_string = "default";
	static int start_bound = 0;
	static int end_bound = 300000000;
	static bool region_defined = false;
	static std::string refdir;
	static int strlen = 100;
	static int klen = 16;

  static IndelGenotyping::FilterMinima filter_minima;

};

/////////////// structures ///////////////////
static const char* shortopts = "ho:i:v:c:e:s:n:d:l:k:";
static const struct option longopts[] = {
	{ "help",       no_argument, NULL, 'h' },
  { "vcf-file",    no_argument, NULL, 'v' },
  { "bam-file",    no_argument, NULL, 'i' },
  { "technology",  no_argument, NULL, 'e' },
  { "chr-choice",  no_argument, NULL, 'c' },
  { "bam-file2",   no_argument, NULL, 's' },
  { "id_string",   no_argument, NULL, 'n' },
  { "refdir",  no_argument, NULL, 'd' },
  { "strlen",  no_argument, NULL, 'l' },
  { "klen", no_argument, NULL, 'k'},
};
/*
      int Samsum_lt = 2000;
      int Minimum_mapq = 20;                 // Minimum read map quality
      int Minimum_mapq_hic = 20;             // Minimum read map quality for HiC technology
      int Minimum_baseq = 8;                 // Minimum read base quality
      int Minimum_baseq_hic = 20;            // Minimum read base quality HiC technology
      int Minimum_baseq_tenx = 20;           // Minimum read base quality 10X technology


*/
///////////////////////////////////////////////////////////////////////////////////////////////////////
static const char *EXTRACT_USAGE_MESSAGE =
"Usage: linker extract_incl_vcf [OPTION] -v /path/to/vcffile.vcf -i /path/to/bamfile.bam -d /path/to/refgenome/directory\n\n"
"\n"
"  Options\n"
"  -i,      input bamfile path \n"
"  -v,      input vcffile path  \n"
"  -e,      long read tech ( tenx, pacbio, nanopore, illumina, hic ) \n"
"  -c,      chromosome name ( chr4 ) or contig name depending on bam \n"
"  -s,      second bam file path \n"
"  -n,      id string for output files \n"
"  -d,      directory containing the reference genome fasta files \n"
"  -l,      length for reference and variant sequences extracted \n"
"  -k,      kmer length for exact match between variant sequence and extracted read \n"
"\n";

///////////////////////////////////////////////////////////////////////////////////////////////////////
static void parse_region_string( std::string samtools_region ) {
	std::string chr_delimiter=":";
	std::string pos_delimiter="-";
	std::size_t index = samtools_region.find(chr_delimiter);
	if ( index != std::string::npos ) {
		std::string cid = samtools_region.substr(0,index);
		std::string gpos = samtools_region.substr(index+1);
		gpos.erase(std::remove(gpos.begin(),gpos.end(),','),gpos.end());
		std::size_t sub_index = gpos.find(pos_delimiter);
        	if ( sub_index != std::string::npos ) {
			opt::chr_choice = cid;
                	std::string start_string = gpos.substr(0,sub_index);
                	std::string end_string = gpos.substr(sub_index+1);
			int start_int = std::stoi(start_string);
			int end_int = std::stoi(end_string);
			if ( end_int > start_int ) {
				opt::chr_choice = cid;
				opt::start_bound = start_int;
				opt::end_bound = end_int;
				opt::region_defined = true;
			}
			else { std::cerr << " region string problem  \n" << EXTRACT_USAGE_MESSAGE; exit(1); }
		}
		else { std::cerr << " region string problem  \n" << EXTRACT_USAGE_MESSAGE; exit(1); }
	}
};
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void parse_extract_hash_options( int argc, char** argv ) {
	bool die = false; //bool vcf_load = false; //bool cov_load = false;
	//if(argc <= 2) { die = true; }
	if (std::string(argv[1]) == "help" || std::string(argv[1]) == "--help") {
		std::cerr << "\n" << EXTRACT_USAGE_MESSAGE;
		exit(1);
	}
        for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
                std::istringstream arg(optarg != NULL ? optarg : "");
                switch (c) {
		case 'h': die = true; break;
                case 'i': arg >> opt::input_bam_file; break;
                case 'v': opt::vcf_load = true; arg >> opt::input_vcf_file; break;
                case 'e': arg >> opt::technology; break;
                case 'c': arg >> opt::chr_choice; break;
                case 's': opt::multiple_bams = true; arg >> opt::second_input_bam_file; break;
                case 'n': arg >> opt::id_string; break;
                case 'd': arg >> opt::refdir; break;
                case 'l': arg >> opt::strlen;
                case 'k': arg >> opt::klen; break;
                }
        }

	parse_region_string( opt::chr_choice );
  std::cout << std::endl;
  std::cout << "############### running link phaser ############### " << std::endl;
  std::cout << "== chromosome                 === " << opt::chr_choice << std::endl;
  std::cout << "== technology                 === " << opt::technology << std::endl;
  std::cout << "== input bam                  === " << opt::input_bam_file << std::endl;
  std::cout << "== reference genome directory === " << opt::refdir << std::endl;
  std::cout << "== kmer length                === " << opt::klen << std::endl;
  if(opt::vcf_load) { std::cout << "== input vcf  === " << opt::input_vcf_file << std::endl; }
  std::cout << "== id string  === " << opt::id_string << std::endl;
  if(opt::multiple_bams) { std::cout << "== second bam === " << opt::second_input_bam_file << std::endl; }
	if ( opt::region_defined ) {
		std::cout << "== start pos  === " << opt::start_bound << std::endl;
		std::cout << "== end pos    === " << opt::end_bound << std::endl;
	}
  std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx " << std::endl;
  std::cout << std::endl;
};

//**************************************************************************
//**************************************************************************

int run_extract_hash_incl_indels(int argc, char *argv[]){

  parse_extract_hash_options(argc, argv);
  std::string output_link_file = "./output/graph_variant_" + opt::id_string + "_" + opt::technology + "_" + opt::chr_choice + ".dat";
  std::string output_hash_file = "./output/graph_hash_" + opt::id_string + "_" + opt::technology + "_" + opt::chr_choice + ".dat";
  std::cout << "- output variant link file: " << output_link_file << std::endl;
  std::cout << "- output hash link file: " << output_hash_file << std::endl;
  std::cout << std::endl;

  read_graph rgraph;
  VariantSite::variant_graph vgraph;

	bool paired;



        //#################### start of code ##########################

  //*****************************
  // Load VCF records from file
  //*****************************
  std::map<std::string,int> chr_str_map;
  IndelGenotyping::contig_name_map(opt::input_bam_file,chr_str_map);
  int chromosome = chr_str_map[opt::chr_choice];


  std::vector<VariantSite::VCFRecord> vvec;

	if( opt::vcf_load ) IndelGenotyping::LoadVCF( opt::input_vcf_file, chromosome, opt::chr_choice, vvec );

	if (opt::region_defined) {
		IndelGenotyping::subset_het_sites( vvec, opt::start_bound, opt::end_bound );
	}

  if (opt::technology == "pacbio" || opt::technology == "tenx" || opt::technology == "illumina" || opt::technology == "nanopore" || opt::technology == "hic") {
    IndelGenotyping::ConnectUpVariantsBamPileup( opt::refdir, vvec, opt::input_bam_file,
                                                 vgraph, rgraph, opt::technology, opt::filter_minima,
                                                 opt::strlen, opt::klen );
    if(opt::multiple_bams) {
      IndelGenotyping::ConnectUpVariantsBamPileup( opt::refdir, vvec, opt::second_input_bam_file,
                                                   vgraph, rgraph, opt::technology, opt::filter_minima,
                                                   opt::strlen, opt::klen );
    }
  }
  else { cout << "error: not a valid technology choice [pacbio,tenx,illumina,nanopore,hic] " << std::endl; return(1); }

  write_variant_link_list(vgraph,rgraph,output_link_file,opt::chr_choice);
  write_hash_link_list(vgraph,rgraph,output_hash_file,opt::chr_choice);
	////////

/*
  BamTools::BamAlignment al;
  int varpos = 65;

  al.Position = 55;  //            |Variant (should be at 15 on read)
  al.AlignedBases =     "AAAAAACTGCAAGCGGGGGCCCCCCC";
  al.QueryBases =  "TTTTTAAAAAACTGCAAGCGGGGGCCCCCCC";
  std:vector<BamTools::CigarOp> cigar = { {'S',5}, {'M',19}, {'S',7}  };
  al.CigarData = cigar;

  std::cout << al.Position << " " << al.QueryBases << " ";
  for(auto x : al.CigarData) std::cout << x.Length << x.Type;
  std::cout << '\n';

  std::cout << IndelGenotyping::VariantPositionOnRead(al, varpos) << '\n';
*/
/*              //                  |
  std::string x = "ATTCGTAGCCGTGATGGCCAATGACGTCCA";
  std::string y = "ATTCGTAGCCGTGATGGCCAATGACGTCCA";
  int k = 6;
  int ksize = IndelGenotyping::GetKmerSize(x, y, 15, k);
  std::cout <<  ksize << " " << x.size() << '\n';
  std::cout << x.substr(15-(ksize/2), ksize) << '\n';
  std::cout << y.substr(15-(ksize/2), ksize) << '\n';
*/
  return(0);
}
