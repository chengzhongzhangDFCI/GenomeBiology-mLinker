#ifndef READ_BAM_INCL_INDELS_H
#define READ_BAM_INCL_INDELS_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

//////////////// linker include ////////////////
#include "read_tree.h"
#include "variant_site_incl_indels.h"
//#include "coord_dict.h"

//////////////// bamtools //////////////////////
#include "api/BamWriter.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "api/BamAlgorithms.h"
#include "api/BamMultiReader.h"
#include "api/BamAux.h"


//////////////// htslib ////////////////////////
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>

#include <htslib/faidx.h>

extern "C" {
  #include "htslib/vcf.h"
  #include "htslib/vcfutils.h"
  #include "htslib/synced_bcf_reader.h"
  #include "htslib/faidx.h"
}

namespace IndelGenotyping{

  struct FilterMinima{
      int Samsum_lt = 2000;
      int Minimum_mapq = 20;                 // Minimum read map quality
      int Minimum_mapq_hic = 20;             // Minimum read map quality for HiC technology
      int Minimum_baseq = 8;                 // Minimum read base quality
      int Minimum_baseq_hic = 20;            // Minimum read base quality HiC technology
      int Minimum_baseq_tenx = 20;           // Minimum read base quality 10X technology
  };

  class ReadKmer{
    public:
      ReadKmer();
      ~ReadKmer();

      std::string Seq;
      std::string Qualities;
  };

  void contig_name_map( std::string inputFilename, std::map<std::string,int> &chr_str_map );

  std::string return_index_string( std::string filename ); // From Rick's code

  static void open_bam_file( BamTools::BamReader &reader, std::string inputFilename ); // From Rick's code

  void subset_het_sites( std::vector<VariantSite::VCFRecord> &t_vvec, int t_start_bound, int t_end_bound );

  std::string technology_hash( std::string tech, BamTools::BamAlignment al );

  std::vector<std::string> split_string( std::string teststring );

  bool get_tag_value( const std::string &t_tag, VariantSite::variant_graph &t_variant_cnx );

  /////////////////////////////////////////////
  int QualityCharToInt( const char t_c );

  char Complement(const char t_c);

  std::string RevComplement(const std::string &t_seq);

  bool PassedBaseQFilter( ReadKmer &t_kmer, std::string &t_technology, FilterMinima &t_filter_minima );

  bool PassedMapQFilter( BamTools::BamAlignment &t_al, int t_samsum_lt = 0, int t_minimum_mapq = 0);

  bool PassedMapQFilterHiC( BamTools::BamAlignment &t_al, int t_minimum_mapq_hic = 0 );

  bool PassedAlignmentFilter( BamTools::BamAlignment &t_al, std::string &t_technology, FilterMinima &t_filter_minima );

  void split(const std::string& t_s, char t_delim, std::vector < std::string >& t_v);

  int LoadVCF( std::string t_vcfname, int t_ChrID, std::string t_Chr, std::vector<VariantSite::VCFRecord> &t_vcfs );

  std::size_t MaxStringSize( const std::vector<std::string> &t_v );

  std::string Capitalize(const std::string &t_s);

  void GenAltRef( VariantSite::VCFRecord &t_vcf, faidx_t *t_faidx, int t_strlen = 0);

  int GetKmerSize( const std::string &t_seq1, const std::string &t_seq2, const int &t_middle_point, const int &t_k );

  void FindReads( VariantSite::VCFRecord &t_vcf, BamTools::BamReader &t_reader, VariantSite::variant_graph &t_variant_cnx, read_graph &t_rgraph,
                  std::string t_technology, FilterMinima &t_filter_minima, int t_k = 0 );

  int HammingDistance(const std::string &t_s1, const std::string &t_s2);

  int VariantPositionOnRead( BamTools::BamAlignment &t_al, const int t_varpos );

  bool DoesReadSupportVariant( const std::string &t_v, const std::string &t_seq, const std::string &t_seq_rc, const std::string &t_qualities,
                               const std::string &t_qualities_rc, const int &t_ssize,
                               ReadKmer &t_kmer, const int &t_vpos_onread, int t_x = 0, int t_k = 0, int t_d = 0 );

  void UpdateGraphs( const int t_pos, const std::string t_vbases, const std::string t_rbases, VariantSite::variant_graph &t_variant_cnx,
                     read_graph &t_rgraph, std::string t_technology, BamTools::BamAlignment &t_al, const  std::string t_kmer, bool t_is_variant  );

  int ConnectUpVariantsBamPileup( const std::string t_refdir,
          std::vector<VariantSite::VCFRecord> &t_vcfs, const std::string t_bamname,
          VariantSite::variant_graph &t_variant_cnx, read_graph &t_rgraph,
          std::string t_technology, FilterMinima &t_filter_minima, int t_strlen = 0, int t_k = 0);
}
#endif
