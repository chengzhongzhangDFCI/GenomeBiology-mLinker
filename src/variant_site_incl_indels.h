#ifndef VARIANT_SITE_INCL_INDELS_H
#define VARIANT_SITE_INCL_INDELS_H

//////////////// c++ include //////////////////
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <iostream>

namespace VariantSite {

  class VCFRecord{
    public:
      VCFRecord();
      ~VCFRecord();

      // Class functions
      void SetVCFRecord(const std::string t_chrom , const int t_chromID, const int t_chromlength,
                        const int t_pos, const std::string t_ref, const std::vector<std::string> t_valt,
                        const float t_qual );
      void SetChrom( const std::string t_chrom );
      void SetChromID( const int t_chromID );
      void SetChromLength( const int t_chromlength );
      void SetPos( const int t_pos );
      void SetRef( const std::string t_ref );
      void SetValt( const std::string t_alt );
      void SetQual( const float t_qual );
      void SetBounded( const bool t_bounded );
      void SetRefFragment( const std::string t_refragment );
      void SetAltFragment( const std::string t_Alt, const std::string t_altfragment );

      std::string GetChrom();
      int GetChromID();
      int GetChromLength();
      int GetPos();
      std::string GetRef();
      std::string GetAlt();
      std::vector<std::string> GetValt();
      float GetQual();
      bool GetBounded();
      std::string GetRefFragment();
      std::map<std::string, std::string> GetAltFragments();
      std::map<std::string, int> GetKmerSize();
      void FindKmerSize( int t_k = 0 );
      void Print();

    private:
    // Class members
      std::string Chrom;
      int ChromID = 0;
      int ChromLength = 0;
      int Pos =0;
      std::string Ref;
      std::vector<std::string> Valt;
      float Qual = 0.;
      bool Bounded = true;
      std::string RefFragment;
      std::map<std::string, std::string> AltFragments;
      std::map<std::string, int> KmerSize;

  };

  /////////////// structures ///////////////////
  struct vcf_entry {
      bool bounded;
      int pos,chromosome_id;
      std::string ref_base,var_base;
  };

///////////////////////////////////////////////
  class variant_node {
    public:

      bool var = false;  // Whether the variant site has the varian (true) or the reference (false) bases
      int pos = 0;   // Variant position
      std::string var_base;                                // Variant base/s at position pos
      std::string ref_base;                                // Reference base/s at position pos
      std::string variant_id;                              // Variant identifier
      std::vector<std::string> connected_readnames;        // List of reads connected to the node
      std::vector<std::string> connected_reads_long_form;  // List of connected reads with long labels

      // Member functions //
      //////////////////////
      //void set_values(int, bool, std::string ,std::string );
      void set_values( const int t_position, const bool t_variant, const std::string t_variant_base, const std::string t_reference_base );
      void add_connected_read( std::string, std::string );
  };

  //////////////// definitions //////////////////
  typedef std::unordered_map<std::string, variant_node> variant_graph;
  typedef std::vector<vcf_entry> vcf_vector;

}
#endif  // VARIANT_SITE_H
