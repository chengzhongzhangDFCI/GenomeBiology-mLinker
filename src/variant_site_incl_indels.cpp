#include "variant_site_incl_indels.h"

//************************
// VCFRecord Class
//************************
VariantSite::VCFRecord::VCFRecord() {};
VariantSite::VCFRecord::~VCFRecord() {};

void VariantSite::VCFRecord::SetVCFRecord(const std::string t_chrom , const int t_chromID, const int t_chromlength,
                        const int t_pos, const std::string t_ref, const std::vector<std::string> t_valt,
                        const float t_qual ){


 Chrom = t_chrom;
 ChromID = t_chromID;
 ChromLength = t_chromlength;
 Pos = t_pos;
 Ref = t_ref;
 Valt = t_valt;
 Qual = t_qual;
}

void VariantSite::VCFRecord::SetChrom( const std::string t_chrom ){
  Chrom = t_chrom;
}
void VariantSite::VCFRecord::SetChromID( const int t_chromID ){
  ChromID = t_chromID;
}
void VariantSite::VCFRecord::SetChromLength( const int t_chromlength ){
  ChromLength = t_chromlength;
}
void VariantSite::VCFRecord::SetPos( const int t_pos ){
  Pos = t_pos;
}
void VariantSite::VCFRecord::SetRef( const std::string t_ref ){
  Ref = t_ref;
}
void VariantSite::VCFRecord::SetValt( const std::string t_alt ){
  Valt.push_back(t_alt);
}
void VariantSite::VCFRecord::SetQual( const float t_qual ){
  Qual = t_qual;
}
void VariantSite::VCFRecord::SetBounded( const bool t_Bounded ){
  Bounded = t_Bounded;
}
void VariantSite::VCFRecord::SetRefFragment( const std::string t_refragment ){
  RefFragment = t_refragment;
}
void VariantSite::VCFRecord::SetAltFragment( const std::string t_Alt, const std::string t_altfragment ){
  AltFragments[t_Alt] = t_altfragment;
}

std::string VariantSite::VCFRecord::GetChrom(){ return(Chrom); }
int VariantSite::VCFRecord::GetChromID(){ return(ChromID); }
int VariantSite::VCFRecord::GetChromLength(){ return(ChromLength); }
int VariantSite::VCFRecord::GetPos(){ return(Pos); }
std::string VariantSite::VCFRecord::GetRef(){ return(Ref); }
std::vector<std::string> VariantSite::VCFRecord::GetValt(){ return(Valt); }
float VariantSite::VCFRecord::GetQual(){ return(Qual); }
bool VariantSite::VCFRecord::GetBounded(){ return(Bounded); }
std::string VariantSite::VCFRecord::GetRefFragment(){ return(RefFragment); }
std::map<std::string, std::string> VariantSite::VCFRecord::GetAltFragments(){ return(AltFragments); }
std::map<std::string, int> VariantSite::VCFRecord::GetKmerSize() { return(KmerSize); }
void VariantSite::VCFRecord::FindKmerSize( int t_k ){

  std::size_t middle_point = RefFragment.size()/2;
  KmerSize["Ref"] = 0; // Kmer size for the reference
  for ( auto v : Valt ){
    KmerSize[v] = 0; // Initialize kmer size for variant v
    int ksize = t_k;
    while( ksize <= 2*middle_point ){
      std::string patt1 = RefFragment.substr(middle_point-(ksize/2), ksize);
      std::string patt2 = AltFragments[v].substr(middle_point-(ksize/2), ksize);
      if ( patt1 != patt2 ){
        KmerSize[v] = ksize;
        if ( ksize > KmerSize["Ref"] ) KmerSize["Ref"] = ksize;
        break;
      }
      else ksize++;
    }
  }
}

void VariantSite::VCFRecord::Print(){
  std::cout << "#CHROM\tPOS\tREF\tALT\tQUAL" << '\n';
  std::cout << Chrom << '\t' << Pos << '\t' << Ref << '\t';
  for(auto x : Valt) std::cout << x << ", ";
  std::cout << '\t' << Qual << '\n';
}


void VariantSite::variant_node::set_values( const int t_position, const bool t_variant, const std::string t_variant_base, const std::string t_reference_base ) {
    pos = t_position;                 // int
    var = t_variant;                  // bool
    var_base = t_variant_base;        // std::string
    ref_base = t_reference_base;      // std::string
    connected_reads_long_form.reserve(100);
    /*
    base_dict['A'] = 0;
    base_dict['T'] = 0;
    base_dict['C'] = 0;
    base_dict['G'] = 0;
    base_dict['D'] = 0;
    base_dict['I'] = 0;
    */
    //total_bases = 0;
    //filter = false;
    if (t_variant == true) { variant_id = std::to_string(t_position) + "_" + t_reference_base + "_" + t_variant_base; }
    else {                 variant_id = std::to_string(t_position) + "_" + t_reference_base + "_" + t_reference_base; }
}

/*
void VariantSite::variant_node::add_base_to_dict( std::string base, std::string readhash_long ) {
    if( base_dict.find(base) != base_dict.end() )base_dict[base] += 1;
    else base_dict[base] = 0;
    total_bases += 1;
    if (readhash_long != "nohash" ) {
    	base_dict_tags[base].push_back(readhash_long) ;
    	base_dict_set[base].insert(readhash_long) ;
    }
}


void VariantSite::variant_node::count_connections() {
    num_hets = connections.size();
    num_reads = connected_reads_long_form.size();
}

void VariantSite::variant_node::add_connection( std::string hashname ) { connections[hashname] += 1; };
*/
void VariantSite::variant_node::add_connected_read( std::string readhash_long, std::string readname ) {
    connected_reads_long_form.push_back(readhash_long);
    connected_readnames.push_back(readname);
}
