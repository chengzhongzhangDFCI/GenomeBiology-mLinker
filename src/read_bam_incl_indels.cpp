/*
*****************************************************
THIS FUNCTION GENERATES A LIST OF 500-bp REFERENCE
SEQUENCES THAT CONTIAIN THE VARIANTS REPORTED IN AN
INPUT VCF FILE
*****************************************************
*/

#include "read_bam_incl_indels.h"

IndelGenotyping::ReadKmer::ReadKmer() {};
IndelGenotyping::ReadKmer::~ReadKmer() {};

//***********************************************************
/////// bamtools prep /////////////
//***********************************************************
void IndelGenotyping::contig_name_map( std::string inputFilename, std::map<std::string,int> &chr_str_map )  {
  BamTools::BamReader reader;
  open_bam_file(reader,inputFilename);
  BamTools::RefVector refnames = reader.GetReferenceData();
  int bi = 0;  for (auto& it : refnames) { chr_str_map[it.RefName] = bi; bi++; }
}
//************************************************************
//************************************************************
std::string IndelGenotyping::return_index_string( std::string filename ) {
  std::string index_filename;
  index_filename = filename.substr(0,filename.size()-1);
  index_filename += "i";
  return index_filename;
};
//*****************************************************************************
//*****************************************************************************
static void IndelGenotyping::open_bam_file( BamTools::BamReader &reader, std::string inputFilename ) {
  std::string indexFilename;
  indexFilename = return_index_string( inputFilename );
  if (!reader.Open(inputFilename)) { std::cerr << "Could not open input BAM file." << std::endl; return; }
  if (!reader.OpenIndex(indexFilename)) {
    indexFilename = inputFilename + ".bai";
    std::cout << "looking for index file " << indexFilename << std::endl;
    if (!reader.OpenIndex(indexFilename)) { std::cerr << "Could not find BAM index file." << std::endl; return; }
  }
  reader.Open(inputFilename);
  reader.OpenIndex(indexFilename);        //BamTools::BamRegion BamRegion;
	return;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void IndelGenotyping::subset_het_sites( std::vector<VariantSite::VCFRecord> &t_vvec, int t_start_bound, int t_end_bound ) {
	for (int i=0; i < t_vvec.size(); i++) {
		if ((t_vvec[i].GetPos() <= t_start_bound) || (t_vvec[i].GetPos() >= t_end_bound)) {
			t_vvec[i].SetBounded(false);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string IndelGenotyping::technology_hash( std::string tech, BamTools::BamAlignment al ) {
  std::string readhash_long;

  //************************************************************
  // The technology hash depends on the sequencing technology
  // As more technology evols more options (if-statements) are
  // added
  //************************************************************

  if ( tech == "pacbio" || tech == "nanopore" ) { //|| tech == "illumina" || tech == "hic"
    std::string readname,start_pos,end_pos;
    std::vector<std::string> splitted_string;
    splitted_string.reserve(10);
    splitted_string = IndelGenotyping::split_string(al.Name);
    readname = splitted_string.back();          //readname should
    start_pos = std::to_string(al.Position);
    end_pos = std::to_string(al.GetEndPosition());
    readhash_long = readname + "_" + start_pos + "_" + end_pos; // readname should be bxn
  }
	else if ( tech == "illumina" || tech == "hic") {
		std::string readname;
		std::vector<std::string> splitted_string;
		splitted_string.reserve(10);
    splitted_string = IndelGenotyping::split_string(al.Name);
    readname = splitted_string.back();
		readhash_long = readname;
	}
  else if ( tech == "tenx" ) {
    std::string bxtag = "BX";
    std::string return_bx_string;
    al.GetTag(bxtag,return_bx_string);
    readhash_long = return_bx_string;
  }
  return readhash_long;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<std::string> IndelGenotyping::split_string( std::string teststring ) {
        std::stringstream test(teststring);
        std::string segment;
        std::vector<std::string> seglist;
        while(std::getline(test,segment,'/')) { seglist.push_back(segment); }
        return seglist;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
bool IndelGenotyping::get_tag_value( const std::string &t_tag, VariantSite::variant_graph &t_variant_cnx ) {

//######################################################
// get_tag_value() checks if t_tag is in t_variant_cnx
// #####################################################
        auto check = t_variant_cnx.find(t_tag);
        if (check == t_variant_cnx.end()) return false;
        return true;
};

//////////////////////////////////////////////////////////////////////////////////
int IndelGenotyping::QualityCharToInt( const char t_ch ){
  int xord = t_ch;
  int qscore = (xord-33);
  return(qscore);
}

//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////

char IndelGenotyping::Complement(const char t_c)
{
 if(t_c == 'A') return('T');
 if(t_c == 'T') return('A');
 if(t_c == 'C') return('G');
 if(t_c == 'G') return('C');
 return('N');
}

std::string IndelGenotyping::RevComplement(const std::string &t_seq)
{
 std::string revc;
 for(auto c : t_seq)
 {
  revc.insert(revc.begin(),IndelGenotyping::Complement(c));
 }
 return(revc);
}

/////////////////////////////////////////////////////////////////////////////
bool IndelGenotyping::PassedBaseQFilter( IndelGenotyping::ReadKmer &t_kmer, std::string &t_technology, IndelGenotyping::FilterMinima &t_filter_minima )
{
  // If t_kmer is empty return false
  //.................................
  if( t_kmer.Seq.empty() ) return(false);

  // Compute the median of the base qualities in t_kmer
  //....................................................
  std::vector<int> vqualities;
  int n(0);
  for( auto ch : t_kmer.Qualities ){
    vqualities.push_back(IndelGenotyping::QualityCharToInt(ch));
    n++;
  }
  std::sort(vqualities.begin(), vqualities.end());

  double bqual(0.);
  if( (n%2) == 0 ) bqual = (vqualities[(n/2)-1] + vqualities[n/2])/2.;
  else bqual = (double)vqualities[n/2];
  //..................................................

  if ( t_technology == "hic" ) {
    if ( bqual >= t_filter_minima.Minimum_baseq_hic ) return(true);
  }
  else if ( t_technology == "tenx" ) {
    if ( bqual >= t_filter_minima.Minimum_baseq_tenx ) return(true);
  }
  else {
    if ( bqual >= t_filter_minima.Minimum_baseq ) return(true);
  }
  return(false);

}

////////////////////////////////////////////////////////////////////////////////////
bool IndelGenotyping::PassedMapQFilter( BamTools::BamAlignment &t_al, const int t_samsum_lt, const int t_minimum_mapq ) {
  if (t_al.IsPrimaryAlignment()) {
    if (t_al.AlignmentFlag < t_samsum_lt) {
      if (!t_al.IsDuplicate()) {
        if (t_al.MapQuality > t_minimum_mapq) {
          return(true);
        }
      }
    }
  }
  return(false);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
bool IndelGenotyping::PassedMapQFilterHiC( BamTools::BamAlignment &t_al, const int t_minimum_mapq_hic ) {
        bool return_bool = false;
  if (!t_al.IsDuplicate()) {
    if (t_al.MapQuality > t_minimum_mapq_hic) {
      return(true);
		}
	}
  return(false);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

bool IndelGenotyping::PassedAlignmentFilter( BamTools::BamAlignment &t_al, std::string &t_technology, IndelGenotyping::FilterMinima &t_filter_minima )
{
  if( t_technology == "highc" ) return( IndelGenotyping::PassedMapQFilterHiC( t_al, t_filter_minima.Minimum_mapq_hic ) );
  else  return(IndelGenotyping::PassedMapQFilter( t_al, t_filter_minima.Samsum_lt, t_filter_minima.Minimum_mapq ) );
}


//***********************************************************************************
// Split a string t_s at the locations of character t_delim. Store sustrings in t_v
//***********************************************************************************
void IndelGenotyping::split(const std::string& t_s, char t_delim, std::vector < std::string >& t_v)
{
    std::size_t i = 0;
    std::size_t pos = t_s.find(t_delim);

    if(pos == std::string::npos) t_v.push_back(t_s);
    else
    {
     while (pos != std::string::npos) {
      t_v.push_back(t_s.substr(i, pos-i));
      i = ++pos;
      pos = t_s.find(t_delim, pos);

      if (pos == std::string::npos)
         t_v.push_back(t_s.substr(i, t_s.length()));
     }
    }
}


int IndelGenotyping::LoadVCF( std::string t_vcfname, int t_ChrID, std::string t_Chr,
                              std::vector<VariantSite::VCFRecord> &t_vcfs )
{
  bcf_srs_t *sr = bcf_sr_init();
  bcf_sr_add_reader(sr,t_vcfname.c_str());
  bcf1_t *line;
  int i = 0;
  while(bcf_sr_next_line(sr)) {
    line = bcf_sr_get_line(sr,0);
    if (line->rid == t_ChrID) {
			std::string chr = t_Chr;
			int chrID = line->rid;
			int chrlength = line->rlen;
 			int pos = line->pos;

			std::vector<std::string> valt;
			char **x = line->d.allele;
      std::string refallele = std::string(x[0]);
			for(auto i = 1; i < line->n_allele; ++i )valt.push_back(std::string(x[i]));
			for(auto i = 0; i < valt.size(); ++i){
        if( ( valt[i].size() == 1 )&&( (int)valt[i][0] == 42 ) ) valt[i] = "D"; // Substitute the unprintable character '*' by 'D' (deletion)
			}
      int qual = line->qual;

      VariantSite::VCFRecord vcf;
      vcf.SetVCFRecord(chr, chrID, chrlength, pos, refallele, valt, qual);
      t_vcfs.push_back(vcf);

		}
  }

  return(0);
}

//********************************
//********************************
std::size_t IndelGenotyping::MaxStringSize( const std::vector<std::string> &t_v ){

  if(t_v.size() == 1) return(t_v[0].size());

  std::size_t smax = t_v[0].size();
  for( auto x : t_v ){
    std::size_t s = x.size();
    if(s > smax) smax = s;
  }
  return(smax);
}
//*********************************************
//*********************************************
std::string IndelGenotyping::Capitalize(const std::string &t_s){

  std::string r;

  std::map<char, char> table;

  table['A'] = 'A';
  table['C'] = 'C';
  table['G'] = 'G';
  table['T'] = 'T';
  table['D'] = 'D';

  table['a'] = 'A';
  table['c'] = 'C';
  table['g'] = 'G';
  table['t'] = 'T';
  table['d'] = 'D';

  for(auto x : t_s) r.push_back(table[x]);

  return(r);

}

//********************************************
//********************************************
void IndelGenotyping::GenAltRef( VariantSite::VCFRecord &t_vcf,
                                  faidx_t *t_faidx, int t_strlen ){

  int p1 = t_vcf.GetPos() + 1 - t_strlen, p2 = t_vcf.GetPos() + 1;

  std::string rseq = t_vcf.GetRef();
  std::size_t m = IndelGenotyping::MaxStringSize(t_vcf.GetValt()), rsize = rseq.size();

  std::string sreg1 = t_vcf.GetChrom() + ':' + std::to_string(p1) + '-' + std::to_string(p1 + t_strlen - 1 );
  std::string sreg2 = t_vcf.GetChrom() + ':' + std::to_string(p2) + '-' + std::to_string(p2 + t_strlen + m + rsize );

  int* len1 = &t_strlen;
  int  ilen2 = t_strlen + m + rsize;
  int* len2 = &ilen2;
  char *s01, *s02;

  s01 = fai_fetch( t_faidx, sreg1.c_str(), len1 );
  s02 = fai_fetch( t_faidx, sreg2.c_str(), len2 );

  std::string s1(s01), s2;
  if( s02 != NULL ) s2 = std::string(s02);
  free(s02);
  free(s01);

  t_vcf.SetRefFragment( IndelGenotyping::Capitalize(s1+s2.substr(0,t_strlen)) ); // Store the reference fragment sequence

  if( !s2.empty() ){
    for(auto x : t_vcf.GetValt()){
      std::string sa;
      std::size_t xsize = x.size();
      if(x == "D") sa = s1 + s2.substr(rsize, t_strlen);
      else{
        if( xsize <= t_strlen ) sa = s1 + x + s2.substr(rsize, t_strlen - xsize);
        else sa = s1 + x.substr(0,t_strlen);
      }
      t_vcf.SetAltFragment( IndelGenotyping::Capitalize(x), IndelGenotyping::Capitalize(sa)); // Store the Alternative fragment sequence
    }
  }
  else{
    std::cout << "Region is outside genome\n";
  }
}

//*************************************
//***************************************
int IndelGenotyping::HammingDistance(const std::string &t_s1, const std::string &t_s2)
{
 int dis(0);
 std::size_t size1 = t_s1.size(), size2 = t_s2.size();
 try
 {
  if(size2 != size1){
   throw std::runtime_error("HammingDistance: String lengths are different\n");
  }
 }
 catch (std::runtime_error err)
 {
  std::cout << err.what();
  return(-1);
 }
 for(auto i = 0; i < size1; ++i)
 {
  if(t_s2[i] != t_s1[i]) dis ++;
 }
 return(dis);

}

//******************************************************
//******************************************************
int IndelGenotyping::VariantPositionOnRead( BamTools::BamAlignment &t_al, const int t_varpos )
{
  int q(0), pos(t_al.Position), flag(0);
  //std::cout << "\npos = " << pos << ", q = " << q << '\n';
  for( auto &cg : t_al.CigarData ){
    if( ( (cg.Type == 'S')||(cg.Type == 'H') )&&(flag == 0) ){
      q += cg.Length;
      flag = 1;
    }
    else{
      if ( (pos + cg.Length >= t_varpos)&&(cg.Type != 'I') ){
        if( cg.Type != 'D' ) q += (t_varpos - pos);
        //std::cout << cg.Length << cg.Type << "Return q = " << q << " at pos = " << pos + (t_varpos - pos) << '\n';
        return(q);
      }
      else{
        if( cg.Type != 'I') pos += cg.Length;
        if( cg.Type != 'D' ) q += cg.Length;
      }
    }
    //std::cout << cg.Length << cg.Type << ": pos = " << pos << ", q = " << q << '\n';
  }
  std::cout << "Error: Read " << t_al.Name << " does not overlap with the variant at positio " << t_varpos << '\n';
  return(q);
}

//************************************************
//************************************************
//bool IndelGenotyping::DoesReadSupportVariant( const std::string &t_v, const std::string &t_seq, const std::string &t_qualities,
//                                              IndelGenotyping::ReadKmer &t_kmer, int &t_flag, int t_x, int t_k, int t_d )
bool IndelGenotyping::DoesReadSupportVariant( const std::string &t_v, const std::string &t_seq, const std::string &t_seq_rc,
                                              const std::string &t_qualities, const std::string &t_qualities_rc, const int &t_ssize,
                                              IndelGenotyping::ReadKmer &t_kmer, const int &t_vpos_onread, int t_x, int t_k, int t_d )
{

  // Extract from t_v the k-mer to find in the read with sequence t_seq
  ///////////////////////////////////////////////////////////////////////
  //std::string seq_rc = IndelGenotyping::RevComplement(seq); // Generate the reverse complement of t_seq
  //std::string qualities_rc;

  // See if read supports kmer at the variant site
  //////////////////////////////////////////////////
  int window = 0; //used 2 before
  std::string pattern = t_v.substr(t_x-(t_k/2), t_k); //t_v.substr(x-(t_k/2), t_k);
  int wa = std::min( (t_k/2) + window, t_vpos_onread);
  int imax = std::min((t_vpos_onread-(t_k/2)+window), t_ssize-t_k);
  for(auto i = t_vpos_onread-wa; i <= imax; ++i )
  {
    if( IndelGenotyping::HammingDistance(t_seq.substr(i, t_k), pattern) <= t_d ){
      t_kmer.Seq = pattern;
      t_kmer.Qualities = t_qualities.substr(i, t_k);
      //t_flag = 1;
      //std::cout << "\nRead k-mer: " << t_kmer.Seq << '\n';
      //std::cout << "Var. k-mer: " << pattern << '\n';
      return(true);
    }

    if( IndelGenotyping::HammingDistance(t_seq_rc.substr(i, t_k), pattern) <= t_d ){
      t_kmer.Seq = pattern;
      t_kmer.Qualities = t_qualities_rc.substr(i, t_k);
      //t_flag = 1;
      //std::cout << "\nRead k-mer (reversed): " << t_kmer.Seq << '\n';
      //std::cout << "Var. k-mer             : " << pattern << '\n';
      return(true);
    }
  }

  return(false);
}

//************************************************
//************************************************
void IndelGenotyping::UpdateGraphs( const int t_pos, const std::string t_vbases, const std::string t_rbases,
                                    VariantSite::variant_graph &t_variant_cnx, read_graph &t_rgraph,
                                    std::string t_technology, BamTools::BamAlignment &t_al,
                                    const std::string t_kmer, bool t_is_variant )
{

  //################################################################################
  // UpdateGraphs():
  //################
  // Updates the variant graph object t_variant_cnx and
  // the read graph object t_rgraph
  //................................................................................
  //................................................................................
  // t_pos        : Variant position
  // t_vbases     : Variant base/s at t_pos
  // t_rbases     : Reference base/s at t_pos
  // technology   : Sequencing technology used
  // t_al         : Alignment record of read overlapping variant site
  // t_kmer       : kmer to search in the read (either reference or variatn kmer)
  // t_is_variant : true if read supports the variant kmer
  //################################################################################
   std::string readhash_long = IndelGenotyping::technology_hash(t_technology, t_al);

  if ( readhash_long.size() > 1 ) {

    std::string variant_id = std::to_string(t_pos) + "_" + t_rbases + "_" + t_vbases;

    bool present = get_tag_value(variant_id, t_variant_cnx);   // Check if variant_id is in t_variant_cnx already

    if(!present) { // If variant_id is not in t_variant_cnx, then create its corresponding node, just add connection otherwise
      VariantSite::variant_node v_node;
      t_variant_cnx[variant_id] = v_node;
      t_variant_cnx[variant_id].set_values(t_pos, t_is_variant, t_vbases, t_rbases);
    }

    std::string het_hash = std::to_string(t_pos) + "_" + t_rbases + "_" + t_kmer;
    if ( t_rgraph.find(readhash_long) != t_rgraph.end() ) t_rgraph[readhash_long].add_connection(het_hash, t_al.Name);
    else {
      t_rgraph[readhash_long].set_name(readhash_long);
      t_rgraph[readhash_long].add_connection(het_hash, t_al.Name);
    }
    t_variant_cnx[variant_id].add_connected_read(readhash_long,t_al.Name);

  }

}
//*****************************************************
//*****************************************************
int IndelGenotyping::GetKmerSize( const std::string &t_seq1, const std::string &t_seq2, const int &t_middle_point, const int &t_k )
{
  int ksize(t_k);
  while( ksize <= 2*t_middle_point ){ // ksize <= t_seq1/2 size
    std::string patt1 = t_seq1.substr(t_middle_point-(ksize/2), ksize);
    std::string patt2 = t_seq2.substr(t_middle_point-(ksize/2), ksize);
    if ( patt1 != patt2 ) return(ksize);
    else ksize++;
  }
  return(ksize);
}
//******************************************************
//******************************************************
void IndelGenotyping::FindReads( VariantSite::VCFRecord &t_vcf, BamTools::BamReader &t_reader,
                                 VariantSite::variant_graph &t_variant_cnx, read_graph &t_rgraph,
                                 std::string t_technology, IndelGenotyping::FilterMinima &t_filter_minima,
                                 int t_k )
{
  std::string chr = t_vcf.GetChrom();
  int refID = t_reader.GetReferenceID(chr), pos = t_vcf.GetPos();
  BamTools::RefVector vreferences = t_reader.GetReferenceData();

  int reflength =  vreferences[refID].RefLength, leftboundary(0), rightboundary(0);

  if( pos >= 2) leftboundary = pos - 2;
  else leftboundary = 0;
  if( pos <= (reflength-2) ) rightboundary = pos + 2;
  else rightboundary = reflength;


  BamTools::BamRegion region(refID, leftboundary, refID, rightboundary);
  t_reader.SetRegion(region);

  ////// Get the extracted reference and variant fragments  /////////
  ///////////////////////////////////////////////////////////////////
  std::map<std::string, std::string> varseqvec = t_vcf.GetAltFragments();
  std::string rseq = t_vcf.GetRefFragment();
  std::string rbases = t_vcf.GetRef();
  ///////////////////////////////////////////////////////////////////

  std::size_t p = rseq.size() / 2;
  BamTools::BamAlignment al;
  while (t_reader.GetNextAlignmentCore(al)) { // Loop over reads
    if ( (al.Position < pos) && (al.GetEndPosition() > pos) ) { // Only collect reads that overlap the variant
      al.BuildCharData();
      if ( IndelGenotyping::PassedAlignmentFilter( al, t_technology, t_filter_minima ) )
      {
        int seq_size(0); // QueryBases size
        std::string seq_rc, qualities_rc; // Reverse Complement of QueryBases and qualities
        for( auto &ch : al.Qualities ){
          qualities_rc.insert(qualities_rc.begin(),ch); // Generate reverse of t_qualities
          seq_rc.insert(seq_rc.begin(),IndelGenotyping::Complement(al.QueryBases[seq_size])); // Revese complement seq
          seq_size ++;
        }
        int varpos_onread = IndelGenotyping::VariantPositionOnRead( al, pos ); // Get position of variant site on read
        bool does_support = false;
        bool passed_baseq = false;
        IndelGenotyping::ReadKmer kmer;
        int kmax(t_k);
        for(auto &v : varseqvec ) { // Loop over variant fragments
          int ksize = IndelGenotyping::GetKmerSize( rseq, v.second, p, t_k );
          if(ksize > kmax) kmax = ksize;
          if( ksize <= 2*p ){ // If the kmer size is not larger than rseq (fragment) size
            does_support = IndelGenotyping::DoesReadSupportVariant( v.second, al.QueryBases, seq_rc,
                                                                    al.Qualities, qualities_rc, seq_size, kmer, varpos_onread, p, ksize, 0 );
            passed_baseq = IndelGenotyping::PassedBaseQFilter( kmer, t_technology, t_filter_minima );
          }
          if( does_support && passed_baseq ){
            t_vcf.Print();
            IndelGenotyping::UpdateGraphs( pos, v.first, rbases, t_variant_cnx, t_rgraph, t_technology, al, kmer.Seq, true );

            std::cout << "Read " << al.Name << " at " << al.Position << " supports variant " << v.first
                      << " at variant position " << varpos_onread << " on read" << '\n';
            std::cout << "Ref. fragment : " << rseq << '\n';
            std::cout << "Var. fragment : " << v.second << '\n';
            std::cout << "Read k-mer: " << kmer.Seq << '\n';
            std::cout << "Ref. k-mer: " << rseq.substr(p-(ksize/2),ksize) << '\n';
            std::cout << "Var. k-mer: " << v.second.substr(p-(ksize/2),ksize)  << '\n';

            does_support = false; // Reset them to false
            passed_baseq = false;
          }
        }
        // For the reference kmer
        if ( kmax <= 2*p ){
          does_support = IndelGenotyping::DoesReadSupportVariant( rseq, al.QueryBases, seq_rc, al.Qualities, qualities_rc, seq_size,
                                                                  kmer, varpos_onread, p, kmax, 0 );
          passed_baseq = IndelGenotyping::PassedBaseQFilter( kmer, t_technology, t_filter_minima );
        }
        if( does_support && passed_baseq ){
          t_vcf.Print();
          IndelGenotyping::UpdateGraphs( pos, rbases, rbases, t_variant_cnx, t_rgraph, t_technology, al, kmer.Seq, false );

          std::cout << "Read " << al.Name << " at " << al.Position << " supports reference " << rbases
                    << " at variant position " << varpos_onread << " on read" << '\n';
          std::cout << "Ref. fragment: " << rseq << '\n';
          std::cout << "Read k-mer: " << kmer.Seq << '\n';
          std::cout << "Ref. k-mer: " << rseq.substr(p-(kmax/2),kmax) << '\n';
        }
      }
    }
  }
}

//********************************************
//********************************************
int IndelGenotyping::ConnectUpVariantsBamPileup( const std::string t_refdir,
          std::vector<VariantSite::VCFRecord> &t_vcfs, const std::string t_bamname,
          VariantSite::variant_graph &t_variant_cnx, read_graph &t_rgraph,
          std::string t_technology, IndelGenotyping::FilterMinima &t_filter_minima, int t_strlen, int t_k )
{

  BamTools::BamReader reader;
	open_bam_file(reader,t_bamname);
  BamTools::BamAlignment al;

  //********************************************************
  // Open the fasta for the chromosome in first vcf record
  //********************************************************

  /// Get the index of first record ///
  /////////////////////////////////////
  int i0(0);
  for( auto v : t_vcfs ){
    if( v.GetBounded() ) break;
    i0++;
  }
  /////////////////////////////////////
  std::string chr = t_vcfs[i0].GetChrom();

  std::string faname = t_refdir + chr + ".fa";

  // Load index file
  //.................
  faidx_t *faidx;
  faidx =  fai_load( faname.c_str() );

  // Generate the alternative reference sequences
  // for the first vcf record
  //..............................................
  IndelGenotyping::GenAltRef( t_vcfs[i0], faidx, t_strlen/2 );
  IndelGenotyping::FindReads( t_vcfs[i0], reader, t_variant_cnx, t_rgraph, t_technology, t_filter_minima, t_k );

  for( auto i = i0+1; i < t_vcfs.size(); ++i ){

    if( t_vcfs[i].GetBounded() ){
      if(t_vcfs[i].GetChrom() != chr){
        fai_destroy(faidx);
        faidx = NULL;
        chr = t_vcfs[i].GetChrom();
        faname = t_refdir + chr + ".fa";

        faidx =  fai_load( faname.c_str() );
      }
      IndelGenotyping::GenAltRef( t_vcfs[i], faidx, t_strlen/2 );
      IndelGenotyping::FindReads( t_vcfs[i], reader, t_variant_cnx, t_rgraph, t_technology, t_filter_minima, t_k );
    }
  }

  fai_destroy(faidx);
  faidx = NULL;

  reader.Close();

  return(0);
}

