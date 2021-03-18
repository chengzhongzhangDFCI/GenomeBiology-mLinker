#ifndef WRITE_LINKER_OUTPUT_INCL_INDELS_H
#define WRITE_LINKER_OUTPUT_INCL_INDELS_H

//////////////// c++ include //////////////////
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <iostream>
#include <fstream>
#include <set>
using namespace std;

//////////////// linker include ////////////////
#include "variant_site_incl_indels.h"
#include "read_tree.h"

/////////////// functions /////////////////////
void write_variant_link_list( std::unordered_map<std::string, VariantSite::variant_node> &var_dict, read_graph &rgraph, std::string outputFile, std::string chr_choice );
void write_hash_link_list( std::unordered_map<std::string, VariantSite::variant_node> &var_dict, read_graph &rgraph, std::string outputFile, std::string chr_choice );

#endif  // WRITE_LINKER_OUTPUT_H
