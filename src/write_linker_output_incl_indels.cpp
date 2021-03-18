#include "write_linker_output_incl_indels.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_variant_link_list( std::unordered_map<std::string, VariantSite::variant_node>& var_dict, read_graph& rgraph, std::string outputFile, std::string chr_choice ) {
        ofstream ofile; ofile.open(outputFile);
        for (auto& it : var_dict) {
                std::string variant_string = it.first;
                std::vector<std::string> hash_connections = it.second.connected_reads_long_form;
                std::vector<std::string> hash_readnames = it.second.connected_readnames;
                for (int j = 0; j < hash_connections.size(); j++) {
                        ofile << variant_string << "\t" << chr_choice << "\t" << it.second.pos << "\t" << hash_connections[j] << "\t" << hash_readnames[j] << endl;
                }
        }
        ofile.close();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_hash_link_list( std::unordered_map<std::string, VariantSite::variant_node>& var_dict, read_graph& rgraph, std::string outputFile, std::string chr_choice ) {
        ofstream ofile; ofile.open(outputFile);
        for (auto& it : rgraph) {
        	std::vector<std::string> het_connections = rgraph[it.first].connected_strings;
        	std::vector<std::string> het_readnames = rgraph[it.first].connected_readnames;
        	for (int j = 0; j < het_connections.size(); j++) {
                	ofile << it.first << "\t" << chr_choice << "\t" << het_connections[j] << "\t" << het_readnames[j] << endl;
                }
        }
        ofile.close();
};

