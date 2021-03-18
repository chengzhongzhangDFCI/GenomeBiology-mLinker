#ifndef RUN_EXTRACT_HASH_INCL_INDELS_H
#define RUN_EXTRACT_HASH_INCL_INDELS_H

////////////////////////////////////////////////
#include "variant_site_incl_indels.h"
#include "read_bam_incl_indels.h"
#include "write_linker_output_incl_indels.h"

static void parse_region_string( std::string samtools_region );
static void parse_extract_hash_options( int argc, char** argv );
int run_extract_hash_incl_indels(int argc, char *argv[]);

#endif
