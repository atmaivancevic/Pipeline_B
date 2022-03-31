#!/bin/bash

# Append species name to all sequence headers
#
# Example usage:
# ELEMENT=<repeat_type> SPECIES=<species_name> bash 1d_append_name_to_headers.sh

cat "$SPECIES"_"$ELEMENT"_combined.fasta | sed 's/^>/>'$SPECIES'_/g' > "$SPECIES"_"$ELEMENT"_final.fasta
