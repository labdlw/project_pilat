#!/bin/bash

# Create base folders
mkdir -p wet-meadow-connectivity/{data/{metadata,genetic/{filtered_vcf,structure_input},sdm/{final_habitat_maps},connectivity/{resistance_maps,graphab_outputs}},scripts/{radseq_pipeline,population_genetics,sdm_modeling,connectivity_analysis,local_adaptation},results/{figures,tables},manuscript,docs}

# Create starter files
touch wet-meadow-connectivity/README.md
touch wet-meadow-connectivity/LICENSE
touch wet-meadow-connectivity/.gitignore
touch wet-meadow-connectivity/data/metadata/sample_metadata.csv
touch wet-meadow-connectivity/manuscript/draft.docx
touch wet-meadow-connectivity/docs/references.bib
