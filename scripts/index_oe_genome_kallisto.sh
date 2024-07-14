#!/bin/bash

index_path='/home/jdubos2/mtstp/data/oe_genome/kallisto_index.idx'
genes_path='/home/jdubos2/mtstp/data/oe_genome/cds.fasta'

kallisto index -i $index_path $genes_path