#!/bin/bash

blastp -query ../31899.10/31899.10.PATRIC.faa -db Cbes_ncbi -out Cbes_blastp_out.tsv -evalue 0.0000001 -outfmt 6 
