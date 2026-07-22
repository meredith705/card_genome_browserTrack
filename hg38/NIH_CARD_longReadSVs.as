table svBed
"Structural variants from long-read sequencing"
(
string chrom;      "Chromosome"
uint   chromStart; "Start position"
uint   chromEnd;   "End position"
string name;       "SV ID"
uint   score;      "Score (unused)"
char[1] strand;    "Strand"
uint   thickStart; "Thick start"
uint   thickEnd;   "Thick end"
uint   reserved;   "RGB color"
string svType;     "Structural variant type (DEL/DUP/INS/INV/BND)"
string svLen;      "Structural variant length"
string alleleFreq; "Allele frequency"
uint   alleleCount; "Allele Count for this variant"
uint   nabecAlleleCount;   "NABEC allele count"
uint   hbccAlleleCount;    "HBCC allele count"
)
