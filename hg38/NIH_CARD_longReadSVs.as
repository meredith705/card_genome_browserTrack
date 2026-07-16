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
uint   carrierCount; "Number of samples carrying this variant"
uint   nabecCount;   "Number of NABEC carrier samples"
uint   hbccCount;    "Number of HBCC carrier samples"
)
