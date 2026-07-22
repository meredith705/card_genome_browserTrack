#!/usr/bin/env bash
# ============================================================================
# build_sv_hub.sh
#
# Extracts structural variants (with allele frequency) from a long-read SV
# VCF and builds a complete UCSC Genome Browser track hub (bigBed, colored
# by SV type) ready to push to GitHub and connect in UCSC.
#
# Requirements (install via conda/mamba/bioconda):
#   bcftools, bedtools (optional), UCSC tools: bedToBigBed, fetchChromSizes
#
#   conda install -c bioconda bcftools ucsc-bedtobigbed ucsc-fetchchromsizes
#
# Usage:
#   ./build_sv_hub.sh your_svs.vcf.gz my_sv_hub/
# ============================================================================
set -euo pipefail

VCF="$1"                  # input SV VCF (bgzipped + tabix-indexed, or plain)
OUTDIR="${2:-sv_hub}"     # output directory for the hub
GENOME="hg38"
HUB_SHORT_LABEL="SVs from the NIH CARD Long-read SV catalogue"
HUB_LONG_LABEL="NIH Center for Alzheimer's and Related Dementias Long-read Structural Variants (205 European Ancestry Samples 146 African American Ancestry Samples)"
TRACK_NAME="NIH_CARD_longReadSVs"

mkdir -p "$OUTDIR/$GENOME"
cd "$OUTDIR"

echo "[1/7] Checking VCF INFO fields for SVTYPE / END / SVLEN / AF ..."
bcftools view -h "../$VCF" | grep "^##INFO" || true
echo "  -> confirm SVTYPE, END, SVLEN, AF (or AC/AN) exist above."
echo "     If AF is missing, add: bcftools +fill-tags in.vcf.gz -- -t AF"

echo "[2/7] Extracting SV records to BED9+6 (SVTYPE, SVLEN, AF, carrier counts by cohort) ..."
# Notes:
#  - For INS, END == POS+1 (a single-base insertion point) since most
#    long-read callers (Sniffles2/cuteSV/pbsv) report END==POS for insertions.
#  - Score column (col5) left at 0; itemRgb (col9) encodes SVTYPE color.
#  - thickStart/thickEnd (cols7-8) set equal to chromStart/chromEnd (required
#    fields for bigBed even if unused for rendering).
#
#  Carrier counting: rather than parsing per-sample GT strings in awk (slow --
#  effectively an interpreted loop over every sample for every SV), we let
#  bcftools' fill-tags plugin do it in C. Giving it a sample->population
#  mapping produces per-population AC_Het_<POP>/AC_Hom_<POP> tags in a single
#  pass over the VCF; carrierCount for a group = AC_Het + AC_Hom for that
#  group (each het genotype = 1 carrier, each hom-alt genotype = 1 carrier).

echo "  building sample -> cohort population map ..."
bcftools query -l "../$VCF" | awk '{
  if ($0 ~ /^NABEC/)      print $0"\tNABEC";
  else if ($0 ~ /^HBCC/)  print $0"\tHBCC";
  else                    print $0 > "/dev/stderr";
}' > sample_pop.txt 2> unclassified_samples.txt

if [ -s unclassified_samples.txt ]; then
  echo "  WARNING: $(wc -l < unclassified_samples.txt) sample(s) matched neither NABEC nor HBCC and were excluded from carrier counting:"
  cat unclassified_samples.txt | sed 's/^/    /'
fi
rm -f unclassified_samples.txt

echo "  computing per-cohort carrier counts with bcftools +fill-tags (single pass) ..."
bcftools +fill-tags "../$VCF" -Ob -o tagged.bcf -- -S sample_pop.txt -t AC_Het,AC_Hom
bcftools index tagged.bcf

bcftools query -f \
  '%CHROM\t%POS0\t%END\t%ID\t0\t.\t%POS0\t%END\t0,0,0\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/AF\t%INFO/AC_Het_NABEC\t%INFO/AC_Hom_NABEC\t%INFO/AC_Het_HBCC\t%INFO/AC_Hom_HBCC\n' \
  tagged.bcf > sv_raw.tsv

# rm -f tagged.bcf tagged.bcf.csi sample_pop.txt

# Fix zero-length INS intervals (END<=START) to be 1bp wide, assign color,
# and sum the pre-computed per-population het+hom counts.
awk -F'\t' 'BEGIN{OFS="\t"}
function num(x) { return (x == "." || x == "") ? 0 : x+0 }
{
  start=$2; end=$3;
  if (end <= start) end = start + 1;   # ensure viewable interval for INS/BND
  svtype=$10;
  color="0,0,0";
  if (svtype ~ /DEL/) color="217,95,2";     # red-orange
  else if (svtype ~ /DUP/) color="27,158,119"; # green
  else if (svtype ~ /INV/) color="117,112,179"; # purple
  else if (svtype ~ /INS/) color="230,171,2";  # orange
  else if (svtype ~ /BND/) color="102,102,102"; # grey
  svlen = $11; gsub("-", "", svlen);   # SVLEN is often negative for DELs; strip the sign
  name = svtype"-"svlen;

  nabec = num($13) + num($14);
  hbcc  = num($15) + num($16);
  other = num($17) + num($18);
  total = nabec + hbcc + other;

  print $1, start, end, name, 0, ".", start, end, color, svtype, $11, $12, total, nabec, hbcc;
}' sv_raw.tsv | sort -k1,1 -k2,2n > "$GENOME/${TRACK_NAME}.bed"

echo "[3/7] Fetching chrom.sizes for $GENOME ..."
fetchChromSizes "$GENOME" > "$GENOME/${GENOME}.chrom.sizes"

echo "[4/7] Building autoSql (.as) schema for extra fields ..."
cat > "$GENOME/${TRACK_NAME}.as" << 'EOF'
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
EOF

echo "[5/7] Converting BED -> bigBed ..."
bedToBigBed -as="$GENOME/${TRACK_NAME}.as" -type=bed9+6 \
  "$GENOME/${TRACK_NAME}.bed" "$GENOME/${GENOME}.chrom.sizes" \
  "$GENOME/${TRACK_NAME}.bb"

rm -f sv_raw.tsv "$GENOME/${TRACK_NAME}.bed"

echo "[6/7] Writing track description HTML page ..."
# UCSC track hubs render a description page (like hgTrackUi) from an .html
# file referenced via the "html" line in trackDb.txt. UCSC supplies the page
# chrome/navigation automatically -- this file only needs the body content,
# following UCSC's standard section conventions.
cat > "$GENOME/${TRACK_NAME}.html" << EOF
<h2>Description</h2>
<p>
This track shows structural variants (SVs) identified by Oxford Nanopore long-read sequencing of post-mortem brain 
tissue (prefrontal cortex) from 351 individuals, generated by the NIH Center for Alzheimer's and Related Dementias 
(NIH CARD) Long-Read Initiative. Structural variants are genomic rearrangements larger than about 50 bp, such as 
deletions, insertions, inversions and duplications; because they alter or move large stretches of DNA at once they 
can have outsized effects on gene dosage, gene regulation and DNA methylation compared with single-nucleotide changes.
</p>
<p>
The cohort combines two studies: 205 samples of European ancestry from the North American Brain Expression Consortium 
(NABEC, dbGaP phs001300) and 146 samples of African and African-admixed ancestry from the NIMH Human Brain Collection 
Core (HBCC, dbGaP phs000979). The track contains more than 228,000 SVs called against GRCh38: about 127,000 insertions, 
102,000 deletions, 431 inversions and one tandem duplication. Each record carries the number of carrier samples overall 
and split by cohort (NABEC and HBCC), together with the allele frequency reported by the source project.

</p>

<h2>Display Conventions and Configuration</h2>
<p>
Items are colored by SV type:
</p>
<ul>
<li><b><span style="color:rgb(217,95,2)">DEL</span></b> &mdash; deletion</li>
<li><b><span style="color:rgb(27,158,119)">DUP</span></b> &mdash; duplication</li>
<li><b><span style="color:rgb(117,112,179)">INV</span></b> &mdash; inversion</li>
<li><b><span style="color:rgb(230,171,2)">INS</span></b> &mdash; insertion (shown as a 1 bp feature at the insertion point)</li>
<li><b><span style="color:rgb(102,102,102)">BND</span></b> &mdash; breakend / translocation</li>
</ul>
<p>
The item name is <code>SVTYPE-SVLEN</code> (e.g. <code>DEL-1200</code>).
Clicking an item shows its details, including:
</p>
<ul>
<li><b>svType</b> / <b>svLen</b> &mdash; structural variant type and length, from the source VCF INFO fields</li>
<li><b>alleleFreq</b> &mdash; allele frequency, from the source VCF INFO/AF field</li>
<li><b>carrierCount</b> &mdash; number of samples with a non-reference, non-missing genotype for this variant</li>
<li><b>nabecCount</b> / <b>hbccCount</b> &mdash; carrier counts restricted to samples whose ID starts with <code>NABEC</code> or <code>HBCC</code>, respectively</li>
</ul>

<h2>Methods</h2>
<p>
NABEC samples were sequenced with ONT R9.4.1 and HBCC samples were sequnced with R10.4.1 PromethION flow cells with median N50 of 27 Kbp and at an average depth of ∼40x genome coverage. 
SVs were called from minimap2 read alignments with Sniffles2 v.2.3 as well as with assembly alignments using Hapdiff. Assemblies were
produced by Shasta v0.11.1 and phased using HapDup v.0.12. Assembly SVs were then merged across samples with Truvari while read SVs were merged 
aross samples with Sniffles2, finally Truvari was again used to merge the read and assembly SVs together as well as across the 
2 cohorts. All of these data were produced on AnVIL Terra using Nanopore Analysis Pipeline (NAPU) workflows which can be found <a href="https://dockstore.org/organizations/NIHCARD" target="_blank">here</a>, 
also see this <a href="https://www.nature.com/articles/s41592-023-01993-x" target="_blank">publication</a>  for more information about the NAPU workflow.
SVs were then converted to bigBed for display with a custom pipeline
(<a href="https://github.com/meredith705/card_genome_browserTrack" target="_blank">source scripts</a>).
</p>

<h2>Data Access</h2>
<p>
The bigBed file underlying this track can be downloaded directly from the
hub's GitHub repository, or queried by region with the UCSC
<a href="https://genome.ucsc.edu/goldenPath/help/api.html" target="_blank">REST API</a>
or command-line <code>bigBedToBed</code> utility. The calls in VCF format is available on 
AnVIL for controlled access through dgGaP: NABEC:phs001300.v6.p1 (sub study phsID: phs003181.v2.p1) and HBCC:phs000979.v4.p2.
</p>

<h2>Credits</h2>
<p>
Thanks to the North American Brain Expression Consortium (NABEC), the NIMH Human Brain Collection Core (HBCC), 
the Banner Sun Health Research Institute Brain and Body Donation Program, and the NIH CARD Long-Read Initiative for 
generating and sharing this dataset, and to Melissa Meredith for preparing the browser track. This work was supported 
by the Intramural Research Programs of the NIA, NINDS, NCI, NHGRI and NIMH, and used the NIH STRIDES Initiative and 
the NIH HPC Biowulf cluster.
</p>

<h2>References</h2>
<p>
K. J. Billingsley et al. <a href="https://www.biorxiv.org/content/10.1101/2024.12.16.628723v1.full" target="_blank">Long-read sequencing of hundreds of diverse brains provides insight into the impact of structural variation on 
gene expression and DNA methylation</a>, bioRxiv (2024)p. 2024.12.16.628723.
<p>
Kolmogorov, M. et al.  <a href="https://www.nature.com/articles/s41592-023-01993-x" target="_blank">Scalable Nanopore sequencing of human genomes provides a comprehensive view of haplotype-resolved variation and 
methylation</a>. Nat. Methods 20, 1483–1492 (2023).


<p>
EOF

echo "[7/7] Writing hub.txt, genomes.txt, trackDb.txt ..."

cat > hub.txt << EOF
hub svHub
shortLabel ${HUB_SHORT_LABEL}
longLabel ${HUB_LONG_LABEL}
genomesFile genomes.txt
email melissa@datatecnica.com
EOF

cat > genomes.txt << EOF
genome ${GENOME}
trackDb ${GENOME}/trackDb.txt
EOF

cat > "$GENOME/trackDb.txt" << EOF
track ${TRACK_NAME}
bigDataUrl ${TRACK_NAME}.bb
shortLabel SVs from the NIH CARD Long-read SV catalogue
longLabel NIH Center for Alzheimer's and Related Dementias Long-read Structural Variants (205 European Ancestry Samples 146 African American Ancestry Samples)
type bigBed 9 +
itemRgb on
visibility pack
maxItems 100000
html ${TRACK_NAME}
EOF

echo ""
echo "Done. Hub built in: $(pwd)"
echo ""
echo "Next steps:"
echo "  1. git init, add these files, push to a public GitHub repo"
echo "  2. Your hub URL will be:"
echo "     https://raw.githubusercontent.com/<user>/<repo>/main/hub.txt"
echo "  3. In UCSC: My Data > Track Hubs > My Hubs > paste that URL > Add Hub"
echo "  4. In UCSC: My Data > Sessions > Save Session (check 'include hub tracks')"
echo "     to get your permanent session link"