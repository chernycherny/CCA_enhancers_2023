


#########################################################################################################
########## Using GRCh37.p13 genome, and refseq curated annotations (Apr 2018) ###########################
#########################################################################################################
cd ~/CCA/genome
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz

# Generate chromosome sizes file
samtools faidx GRCh37.p13.genome.fa
cut -f1,2 GRCh37.p13.genome.fa.fai > GRCh37.p13.genome.chrom.sizes



## refseq_curated annotations downloaded from UCSC table browser in Apr 2018: track NCBI refseq, table refseq curated
## Convert to gtf: see http://genomewiki.ucsc.edu/index.php/Genes_in_gtf_or_gff_format, https://bioinformatics.stackexchange.com/questions/2548/hg38-gtf-file-with-refseq-annotations
# remove first column and first line
cut -f 2- genes_refseq_curated.txt > genes_refseq_curated.removeCol1.txt.tmp
tail -n +2 genes_refseq_curated.removeCol1.txt.tmp > genes_refseq_curated.removeCol1.txt
rm genes_refseq_curated.removeCol1.txt.tmp
# convert genepred to gtf format
~/UCSC_tools/genePredToGtf file genes_refseq_curated.removeCol1.txt genes_refseq_curated.gtf
rm genes_refseq_curated.removeCol1.txt

# keep only chr1-22, x, y, m
egrep "^chr[0-9XYM]+\s" genes_refseq_curated.gtf > genes_refseq_curated.filt.gtf
mv genes_refseq_curated.filt.gtf genes_refseq_curated.gtf

# generate list of gene lengths
~/gtftools/GTFtools_0.6.5/gtftools.py -l genes_refseq_curated.gene_lengths.bed genes_refseq_curated.gtf
~/gtftools/GTFtools_0.6.5/gtftools.py -r genes_refseq_curated.transcript_lengths.bed genes_refseq_curated.gtf




# Generate transcript_id to gene_id, gene_name mapping file
echo -e "transcript_id\tgene_id\tgene_name" > genes_refseq_curated.transcript_to_gene_ids.mapping
awk '$3 == "transcript" {print $0} ' < genes_refseq_curated.gtf | \
        cut -f 9 | \
        awk '{if ($1!="gene_id" || $3!="transcript_id" || $5!="gene_name") {print "format problem", $1, $3, $5} else {print $4, $2, $6}}' | \
        tr -d '";' | \
        tr ' ' '\t' \
>> genes_refseq_curated.transcript_to_gene_ids.mapping


# build STAR reference
cd ~/CCA/genome
mkdir STAR_GRCh37.p13_refseq
qsub -pe smp 8 ~/STAR/STAR-master/bin/Linux_x86_64/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ~/CCA/genome/STAR_GRCh37.p13_refseq --genomeFastaFiles ~/CCA/genome/GRCh37.p13.genome.fa --sjdbGTFfile ~/CCA/genome/genes_refseq_curated.gtf --sjdbOverhang 149


# build RSEM reference
mkdir RSEM_GRCh37.p13_refseq
qsub rsem-prepare-reference --gtf ~/CCA/genome/genes_refseq_curated.gtf ~/CCA/genome/GRCh37.p13.genome.fa ~/CCA/genome/RSEM_GRCh37.p13_refseq/GRCh37.p13_refseq



# build REF_FLAT for picardmetrics
cd ~/CCA/genome
picardmetrics refFlat genes_refseq_curated.gtf




####### Build blacklist meta data. rRNA, snoRNA, MT, Metazoa_SRP, RN7 ##########
# First, generate blacklisted transcripts using gencode gtf. Then convert that to refseq IDs
# Intervals for rRNA transcripts.
grep 'gene_type "rRNA"' gencode.v19.annotation.gtf | \
        awk '$3 == "transcript"' | \
        cut -f9 | \
        awk '{print $4, $2, $10}' | \
        tr -d '";' | \
        tr ' ' '\t' \
>> GRCh37.p13.blacklist.meta

# Intervals for snoRNA transcripts.
grep 'gene_type "snoRNA"' gencode.v19.annotation.gtf | \
        awk '$3 == "transcript"' | \
        cut -f9 | \
        awk '{print $4, $2, $10}' | \
        tr -d '";' | \
        tr ' ' '\t' \
>> GRCh37.p13.blacklist.meta

# Intervals for MTgenes transcripts.
grep 'chrM' gencode.v19.annotation.gtf | \
        awk '$3 == "transcript"' | \
        cut -f9 | \
        awk '{print $4, $2, $10}' | \
        tr -d '";' | \
        tr ' ' '\t' \
>> GRCh37.p13.blacklist.meta

# Intervals for Metazoa_SRP transcripts.
grep 'gene_name "Metazoa_SRP"' gencode.v19.annotation.gtf | \
        awk '$3 == "transcript"' | \
        cut -f9 | \
        awk '{print $4, $2, $10}' | \
        tr -d '";' | \
        tr ' ' '\t' \
>> GRCh37.p13.blacklist.meta

# Intervals for RN7 class transcripts.
grep 'gene_name "RN7' gencode.v19.annotation.gtf | \
        awk '$3 == "transcript"' | \
        cut -f9 | \
        awk '{print $4, $2, $10}' | \
        tr -d '";' | \
        tr ' ' '\t' \
>> GRCh37.p13.blacklist.meta

# convert to refseq IDs
cut -f 1 GRCh37.p13.blacklist.meta | cut -d . -f 1 > tmp_GRCh37.p13.blacklist.meta.transcriptIDs
# use https://biodbnet-abcc.ncifcrf.gov/db/db2db.php to convert ensembl transcript IDs to refseq mRNA/noncoding-RNA IDs
cat tmp_GRCh37.p13.blacklist.meta.transcriptIDs.mapped | tail -n +2 | cut -f 2- | tr ';' ' ' | tr '-' ' ' | awk '{for(i=1; i<=NF; i++){print $i}}' > tmp_GRCh37.p13.blacklist.meta.transcriptIDs.mapped.format
# run in R
blacklist_transcripts = unique(read.table("tmp_GRCh37.p13.blacklist.meta.transcriptIDs.mapped.format", stringsAsFactors=F, header=F)[,1]) # 403
transcript_gene_map = read.table("genes_refseq_curated.transcript_to_gene_ids.mapping", stringsAsFactors=F, header=T) # 58384
transcript_gene_map$transcript_id2 = sapply(strsplit(transcript_gene_map$transcript_id, split=".", fixed=T), function(x){x[1]})
blacklist_refseq = transcript_gene_map[transcript_gene_map$transcript_id2 %in% blacklist_transcripts,c(1:3)] # 360
# Warning: using this mapping to create blacklist includes a few wrong genes (eg PKIB), and misses many genes (RN7, many SNORA, SNORD, etc)
# Thus, manually curate my own list, based on the mapping results!
# Manual curation: RN7S, RNA5S, SCARNA, SNORA, SNORD, RNU
blacklist_refseq = transcript_gene_map[grep("^RN7S", transcript_gene_map$gene_name), c(1:3)]
blacklist_refseq = rbind(blacklist_refseq, transcript_gene_map[grep("^RNA5S", transcript_gene_map$gene_name), c(1:3)])
blacklist_refseq = rbind(blacklist_refseq, transcript_gene_map[grep("^SCARNA\\d", transcript_gene_map$gene_name), c(1:3)])
blacklist_refseq = rbind(blacklist_refseq, transcript_gene_map[grep("^SNORA\\d", transcript_gene_map$gene_name), c(1:3)])
blacklist_refseq = rbind(blacklist_refseq, transcript_gene_map[grep("^SNORD\\d", transcript_gene_map$gene_name), c(1:3)])
blacklist_refseq = rbind(blacklist_refseq, transcript_gene_map[grep("^RNU\\d", transcript_gene_map$gene_name), c(1:3)]) # 501
write.table(blacklist_refseq, file="genes_refseq_curated.blacklist.meta", row.names=F, col.names=F, sep="\t", quote=F)
