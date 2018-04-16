#global settings
shell.prefix("set -o pipefail; ")
shell.prefix("set -e; ")
shell.prefix("set -u; ")
localrules: all
configfile: "config.json"

rule all:
        input:
                expand("03.ordered_mapping_bowtie2/{library}/00.hg38/{library}.hg38.sam", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/01.miRNA/{library}.miRNA.sam", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/02.piRNA/{library}.piRNA.sam", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/03.Y_RNA/{library}.Y_RNA.sam", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/04.snRNA/{library}.snRNA.sam", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/05.srpRNA/{library}.srpRNA.sam", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/06.tRNA/{library}.tRNA.sam", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/07.other_lncRNA/{library}.other_lncRNA.sam",library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/08.mRNA/{library}.mRNA.sam",library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/01.miRNA/{library}.miRNA.sam", library=config["Library"]),
                expand("05.counts_HTSeq_bowtie2/{library}/01.miRNA/{library}.miRNA.counts", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/02.piRNA/{library}.piRNA.sam", library=config["Library"]),
                expand("05.counts_HTSeq_bowtie2/{library}/02.piRNA/{library}.piRNA.counts", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/03.Y_RNA/{library}.Y_RNA.sam", library=config["Library"]),
                expand("05.counts_HTSeq_bowtie2/{library}/03.Y_RNA/{library}.Y_RNA.counts", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/04.snRNA/{library}.snRNA.sam", library=config["Library"]),
                expand("05.counts_HTSeq_bowtie2/{library}/04.snRNA/{library}.snRNA.counts", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/05.srpRNA/{library}.srpRNA.sam", library=config["Library"]),
                expand("05.counts_HTSeq_bowtie2/{library}/05.srpRNA/{library}.srpRNA.counts", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/06.tRNA/{library}.tRNA.sam", library=config["Library"]),
                expand("05.counts_HTSeq_bowtie2/{library}/06.tRNA/{library}.tRNA.counts", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/07.other_lncRNA/{library}.other_lncRNA.sam", library=config["Library"]),
                expand("05.counts_HTSeq_bowtie2/{library}/07.other_lncRNA/{library}.other_lncRNA_withStrand.counts", library=config["Library"]),
                expand("05.counts_HTSeq_bowtie2/{library}/07.other_lncRNA/{library}.other_lncRNA_noStrand.counts", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/08.mRNA/{library}.mRNA.sam", library=config["Library"]),
                expand("05.counts_HTSeq_bowtie2/{library}/08.mRNA/{library}.mRNA.counts", library=config["Library"])


rule fetch_mapped_miRNA_sam:
        input:
                miRNA_sam = "03.ordered_mapping_bowtie2/{library}/01.miRNA/{library}.miRNA.sam"
        output:
                miRNA_fastq = "03.ordered_mapping_bowtie2/{library}/01.miRNA/{library}.miRNA.mapped.fastq",
                miRNA_2pass_sam = "03.ordered_mapping_bowtie2/{library}/01.miRNA/{library}.miRNA_2pass_mapped.sam",
                miRNA_counts = "05.counts_HTSeq_bowtie2/{library}/01.miRNA/{library}.miRNA.counts"
        params:
                hg38_index = config["References"]["bowtie2_hg38_genome_index_dir"],
                miRNA_GTF = config["References"]["miRNA_GTF"],
                cpu = config["Remove_rRNA"]["cpu"],
                jobname = "{library}.fetch_mapped_miRNA_sam"
        shell:
                """
                samtools view -bhF 4 {input.miRNA_sam} | bedtools bamtofastq -i - -fq {output.miRNA_fastq}
                bowtie2 -p {params.cpu} --sensitive-local -x {params.hg38_index} {output.miRNA_fastq} -S {output.miRNA_2pass_sam}
                htseq-count -m intersection-strict --idattr=Name --type=miRNA_primary_transcript {output.miRNA_2pass_sam} {params.miRNA_GTF} > {output.miRNA_counts}
                """

rule fetch_mapped_piRNA_sam:
        input:
                piRNA_sam = "03.ordered_mapping_bowtie2/{library}/02.piRNA/{library}.piRNA.sam"
        output:
                piRNA_fastq = "03.ordered_mapping_bowtie2/{library}/02.piRNA/{library}.piRNA.mapped.fastq",
                piRNA_2pass_sam = "03.ordered_mapping_bowtie2/{library}/02.piRNA/{library}.piRNA_2pass_mapped.sam",
                piRNA_counts = "05.counts_HTSeq_bowtie2/{library}/02.piRNA/{library}.piRNA.counts"
        params:
                hg38_index = config["References"]["bowtie2_hg38_genome_index_dir"],
                piRNA_GTF = config["References"]["piRNA_GTF"],
                cpu = config["Remove_rRNA"]["cpu"],
                jobname = "{library}.fetch_mapped_piRNA_sam"
        shell:
                """
                samtools view -bhF 4 {input.piRNA_sam} | bedtools bamtofastq -i - -fq {output.piRNA_fastq}
                bowtie2 -p {params.cpu} --sensitive-local -x {params.hg38_index} {output.piRNA_fastq} -S {output.piRNA_2pass_sam}
                htseq-count -m intersection-strict --idattr=transcript_id --type=exon {output.piRNA_2pass_sam} {params.piRNA_GTF} > {output.piRNA_counts}
                """

rule fetch_mapped_Y_RNA_sam:
        input:
                Y_RNA_sam = "03.ordered_mapping_bowtie2/{library}/03.Y_RNA/{library}.Y_RNA.sam"
        output:
                Y_RNA_fastq = "03.ordered_mapping_bowtie2/{library}/03.Y_RNA/{library}.Y_RNA.mapped.fastq",
                Y_RNA_2pass_sam = "03.ordered_mapping_bowtie2/{library}/03.Y_RNA/{library}.Y_RNA_2pass_mapped.sam",
                Y_RNA_counts = "05.counts_HTSeq_bowtie2/{library}/03.Y_RNA/{library}.Y_RNA.counts"
        params:
                hg38_index = config["References"]["bowtie2_hg38_genome_index_dir"],
                Y_RNA_GTF = config["References"]["Y_RNA_GTF"],
                cpu = config["Remove_rRNA"]["cpu"],
                jobname = "{library}.fetch_mapped_Y_RNA_sam"
        shell:
                """
                samtools view -bhF 4 {input.Y_RNA_sam} | bedtools bamtofastq -i - -fq {output.Y_RNA_fastq}
                bowtie2 -p {params.cpu} --sensitive-local -x {params.hg38_index} {output.Y_RNA_fastq} -S {output.Y_RNA_2pass_sam}
                htseq-count -m intersection-strict --idattr=transcript_name --type=exon {output.Y_RNA_2pass_sam} {params.Y_RNA_GTF} > {output.Y_RNA_counts}
                """

rule fetch_mapped_snRNA_sam:
        input:
                snRNA_sam = "03.ordered_mapping_bowtie2/{library}/04.snRNA/{library}.snRNA.sam"
        output:
                snRNA_fastq = "03.ordered_mapping_bowtie2/{library}/04.snRNA/{library}.snRNA.mapped.fastq",
                snRNA_2pass_sam = "03.ordered_mapping_bowtie2/{library}/04.snRNA/{library}.snRNA_2pass_mapped.sam",
                snRNA_counts = "05.counts_HTSeq_bowtie2/{library}/04.snRNA/{library}.snRNA.counts"
        params:
                hg38_index = config["References"]["bowtie2_hg38_genome_index_dir"],
                snRNA_GTF = config["References"]["snRNA_GTF"],
                cpu = config["Remove_rRNA"]["cpu"],
                jobname = "{library}.fetch_mapped_snRNA_sam"
        shell:
                """
                samtools view -bhF 4 {input.snRNA_sam} | bedtools bamtofastq -i - -fq {output.snRNA_fastq}
                bowtie2 -p {params.cpu} --sensitive-local -x {params.hg38_index} {output.snRNA_fastq} -S {output.snRNA_2pass_sam}
                htseq-count -m intersection-strict --idattr=transcript_name --type=snRNA {output.snRNA_2pass_sam} {params.snRNA_GTF} > {output.snRNA_counts}
                """

rule fetch_mapped_srpRNA_sam:
        input:
                srpRNA_sam = "03.ordered_mapping_bowtie2/{library}/05.srpRNA/{library}.srpRNA.sam"
        output:
                srpRNA_fastq = "03.ordered_mapping_bowtie2/{library}/05.srpRNA/{library}.srpRNA.mapped.fastq",
                srpRNA_2pass_sam = "03.ordered_mapping_bowtie2/{library}/05.srpRNA/{library}.srpRNA_2pass_mapped.sam",
                srpRNA_counts = "05.counts_HTSeq_bowtie2/{library}/05.srpRNA/{library}.srpRNA.counts"
        params:
                hg38_index = config["References"]["bowtie2_hg38_genome_index_dir"],
                srpRNA_GTF = config["References"]["srpRNA_GTF"],
                cpu = config["Remove_rRNA"]["cpu"],
                jobname = "{library}.fetch_mapped_srpRNA_sam"
        shell:
                """
                samtools view -bhF 4 {input.srpRNA_sam} | bedtools bamtofastq -i - -fq {output.srpRNA_fastq}
                bowtie2 -p {params.cpu} --sensitive-local -x {params.hg38_index} {output.srpRNA_fastq} -S {output.srpRNA_2pass_sam}
                htseq-count -m intersection-strict --idattr=transcript_name --type=exon {output.srpRNA_2pass_sam} {params.srpRNA_GTF} > {output.srpRNA_counts}
                """

rule fetch_mapped_tRNA_sam:
        input:
                tRNA_sam = "03.ordered_mapping_bowtie2/{library}/06.tRNA/{library}.tRNA.sam"
        output:
                tRNA_fastq = "03.ordered_mapping_bowtie2/{library}/06.tRNA/{library}.tRNA.mapped.fastq",
                tRNA_2pass_sam = "03.ordered_mapping_bowtie2/{library}/06.tRNA/{library}.tRNA_2pass_mapped.sam",
                tRNA_counts = "05.counts_HTSeq_bowtie2/{library}/06.tRNA/{library}.tRNA.counts"
        params:
                hg38_index = config["References"]["bowtie2_hg38_genome_index_dir"],
                tRNA_GTF = config["References"]["tRNA_GTF"],
                cpu = config["Remove_rRNA"]["cpu"],
                jobname = "{library}.fetch_mapped_tRNA_sam"
        shell:
                """
                samtools view -bhF 4 {input.tRNA_sam} | bedtools bamtofastq -i - -fq {output.tRNA_fastq}
                bowtie2 -p {params.cpu} --sensitive-local -x {params.hg38_index} {output.tRNA_fastq} -S {output.tRNA_2pass_sam}
                htseq-count -m intersection-strict --idattr=transcript_name --type=tRNA {output.tRNA_2pass_sam} {params.tRNA_GTF} > {output.tRNA_counts}
                """

rule fetch_mapped_other_lncRNA_sam:
        input:
                other_lncRNA_sam = "03.ordered_mapping_bowtie2/{library}/07.other_lncRNA/{library}.other_lncRNA.sam"
        output:
                other_lncRNA_fastq = "03.ordered_mapping_bowtie2/{library}/07.other_lncRNA/{library}.other_lncRNA.mapped.fastq",
                other_lncRNA_2pass_sam = "03.ordered_mapping_bowtie2/{library}/07.other_lncRNA/{library}.other_lncRNA_2pass_mapped.sam",
                other_lncRNA_withStrand_counts = "05.counts_HTSeq_bowtie2/{library}/07.other_lncRNA/{library}.other_lncRNA_withStrand.counts",
                other_lncRNA_noStrand_counts = "05.counts_HTSeq_bowtie2/{library}/07.other_lncRNA/{library}.other_lncRNA_noStrand.counts"
        params:
                hg38_index = config["References"]["bowtie2_hg38_genome_index_dir"],
                other_lncRNA_All_withStrand_GTF = config["References"]["other_lncRNA_All_withStrand_GTF"],
                other_lncRNA_All_noStrand_GTF = config["References"]["other_lncRNA_All_noStrand_GTF"],
                cpu = config["Remove_rRNA"]["cpu"],
                jobname = "{library}.fetch_mapped_other_lncRNA_sam"
        shell:
                """
                samtools view -bhF 4 {input.other_lncRNA_sam} | bedtools bamtofastq -i - -fq {output.other_lncRNA_fastq}
                bowtie2 -p {params.cpu} --sensitive-local -x {params.hg38_index} {output.other_lncRNA_fastq} -S {output.other_lncRNA_2pass_sam}
                htseq-count -m intersection-strict --idattr=transcript_id --type=exon {output.other_lncRNA_2pass_sam} {params.other_lncRNA_All_withStrand_GTF} > {output.other_lncRNA_withStrand_counts}
                htseq-count -s no -m intersection-strict --idattr=transcript_id --type=exon {output.other_lncRNA_2pass_sam} {params.other_lncRNA_All_noStrand_GTF} > {output.other_lncRNA_noStrand_counts}
                """

rule fetch_mapped_mRNA_sam:
        input:
                mRNA_sam = "03.ordered_mapping_bowtie2/{library}/08.mRNA/{library}.mRNA.sam"
        output:
                mRNA_fastq = "03.ordered_mapping_bowtie2/{library}/08.mRNA/{library}.mRNA.mapped.fastq",
                mRNA_2pass_sam = "03.ordered_mapping_bowtie2/{library}/08.mRNA/{library}.mRNA_2pass_mapped.sam",
                mRNA_counts = "05.counts_HTSeq_bowtie2/{library}/08.mRNA/{library}.mRNA.counts"
        params:
                hg38_index = config["References"]["bowtie2_hg38_genome_index_dir"],
                mRNA_GTF = config["References"]["mRNA_GTF"],
                cpu = config["Remove_rRNA"]["cpu"],
                jobname = "{library}.fetch_mapped_mRNA_sam"
        shell:
                """
                samtools view -bhF 4 {input.mRNA_sam} | bedtools bamtofastq -i - -fq {output.mRNA_fastq}
                bowtie2 -p {params.cpu} --sensitive-local -x {params.hg38_index} {output.mRNA_fastq} -S {output.mRNA_2pass_sam}
                htseq-count -m intersection-strict --idattr=transcript_name --type=exon {output.mRNA_2pass_sam} {params.mRNA_GTF} > {output.mRNA_counts}
                """
