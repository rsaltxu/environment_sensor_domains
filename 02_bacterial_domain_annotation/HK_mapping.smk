"""
Snakefile to map processed shotgun metagenomic reads to ORFs from PRODIGAL metagenome annotations, per sample
"""

samples=[]
with open('/data/sample_list.txt') as sample_list_file:
    for line in sample_list_file:
        samples.append(line[:-1])

rule all:
    input:
        expand('/data/diamond_alignments/{sample}_read_align.tsv',sample=samples)

rule concat_reads:
    output:
        temp('/data/concat_reads/{sample}_collapsed.fq.gz')
    shell:
        """
        cat /data/reads/{wildcards.sample}_*_trimmed_R1_paired.fq.gz \
        /data/reads/{wildcards.sample}_*_trimmed_R2_paired.fq.gz \
         > /data/concat_reads/{wildcards.sample}_collapsed.fq.gz
        """

rule make_diamond_db:
    input:
        '/australia/Annotation/Prodigal/MEGAHIT/sample_{sample}/MEGAHIT-sample_{sample}_prodigal.faa.gz'
    output:
        temp('/data/diamond_dbs/proteome_{sample}.dmnd')
    shell:
        """diamond makedb --in {input} -d {output}"""

rule perform_alignment:
    input:
        collapsed_reads='/data/concat_reads/{sample}_collapsed.fq.gz',
        database='/data/diamond_dbs/proteome_{sample}.dmnd'
    output:
        '/data/diamond_alignments/{sample}_read_align.tsv'
    threads:
        18
    shell:
        """diamond blastx -f 6 --very-sensitive \
        -d {input.database} \
        -q {input.collapsed_reads} \
        -o {output} \
        --quiet -k 1 --id 80 -p {threads}
        """
