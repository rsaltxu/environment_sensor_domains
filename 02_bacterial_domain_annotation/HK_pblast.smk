"""
Snakefile to map processed ORFs from PRODIGAL metagenome annotations to environmental sensor domain DIAMOND database, per sample
"""
with open('/data/sample_list.txt') as sample_list_file:
    for line in sample_list_file:
        samples.append(line[:-1])

rule all:
    input:
        expand('/data/diamond_pblast/normal/blastp_sample_{sample}.tsv',sample=samples)


rule db_blastp_normal:
    input:
        proteome='/australia/Annotation/Prodigal/MEGAHIT/sample_{sample}/MEGAHIT-sample_{sample}_prodigal.faa.gz'
    output:
        output_file='/data/diamond_pblast/normal/blastp_sample_{sample}.tsv'
    threads:
        18
    shell:
        """diamond blastp -f 6 --very-sensitive -v -d /home/roman.sanz/data/diamond_dbs/HK_db_diamond.dmnd \
        -q {input.proteome} \
        -o {output.output_file} \
        --quiet -k 1 -p {threads}
        """
