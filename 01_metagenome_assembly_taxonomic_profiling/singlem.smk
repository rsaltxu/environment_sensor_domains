samples=[]
with open('/data/sample_list.txt') as sample_list_file:
    for line in sample_list_file:
        samples.append(line[:-1])

rule all:
    input:
        expand('/data/singlem/australia/{sample}_profile_with_extras.tsv',sample=samples)

rule singlem:
    output:
        output='/data/singlem/australia/{sample}_profile.tsv'
    threads:
        18
    conda:
        "/miniforge3/envs/singlem"
    shell:
        '''
        singlem pipe -1 /HK_annotation_base/trimmed/{wildcards.sample}*_trimmed_R1_paired.fq.gz \
        -2 /HK_annotation_base/trimmed/{wildcards.sample}*_trimmed_R2_paired.fq.gz \
        -p /data/singlem/australia/{wildcards.sample}_profile.tsv --threads {threads}
        '''
rule summarise:
    output:
        output='/data/singlem/australia/{sample}_profile_with_extras.tsv'
    threads:
        1
    conda:
        "/miniforge3/envs/singlem"
    shell:
        '''
        singlem summarise --input-taxonomic-profile /data/singlem/australia/{wildcards.sample}_profile.tsv \
        --output-taxonomic-profile-with-extras /data/singlem/australia/{wildcards.sample}_profile_with_extras.tsv
        '''