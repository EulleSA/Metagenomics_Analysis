from os.path import join
import git
from pathlib import Path

# Download eggnog-mapper via git

dir = Path("eggnog-mapper")

if dir.exists():
    pass
else:
    git.Git("./").clone("https://github.com/jhcepas/eggnog-mapper.git")

phredQuality = "20"

#funcionaldb = "databases/funcionaldb"

EGGNOG_OUTPUT = ["{eggnog_output}".format(eggnog_output=eggnog_output) for eggnog_output in ["seed_orthologs","annotations"]]

IDs =  glob_wildcards("data_raw/{sample}.fastq")

print(IDs.sample[0])

rule all:
    input:
        expand("{sample}_preprocessing.fastq",sample=IDs.sample[0]),
	["results/eggnog-mapper/{sample}.emapper.{eggnog_out}".format(sample=IDs.sample[0],eggnog_out=eggnog_output) for eggnog_output in EGGNOG_OUTPUT],
        "eggnog-mapper/data/eggnog.db",
        "eggnog-mapper/data/eggnog_proteins.dmnd",
        ["results/{sample}_count_pathways".format(sample=IDs.sample[0])],
        ["filtered_fastq/{sample}_filtred.fastq".format(sample=IDs.sample[0])]    
 
rule qualityControl:    
    input:
        "data_raw/{wildcards.sample}.fastq"
    output:
        preprocessing = "{wildcards.sample}_preprocessing.fastq",
        html_report = "{wildcards.sample}_report.html",
        json_report = "{wildcards.sample}_report.json"
	
    threads: workflow.cores
	
    shell:
        "fastp -i {input} -o {output.preprocessing} -q {phredQuality} -h {output.html_report} -j {output.json_report}"
 
rule download_funcionaldb:    
    output:
        "eggnog-mapper/data/eggnog.db",
        "eggnog-mapper/data/eggnog_proteins.dmnd"
    shell:
        "python2 eggnog-mapper/download_eggnog_data.py" 

rule functional_annot:    
    input:
        ["{sample}_preprocessing.fastq".format(sample=IDs.sample[0])] 
 
    output:
        seed_orthologs = ["results/eggnog-mapper/{sample}.emapper.seed_orthologs".format(sample=IDs.sample[0])],
        annotation = ["results/eggnog-mapper{sample}.emapper.annotations".format(sample=IDs.sample[0])], 
        file_name = ["{sample}".format(sample=IDs.sample[0])]
    threads: workflow.cores
    
    shell:
        "python2 eggnog-mapper/emapper.py -i {input} -o {output.file_name} -m diamond --translate --cpu {threads}"

rule treating_eggnog_output:
    input:
        annotation ="results/eggnog-mapper/{sample}.emapper.annotations".format(sample=IDs.sample[0]),
        pipeline = "scripts/r/pipeline.R"
    output:
        ["results/{sample}_count_pathways".format(sample=IDs.sample[0])],
        ["results/{sample}_reads_ids.txt".format(sample=IDs.sample[0])]
    shell:
        "Rscript --vanilla {input.pipeline} {input.annotation}" 

rule metagenoma_filtered:
    input:
        fastq_processing = ["{sample}_preprocessing.fastq".format(sample=IDs.sample[0])],
        list_ids = ["results/{sample}_reads_ids.txt".format(sample=IDs.sample[0])]
    output:
        ["filtered_fastq/{sample}_filtred.fastq".format(sample=IDs.sample[0])]
    threads:workflow.cores
    shell:
        "seqtk subseq {input.fastq_processing} {input.list_ids} > {output}"

#rule download-nrdb:
#    input:
#    output:
#    threads: workflow.cores
#    shell:
    

#rule alignment:
#    input:

#    output:

#    threads:workflow.cores
#    shell:
#        "diamond blastx {threads}"

#rule download_basta:

#rule basta:

#rule pythonnnn:
