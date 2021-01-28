from os.path import join
import git
from pathlib import Path

# Download eggnog-mapper via git

dir = Path("./tools/eggnog-mapper")
#dir2 = Path("./tools/BASTA")

if dir.exists():
    pass
else:
#    git.Git("./tools/").clone("https://github.com/timkahlke/BASTA.git")
    git.Git("./tools/").clone("https://github.com/jhcepas/eggnog-mapper.git")

phredQuality = "20"

#funcionaldb = "databases/funcionaldb"

EGGNOG_OUTPUT = ["{eggnog_output}".format(eggnog_output=eggnog_output) for eggnog_output in ["seed_orthologs","annotations"]]
TAXONOMIC_DB = "databases/taxonomic_db"
results_eggnogmap = "results/eggnog-mapper/"

IDs =  glob_wildcards("data_raw/{sample}.fastq")
file_name = IDs.sample[0].replace("_R1","")

print(file_name)

ruleorder: qualityControl_paired_end > qualityControl_single_end
ruleorder: remove_duplicates_paired > remove_duplicates_single

#expand("{sample}_preprocessing.fastq",sample=IDs.sample[0].replace("_R1","")),

rule all:
    input:
        ["results/eggnog-mapper/{sample}.emapper.{eggnog_out}".format(sample=IDs.sample[0].replace("_R1",""),eggnog_out=eggnog_output) for eggnog_output in EGGNOG_OUTPUT],
        "eggnog-mapper/data/eggnog.db",
        "eggnog-mapper/data/eggnog_proteins.dmnd",
        ["results/{sample}_count_pathways".format(sample=IDs.sample[0].replace("_R1",""))],
        ["filtered_fastq/{sample}_filtred.fastq".format(sample=IDs.sample[0].replace("_R1",""))],
        join(TAXONOMIC_DB,"nr0924.dmnd"),
        ["results/alignment/{sample}.m8".format(sample=IDs.sample[0].replace("_R1",""))],
        ["results/taxonomic/{sample}_best.txt".format(sample=IDs.sample[0].replace("_R1",""))]

rule qualityControl_single_end:
    input:
        ["data_raw/{sample}.fastq".format(sample=IDs.sample[0].replace("_R1",""))]
    output:
        #preprocessing = "{wildcards.sample}_preprocessing.fastq",
        #html_report = "{wildcards.sample}_report.html",
        #json_report = "{wildcards.sample}_report.json"
        preprocessing = ["data_raw/{sample}_preprocessing.fastq".format(sample=IDs.sample[0].replace("_R1",""))],
        html_report = ["data_raw/{sample}_report.html".format(sample=IDs.sample[0].replace("_R1",""))],
        json_report = ["data_raw/{sample}_report.json".format(sample=IDs.sample[0].replace("_R1",""))]
    threads: workflow.cores
    shell:
        "fastp -i {input} -o {output.preprocessing} -q {phredQuality} -h {output.html_report} -j {output.json_report}"

rule qualityControl_paired_end:
    input:
        forward = ["data_raw/{sample}_R1.fastq".format(sample=IDs.sample[0].replace("_R1",""))],
        reversee = ["data_raw/{sample}_R2.fastq".format(sample=IDs.sample[1].replace("_R2",""))]
    output:
        merged = ["qualityControl/{sample}_merged.fastq".format(sample=IDs.sample[0].replace("_R1",""))],
        forward_unmerged = ["qualityControl/{sample}_1_unmerged".format(sample=IDs.sample[0].replace("_R1",""))],
        reverse_unmerged = ["qualityControl/{sample}_2_unmerged".format(sample=IDs.sample[0].replace("_R1",""))],
        html_report = ["qualityControl/{sample}_report.html".format(sample=IDs.sample[0].replace("_R1",""))],
        json_report = ["qualityControl/{sample}_report.json".format(sample=IDs.sample[0].replace("_R1",""))]
    threads: workflow.cores
    shell:
        "fastp -i {input.forward} -I {input.reversee} -o {output.forward_unmerged} -O {output.reverse_unmerged} -q {phredQuality} --detect_adapter_for_pe -h {output.html_report} -j {output.json_report} -m --merged_out {output.merged} -w {threads}"

rule remove_duplicates_single:
    input:
        ["qualityControl/{sample}_preprocessing.fastq".format(sample=IDs.sample[0].replace("_R1",""))]
    output:
        ["qualityControl/collapse/{sample}_collapse.fasta".format(sample=IDs.sample[0].replace("_R1",""))]
    shell:
        "fastx_collapser -i {input} -o {output}"

rule remove_duplicates_paired:
    input:
        ["qualityControl/{sample}_merged.fastq".format(sample=IDs.sample[0].replace("_R1",""))]
    output:
        ["qualityControl/collapse/{sample}_collapse.fasta".format(sample=IDs.sample[0].replace("_R1",""))]
    shell:
        "fastx_collapser -i {input} -o {output}"
    
rule download_funcionaldb:
    output:
        "eggnog-mapper/data/eggnog.db",
        "eggnog-mapper/data/eggnog_proteins.dmnd"
    shell:
        "python2 eggnog-mapper/download_eggnog_data.py"

rule functional_annot:
    input:
        reads = ["qualityControl/collapse/{sample}_collapse.fasta".format(sample=IDs.sample[0].replace("_R1",""))]
        #sample = ["data_raw/{sample}_R1.fastq".format(sample=IDs.sample[0].replace("_R1",""))]
    output:
        seed_orthologs = ["results/eggnog-mapper/{sample}.emapper.seed_orthologs".format(sample=IDs.sample[0].replace("_R1",""))],
        annotation = ["results/eggnog-mapper/{sample}.emapper.annotations".format(sample=IDs.sample[0].replace("_R1",""))],
        #file_name = ["{results_eggnogmap}/{sample}".format(results_eggnogmap=results_eggnogmap,sample=IDs.sample[0].replace("_R1",""))]
    threads: workflow.cores
    shell:
        "python2 eggnog-mapper/emapper.py -i {input.reads} -o {file_name} --output_dir {results_eggnogmap} -m diamond --translate --no_file_comments --cpu {threads}"

rule treating_eggnog_output:
    input:
        annotation =["results/eggnog-mapper/{sample}.emapper.annotations".format(sample=IDs.sample[0].replace("_R1",""))],
        pipeline = "scripts/r/pipeline.R"
    output:
        ["results/{sample}_count_pathways".format(sample=IDs.sample[0].replace("_R1",""))],
        ["results/{sample}_reads_ids.txt".format(sample=IDs.sample[0].replace("_R1",""))]
    shell:
        "Rscript --vanilla {input.pipeline} {input.annotation}"

rule metagenoma_filtered:
    input:
        fastq_processing = ["qualityControl/collapse/{sample}_collapse.fasta".format(sample=IDs.sample[0].replace("_R1",""))],
        list_ids = ["results/{sample}_reads_ids.txt".format(sample=IDs.sample[0].replace("_R1",""))]
    output:
        ["filtered_fastq/{sample}_filtred.fastq".format(sample=IDs.sample[0].replace("_R1",""))]
    threads:workflow.cores
    shell:
        "seqtk subseq {input.fastq_processing} {input.list_ids} > {output}"

#rule diamondMakeDB_Taxonomic:
#    output:
#        index = join(TAXONOMIC_DB,"nr0924.dmnd")
#    threads: workflow.cores
#    shell: "wget -P {TAXONOMIC_DB} \
#        ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz \
#        && pigz -d nr.gz -p {threads} && diamond makedb --in nr -d {output.index} --threads {threads}"

rule alignment_metagenome_filtered:
    input:
        index = join(TAXONOMIC_DB,"nr0924.dmnd"),
        reads = ["filtered_fastq/{sample}_filtred.fastq".format(sample=IDs.sample[0].replace("_R1",""))]
    output:
        matches = ["results/alignment/{sample}.m8".format(sample=IDs.sample[0].replace("_R1",""))],
        unaligned = ["results/alignment/{sample}_unaligned.fasta".format(sample=IDs.sample[0].replace("_R1",""))]
    threads:workflow.cores
    shell:
        "diamond blastx -d {input.index} -q {input.reads} -o {output.matches} --top 3 --more-sensitive --un {output.unaligned} --threads {threads} "

#rule download_bastadb:
#    output:
#        join(TAXONOMIC_DB,"complete_taxa.db"),
#        join(TAXONOMIC_DB,"prot_mapping.db")
    #conda: env
#    shell: "../miniconda3/envs/cactus/bin/basta taxonomy -d {TAXONOMIC_DB} \
#        && ../miniconda3/envs/cactus/bin/basta download prot -d {TAXONOMIC_DB}"

rule taxonomicClassification:
    input:
        matches = ["results/alignment/{sample}.m8".format(sample=IDs.sample[0].replace("_R1",""))],
        leveldb = {TAXONOMIC_DB}
    output:
        out = ["results/taxonomic/{sample}_tax.txt".format(sample=IDs.sample[0].replace("_R1",""))],
        lca = ["results/taxonomic/{sample}_lca.html".format(sample=IDs.sample[0].replace("_R1",""))],
        best = ["results/taxonomic/{sample}_best.html".format(sample=IDs.sample[0].replace("_R1",""))],
        
        lca_txt = ["results/taxonomic/{sample}_lca.txt".format(sample=IDs.sample[0].replace("_R1",""))],
        best_txt = ["results/taxonomic/{sample}_best.txt".format(sample=IDs.sample[0].replace("_R1",""))]
    #conda: env
    threads: workflow.cores
    shell: "basta sequence {input.matches} {output.out} prot -d {input.leveldb} -l 1 -m 1 -b True \
        && awk -F \"\t\" '{{print $1\"\t\"$2}}' {output.out} > {output.lca_txt} \
        && awk -F \"\t\" '{{print $1\"\t\"$3}}' {output.out} > {output.best_txt} \
        && basta2krona.py {output.lca_txt} {output.lca} \
        && basta2krona.py {output.best_txt} {output.best}"

#rule basta:

#rule pythonnnn:
