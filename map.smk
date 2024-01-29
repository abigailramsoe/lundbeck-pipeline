
import pandas as pd
import glob

## --------------------------------------------------------------------------------
## config parameters

REF = config["ref"]
seed = config["seed"]
suffix = config["suffix"]
a1 = config["a1"]
a2 = config["a2"]
units = config["units"]
out = config["outfol"]


## --------------------------------------------------------------------------------
## helper functions

def make_rg_SE(wildcards):
   sm = str(wildcards.lib).split("-")[0][0:6]
   lb = str(wildcards.lib).split("-")[3]
   id = str(wildcards.lib) + "-collapsed"
   txt = "@RG\\tID:%s\\tSM:%s\\tCN:CGG\\tPL:ILLUMINA\\tLB:%s\\tDS:seqcenter@sund.ku.dk" % (id, sm, lb)
   return txt

def make_rg_PE(wildcards):
   sm = str(wildcards.lib).split("-")[0][0:6]
   lb = str(wildcards.lib).split("-")[3]
   id = str(wildcards.lib) + "-paired"
   txt = "@RG\\tID:%s\\tSM:%s\\tCN:CGG\\tPL:ILLUMINA\\tLB:%s\\tDS:seqcenter@sund.ku.dk" % (id, sm, lb)
   return txt

def get_new_name(wildcards):
    cram = [x.replace("-","-")+suffix for x in wildcards.libs]
    return cram


def get_new_name(wildcards):
    old = new_name_dict[wildcards.cram]
    return old

def get_mem_mb(wildcards, threads):
    return threads * 150

## --------------------------------------------------------------------------------
## dirs

TRIM_TYPE = ["collapsed.gz", "pair1.truncated.gz", "pair2.truncated.gz"]

## --------------------------------------------------------------------------------
## pool level parameters

fqdir="./"


## --------------------------------------------------------------------------------
## output lists

unit_df = pd.read_table(config["units"], comment="#").set_index(["fastq"])
print(unit_df)
fastqs = [x.split("/")[-1] for x in glob.glob(fqdir + "*.fastq.gz")]
fastqs_noext = list(set(["_".join(x.replace(".fastq.gz", "").split("_")[0:3]) for x in fastqs]))
fastqs_trim = [x + ".paired" for x in fastqs_noext] + [x + ".collapsed" for x in fastqs_noext]
fastqs_libs = list(set([x.split("_")[0] for x in fastqs_noext]))
cram_names = [x.replace("-","_")+suffix for x in fastqs_libs]

new_name_dict = {}
for x in fastqs_libs:
    new = x.replace("-","_")+suffix
    new_name_dict[new] = out + "crams/" + x + ".cram"

## --------------------------------------------------------------------------------
## all

rule all:
    input:
            expand(out + "crams/{libs}{suffix}.cram", libs = fastqs_libs, suffix=suffix),


## --------------------------------------------------------------------------------
## wildcard constraints


wildcard_constraints:
    fq_fullname = '|'.join(fastqs_noext),
    fq_with_trim = '|'.join(fastqs_trim),
    libs = '|'.join(fastqs_libs),
    cram = '|'.join(new_name_dict.keys())


## --------------------------------------------------------------------------------
##  rules

rule trim:
    input:
        R1 = fqdir + "{fq_fullname}_R1_001.fastq.gz",
        R2 = fqdir + "{fq_fullname}_R2_001.fastq.gz"
    output:
        r1 = out + "trimmed/{fq_fullname}."+TRIM_TYPE[1],
        r2 = out + "trimmed/{fq_fullname}."+TRIM_TYPE[2],
        c = out + "trimmed/{fq_fullname}."+TRIM_TYPE[0],
        s = out + "trimmed/{fq_fullname}.settings",
    params:
        out = out + "trimmed/{fq_fullname}"
    threads: 6
    shell:
        "/maps/projects/lundbeck/apps/adapterremoval-2.3.2/build/AdapterRemoval --threads {threads} --file1 {input.R1} --file2 {input.R2} --minlength 30 --adapter1 {a1} --adapter2 {a2} --collapse-conservatively --basename {params.out} --gzip "


rule aln:
    input:
         out + "trimmed/{fq_fullname}.{trim}"
    output:
         out + "map/{fq_fullname}.{trim}.sai"
    threads: 12
    shell:
        "bwa aln -l {seed} -t {threads} {REF} {input} > {output}"


rule samse:
    input:
        fq1 = out + "trimmed/{lib}.collapsed.gz",
        sai1 = out + "map/{lib}.collapsed.gz.sai",
    output:
        out + "sams/{lib}.collapsed.sam"
    params: make_rg_SE

    shell:
        "bwa samse -r \"{params}\" {REF} {input.sai1} {input.fq1} > {output}"


rule sampe:
    input:
        fq1 = out + "trimmed/{lib}.pair1.truncated.gz",
        sai1 = out + "map/{lib}.pair1.truncated.gz.sai",
        fq2 = out + "trimmed/{lib}.pair2.truncated.gz",
        sai2 = out + "map/{lib}.pair2.truncated.gz.sai",
    output:
        out + "sams/{lib}.paired.sam"
    params: make_rg_PE

    shell:
        "bwa sampe -r \"{params}\" {REF} {input.sai1} {input.sai2} {input.fq1} {input.fq2} > {output}"

rule sam2bam:
    input:
        out + "sams/{fq_with_trim}.sam"
    output:
        out + "bams/{fq_with_trim}.bam"
    threads: 6

    shell:
        "samtools sort {input} -@{threads} -m4G > {output}"

rule mapped:
    input:
        out + "bams/{fq_with_trim}.bam"
    output:
        out + "bams/{fq_with_trim}.mapped.bam"
    threads: 6

    shell:
        "samtools view -F 4 -b {input} -@{threads} > {output}"


rule merge:
    input:
        lambda wildcards: [out+"bams/"+x+".mapped.bam" for x in fastqs_trim if wildcards.libs in x],
    output:
        out + "bams/{libs}.bam"
    threads: 6

    shell:
        "samtools merge {output} {input} -@{threads}"

rule md:
    input:
        out + "bams/{libs}.bam"
    output:
        out + "bams/{libs}.md.bam"

    params:
        metrics = out + "bams/{libs}.metrics",
        tmpdir = out + "bams/"
    threads: 6
    shell:
         "java -Djava.io.tmpdir={params.tmpdir} -XX:ParallelGCThreads={threads} -Xmx2g -jar /maps/projects/lundbeck/apps/picard/build/libs/picard.jar MarkDuplicates OPTICAL_DUPLICATE_PIXEL_DISTANCE=12000 I={input} o={output} REMOVE_DUPLICATES=false METRICS_FILE={params.metrics} TAGGING_POLICY=All VALIDATION_STRINGENCY=LENIENT"

rule cram:
    input:
        out + "bams/{libs}.md.bam"
    output:
        out + "crams/{libs}{suffix}.cram"

    threads: 6
    shell:
        "samtools view {input} -F 4 -C -T {REF} -o {output} -@{threads}"
