
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
    sm = unit_df.loc[(wildcards.lib_lane), "sm"]
    lb = unit_df.loc[(wildcards.lib_lane), "lb"]
    id = unit_df.loc[(wildcards.lib_lane), "id"] + "-collapsed"
    txt = "@RG\\tID:%s\\tSM:%s\\tCN:CGG\\tPL:ILLUMINA\\tLB:%s\\tDS:seqcenter@sund.ku.dk" % (id, sm, lb)
    return txt

def make_rg_PE(wildcards):
    sm = unit_df.loc[(wildcards.lib_lane), "sm"]
    lb = unit_df.loc[(wildcards.lib_lane), "lb"]
    id = unit_df.loc[(wildcards.lib_lane), "id"] + "-paired"
    txt = "@RG\\tID:%s\\tSM:%s\\tCN:CGG\\tPL:ILLUMINA\\tLB:%s\\tDS:seqcenter@sund.ku.dk" % (id, sm, lb)
    return txt

def get_fastq(wildcards):
    return [unit_df.loc[(wildcards.lib_lane), "R1"],unit_df.loc[(wildcards.lib_lane), "R2"]]

## --------------------------------------------------------------------------------
## output lists

unit_df = pd.read_table(config["units"], comment="#").set_index(["idlane"])
TRIM_TYPE = ["collapsed.gz", "pair1.truncated.gz", "pair2.truncated.gz"]



## --------------------------------------------------------------------------------
## all

rule all:
    input:
        expand(out + "bams/{lib}.md.bam", lib=list(set(unit_df["id"])))



## --------------------------------------------------------------------------------
## wildcard constraints


wildcard_constraints:
    lib_lane = '|'.join(list(set(unit_df.index))),
    lib = '|'.join(list(set(unit_df["id"]))),
    sepe = '|'.join(["paired","collapsed"]),
    trim = '|'.join(TRIM_TYPE)


## --------------------------------------------------------------------------------
##  rules

rule trim:
    input:
        get_fastq
    output:
        c = out + "trimmed/{lib_lane}."+TRIM_TYPE[0],
        r1 = out + "trimmed/{lib_lane}."+TRIM_TYPE[1],
        r2 = out + "trimmed/{lib_lane}."+TRIM_TYPE[2],
        s = out + "trimmed/{lib_lane}.settings",
    params:
        out = out + "trimmed/{lib_lane}"
    threads: 6
    shell:
        "AdapterRemoval --threads {threads} --file1 {input[0]} --file2 {input[1]} --minlength 30 --adapter1 {a1} --adapter2 {a2} --collapse-conservatively --basename {params.out} --gzip "


rule aln:
    input:
         out + "trimmed/{lib_lane}.{trim}"
    output:
         temp(out + "tmp/{lib_lane}.{trim}.sai")
    threads: 12
    shell:
        "bwa aln -l {seed} -t {threads} {REF} {input} > {output}"


rule samse:
    input:
        fq1 = out + "trimmed/{lib_lane}.collapsed.gz",
        sai1 = out + "tmp/{lib_lane}.collapsed.gz.sai",
    output:
        temp(out + "tmp/{lib_lane}.collapsed.sam")
    params: make_rg_SE
    shell:
        "bwa samse -r \"{params}\" {REF} {input.sai1} {input.fq1} > {output}"


rule sampe:
    input:
        fq1 = out + "trimmed/{lib_lane}.pair1.truncated.gz",
        sai1 = out + "tmp/{lib_lane}.pair1.truncated.gz.sai",
        fq2 = out + "trimmed/{lib_lane}.pair2.truncated.gz",
        sai2 = out + "tmp/{lib_lane}.pair2.truncated.gz.sai",
    output:
        temp(out + "tmp/{lib_lane}.paired.sam")
    params: make_rg_PE

    shell:
        "bwa sampe -r \"{params}\" {REF} {input.sai1} {input.sai2} {input.fq1} {input.fq2} > {output}"

rule sam2bam:
    input:
        out + "tmp/{lib_lane}.{sepe}.sam"
    output:
        temp(out + "tmp/{lib_lane}.{sepe}.bam")
    threads: 6

    shell:
        "samtools sort {input} -@{threads} -m4G > {output}"

rule mapped:
    input:
        out + "tmp/{lib_lane}.{sepe}.bam"
    output:
        temp(out + "tmp/{lib_lane}.{sepe}.mapped.bam")
    threads: 6

    shell:
        "samtools view -F 4 -b {input} -@{threads} > {output}"


rule merge:
    input:
        lambda wildcards: [out+"tmp/"+x+"."+y+".mapped.bam" for y in ["paired", "collapsed"] for x in list(set(unit_df.index)) if wildcards.lib in x],
    output:
        temp(out + "tmp/{lib}.merged.bam")
    threads: 6

    shell:
        "samtools merge {output} {input} -@{threads}"

rule md:
    input:
        out + "tmp/{lib}.merged.bam"
    output:
        out + "bams/{lib}.md.bam"
    params:
        metrics = out + "bams/{lib}.metrics",
        tmpdir = out + "bams/"
    threads: 6
    shell:
         "java -Djava.io.tmpdir={params.tmpdir} -XX:ParallelGCThreads={threads} -Xmx2g -jar /maps/projects/lundbeck/apps/picard/build/libs/picard.jar MarkDuplicates OPTICAL_DUPLICATE_PIXEL_DISTANCE=12000 I={input} o={output} REMOVE_DUPLICATES=false METRICS_FILE={params.metrics} TAGGING_POLICY=All VALIDATION_STRINGENCY=LENIENT"
