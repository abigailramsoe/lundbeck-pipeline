
import pandas as pd

## --------------------------------------------------------------------------------
## global parameters from config file

REF = config["ref"]
units = config["units"]
dir = config["outfol"]
angsd = config["angsd_dir"]
bam2prof = config["bam2prof_exec"]
superduper = config["superduper_exec"]

## --------------------------------------------------------------------------------
## helpers

CHROMS = range(1, 23)

unit_df = pd.read_table(config["units"]).set_index(["sampleId"])
all_bams = unit_df["path"]

## --------------------------------------------------------------------------------
## functions

def get_bam(wildcards):
    return unit_df.loc[(wildcards.id), "path"]


## --------------------------------------------------------------------------------
## constraints

wildcard_constraints:
    full = '|'.join(all_bams),
    #lib = '|'.join(SAMPLES),

## --------------------------------------------------------------------------------
## targets

rule all:
    input:
        expand(dir + "superduper/{id}.dupstat.txt", id = unit_df.index),
        expand(dir + "haplo/{id}.haplo", id = unit_df.index),
        expand(dir + "depth/{id}.depth", id = unit_df.index),
        expand(dir + "sex/{id}.sex", id = unit_df.index),
        expand(dir + "angsdX/{id}.res", id = unit_df.index),
        expand(dir + "damage/{id}.prof", id = unit_df.index),
        expand(expand(dir + "contamix/{{id}}_{cov}x_dif{dif}.mt.summary.txt", zip, cov = [1, 5], dif = [0.5, 0.7]), zip, id = unit_df.index),


onsuccess:
    print("Workflow finished, no error")
    #tt, sdate = config["index"].split("/")[-1].split(".")[0].split("-")
    #csv_name = "%s-%s.csv" % (sdate, tt)
    #shellmsg = 'bash scripts/csv_wrapper.sh '+config["index"]+' '+config["pipe"] + '> '+dir+'/stats-per-run/csv/'+ csv_name
    #shell(shellmsg)
    #shellmsg = 'mail -s "%s stats FINISHED" abigail@palaeome.org < {log}' % RUN
    #shell(shellmsg)



onerror:
    print("An error occurred")
    #shellmsg = 'mail -s "%s stats ERROR" abigail@palaeome.org < {log}' % RUN
    #shell(shellmsg)
## --------------------------------------------------------------------------------
## rules




rule haplo:
    input:
        get_bam
    output:
        dir + "haplo/{id}.haplo"
    threads: 1
    shell:
        "bash analyses/haplogrep.sh -b {input} -f {REF} -o {output}"

rule depth:
    input:
        get_bam
    output:
        dir + "depth/{id}.depth"
    threads: 4
    shell:
        "bash analyses/depth.sh -b {input} -t {threads} -o {output}"


rule sex:
    input:
        get_bam
    output:
        dir + "sex/{id}.sex"
    threads: 4
    shell:
        "bash analyses/sex.sh -b {input} -f {REF} -o {output}"


rule angsdX:
    input:
        get_bam
    output:
        dir + "angsdX/{id}.res"
    threads: 1
    shell:
        "bash analyses/angsdX.sh -b {input} -e {angsd} -f {REF} -o {output}"


rule dmg:
    input:
        get_bam
    output:
        dir + "damage/{id}.prof"
    threads: 1
    shell:
        "bash analyses/damage.sh -b {input} -e {bam2prof} -o {output}"



rule contamix:
    input:
        get_bam
    output:
        dir + "contamix/{id}_{cov}x_dif{dif}.mt.summary.txt"
    params:
        dir + "contamix/{id}_{cov}x_dif{dif}"
    threads: 1
    shell:
        "bash analyses/contamix.sh -b {input} -f {REF} -o {params} -c {wildcards.cov} -d {wildcards.dif}"


rule superduper:
    input:
        get_bam
    output:
        dir + "superduper/{id}.dupstat.txt",
    params:
        dir + "superduper/{id}",
    threads: 4
    shell:
        "bash analyses/superduper.sh -b {input} -o {params} -t {threads} -e {superduper} -f {REF}"
