
import pandas as pd
import sys

## --------------------------------------------------------------------------------
## global parameters from config file

REF = config["ref"]
units = config["units"]
dir = config["outfol"]
angsd = config["angsd_dir"]
bam2prof = config["bam2prof_exec"]
lib_complexity = config["lib_complexity"]

## --------------------------------------------------------------------------------
## functions

def get_bam(wildcards):
    return unit_df.loc[(wildcards.id), "path"]

def ar_dir_error():
    print("Could not find AdapterRemoval settings directory automatically from units")
    print("Give up on getting library projections (library_complexity: no)")
    exit(1)


## --------------------------------------------------------------------------------
## helpers

unit_df = pd.read_table(config["units"]).set_index(["sampleId"])

ar_settings_dir = ""
if lib_complexity:
    import os
    dirs = [os.path.dirname(x) for x in list(unit_df["path"])]
    if len(set(dirs)) != 1: ar_dir_error()
    ar_dir = dirs[0]
    if not os.path.isdir(ar_dir): ar_dir_error()
    ar_settings_dir = "/".join(ar_dir.split("/")[0:-1]) + "/stats/"
    print("Doing library complexity analysis")
    superduper = config["superduper_exec"]

script_path = sys.path[0] + "/"





## --------------------------------------------------------------------------------
## targets

simple_res = expand(dir + "flagstat/{id}.txt", id = unit_df.index) + \
            expand(dir + "haplo/{id}.haplo", id = unit_df.index) + \
            expand(dir + "depth/{id}.depth", id = unit_df.index) + \
            expand(dir + "sex/{id}.sex", id = unit_df.index) + \
            expand(dir + "angsdX/{id}.res", id = unit_df.index) + \
            expand(dir + "damage/{id}.prof", id = unit_df.index) + \
            expand(expand(dir + "contamix/{{id}}_{cov}x_dif{dif}.mt.summary.txt", zip, cov = [1, 5], dif = [0.5, 0.7]), zip, id = unit_df.index)

complexity_res = [dir + "complexity/ar.full.txt"] + \
                expand(expand(dir + "complexity/{{id}}_DEPTH-{depth}.complexity.out", depth = [0.1,0.5,0.7,1,2,4,8,12]), zip, id = unit_df.index) + \
                expand(dir + "superduper/{id}.dupstat.txt", id = unit_df.index)

all_res = simple_res
if lib_complexity:
    all_res += complexity_res

rule all:
    input:
        all_res,
        expand(dir + "results.csv"),

onsuccess:
    print("Workflow finished, no error")

onerror:
    print("An error occurred")

## --------------------------------------------------------------------------------
## rules

rule flagstat:
    input:
        get_bam
    output:
        dir + "flagstat/{id}.txt"
    shell:
        "samtools flagstat {input} > {output}"


rule concat_ar:
    input:
        ar_settings_dir
    output:
        dir + "complexity/ar.full.txt",
        dir + "complexity/totreads.txt",
        dir + "complexity/discm1.txt",
    params:
        dir + "complexity/"
    shell:
        "bash {script_path}/helpers/parse_AR.sh {input} {params} {units}"

rule results:
    input:
        all_res
    output:
        dir + "results.csv",
    threads: 1
    shell:
        "bash {script_path}/helpers/csv_wrapper.sh {units} {dir} {lib_complexity} > {output}"

rule complex:
    input:
        totreads = dir + "complexity/totreads.txt",
        dupstat = dir + "superduper/{id}.dupstat.txt",
        table = dir + "superduper/{id}.table.txt",
    output:
        dir + "complexity/{id}_DEPTH-{depth}.complexity.out"
    params:
        dir + "complexity/{id}"
    threads: 1
    shell:
        "bash {script_path}/analyses/complex.sh -t {input.totreads} -a {input.table} -u {input.dupstat} -o {params} -d {wildcards.depth}"


rule haplo:
    input:
        get_bam
    output:
        dir + "haplo/{id}.haplo"
    threads: 1
    shell:
        "bash {script_path}/analyses/haplogrep.sh -b {input} -f {REF} -o {output}"

rule depth:
    input:
        get_bam
    output:
        dir + "depth/{id}.depth"
    threads: 4
    shell:
        "bash {script_path}/analyses/depth.sh -b {input} -t {threads} -o {output}"


rule sex:
    input:
        get_bam
    output:
        dir + "sex/{id}.sex"
    threads: 4
    shell:
        "bash {script_path}/analyses/sex.sh -b {input} -f {REF} -o {output}"


rule angsdX:
    input:
        get_bam
    output:
        dir + "angsdX/{id}.res"
    threads: 1
    shell:
        "bash {script_path}/analyses/angsdX.sh -b {input} -e {angsd} -f {REF} -o {output}"


rule dmg:
    input:
        get_bam
    output:
        dir + "damage/{id}.prof"
    threads: 1
    shell:
        "bash {script_path}/analyses/damage.sh -b {input} -e {bam2prof} -o {output}"


rule contamix:
    input:
        get_bam
    output:
        dir + "contamix/{id}_{cov}x_dif{dif}.mt.summary.txt"
    params:
        dir + "contamix/{id}_{cov}x_dif{dif}"
    threads: 1
    shell:
        "bash {script_path}/analyses/contamix.sh -b {input} -f {REF} -o {params} -c {wildcards.cov} -d {wildcards.dif}"


rule superduper:
    input:
        get_bam
    output:
        dir + "superduper/{id}.dupstat.txt",
        dir + "superduper/{id}.table.txt",
    params:
        dir + "superduper/{id}",
    threads: 4
    shell:
        "bash {script_path}/analyses/superduper.sh -b {input} -o {params} -t {threads} -e {superduper} -f {REF}"
