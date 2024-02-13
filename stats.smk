
import pandas as pd

## --------------------------------------------------------------------------------
## global parameters from config file

REF = config["ref"]
units = config["units"]
dir = config["outfol"]
angsd = config["angsd_dir"]

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
        #expand(dir + "results/{cgg}/superduper/{lib}.dupstat.txt", zip, lib = SAMPLES, cgg = CGGs),
        #expand(dir + "results/{cgg}/superduper/{lib}.table.txt", zip, lib = SAMPLES, cgg = CGGs),
        expand(dir + "haplo/{id}.haplo", id = unit_df.index),
        expand(dir + "depth/{id}.depth", id = unit_df.index),
        expand(dir + "sex/{id}.sex", id = unit_df.index),
        expand(dir + "angsdX/{id}.res", id = unit_df.index),

        #expand(dir + "results/{cgg}/contX/{lib}.res", zip, lib = SAMPLES, cgg = CGGs),
        #expand(dir + "results/{cgg}/metadamage/{lib}.prof", zip, lib = SAMPLES, cgg = CGGs),
        #expand(expand(dir + "results/{{cgg}}/contMT_{cov}xdif{dif}/{{lib}}_mt.summary.txt", zip, cov = [1, 5], dif = [5, 7]), zip, lib = SAMPLES, cgg = CGGs), # fix this s


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
        "bash analyses/angsdX.sh -b {input} -a {angsd} -f {REF} -o {output}"


rule dmg:
    input:
        cram = dir + "{cgg}/{lib}.cram",
        idx = dir + "{cgg}/{lib}.cram.crai"
    output:
        dir + "results/{cgg}/metadamage/{lib}.prof"
    threads: 1
    shell:
        "mkdir -p $(dirname {output}); bash /projects/lundbeck/apps/scripts/dmg.sh {input.cram} {output}"



rule contamix:
    input:
        cram = dir + "{cgg}/{lib}.cram",
        idx = dir + "{cgg}/{lib}.cram.crai"
    output:
        dir + "results/{cgg}/contMT_{cov}xdif{dif}/{lib}_mt.summary.txt"
    params:
        dir + "results/{cgg}/contMT_{cov}xdif{dif}/"
    threads: 3
    shell:
        "mkdir -p $(dirname {output}); bash /projects/lundbeck/apps/scripts/contamix.sh {input.cram} {wildcards.cov} {wildcards.dif} {params} {PIPE}"


rule superduper:
    input:
        cram = dir + "{cgg}/{lib}.cram",
        idx = dir + "{cgg}/{lib}.cram.crai"
    output:
        dupstat = dir + "results/{cgg}/superduper/{lib}.dupstat.txt",
        table = dir + "results/{cgg}/superduper/{lib}.table.txt"
    params:
        dir + "results/{cgg}/superduper/{lib}"
    threads: 4
    shell:
        "mkdir -p $(dirname {output.table}); bash /projects/lundbeck/apps/scripts/superduper.sh {input.cram} {PIPE} {params} {threads}; touch {output.dupstat}; touch {output.table}"
