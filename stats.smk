
import pandas as pd

## --------------------------------------------------------------------------------
## global parameters from config file

REF = config["ref"]
units = config["units"]
dir = config["outfol"]


## --------------------------------------------------------------------------------
## helpers

CHROMS = range(1, 23)

unit_df = pd.read_table(config["units"]).set_index(["sampleId"])
print(unit_df)
all_bams = unit_df["path"]
#unit_df["CGG"] = [x.split("/")[4] for x in unit_df["Cram"]]
#unit_df["basename"] = [x.split("/")[-1].split(".")[0] for x in unit_df["Cram"]]
#SAMPLES = unit_df["basename"].tolist()
#CGGs = ["CGG" + x.split("_")[0][0:6] for x in unit_df["basename"].tolist()]
#RUN = unit_df["Run"].tolist()[0]
#unit_df = unit_df.set_index(["basename"])
print('|'.join(all_bams))
def get_extension(lib_name):
    if lib_name.endswith('.cram'):
        return 'crai'
    else:
        return 'bai'

## --------------------------------------------------------------------------------
## functions

def get_bam(wildcards):
    return unit_df.loc[(wildcards.id), "path"]


## --------------------------------------------------------------------------------
## constraints

wildcard_constraints:
    full = '|'.join(all_bams),
    ext = "|".join(["crai", "bai"])
    #lib = '|'.join(SAMPLES),

## --------------------------------------------------------------------------------
## targets

rule all:
    input:
        #expand("{full}.{ext}", full=all_bams, ext=[get_extension(x) for x in all_bams]),
        #expand(dir + "results/{cgg}/depth/{lib}.depth", zip, lib = SAMPLES, cgg = CGGs),
        #expand(dir + "results/{cgg}/superduper/{lib}.dupstat.txt", zip, lib = SAMPLES, cgg = CGGs),
        #expand(dir + "results/{cgg}/superduper/{lib}.table.txt", zip, lib = SAMPLES, cgg = CGGs),
        expand(dir + "results/haplo/{id}.haplo", id = unit_df.index),
        #expand(dir + "results/{cgg}/sex/{lib}.sex", zip, lib = SAMPLES, cgg = CGGs),
        #expand(dir + "results/{cgg}/contX/{lib}.res", zip, lib = SAMPLES, cgg = CGGs),
        #expand(dir + "results/{cgg}/metadamage/{lib}.prof", zip, lib = SAMPLES, cgg = CGGs),
        #expand(expand(dir + "results/{{cgg}}/contMT_{cov}xdif{dif}/{{lib}}_mt.summary.txt", zip, cov = [1, 5], dif = [5, 7]), zip, lib = SAMPLES, cgg = CGGs), # fix this s
        #expand(expand(dir + "results/{{cgg}}/complex/{{lib}}_DEPTH-{depth}.complexity.out", depth = [0.1,0.5,0.7,1,2,4,8,12]), zip, lib = SAMPLES, cgg = CGGs), # fix this s


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
        dir + "results/haplo/{id}.haplo"
    params:
        basename = dir + "results/haplo/{id}",
    threads: 1
    shell:
        "bash analyses/haplogrep.sh -b {input} -f {REF} -o {params.basename}"


rule complex:
    input:
        cram = dir + "{cgg}/{lib}.cram",
        idx = dir + "{cgg}/{lib}.cram.crai",
        totreads = dir + "stats-per-run/hi/totreads.txt",
        dupstat = dir + "results/{cgg}/superduper/{lib}.dupstat.txt",
        table = dir + "results/{cgg}/superduper/{lib}.table.txt",
    output:
        dir + "results/{cgg}/complex/{lib}_DEPTH-{depth}.complexity.out"
    params:
        dir + "results/{cgg}/complex/"
    threads: 1
    shell:
        "mkdir -p $(dirname {output}); bash /projects/lundbeck/apps/scripts/complex.sh {input.cram} {input.totreads} {input.table} {input.dupstat} {params} {wildcards.depth}"



rule depth:
    input:
        get_bam
    output:
        dir + "results/{cgg}/depth/{lib}.depth"
    threads: 4
    shell:
        "mkdir -p $(dirname {output}); bash /projects/lundbeck/apps/scripts/depth.sh {input.cram} {PIPE} {threads} > {output}"


rule sex:
    input:
        cram = dir + "{cgg}/{lib}.cram",
        idx = dir + "{cgg}/{lib}.cram.crai",
    output:
        dir + "results/{cgg}/sex/{lib}.sex"
    threads: 1
    shell:
        "mkdir -p $(dirname {output}); bash /projects/lundbeck/apps/scripts/sex.sh {input.cram} new > {output}"


rule angsdX:
    input:
        cram = dir + "{cgg}/{lib}.cram",
        idx = dir + "{cgg}/{lib}.cram.crai"
    output:
        dir + "results/{cgg}/contX/{lib}.res"
    threads: 1
    shell:
        "mkdir -p $(dirname {output}); bash /projects/lundbeck/apps/scripts/angsdX.sh {input.cram} {PIPE} {output}"


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