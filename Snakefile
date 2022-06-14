from typing import Any, List
import os

configfile: "config.yaml"

ITER_TOOLS = config['polish_rounds'].copy()
del ITER_TOOLS['flye']
ITER_TOOLS['minimap'] = ITER_TOOLS['racon']
ITER_TOOLS['bwa'] = ITER_TOOLS['pilon']

def racon_input(wildcards)->str:
    global config
    n: int = int(wildcards.n)
    if n == 1:
        return os.path.join(config['out_dir'],
                            wildcards.strain,
                            'flye',
                            'assembly.fasta')
    else:
        return os.path.join(config['out_dir'],
                            wildcards.strain,
                            'racon',
                            f'iteration_{n-1}',
                            'racon_polished.fasta')


def pilon_input(wildcards)->str:
    global config
    n: int = int(wildcards.n)
    if n == 1:
        return os.path.join(config['out_dir'],
                            wildcards.strain,
                            'medaka',
                             'consensus.fasta')
    else:
        return os.path.join(config['out_dir'],
                            wildcards.strain,
                            'pilon',
                            f'iteration_{n-1}',
                            'pilon.fasta')


int_keys: List[int] = [x for x in config['strains'].keys() if isinstance(x, int)]
read_keys: List[str] = ['strains', 'forward', 'reverse']
key: Any
for key in int_keys:
    read_key: str
    for read_key in read_keys:
        config[read_key][str(key)] = config[read_key][key]
        del config[read_key][key]

if not os.path.exists(config['out_dir']):
    os.makedirs(config['out_dir'])
    
if not os.path.exists(config['out_dir']):
    os.makedirs(config['out_dir'])

strain: str
for strain in config['strains'].keys():
    os.makedirs(os.path.join(config['out_dir'], f'{strain}'), exist_ok=True)
    os.makedirs(os.path.join(config['out_dir'], f'{strain}', 'medaka'), exist_ok=True)
    tool: str
    for tool in ITER_TOOLS.keys():
        it: int
        for it in range(1, ITER_TOOLS[tool] + 1):
            os.makedirs(os.path.join(config['out_dir'], f'{strain}', tool, f'iteration_{it}'), exist_ok=True)


rule all:
    input:
#        expand("{outdir}/{strain}/busco/short_summary.{strain}.txt",
        expand("{outdir}/{strain}/cry_processor/raw_full_{strain}.fasta",
               outdir = config['out_dir'],
               strain = config['strains']),
        expand("{outdir}/{strain}/idops/idops_scan.tsv",
               outdir = config['out_dir'],
               strain = config['strains']),
        expand("{outdir}/{strain}/sortpred/predicted_sortase.csv",
               outdir = config['out_dir'],
               strain = config['strains']),
        expand("{outdir}/{strain}/phigaro/pilon.phigaro.tsv",
               outdir = config['out_dir'],
               strain = config['strains']),
        expand("{outdir}/{strain}/antismash/{strain}.gbk",
               outdir = config['out_dir'],
               strain = config['strains']),
        expand("{outdir}/{strain}/deepbgc/README.txt",
               outdir = config['out_dir'],
               strain = config['strains']),
        expand("{outdir}/{strain}/eggnog/{strain}.emapper.annotations",
               outdir = config['out_dir'],
               strain = config['strains'])
    output: touch('.status')


rule assemble_flye:
    input:
        reads = lambda wildcards: config['strains'][wildcards.strain]
    output:
        "{outdir}/{strain}/flye/assembly.fasta"
    params:
        outdir = "{outdir}/{strain}/flye",
        polish_rounds = config['polish_rounds']['flye']
    threads: config['threads']
    shell:
        "mkdir -p {params.outdir} && \
         flye \
          --nano-raw {input.reads} \
          -t {threads} \
          -i {params.polish_rounds} \
          -o {params.outdir}"


rule racon_polish:
    input:
        assembly = racon_input,
        reads = lambda wildcards: config['strains'][wildcards.strain]
    output:
        minimap_out = "{outdir}/{strain}/minimap2/iteration_{n}/aln.sam",
        racon_out = "{outdir}/{strain}/racon/iteration_{n}/racon_polished.fasta"
    params:
        outdir = config['out_dir'],
        minimap_dir = "{outdir}/{strain}/minimap2/iteration_{n}",
        racon_dir = "{outdir}/{strain}/racon/iteration_{n}"
    wildcard_constraints:
        n = f"[1-{ITER_TOOLS['racon']}]"
    threads: config['threads']
    shell:
        "mkdir -p {params.minimap_dir} && \
         mkdir -p {params.racon_dir} && \
         minimap2 \
           -t {threads} \
           -ax \
           map-ont \
           {input.assembly} \
           {input.reads} \
           > {output.minimap_out} && \
         sed \
           -i \
           '/\[/d' \
           {output.minimap_out} && \
         racon \
           -t {threads} \
           {input.reads} \
           {output.minimap_out} \
           {input.assembly} \
           > {output.racon_out} && \
         sed \
           -i \
           '/\[/d' \
           {output.racon_out}"


rule medaka:
    input:
        assembly = "{outdir}/{strain}/racon/iteration_4/racon_polished.fasta",
        reads = lambda wildcards: config['strains'][wildcards.strain]
    output:
        "{outdir}/{strain}/medaka/consensus.fasta"
    params:
        outdir = config['out_dir'],
        medaka_dir = "{outdir}/{strain}/medaka"
    threads: config['threads']
    shell:
        "mkdir -p {params.medaka_dir} && \
         medaka_consensus \
           -i {input.reads} \
           -d {input.assembly} \
           -o {params.medaka_dir} \
           -t {threads}"


rule pilon_polish:
        input:
            ref = pilon_input,
            f = lambda wildcards: config['forward'][wildcards.strain],
            r =  lambda wildcards: config['reverse'][wildcards.strain]
        output:
            bwa = "{outdir}/{strain}/bwa/iteration_{n}/bwa.sam",
            pilon = "{outdir}/{strain}/pilon/iteration_{n}/pilon.fasta"
        params:
            bwa_dir = "{outdir}/{strain}/bwa/iteration_{n}",
            pilon_dir = "{outdir}/{strain}/pilon/iteration_{n}"
        wildcard_constraints:
            n = f"[1-{ITER_TOOLS['pilon']}]"
        threads: config['threads']
        shell:
            "mkdir -p {params.bwa_dir} && \
             mkdir -p {params.pilon_dir} && \
             bwa index \
               -a bwtsw \
               {input.ref} && \
             bwa mem \
               {input.ref} \
               {input.f} \
               {input.r} \
               -t {threads} \
               > {output.bwa} && \
             samtools view \
               -S \
               -b \
               {output.bwa} \
               > {output.bwa}.bam && \
             samtools sort \
               -@ config['threads'] \
               -o {output.bwa}.bam \
               {output.bwa}.bam && \
             samtools index \
               {output.bwa}.bam && \
             pilon \
               --genome {input.ref} \
               --frags {output.bwa}.bam \
               --outdir {params.pilon_dir} && \
             sed -i 's/_pilon//g' {output.pilon}"


rule prokka:
    input:
        contigs = "{outdir}/{strain}/pilon/iteration_4/pilon.fasta",
        model = config['prokka']['prodigal_model'],
        ref_prot = config['prokka']['reference_proteins']
    output:
        "{outdir}/{strain}/prokka/{strain}.faa"
    params:
        out_dir = "{outdir}/{strain}/prokka",
        prefix = lambda wildcards: wildcards.strain
    conda:
        "/home/yura/anaconda3/envs/prokka.yaml"
    threads: config['threads']
    shell:
        "prokka \
          --addgenes \
          --force \
          --cpus {threads} \
          --prefix {params.prefix} \
          --prodigaltf {input.model} \
          --proteins {input.ref_prot} \
          --outdir {params.out_dir} \
          {input.contigs}"


rule cry_processor:
    input:
        "{outdir}/{strain}/prokka/{strain}.faa"
    output:
        "{outdir}/{strain}/cry_processor/raw_full_{strain}.fasta"
    params:
        all_domains = 1 if config['cry_processor']['all_domains'] else 2,
        mode = config['cry_processor']['mode'],
        out_dir = "{outdir}/{strain}/cry_processor"
    threads: config['threads']
    shell:
        "cry_processor.py \
          --annotate \
          --force \
          -fi {input} \
          -pr {params.all_domains} \
          -r {params.mode} \
          -od {params.out_dir}"


rule idops:
    input:
        "{outdir}/{strain}/prokka/{strain}.faa"
    output:
        "{outdir}/{strain}/idops/idops_scan.tsv"
    params:
         out_dir = "{outdir}/{strain}/idops"
    threads: config['threads']
    conda:
       "/home/yura/anaconda3/envs/idops.yaml"
    shell:
        "idops \
          -o {params.out_dir} \
          -t \
          {input}"


rule phigaro:
    input:
        "{outdir}/{strain}/pilon/iteration_4/pilon.fasta"
    output:
        "{outdir}/{strain}/phigaro/pilon.phigaro.tsv"
    params:
         out_dir = "{outdir}/{strain}/phigaro"
    conda:
        "/home/yura/anaconda3/envs/phigaro.yaml"
    threads: config['threads']
    shell:
        "phigaro \
          -e tsv bed gff \
          -f {input} \
          -o {params.out_dir} \
          -p \
          -t {threads} \
          --delete-shorts \
          --not-open"


rule sortpred_prep:
    input:
        "{outdir}/{strain}/prokka/{strain}.faa"
    output:
        "{outdir}/{strain}/prokka/{strain}_short.faa"
    threads:
        config['threads']
    shell:
        "lengthFilter.py \
           -i {input} \
           -o {output}"


rule sortpred:
    input:
        "{outdir}/{strain}/prokka/{strain}_short.faa"
    output:
        "{outdir}/{strain}/sortpred/predicted_sortase.csv"
    params:
        script = "/home/yura/bacillus_pangenomics/SortPred_standalone/server.R",
        out_dir = "{outdir}/{strain}/sortpred",
        model1 = config['sortpred']['model1'],
        model2 = config['sortpred']['model2']
    conda:
        "/home/yura/anaconda3/envs/R_bacillus.yaml"
    shell:
        "Rscript {params.script} \
            --layer1 {params.model1} \
            --layer2 {params.model2} \
            -i {input} \
            -o {params.out_dir}"


rule antismash:
    input:
        "{outdir}/{strain}/pilon/iteration_4/pilon.fasta"
    output:
        "{outdir}/{strain}/antismash/{strain}.gbk"
    params:
        out_dir = "{outdir}/{strain}/antismash",
        features = "{outdir}/{strain}/prokka/{strain}.gff",
        strain = lambda wildcards: wildcards.strain
    threads: config['threads']
    conda:
        "/home/yura/anaconda3/envs/antismash.yaml"
    shell:
        "antismash \
          --cb-knownclusters \
          --fullhmmer \
          --cpus {threads} \
          --output-basename {params.strain} \
          --output-dir {params.out_dir} \
          --genefinding-gff3 {params.features} \
          {input}"


rule deepbgc:
    input:
        "{outdir}/{strain}/pilon/iteration_4/pilon.fasta"
    output:
        "{outdir}/{strain}/deepbgc/README.txt"
    params:
        out_dir = "{outdir}/{strain}/deepbgc",
        score = config['deepbgc_score']
    conda:
        "/home/yura/anaconda3/envs/deepbgc.yaml"
    shell:
        "export DEEPBGC_DOWNLOADS_DIR=`pwd`/deepdata ; \
        mkdir `pwd`/deepdata && \
        deepbgc download ; \
        deepbgc pipeline \
          -o {params.out_dir} \
          -s {params.score} \
          {input}"


rule eggnog:
    input:
        "{outdir}/{strain}/prokka/{strain}.faa"
    output:
        "{outdir}/{strain}/eggnog/{strain}.emapper.annotations"
    params:
        data = config['eggnog']['data'],
        evalue = config['eggnog']['evalue'],
        pident = config['eggnog']['pident'],
        out_dir = "{outdir}/{strain}/eggnog",
        strain = lambda wildcards: wildcards.strain,
        tool = config['eggnog']['tool']
    conda:
        "/home/yura/anaconda3/envs/eggnog.yaml"
    threads: config['threads']
    shell:
        "mkdir {params.data} && \
         download_eggnog_data.py \
            --data_dir {params.data} \
           -P \
           -M \
           -H \
           -d 2 \
           -y \
           -f ; \
         emapper.py \
           --cpu {threads} \
           --override \
           --data_dir {params.data} \
           -i {input} \
           --itype proteins \
           --pident {params.pident} \
           --evalue {params.evalue} \
           --output {params.strain} \
           --output_dir {params.out_dir}"


#rule make_mmseqs2_db:
#    input:
#    output:
#        directory()
#    params:
#    threads: config['threads']
#    conda:
#        "/home/yura/anaconda3/"
#    shell:
#        "if [[ ! -z {output} ]]
#         then
#             <>
#         fi"

#rule mmseqs2:
#    input:
#        annot_file = "{out_dir}/{strain}/prokka/{strain}.faa",
#        annot_db = config['mmseqs2']['annot_db'] if os.path.exists(config['mmseqs2']['annot_db']) else make_mmseqs2_db.output
#    output:
#        "{out_dir}/{strain}/mmseqs"
#    params:
#    conda:
#    threads: config['threads']
#    shell:
#    

#rule busco:
#    input:
#        "{outdir}/{strain}/prokka/{strain}.faa"
#    output:
#        directory("{outdir}/{strain}/busco")
#    params:
#        evalue = config['busco']['evalue'],
#        lineage = config['busco']['lineage'],
#        strain = lambda wildcards: wildcards.strain
#    conda:
#        "/home/yura/anaconda3/envs/busco.yaml"
#    threads: config['threads']
#    shell:
#        "busco \
#          --force \
#          --out_path {output} \
#          -c {threads} \
#          -e {params.evalue} \
#          -i {input} \
#          -l {params.lineage} \
#          -m prot \
#          -o {params.strain}"
