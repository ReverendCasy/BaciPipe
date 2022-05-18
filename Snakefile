from typing import Any, List
import os

configfile: "config_2022.yaml"
ITER_TOOLS = ['bwa', 'minimap', 'pilon', 'racon']
ITERS = list(range(1, 5))

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
    for tool in ITER_TOOLS:
        it: int
        for it in ITERS:
            os.makedirs(os.path.join(config['out_dir'], f'{strain}', tool, f'iteration_{it}'), exist_ok=True)


rule all:
    input:
        expand("{outdir}/{strain}/sortpred/predicted_sortase.csv",
               outdir=config['out_dir'],
               strain=config['strains'])
    output: touch('.status')


rule assemble_flye:
    input:
        reads = lambda wildcards: config['strains'][wildcards.strain]
    output:
        "{outdir}/{strain}/flye/assembly.fasta"
    params:
        outdir = "{outdir}/{strain}/flye",
        polish_rounds = config['polish_rounds']
    threads: 72
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
        n = "[1-4]"
    threads: 72
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
    threads: 72
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
            n = "[1-4]"
        threads: 72
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
               -@ 72 \
               -o {output.bwa}.bam \
               {output.bwa}.bam && \
             samtools index \
               {output.bwa}.bam && \
             pilon \
               --genome {input.ref} \
               --frags {output.bwa}.bam \
               --outdir {params.pilon_dir}"


rule prokka:
    input:
        contigs = "{outdir}/{strain}/pilon/iteration_4/pilon.fasta",
        model = config['prodigal_model'],
        ref_prot = config['reference_proteins']
    output:
        "{outdir}/{strain}/prokka/{strain}.faa"
    params:
        out_dir = "{outdir}/{strain}/prokka",
        prefix = lambda wildcards: wildcards.strain
    conda:
        "/home/yura/anaconda3/envs/prokka.yaml"
    threads: 72
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


rule idops:
    input:
        "{outdir}/{strain}/prokka/{strain}.faa"
    output:
        "{outdir}/{strain}/idops/idops_scan.tsv"
    params:
         out_dir = "{outdir}/{strain}/idops"
    threads: 72
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
    threads: 72
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
        72
    shell:
        "./lengthFilter.py \
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
        model1 = config['sp_model1'],
        model2 = config['sp_model2']
    conda:
        "/home/yura/anaconda3/envs/R_bacillus.yaml"
    shell:
        "Rscript {params.script} \
            --layer1 {params.model1} \
            --layer2 {params.model2} \
            -i {input} \
            -o {params.out_dir}"
