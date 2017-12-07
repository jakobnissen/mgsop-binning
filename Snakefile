configfile: srcdir('config.yaml')

include: 'scripts/parseconfig.snakefile'

workdir: config['binning_workdir']

# Help with avoiding missing file errors if there's latency on network.
shell.prefix('sleep 10; ')

rule all:
    input:        
        # Gather contigs to one file
        #['contigassembler/onehalf.fna', 'contigassembler/otherhalf.fna']
        # Run contigAssembler
        #'contigassembler/merged.fna'
        # Discard small contigs
        #'contigs/contigs.fna'
        # BWA index contigs
        #'contigs/contigs.fna.bwt'
        # BWA mem experiments to contigs
        #expand('mapped/{experiment}.bam, experiment=config['experiments'])
        # Samtools sort
        #expand('mapped/{experiment}.sorted.bam', experiment=config['experiments'])
        # Metabat create depths file
        #'metabat/depths.txt'
        # Metabat create bins
        dynamic('metabat/bins/bin.{binid}.fna')
        
# Concats contigs
rule concat_contigs:
    input: expand(os.path.join(config['assembly_workdir'], 'assembly/{sample}.fna', sample=config['samples']))
    output:
        onehalf = temp('contigs/onehalf.fna'),
        otherhalf = temp('contigs/otherhalf.fna')
    params:
        walltime = 86400,
        ppn = 1,
        mem = '1gb'
    shell: 'cat {input} > {output.onehalf} && touch {output.otherhalf}'
         
# Run contigAssembler
rule contigassembler:
    input: rules.concat_contigs.output
    output: 'contigassembler/merged.fna'
    params:
        walltime = 864000,
        ppn = CORES,
        mem = '95gb',
        contigassembler = srcdir('contigAssembler.pl')
    shell: '{params.contigassembler} -i {input[0]} -j {input[1]} -o {output} -c {params.ppn}'

def filter_contigs_input(wildcards):
    if ASSEMBLE_CONTIGS:
        return rules.contigassembler.output
    else:
        return expand(os.path.join(config['assembly_workdir'], 'assembly/{sample}.fna', sample=config['samples']))

# Discard small contigs and rename them to ensure unique names
rule filter_contigs:
    input: filter_contigs_input
    output: protected('contigs/contigs.fna')
    threads: CORES
    params:
        walltime = 864000,
        ppn = 1,
        mem = '95gb',
        min_contig_length = MIN_CONTIG_LENGTH
    message:
        '\n    input: {input}\n'
        '    output: {output}\n'
        '    minimum gene length: {params.min_contig_length}\n'
        '    script: filtercontigs.py'
    script: 'scripts/filtercontigs.py'    

# BWA index contigs
rule index_contigs:
    input: rules.filter_contigs.output
    output:
        'contigs/contigs.fna.bwt',
        'contigs/contigs.fna.amb',
        'contigs/contigs.fna.ann',
        'contigs/contigs.fna.pac',
        'contigs/contigs.fna.sa'
    log:
        out = 'log/index_contigs/bwa.out.log',
        err = 'log/index_contigs/bwa.err.log',
    params:
        walltime = 864000,
        ppn = 1,
        mem = '15gb'
    shell:
        '{config[bwa_path]} index {input} 2> {log.err} > {log.out}'

# This allows map_to_contigs to map either single end or paired end reads.
def map_input(wildcards):
    if IS_PE[wildcards.experiment]:
        return [os.path.join(config['assembly_workdir'], 'trim/{}.fw.fastq.gz'.format(wildcards.experiment)),
                os.path.join(config['assembly_workdir'], 'trim/{}.rv.fastq.gz'.format(wildcards.experiment))]
    else:
        return os.path.join(config['assembly_workdir'], 'trim/{}.fastq.gz'.format(wildcards.experiment))

def readgroupstring(wildcards):
    for smp, experiments in config['samples'].items():
        if wildcards.experiment in experiments:
            sample = smp
            break
        
    template = '"@RG\\tID:{}\\tSM:{}\\tPL:{}"'
    return template.format(wildcards.experiment, sample, config['platform'])

rule map_to_contigs:
    input:
        contigfile = rules.filter_contigs.output,
        contigfilelist = rules.index_contigs.output,
        fastq = map_input
    output: temp('mapped/{experiment}.bam')
    threads: CORES
    log: 
        bwa = 'log/mapping/{experiment}.log',
        samtools = 'log/mapping/{experiment}.log'
    params:
        walltime = 864000,
        ppn = CORES,
        mem = '90gb',
        readgroup = readgroupstring
    shell:
        '{config[bwa_path]} mem -t {threads} {input.contigfile} {input.fastq} '
        '-R {params.readgroup} 2> {log.bwa} | '
        '{config[samtools_path]} view -F 2560 -b -h -q 20 --threads {threads} - > ' 
        '{output} 2> {log.samtools}'
        
########################## SAMTOOLS SORT ##################################
# Due to a bug in Samtools sort, the memory may far exceed (>260%) specified
# amount. This problem is worse for multi-thread sorting. If a submitted job
# exceeds given memory, it will be killed midway.
# Too little specified memory will make the program crash with sufficiently
# large files due to another bug with headers.

# Therefore, each process is specified to use 10 GB, single thread, but given
# 25 GB by the scheduler (and 7 cores to block all 25 GB).
###########################################################################

rule namesort:
    input: rules.map_to_contigs.output
    output: temp('mapped/{experiment}.namesorted.bam')
    threads: 4
    log: 'log/bam/{experiment}.namesort.log'
    params:
        walltime = 864000,
        ppn = 4,
        mem = '12gb'
    shell:
        '{config[samtools_path]} sort -n '
        '-m 10G --threads 4 -T {input}.tmp '
        '{input} -o {output} 2> {log}'

rule fixmate:
    input: rules.namesort.output
    output: temp('mapped/{experiment}.matefixed.bam')
    threads: 4
    log: 'log/bam/{experiment}.fixmate.log'
    params:
        walltime = 864000,
        ppn = 4,
        mem = '12gb'
    shell: '{config[samtools_path]} fixmate -m -@ 3 {input} {output} 2> {log}'

def possort_input(wildcards):
    if REMOVE_DUPLICATES:
        return rules.fixmate.output
    else:
        return rules.map_to_contigs.output

rule possort:
    input: possort_input
    output: 'mapped/{experiment}.possort.bam'
    log: 'log/bam/{experiment}.possort.log'
    params:
        walltime = 864000,
        ppn = 4,
        mem = '12gb'
    shell:
        '{config[samtools_path]} sort '
        '-m 10G --threads 4 -T {input}.tmp '
        '{input} -o {output} 2> {log}'

rule markdup:
    input: rules.possort.output
    output: temp('mapped/{experiment}.markdup.bam')
    threads: 4
    log: 'log/bam/{experiment}.markdup.log'
    params:
        walltime = 864000,
        ppn = 4,
        mem = '12gb'
    shell: '{config[samtools_path]} markdup -@ 3 {input} {output} 2> {log}'

rule rmpossorted:
    input:
        possort = rules.possort.output,
        markdup = rules.markdup.output
    output: temp(touch('flags/rmpossorted.{experiment}.done'))
    log: 'log/rmpossorted.log'
    params:
        walltime = 86400,
        ppn = 1,
        mem = '100 mb'
    shell: 'rm {input.possort} 2> {log}'

rule dupfilter:
    input:
        markdup = rules.markdup.output,
        rmpossortedflag = rules.rmpossorted.output
    output: 'mapped/{experiment}.filtered.bam'
    log: 'log/bam/{experiment}.filtered.bam'
    params:
        walltime = 864000,
        ppn = 1,
        mem = '4gb'
    shell: '{config[samtools_path]} view -F 3584 -b -h {input.markdup} > {output} 2> {log}'

def mergebam_input(wildcards):
    if REMOVE_DUPLICATES:
        return expand('mapped/{exp}.filtered.bam', exp=config['samples'][sample])
    else:
        return expand('mapped/{exp}.possort.bam', exp=config['samples'][sample])

rule mergebam:
    input: mergebam_input
    output: 'mapped/{sample}.merged.bam'
    threads: 4
    log: 'log/bam/{sample}.merged.log'
    params:
        walltime= 86400,
        ppn = 4,
        mem = '12gb'
    shell: '{config[samtools_path]} merge -@ 3 {output} {input}'

rule rmunmerged:
    input:
        merged = rules.mergebam.output,
        unmerged = mergebam_input
    output: temp(touch('flags/rmunmerged.{sample}.done'))
    log: 'log/rmpossorted.log'
    params:
        walltime = 86400,
        ppn = 1,
        mem = '100 mb'
    shell: 'rm {input.unmerged} 2> {log}'

def depths_input(wildcards):
    inputfiles = list()

    for sample, experiments in config['samples'].items():
        if len(experiments) == 1:
            if REMOVE_DUPLICATES:
                inputfiles.append('mapped/{}.filtered.bam'.format(experiments[0]))
            else:
                inputfiles.append('mapped/{}.possort.bam'.format(experiments[0]))

        else:
            inputfiles.append('mapped/{}.merged.bam'.format(sample))

    return inputfiles 

def depths_flags(wildcards):
    flags = list()
    for sample, experiments in config['samples'].items():
        if len(experiments) > 1: 
            flags.append('flags/rmunmerged.{}.done'.format(sample))

    return flags

rule metabat_create_depths_file:
    input:
        bamfiles = depths_input,
        flags = depths_flags
    output: 'metabat/depths.txt'
    log: 'log/metabat/jgi.log'
    params:
        walltime = 864000,
        ppn = 1,
        mem = '5gb',
        min_contig_length = MIN_CONTIG_LENGTH
    shell: '{config[jgi_path]} --minContigLength {params.min_contig_length} '
           '--outputDepth {output} {input.bamfiles} 2> {log}'

rule metabat_create_bins:
    input:
        depths = rules.metabat_create_depths_file.output,
        contigs = rules.filter_contigs.output,
    output: dynamic('metabat/bins/bin.{binid}.fna')
    log: 'log/metabat/metabat.log'
    threads: CORES
    params:
        min_contig_length = MIN_CONTIG_LENGTH,
        walltime = 864000,
        ppn = CORES,
        mem = '95gb'
    shell:
        '{config[metabat_path]} -i {input.contigs} -a {input.depths} -m {params.min_contig_length} '
        '-o {output} -t {params.ppn}'

def cleanup(directory):
    # Move the numerous "snakejob" files that Computerome can generate.
    files = [f for f in os.listdir(directory) if f.startswith('snakejob.')]
    snakejobdir = os.path.join(directory, 'snakejobs')

    try:
        os.mkdir(snakejobdir)
    except FileExistsError:
        if os.path.isdir(snakejobdir):
            pass
        else:
            print('A non-directory file called "snakejobs" already exists. Cannot clean.')
            return

    for file in files:
        source = os.path.join(directory, file)
        destination = os.path.join(snakejobdir, file)

        os.rename(source, destination)

onsuccess:
    cleanup(os.getcwd())

onerror:
    cleanup(os.getcwd())