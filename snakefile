
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()
GS_PREFIX = "me-inbre-rnaseq-pipelinev2"

configfile: "config.yaml"
SAMPLES=config["samples"]
print(SAMPLES)

rule all:
    input:
        #comment out any file here that you don't want the pipeline to run to
        #raw fastqc + multiqc
        #"qc/raw/multiqc/multiqc_report.html",
        #trimmomatic
        expand(["data/trimmed/{sample}_trimmed_1.fastq"],sample=SAMPLES),
        #trimmed fastqc + multiqc
        #expand(["qc/trimmed/fastqc/{sample}_trimmed_1_fastqc.html"],sample=SAMPLES), 
        #"qc/trimmed/multiqc/multiqc_report.html", 
        #bwa mem + samtools sort + samtools index
        expand(["data/aligned/{sample}.bam"],sample=SAMPLES), 
        expand(["data/sorted/{sample}_sorted.bam"],sample=SAMPLES),
        expand(["data/sorted/{sample}_sorted.bam.bai"],sample=SAMPLES), 
        #samtools flagstat
        expand(["data/sorted/{sample}_stats.txt"],sample=SAMPLES), 
        #build salmon index
        #"data/reference/transcriptome_index/refseq.bin", 
        #run salmon quant
        #expand(["data/salmon/{sample}/quant.sf"],sample=SAMPLES)

rule fastqc_raw:
    input:
       expand(["data/raw_fastqSub/{sample}_1.fastq",
                 "data/raw_fastqSub/{sample}_2.fastq"], sample=SAMPLES)
    output:
        expand(["qc/raw/fastqc/{sample}_1_fastqc.html", "qc/raw/fastqc/{sample}_2_fastqc.html"], sample=SAMPLES)
    conda:
        "envs/fastqc.yaml"
    params:
        outdir="qc/raw/fastqc"
    threads: 4
    shell:
        "fastqc -o {params.outdir} -t {threads} {input}"

rule multiqc_raw: 
    input: 
        expand(["qc/raw/fastqc/{sample}_1_fastqc.html",
                "qc/raw/fastqc/{sample}_2_fastqc.html"], sample=SAMPLES)
    output:
        "qc/raw/multiqc/multiqc_report.html"
    conda: 
        "envs/fastqc.yaml"
    params:
        indir="qc/raw/fastqc"
    shell: 
        "multiqc {params.indir} -o qc/raw/multiqc -f"


rule trimmomatic_pe_fq:
    input:
        r1="data/raw_fastqSub/{sample}_1.fastq",
        r2="data/raw_fastqSub/{sample}_2.fastq"
    output:
        r1="data/trimmed/{sample}_trimmed_1.fastq",
        r2="data/trimmed/{sample}_trimmed_2.fastq",
        # reads where trimming entirely removed the mate
        r1_unpaired="data/trimmed/{sample}_trimmed_1.unpaired.fastq",
        r2_unpaired="data/trimmed/{sample}_trimmed_2.unpaired.fastq"
    params:
        # list of trimmers (see manual)
        trimmer=["ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36"],
        # optional parameters
        extra="",
        compression_level="-9"
    threads: 12
    conda:
        "envs/trimmomatic.yaml"
    shell:
        "trimmomatic PE -threads {threads} "
        "{input.r1} {input.r2} "
        "{output.r1} {output.r1_unpaired} "
        "{output.r2} {output.r2_unpaired} "
        "{params.trimmer}"

rule fastqc_trimmed:
    input:
       expand(["data/trimmed/{sample}_trimmed_1.fastq",
                 "data/trimmed/{sample}_trimmed_2.fastq"], sample=SAMPLES)
    output:
        expand(["qc/trimmed/fastqc/{sample}_trimmed_1_fastqc.html", "qc/trimmed/fastqc/{sample}_trimmed_2_fastqc.html"], sample=SAMPLES)
    conda:
        "envs/fastqc.yaml"
    params:
        outdir="qc/trimmed/fastqc"
    threads: 4
    shell:
        "fastqc -o {params.outdir} -t {threads} {input}"

rule multiqc_trimmed:
    input:
        expand(["data/trimmed/{sample}_trimmed_1.fastq",
                 "data/trimmed/{sample}_trimmed_2.fastq"], sample=SAMPLES)
    output:
        "qc/trimmed/multiqc/multiqc_report.html"
    conda: 
        "envs/fastqc.yaml"
    params:
        indir="qc/trimmed/fastqc"
    shell:
        "multiqc {params.indir} -o qc/trimmed/multiqc -f"
        
rule bwa_mem:
    input:
        reads=["data/trimmed/{sample}_trimmed_1.fastq", "data/trimmed/{sample}_trimmed_2.fastq"]
    output:
        "data/aligned/{sample}.bam"
    conda: 
        "envs/bwa.yaml"
    params:
        index="data/reference/M_chelonae_NZ_CP007220.fasta",
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'", #make sure this gets changed if there are multiple fastq runs per sample
        sort="none",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        sort_extra=""            # Extra args for samtools/picard.
    threads: 12
    wrapper:
        "master/bio/bwa/mem/"

rule samtools_sort:
    input:
        "data/aligned/{sample}.bam"
    output:
        "data/sorted/{sample}_sorted.bam"
    conda:
        "envs/samtools.yaml"
    threads: 8
    shell:
        "samtools sort --threads {threads} -o {output} {input}"

rule samtools_index:
    input:
        "data/sorted/{sample}_sorted.bam"
    output:
        "data/sorted/{sample}_sorted.bam.bai"
    params:
        "" # optional params string
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools index {params} {input} {output}"

rule samtools_flagstat:
    input:
        "data/sorted/{sample}_sorted.bam"
    output:
        "data/sorted/{sample}_stats.txt"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools flagstat {input} > {output}"

rule salmon_index:
    input:
        "data/reference/transcript_SRR13349122.fasta"
    output:
       directory("data/reference/transcriptome_index")
    threads: 4
    conda:
        "envs/salmon.yaml"
    params:
        # optional parameters
        decoys="data/reference/decoys.txt", 
    shell:
        "salmon index -t {input} -i {output} -d {params.decoys} --threads {threads}"

rule salmon_quant_reads:
    input:
        # If you have multiple fastq files for a single sample (e.g. technical replicates)
        # use a list for r1 and r2.
        r1 = "data/trimmed/{sample}_trimmed_1.fastq",
        r2 = "data/trimmed/{sample}_trimmed_2.fastq",
        index = "data/reference/transcriptome_index"
    output:
        quant = 'data/salmon/{sample}/quant.sf',
        lib = 'data/salmon/{sample}/lib_format_counts.json'
    params:
        # optional parameters
        libtype ="A",
        #zip_ext = bz2 # req'd for bz2 files ('bz2'); optional for gz files('gz')
        extra=""
    threads: 2
    conda:
        "envs/salmon.yaml"
    wrapper:
       "master/bio/salmon/quant"

rule trinity:
    input:
        left=["data/trimmed/{sample}_trimmed_1.fastq"],
        right=["data/trimmed/{sample}_trimmed_1.fastq"] 
    output:
        "reference/trinity_out_dir/{sample}.fasta"
    params:
        extra=""
    threads: 32
    wrapper:
        "master/bio/trinity"


