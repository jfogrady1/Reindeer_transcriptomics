
rule all:
    input:
        expand('/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Fastqc/{individual}_{N}_fastqc.html', individual = config["species"], N = (1,2)),
         expand('/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/trimming/{individual}_{N}_trimmed_fastqc.html', individual = config["species"], N = (1,2)),
        '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Fastqc/multiqc_report.html',
         directory("/home/workspace/jogrady/Reindeer_transcriptomics/data/RNA_seq/reindeer_star_genome/"),
         directory("/home/workspace/jogrady/Reindeer_transcriptomics/data/RNA_seq/human_star_genome/"),
         directory("/home/workspace/jogrady/Reindeer_transcriptomics/data/RNA_seq/bovine_star_genome/"),
         expand('/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/trimming/{individual}_{N}_trimmed.fastq.gz', individual = config["species"], N = (1,2)),
         expand('/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/{individual}_Log.out', individual=config["species"]),
         expand('/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/{group}_counts.txt', group = ["Reindeer", "Human", "Bovine"]),
         expand('/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/{group}_count_matrix_clean.txt',group = ["Reindeer", "Human", "Bovine"]),
         plot2 = "/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Globin_mRNA_Bovine_V_Reindeer_V_Human.pdf"


rule fastqc:
    input:
        reads = lambda wildcards: expand(f'/home/workspace/jogrady/Reindeer_transcriptomics/data/RNA_seq/RAW/{config["species"][wildcards.individual]}_{{N}}.fq.gz', N = (1,2))
    output:
        reads = expand('/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Fastqc/{{individual}}_{N}_fastqc.zip', N = (1,2)),
        html = expand('/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Fastqc/{{individual}}_{N}_fastqc.html', N = (1,2))
    threads:
        20
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o /home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Fastqc/'

rule multiqc:
    input:
        reads = expand('/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Fastqc/{individual}_{N}_fastqc.zip', individual = config['species'], N = (1,2))
    output:
        report='/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Fastqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o //home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Fastqc
        """


rule reindeer_index:
    input:
        fa="/home/workspace/jogrady/Reindeer_transcriptomics/data/a3s.fi/rta2_ref/PO1266.assembly.RepeatMasked.fasta",
        gtf= "/home/workspace/jogrady/Reindeer_transcriptomics/data/a3s.fi/rta2_ref/Reindeer_inhouse.gtf"
    output:
        STAR_dir = directory("/home/workspace/jogrady/Reindeer_transcriptomics/data/RNA_seq/reindeer_star_genome/")  # Please provide the output files names from this step.
    threads: 20
    shell:
        '''
        mkdir /home/workspace/jogrady/Reindeer_transcriptomics/data/RNA_seq/reindeer_star_genome/ # STAR cannot make directory
        STAR-2.7.1a --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir /home/workspace/jogrady/Reindeer_transcriptomics/data/RNA_seq/reindeer_star_genome/ \
        --genomeFastaFiles {input.fa} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 100 \
        --limitGenomeGenerateRAM 96276080682
        ''' 

rule human_index:
    input:
        fa="/home/workspace/jogrady/Reindeer_transcriptomics/data/human_data/GCA_000001405.15_GRCh38_full_analysis_set.fna",
        gtf= "/home/workspace/jogrady/Reindeer_transcriptomics/data/human_data/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf"
    output:
        STAR_dir = directory("/home/workspace/jogrady/Reindeer_transcriptomics/data/RNA_seq/human_star_genome/")  # Please provide the output files names from this step.
    threads: 20
    shell:
        '''
        mkdir /home/workspace/jogrady/Reindeer_transcriptomics/data/RNA_seq/human_star_genome/ # STAR cannot make directory
        STAR-2.7.1a --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir /home/workspace/jogrady/Reindeer_transcriptomics/data/RNA_seq/human_star_genome/ \
        --genomeFastaFiles {input.fa} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 74 \
        ''' 

rule bovine_index:
    input:
        fa="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        gtf="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf" 
    output:
        STAR_dir = directory("/home/workspace/jogrady/Reindeer_transcriptomics/data/RNA_seq/bovine_star_genome/")  # Please provide the output files names from this step.
    threads: 20
    shell:
        '''
        mkdir /home/workspace/jogrady/Reindeer_transcriptomics/data/RNA_seq/bovine_star_genome/ # STAR cannot make directory
        STAR-2.7.1a --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir /home/workspace/jogrady/Reindeer_transcriptomics/data/RNA_seq/bovine_star_genome/ \
        --genomeFastaFiles {input.fa} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 149
        ''' 


rule trimming:
    input:
        fastp = "/home/workspace/jogrady/Reindeer_transcriptomics/data/fastp",
        reads = lambda wildcards: expand(f'/home/workspace/jogrady/Reindeer_transcriptomics/data/RNA_seq/RAW/{config["species"][wildcards.individual]}_{{N}}.fq.gz', N=(1,2))

    output:
        trimmed_reads=expand('/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/trimming/{{individual}}_{N}_trimmed.fastq.gz', N = (1,2))
    
    params:
        prefix = "{individual}"
    shell:
        """
        {input.fastp} --thread 10 --in1 {input.reads[0]} --in2 {input.reads[1]} -o {output.trimmed_reads[0]} -O {output.trimmed_reads[1]} --thread 10 \
        --cut_front 1 --cut_tail 1 --length_required 25 --cut_mean_quality 30 --report_title {params.prefix} --json /home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/trimming/{params.prefix}.json --html /home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/trimming/{params.prefix}.html

        #mv *.html /home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/trimming/
        #mv *.json /home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/trimming/
        """

rule fastqc_trimmed:
    input:
        reads = lambda wildcards: expand(f'/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/trimming/{config["species"][wildcards.individual]}_{{N}}_trimmed.fastq.gz', N = (1,2))
    output:
        reads = expand('/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/trimming/{{individual}}_{N}_trimmed_fastqc.zip', N = (1,2)),
        html = expand('/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/trimming/{{individual}}_{N}_trimmed_fastqc.html', N = (1,2))
    threads:
        20
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o /home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/trimming/'


rule Alignment_bovine:
    input:
        genome = rules.bovine_index.output.STAR_dir,
        reads = lambda wildcards: expand(f'/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/trimming/{config["bovine"][wildcards.individual]}_{{N}}_trimmed.fastq.gz', N=(1,2)),
    params:
        prefix = lambda wildcards: f'{config["bovine"][wildcards.individual]}'
    output:
        aligned = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/{individual}_Aligned.sortedByCoord.out.bam',
        finallog = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/{individual}_Log.final.out',
        interlog = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/{individual}_Log.progress.out',
        initiallog = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/{individual}_Log.out'
    threads: 10
    shell:
        '''
        echo {threads}
        STAR-2.7.1a  --genomeLoad LoadAndKeep --genomeDir {input.genome} --runThreadN {threads} \
        --readFilesIn {input.reads[0]} {input.reads[1]} --readFilesCommand gunzip -c \
        --outFileNamePrefix /home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/{params.prefix}_ --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 18000000000
        '''

rule Alignment_reindeer:
    input:
        genome = rules.reindeer_index.output.STAR_dir,
        reads = lambda wildcards: expand(f'/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/trimming/{config["reindeer"][wildcards.individual]}_{{N}}_trimmed.fastq.gz', N=(1,2)),
    params:
        prefix = lambda wildcards: f'{config["reindeer"][wildcards.individual]}'
    output:
        aligned = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/{individual}_Aligned.sortedByCoord.out.bam',
        finallog = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/{individual}_Log.final.out',
        interlog = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/{individual}_Log.progress.out',
        initiallog = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/{individual}_Log.out'
    threads: 8
    shell:
        '''
        STAR-2.7.1a  --genomeLoad LoadAndKeep --genomeDir {input.genome} --runThreadN {threads} \
        --readFilesIn {input.reads[0]} {input.reads[1]} --readFilesCommand gunzip -c \
        --outFileNamePrefix /home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/{params.prefix}_ --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 18000000000
        '''

rule Alignment_human:
    input:
        genome = rules.human_index.output.STAR_dir,
        reads = lambda wildcards: expand(f'/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/trimming/{config["human"][wildcards.individual]}_{{N}}_trimmed.fastq.gz', N=(1,2)),
    params:
        prefix = lambda wildcards: f'{config["human"][wildcards.individual]}'
    output:
        aligned = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/{individual}_Aligned.sortedByCoord.out.bam',
        finallog = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/{individual}_Log.final.out',
        interlog = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/{individual}_Log.progress.out',
        initiallog = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/{individual}_Log.out'
    threads: 10
    shell:
        '''
        STAR-2.7.1a  --genomeLoad LoadAndKeep --genomeDir {input.genome} --runThreadN {threads} \
        --readFilesIn {input.reads[0]} {input.reads[1]} --readFilesCommand gunzip -c \
        --outFileNamePrefix /home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/{params.prefix}_ --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 18000000000
        '''

rule quantification_reindeer:
    input:
        bam = expand('/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/{sample}_Aligned.sortedByCoord.out.bam', sample = config["reindeer"]),
        annotation="/home/workspace/jogrady/Reindeer_transcriptomics/data/a3s.fi/rta2_ref/Reindeer_inhouse.gtf"
    output:
        count_matrix = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Reindeer_counts.txt'
    threads: 15
    shell:
        '''
        # use new version of feature counts
        /home/workspace/jogrady/Reindeer_transcriptomics/bin/subread-2.0.6-Linux-x86_64/bin/featureCounts -a {input.annotation} -o {output.count_matrix} {input.bam} -B -p -C -R BAM -T {threads} -s 0 -t gene -g gene_id
        '''

rule quantification_bovine:
    input:
        bam = expand('/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/{sample}_Aligned.sortedByCoord.out.bam', sample = config["bovine"]),
        annotation="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf"
    output:
        count_matrix = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Bovine_counts.txt'
    threads: 15
    shell:
        '''
        # use new version of feature counts
        /home/workspace/jogrady/Reindeer_transcriptomics/bin/subread-2.0.6-Linux-x86_64/bin/featureCounts -a {input.annotation} -o {output.count_matrix} {input.bam} -B -p -C -R BAM -T {threads} -s 0 -t gene -g gene_id
        '''

rule quantification_human:
    input:
        bam = expand('/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/{sample}_Aligned.sortedByCoord.out.bam', sample = config["human"]),
        annotation="/home/workspace/jogrady/Reindeer_transcriptomics/data/human_data/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf"
    output:
        count_matrix = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Human_counts.txt'
    threads: 15
    shell:
        '''
        # use new version of feature counts
        /home/workspace/jogrady/Reindeer_transcriptomics/bin/subread-2.0.6-Linux-x86_64/bin/featureCounts -a {input.annotation} -o {output.count_matrix} {input.bam} -B -p -C -R BAM -T {threads} -s 0 -t gene -g gene_id
        '''


rule cleanup_Reindeer:
    input:
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/RNA_seq/FC_cleanup.sh",
        count_matrix_in = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Reindeer_counts.txt',
        final_script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/RNA_seq/FC_cleanup.py",
    output:
        #temporary = "/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Reindeer_counts_temp.txt",
        count_matrix = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Reindeer_counts2.txt',
        count_matrix_2 = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Reindeer_counts_clean.txt',
        cleaned = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Reindeer_count_matrix_clean.txt',
        temporary_1 = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Reindeer_counts_copy.txt'
    shell:
        '''
        cat {input.count_matrix_in} > {output.temporary_1} 
        tail -n+2 {output.temporary_1} > {output.count_matrix}
        cut -f 1,7-11 {output.count_matrix} > {output.count_matrix_2} # number of samples, get rid of fields we do not want
        sed -i 's#/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/##g' {output.count_matrix_2}
        sed -i 's/'"_Aligned\.sortedByCoord\.out\.bam"'//g' {output.count_matrix_2} 
        cat {output.count_matrix_2} > {output.cleaned} 
        '''

rule cleanup_Bovine:
    input:
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/RNA_seq/FC_cleanup.sh",
        count_matrix_in = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Bovine_counts.txt',
        final_script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/RNA_seq/FC_cleanup.py",
    output:
        #temporary = "/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Reindeer_counts_temp.txt",
        count_matrix = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Bovine_counts2.txt',
        count_matrix_2 = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Bovine_counts_clean.txt',
        cleaned = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Bovine_count_matrix_clean.txt',
        temporary_1 = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Bovine_counts_copy.txt'
    shell:
        '''
        cat {input.count_matrix_in} > {output.temporary_1} 
        tail -n+2 {output.temporary_1} > {output.count_matrix}
        cut -f 1,7-11 {output.count_matrix} > {output.count_matrix_2} # number of samples, get rid of fields we do not want
        sed -i 's#/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/##g' {output.count_matrix_2}
        sed -i 's/'"_Aligned\.sortedByCoord\.out\.bam"'//g' {output.count_matrix_2} 
        cat {output.count_matrix_2} > {output.cleaned} 
        '''



rule cleanup_Human:
    input:
        script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/RNA_seq/FC_cleanup.sh",
        count_matrix_in = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Human_counts.txt',
        final_script = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/RNA_seq/FC_cleanup.py",
    output:
        #temporary = "/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Reindeer_counts_temp.txt",
        count_matrix = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Human_counts2.txt',
        count_matrix_2 = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Human_counts_clean.txt',
        cleaned = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Human_count_matrix_clean.txt',
        temporary_1 = '/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Human_counts_copy.txt'
    shell:
        '''
        cat {input.count_matrix_in} > {output.temporary_1} 
        tail -n+2 {output.temporary_1} > {output.count_matrix}
        cut -f 1,7-11 {output.count_matrix} > {output.count_matrix_2} # number of samples, get rid of fields we do not want
        sed -i 's#/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Alignment/##g' {output.count_matrix_2}
        sed -i 's/'"_Aligned\.sortedByCoord\.out\.bam"'//g' {output.count_matrix_2} 
        cat {output.count_matrix_2} > {output.cleaned} 
        '''
        
rule TPM_visualisation:
    input:
        bovine_counts = "/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Bovine_count_matrix_clean.txt",
        reindeer_counts = "/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Reindeer_count_matrix_clean.txt",
        human_counts = "/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Quantification/Human_count_matrix_clean.txt",
        bovine_annotation = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf",
        reindeer_annotation = "/home/workspace/jogrady/Reindeer_transcriptomics/data/a3s.fi/rta2_ref/Reindeer_inhouse.gtf",
        human_annotation = "/home/workspace/jogrady/Reindeer_transcriptomics/data/human_data/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf",
        script = "/home/workspace/jogrady/Reindeer_transcriptomics/bin/TPM_Globin_Quant.R"
    output:
        plot1 = "/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Globin_mRNA_Bovine_V_Reindeer.pdf",
        plot2 = "/home/workspace/jogrady/Reindeer_transcriptomics/results/RNA-seq/Globin_mRNA_Bovine_V_Reindeer_V_Human.pdf"

    shell:
        '''
        Rscript {input.script} {input.bovine_counts} {input.reindeer_counts} {input.human_counts} {input.bovine_annotation} {input.reindeer_annotation} {input.human_annotation} {output.plot1} {output.plot2} 
        '''