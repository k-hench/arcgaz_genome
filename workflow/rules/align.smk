'''
# execute
snakemake \
    --use-conda \
    --cores 1 \
    --configfile workflow/config.yml

# with singularity
snakemake \
    --use-conda \
    --use-singularity \
    --singularity-args "-B <path_to_shared_directory>"\
    --cores 1 \
    --configfile workflow/config.yml

# on CeBiTec
snakemake \
  --jobs 10 \
  --latency-wait 30 \
  --use-conda -p \
  --configfile workflow/config.yml \
  --default-resources mem_mb=25600 \
  --cluster '
     qsub \
     -V -cwd \
     -P fair_share \
     -l idle=1 \
     -l si_flag=1 \
     -pe multislot {threads} \
     -l vf={resources.mem_mb}' \
  -R align

# dryrun
snakemake \
    -n \
    --configfile workflow/config.yml
'''
container: "docker://khench/msa_envs:v0.1"

rule align:
    input:
      expand( 'results/psl/{species}-18.psl.gz', species = G_QUERY )
      #"img/alignment.svg"

rule lastdb_index:
    input:
      fasta = 'data/genomes/{ref}.fa.gz'
    output:
      temp('results/genome/{ref}lastdb_index.done')
    params:
      indexBase = 'data/genomes/{ref}',
      refSizeFile = 'results/genome/{ref}.size',
    log:
      'logs/db_idx/{ref}_lastdbIndex.log'
    conda:
      'msa_align'
    shell:
      """
      faSize -detailed {input.fasta} > {params.refSizeFile} 2>{log} && lastdb -R 10 -u YASS -c {params.indexBase} {input.fasta} 2>{log} && touch {output} 2>{log}
      """

rule build_index:
    input:
      str(rules.lastdb_index.output).format( ref = G_REF ),
      fastaFile="data/genomes/{species}.fa.gz"
    output:
      "results/genome/{species}.size"
    params:
      indexBase='data/genomes/{ref}'.format( ref = G_REF ),
      speciesSizeFile='results/genome/{species}.size',
      refNibDir='results/genome/nib',
      refFastaFile="data/genomes/{ref}.fa.gz".format( ref = G_REF ),
      refNib2Bit='results/genome/nib/{ref}.2bit'.format( ref = G_REF ),
    log:
      'logs/db_idx/{species}_index.log'
    conda:
      'msa_ucsc'
    threads: 1
    shell:
      """
      mkdir -p {params.refNibDir} && \
      faToTwoBit {params.refFastaFile} {params.refNib2Bit} && \
      faSize -detailed {input.fastaFile} > {output}
      """

rule align_single_last:
    input:
      str(rules.lastdb_index.output).format( ref = G_REF ),
      fastaFile = "data/genomes/{species}.fa.gz",
      speciesSizeFile = 'results/genome/{species}.size',
    output:
      maf = 'results/maf/{species}.maf.gz'
    params:
      indexBase = 'data/genomes/{ref}'.format( ref = G_REF ),
      lastParams = config[ 'lastParams' ],
      mafBase = 'results/maf/{species}.maf'
    log:
      'logs/align/{species}_align.log'
    conda:
      'msa_align'
    threads: 1
    shell:
      """
      lastal {params.lastParams} {params.indexBase} {input.fastaFile} 2>{log} 1>{params.mafBase}
      gzip {params.mafBase}
      """

rule maf_to_psl:
    input:
      maf = 'results/maf/{species}.maf.gz'
    output:
      psl = 'results/psl/{species}.psl.gz'
    params:
      pslBase = 'results/psl/{species}.psl'
    log:
      'logs/psl/{species}_psl.log'
    conda:
      'msa_align'
    threads: 1
    shell:
      """
      zcat {input.maf} | maf-convert psl 2>{log} 1>{params.pslBase}
      gzip {params.pslBase}
      """

rule slim_psl:
    input: 'results/psl/{species}.psl.gz'
    output: 'results/psl/{species}-18.psl.gz'
    shell:
      """
      zcat {input} | cut -f 1-18  | gzip > {output}
      """


rule plot_alignments:
    input:
      expand( 'results/psl/{species}-18.psl.gz', species = G_QUERY )
    output:
      "img/alignment.svg"
    log:
      'logs/plots/alignment_plot.log'
    shell:
     """
     # Rscript R/plot_alignments.R  2> {log} 1> {log}
     echo 'dummy' > img/alignment.svg
     """