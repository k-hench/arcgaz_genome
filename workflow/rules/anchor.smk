'''
# execute
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
    --configfile workflow/config.yml \
    -R anchor
'''

rule anchor:
    input:
      # expand( 'results/genomes/arcgaz_anchored_{haplotype}.fa.gz', haplotype = HAPLOTYPES)
      expand( 'results/anchoring/arcgaz_dt_{haplotype}_hardmasked/anchored_arcgaz_dt_{haplotype}.fasta.gz', haplotype = HAPLOTYPES )

rule create_anchoring_beds:
    input:
    output:
      expand('results/anchoring/arcgaz_dt_{haplotype}_hardmasked/anchored_arcgaz_dt_{haplotype}.bed', haplotype = HAPLOTYPES)
    log:
      'logs/create_anchoring_bed.log'
    shell: 
      '''
      Rscript R/create_anchoring_beds.R &> {log}
      '''

rule run_allmaps:
    input:
      genome = "data/genomes/arcgaz_dt_{haplotype}_masked.fa.gz",
      bed = 'results/anchoring/arcgaz_dt_{haplotype}_hardmasked/anchored_arcgaz_dt_{haplotype}.bed'
    output:
      anchord = 'results/anchoring/arcgaz_dt_{haplotype}_hardmasked/anchored_arcgaz_dt_{haplotype}.fasta'
    conda: "allmaps"
    log:
      'logs/allmaps/allmaps_{haplotype}.log'
    shell:
      '''
      cd results/anchoring/arcgaz_dt_{wildcards.haplotype}_hardmasked
      pwd &> ../../../{log}
      python -m jcvi.assembly.allmaps \
        path \
        anchored_arcgaz_dt_{wildcards.haplotype}.bed \
        ../../../{input.genome} &> ../../../{log}
      '''

rule zip_anchored:
    input:
      anchord = "results/anchoring/arcgaz_dt_{haplotype}_hardmasked/anchored_arcgaz_dt_{haplotype}.fasta"
    output:
      ziped = "results/anchoring/arcgaz_dt_{haplotype}_hardmasked/anchored_arcgaz_dt_{haplotype}.fasta.gz"
    conda: "popgen_basics"
    shell:
      '''
      bgzip {input.anchord}
      samtools faidx {output.ziped}
      '''