"""
snakemake --configfile workflow/config.yml --rerun-triggers mtime -n -R cov_by_fam

snakemake --jobs 50 \
    --configfile workflow/config.yml \
    --latency-wait 30 \
    -p \
    --default-resources mem_mb=51200 threads=1 \
    --use-singularity \
    --singularity-args "--bind $CDATA" \
    --use-conda \
    --rerun-triggers mtime \
    --cluster '
      sbatch_delay \
        --export=ALL \
        -n {threads} \
        -e logs/{name}.{jobid}.err \
        -o logs/{name}.{jobid}.out \
        --mem={resources.mem_mb}' \
        --jn job_c.{name}.{jobid}.sh \
        -R cov_by_fam
"""


wildcard_constraints:
    fam = "[^_]*",
    mscaf = "[^-]*"
    


fams = { "pho": ["halgry", "lepwed", "mirang", "mirleo", "neosch", "odoros", "phovit"],
         "ota": ["calurs", "eumjub", "arcgaz", "zalcal"]}


fams_c = { "pho": "halgry,lepwed,mirang,mirleo,neosch,odoros,phovit",
           "ota": "calurs,eumjub,arcgaz,zalcal"}

rule cov_by_fam:
    input: expand( "results/neutral_tree/cov/fam/{fam}-{mscaf}.collapsed.bed.gz", fam = ["pho", "ota"], mscaf = SCFS )


rule alignment_coverage_fam:
    input:
      hal = 'results/cactus/{name}.hal'.format(name = P_NAME)
    output:
      wig = "results/neutral_tree/wig/fam/{fam}-{mscaf}.wig.gz"
    params:
      prefix = "results/neutral_tree/wig/fam/{fam}-{mscaf}.wig",
      specs = lambda wc: "{f}".format(f = [wc.fam] )
    container: c_cactus
    shell:
      """
      halAlignmentDepth \
        {input.hal} \
        {REF_SPEC} \
        --targetGenomes {params.specs} \
        --refSequence {wildcards.mscaf} \
        --outWiggle {params.prefix}
      gzip {params.prefix}
      """

# ultimately we want a bed file as mask, so we convert the wig to bed format
rule wig_to_bed_fam:
    input:
      wig = "results/neutral_tree/wig/fam/{fam}-{mscaf}.wig.gz"
    output:
      bed = "results/neutral_tree/cov/raw/fam/{fam}-{mscaf}.bed.gz" 
    conda: "bedops"
    shell:
      """
      zcat {input.wig} | wig2bed | gzip > {output.bed}
      """

# unfortunetely, the original bed is single bp elements,
# so we collapse them into chunks of equal coverage 
rule collapse_cov_bed_fam:
    input:
      bed = "results/neutral_tree/cov/raw/fam/{fam}-{mscaf}.bed.gz"
    output:
      bed = "results/neutral_tree/cov/fam/{fam}-{mscaf}.collapsed.bed.gz"
    log: "logs/collapse_cov_bed_{fam}-{mscaf}.log"
    container: c_conda
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/collapse_bed_coverage.R {input.bed} {output.bed} &>> {log}
      """