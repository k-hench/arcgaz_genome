"""
snakemake --configfile workflow/config.yml --rerun-triggers mtime -n -R convert_hal

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
    sbatch \
      --export ALL \
      -n {threads} \
      -e logs/{name}.{jobid}.err \
      -o logs/{name}.{jobid}.out \
      --mem={resources.mem_mb}' \
      --jn job_c.{name}.{jobid}.sh \
      -R convert_hal
"""

REF_SPEC = "arcgaz"
TIP_SPECS = "calurs,eumjub,halgry,lepwed,mirang,mirleo,neosch,odoros,phovit,zalcal"
MSCAFS = [ "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "x"]

rule convert_hal:
    input: 
      maf = expand("results/pinniped/maf/{name}_{mscaf}.maf", name = P_NAME, mscaf = MSCAFS),
      snps = expand("results/pinniped/snps/{name}_{mscaf}.tsv", name = P_NAME, mscaf = MSCAFS),

rule hal_to_maf:
    input:
      hal = 'results/cactus/{name}.hal'.format(name = P_NAME)
    output:
      maf = "results/pinniped/maf/{name}_{mscaf}.maf"
    params:
      js = "results/cactus/scratch/pinniped_set/tmp/js_{name}_{mscaf}"
    container: c_cactus
    shell:
      """
      mkdir -p {params.js}
      cactus-hal2maf \
        {params.js} \
        {input.hal} \
        {output.maf} \
        --refGenome {REF_SPEC} \
        --refSequence mscaf_a1_{wildcards.mscaf} \
        --chunkSize 1000000 \
        --noAncestors
      
      rm -r {params.js}
      """

rule hal_to_snps:
    input:
      hal = 'results/cactus/{name}.hal'.format(name = P_NAME)
    output:
      tsv = "results/pinniped/snps/{name}_{mscaf}.tsv"
    container: c_cactus
    shell:
      """
      halSnps \
        --refSequence mscaf_a1_{wildcards.mscaf} \
        {input.hal} \
        {REF_SPEC} {TIP_SPECS} \
        --unique \
        --tsv {output.tsv}
      """