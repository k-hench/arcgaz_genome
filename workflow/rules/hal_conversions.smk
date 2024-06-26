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

rule convert_hal:
    input: 
      maf = expand("results/pinniped/maf/{name}_{mscaf}.maf", name = P_NAME, mscaf = MSCAFS),
      snps = expand("results/pinniped/snps/{name}_{mscaf}.tsv", name = P_NAME, mscaf = MSCAFS),

rule hal_to_maf:
    input:
      hal = 'results/cactus/{name}.hal'.format(name = P_NAME)
    output:
      maf = "results/pinniped/maf/{name}_{mscaf}.maf"
    log: "logs/hal_to_maf_{name}_{mscaf}.log"
    params:
      sif = c_cactus,
      js = "results/cactus/scratch/pinniped_set/",
      local_js = "js_{name}_{mscaf}",
      run = "run_{name}_{mscaf}"
    shell:
      """
      readonly CACTUS_IMAGE={params.sif} 
      readonly CACTUS_SCRATCH={params.js}\

      mkdir -p {params.js}{params.run}

      apptainer exec --cleanenv \
        --fakeroot --overlay ${{CACTUS_SCRATCH}} \
        --bind ${{CACTUS_SCRATCH}}/tmp:/tmp,{params.js}{params.run}:/run,$(pwd),{s_bind_paths} \
        --env PYTHONNOUSERSITE=1 \
        {params.sif} \
        cactus-hal2maf \
        /tmp/{params.local_js} \
        {input.hal} \
        {output.maf} \
        --refGenome {REF_SPEC} \
        --refSequence mscaf_a1_{wildcards.mscaf} \
        --dupeMode single \
        --filterGapCausingDupes \
        --chunkSize 1000000 \
        --noAncestors 2> {log}
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
        --noDupes \
        --unique \
        --tsv {output.tsv}
      """