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

rule cov_by_fam:
    input: 
      collapsed_cov = expand( "results/neutral_tree/cov/fam/{fam}-{mscaf}.collapsed.bed.gz", fam = ["pho", "ota"], mscaf = SCFS ),
      cov_by = expand( "results/neutral_tree/cov/by_{by}/{mscaf}.tsv.gz", mscaf = SCFS, by = [ "win", "busco" ] ),
      cov_by_fam = expand( "results/neutral_tree/cov/by_{by}/fam/{fam}-{mscaf}.tsv.gz", fam = ["pho", "ota"], mscaf = SCFS, by = [ "win", "busco" ] ),
      cov_combined_win = expand( "results/neutral_tree/cov/by_win/combined_{mscaf}.tsv.gz", mscaf = SCFS )

wildcard_constraints:
    fam = "[^_]*",
    mscaf = "[^-]*"
    

fams = { "pho": ["halgry", "lepwed", "mirang", "mirleo", "neosch", "odoros", "phovit"],
         "ota": ["calurs", "eumjub", "arcgaz", "zalcal"]}


fams_c = { "pho": "halgry,lepwed,mirang,mirleo,neosch,odoros,phovit",
           "ota": "calurs,eumjub,arcgaz,zalcal"}


rule alignment_coverage_fam:
    input:
      hal = 'results/cactus/{name}.hal'.format(name = P_NAME)
    output:
      wig = "results/neutral_tree/wig/fam/{fam}-{mscaf}.wig.gz"
    params:
      prefix = "results/neutral_tree/wig/fam/{fam}-{mscaf}.wig",
      specs = lambda wc: "{f}".format(f = fams_c[wc.fam] )
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

rule coverage_50k_intersect:
    input:
      win = "data/genomes/arcgaz_anc_h1_w50k_s25k.bed.gz",
      cov = "results/neutral_tree/cov/{mscaf}.collapsed.bed.gz"
    output:
      tsv = "results/neutral_tree/cov/by_win/{mscaf}.tsv.gz"
    params:
      prefix = "results/neutral_tree/cov/by_win/{mscaf}.tsv"
    container: c_conda
    conda: "popgen_basics"
    shell:
      """
      echo -e "chr\tstart\tend\twin_idx\tc_start\tc_end\tcov" > {params.prefix}

      bedtools intersect \
          -a {input.win} \
          -b {input.cov} \
          -wa -wb | \
          cut -f 1,2,3,4,6,7,8 >> {params.prefix}
      
      gzip {params.prefix}
      """

rule coverage_50k_intersect_fam:
    input:
      win = "data/genomes/arcgaz_anc_h1_w50k_s25k.bed.gz",
      cov = "results/neutral_tree/cov/fam/{fam}-{mscaf}.collapsed.bed.gz"
    output:
      tsv = "results/neutral_tree/cov/by_win/fam/{fam}-{mscaf}.tsv.gz"
    params:
      prefix = "results/neutral_tree/cov/by_win/fam/{fam}-{mscaf}.tsv"
    container: c_conda
    conda: "popgen_basics"
    shell:
      """
      echo -e "chr\tstart\tend\twin_idx\tc_start\tc_end\tcov" > {params.prefix}

      bedtools intersect \
          -a {input.win} \
          -b {input.cov} \
          -wa -wb | \
          cut -f 1,2,3,4,6,7,8 >> {params.prefix}
      
      gzip {params.prefix}
      """

rule comine_win_coverage:
    input:
      a  = "results/neutral_tree/cov/by_win/{mscaf}.tsv.gz",
      ota = "results/neutral_tree/cov/by_win/fam/ota-{mscaf}.tsv.gz",
      pho = "results/neutral_tree/cov/by_win/fam/pho-{mscaf}.tsv.gz"
    output:
      tsv = "results/neutral_tree/cov/by_win/combined/combined_{mscaf}.tsv.gz"
    log: "logs/comine_win_coverage_{mscaf}.log"
    container: c_conda
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/combine_alignment_coverage.R {wildcards.mscaf} &>> {log}
      """

rule coverage_busco_intersect:
    input:
      win = "results/pinniped/complete_buscos.bed.gz",
      cov = "results/neutral_tree/cov/{mscaf}.collapsed.bed.gz"
    output:
      tsv = "results/neutral_tree/cov/by_busco/{mscaf}.tsv.gz"
    params:
      prefix = "results/neutral_tree/cov/by_busco/{mscaf}.tsv"
    container: c_conda
    conda: "popgen_basics"
    shell:
      """
      echo -e "chr\tstart\tend\tbusco_id\tc_start\tc_end\tcov" > {params.prefix}

      bedtools intersect \
          -a {input.win} \
          -b {input.cov} \
          -wa -wb | \
          cut -f 1,2,3,4,6,7,8 >> {params.prefix}
      
      gzip {params.prefix}
      """

rule coverage_busco_intersect_fam:
    input:
      win = "results/pinniped/complete_buscos.bed.gz",
      cov = "results/neutral_tree/cov/fam/{fam}-{mscaf}.collapsed.bed.gz"
    output:
      tsv = "results/neutral_tree/cov/by_busco/fam/{fam}-{mscaf}.tsv.gz"
    params:
      prefix = "results/neutral_tree/cov/by_busco/fam/{fam}-{mscaf}.tsv"
    container: c_conda
    conda: "popgen_basics"
    shell:
      """
      echo -e "chr\tstart\tend\tbusco_id\tc_start\tc_end\tcov" > {params.prefix}

      bedtools intersect \
          -a {input.win} \
          -b {input.cov} \
          -wa -wb | \
          cut -f 1,2,3,4,6,7,8 >> {params.prefix}
      
      gzip {params.prefix}
      """
