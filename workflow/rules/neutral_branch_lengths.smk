"""
snakemake --configfile workflow/config.yml --rerun-triggers mtime -n -R create_neutral_tree

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
      -R create_neutral_tree
"""

WIN_SIZE = 5000
WIN_N = 1000
SCFS = expand( "mscaf_a1_{scf}", scf = MSCAFS)

rule create_neutral_tree:
    input:
      win_bed = expand( "results/neutral_tree/win/windows_{mscaf}.bed.gz", mscaf = SCFS )
      #,
      #multifa = "results/neutral_tree/multifa/combined_windows.fa",
      #rooted_tree = "results/neutral_tree/rerooted.tree",
      #gerp = expand( "results/maf/{mscaf}.maf.rates", mscaf = SCFS )
      #winmaf = "results/neutral_tree/win/windows.maf.gz"
      #stats = expand( "results/neutral_tree/stats/{mscaf}_statistics.csv", mscaf = SCFS )


# we need to determine what part of the genome is covered 
# by a the alignment of all other (tip) species
# first step for this is to create a wig file from the hal
rule alignment_coverage:
    input:
      hal = 'results/cactus/{name}.hal'.format(name = P_NAME)
    output:
      wig = "results/neutral_tree/wig/{mscaf}.wig.gz"
    params:
      prefix = "results/neutral_tree/wig/{mscaf}.wig"
    container: c_cactus
    shell:
      """
      halAlignmentDepth \
        {input.hal} \
        {REF_SPEC} \
        --noAncestors \
        --refSequence {wildcards.mscaf} \
        --outWiggle {params.prefix}
      gzip {params.prefix}
      """

# ultimately we want a bed file as mask, so we convert the wig to bed format
rule wig_to_bed:
    input:
      wig = "results/neutral_tree/wig/{mscaf}.wig.gz"
    output:
      bed = "results/neutral_tree/cov/raw/{mscaf}.bed.gz" 
    conda: "bedops"
    shell:
      """
      zcat {input.wig} | wig2bed | gzip > {output.bed}
      """

# unfortunetely, the original bed is single bp elements,
# so we collapse them into chunks of equal coverage 
rule collapse_cov_bed:
    input:
      bed = "results/neutral_tree/cov/raw/{mscaf}.bed.gz"
    output:
      bed = "results/neutral_tree/cov/{mscaf}.collapsed.bed.gz"
    log: "logs/collapse_cov_bed_{mscaf}.log"
    conda: "r_tidy"
    shell:
      """
      (set -o posix ; set) > {output}
      echo "===========================" >> {output}
      hostname >> {output}
      echo "-------" >> {output}
      Rscript --vanilla R/collapse_bed_coverage.R {input.bed} {output.bed} &>> {log}
      """

# now we filter the coverage bed to a minimum coverage
# and attach a fourth (dummy) column to match the bedgraph format
# for maffilter
rule filter_maf_coverage:
   input:
      bed = "results/neutral_tree/cov/{mscaf}.collapsed.bed.gz"
   output:
     bed = "results/neutral_tree/cov/filtered/{mscaf}.bed.gz"
   params:
     min_cov = 10
   shell:
     """
     zcat {input.bed} | \
       awk '$4>{params.min_cov}{{print $0}}' | \
       gzip > {output.bed}
     """

# a second mask is created from the genome annotation
# (we want to exclude all CDS)
rule create_cds_mask:
    input:
      gff = "data/genomes/arcgaz_anc_h1.annotation.gff.gz"
    output:
      bed = "results/neutral_tree/masks/cds_{mscaf}.bed.gz"
    shell:
      """
      zgrep "CDS" {input.gff} | \
        grep {wildcards.mscaf} | \
        cut -f 1,4,5 | \
        gzip > {output.bed}
      """

rule create_windows:
    input:
      bed_cov = expand( "results/neutral_tree/cov/filtered/{mscaf}.bed.gz", mscaf = SCFS ),
      bed_cds = expand( "results/neutral_tree/masks/cds_{mscaf}.bed.gz", mscaf = SCFS ),
      genome = "data/genomes/arcgaz_anc_h1.genome"
    output:
      bed_cov = "results/neutral_tree/cov/filtered/cov.bed.gz",
      bed_cds = "results/neutral_tree/masks/cds.bed.gz",
      win_proto = "results/neutral_tree/win/proto.bed.gz",
      bed_win = "results/neutral_tree/win/windows.bed.gz",
      win_n_scaf = "results/neutral_tree/win/win_n_scaf.txt"
    params:
      win_size = WIN_SIZE,
      win_n = WIN_N,
      win_proto_prefix = "results/neutral_tree/win/proto.bed",
      win_seed = 42
    conda: "popgen_basics"
    log: "logs/win.log"
    shell:
      """
      zcat {input.bed_cov} | gzip > {output.bed_cov}
      zcat {input.bed_cds} | gzip > {output.bed_cds}

      mkdir -p results/neutral_tree/win/
      for k in $(seq 1 {params.win_n}); do echo -e "{SCFS[0]}\t0\t{params.win_size}" >> {params.win_proto_prefix}; done
      gzip {params.win_proto_prefix}

      bedtools shuffle \
        -i {output.win_proto} \
        -g {input.genome} \
        -incl {output.bed_cov} \
        -excl {output.bed_cds} \
        -seed {params.win_seed} \
        -noOverlapping 2>> {log} | \
        sort -k 1,1 -k2,2n | \
        gzip > {output.bed_win}
      
      for k in {SCFS}; do WN=$(zgrep ${{k}} {output.bed_win} | wc -l); echo -e "${{k}}\t${{WN}}" >> {output.win_n_scaf}; done 
      """

rule windows_by_scaffold:
    input:
      bed_win = "results/neutral_tree/win/windows.bed.gz"
    output:
      bed_win = "results/neutral_tree/win/windows_{mscaf}.bed.gz"
    params:
      bed_prefix = "results/neutral_tree/win/windows_{mscaf}.bed"
    shell:
      """
      echo 'track type=bedGraph name="BedGraph Format" description="BedGraph format"' > {params.bed_prefix}
      zgrep {wildcards.mscaf} {input.bed_win} | awk '{{print $0"\t0"}}' >> {params.bed_prefix}
      gzip {params.bed_prefix}
      """

