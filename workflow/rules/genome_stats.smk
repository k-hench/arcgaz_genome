"""
snakemake -n --configfile workflow/config.yml -R fa_stats

snakemake \
  --jobs 10 \
  --latency-wait 30 \
  --use-conda -p \
  --configfile workflow/config.yml \
  --default-resources mem_mb=25600 disk_mb=5000 \
  --cluster '
     qsub \
     -V -cwd \
     -P fair_share \
     -l idle=1 \
     -l si_flag=1 \
     -pe multislot {threads} \
     -l vf={resources.mem_mb}' \
  -R fa_stats
"""

rule fa_stats:
  input: 
    expand("results/fa_stats/{gn}_gaps.bed.gz", gn = U_ALL),
    expand("results/fa_stats/{gn}_stats_seqkit.tsv", gn = U_ALL),
    expand("results/fa_stats/{gn}_stats_a.json", gn = U_ALL),
    expand("results/fa_stats/{gn}_stats_pb.txt", gn = U_ALL),
    expand("results/fa_stats/{gn}_stats_bp.yml", gn = U_ALL)

rule gap_bed:
  input: "data/genomes/{genome}.fa.gz"
  output: "results/fa_stats/{genome}_gaps.bed.gz"
  conda: "fa_stats"
  script:
    """
    python py/generate_masked_ranges.py {input} | gzip > {output}
    """

rule seqkit_stats:
  input: "data/genomes/{genome}.fa.gz"
  output: "results/fa_stats/{genome}_stats_seqkit.tsv"
  conda: "fa_stats"
  script:
    """
    seqkit stats -a {input} > {output}
    """

rule assembly_stats:
  input: "data/genomes/{genome}.fa.gz"
  output: "results/fa_stats/{genome}_stats_a.json"
  conda: "fa_stats"
  script:
    """
    zcat {input} | assembly_stats /dev/stdin > {output}
    """

rule pbjelly_stats:
  input: "data/genomes/{genome}.fa.gz"
  output: "results/fa_stats/{genome}_stats_pb.txt"
  conda: "pb_jelly_utils"
  script:
    """
    zcat {input} | python py/summarizeAssembly.py /dev/stdin > {output}
    """

rule bp_stats:
  input: "data/genomes/{genome}.fa.gz"
  output: "results/fa_stats/{genome}_stats_bp.yml"
  conda: "fa_stats"
  script:
    """
    python py/fa_bp_summary.py {input} {output}
    """