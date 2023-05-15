"""
snakemake -n --configfile workflow/config.yml -R fa_stats

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
  log: 'logs/fa_stats/{genome}_bed.log'
  conda: "fa_stats"
  shell:
    """
    python py/generate_masked_ranges.py {input} | gzip > {output} 2>{log}
    """

rule seqkit_stats:
  input: "data/genomes/{genome}.fa.gz"
  output: "results/fa_stats/{genome}_stats_seqkit.tsv"
  log: 'logs/fa_stats/{genome}_seqkit.log'
  conda: "fa_stats"
  shell:
    """
    seqkit stats -a {input} > {output} 2>{log}
    """

rule assembly_stats:
  input: "data/genomes/{genome}.fa.gz"
  output: "results/fa_stats/{genome}_stats_a.json"
  log: 'logs/fa_stats/{genome}_astats.log'
  conda: "fa_stats"
  shell:
    """
    zcat {input} | assembly_stats /dev/stdin > {output} 2>{log}
    """

rule pbjelly_stats:
  input: "data/genomes/{genome}.fa.gz"
  output: "results/fa_stats/{genome}_stats_pb.txt"
  log: 'logs/fa_stats/{genome}_pbstats.log'
  conda: "pb_jelly_utils"
  shell:
    """
    zcat {input} | python py/summarizeAssembly.py /dev/stdin > {output} 2>{log}
    """

rule bp_stats:
  input: "data/genomes/{genome}.fa.gz"
  output: "results/fa_stats/{genome}_stats_bp.yml"
  log: 'logs/fa_stats/{genome}_bpstats.log'
  conda: "fa_stats"
  shell:
    """
    python py/fa_bp_summary.py {input} {output} 2>{log}
    """