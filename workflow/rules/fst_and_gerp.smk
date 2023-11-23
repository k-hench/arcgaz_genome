"""
snakemake --configfile workflow/config.yml --rerun-triggers mtime -n -R win_and_busco

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
        --export=ALL \
        -n {threads} \
        -e logs/{name}.{jobid}.err \
        -o logs/{name}.{jobid}.out \
        --mem={resources.mem_mb}' \
        --jn job_c.{name}.{jobid}.sh \
        -R win_and_busco
"""

rule win_and_busco:
  input:
    win_summaries = expand( "results/pinniped/{stat}/tsv/{stat}_snps_{mscaf}_summary.tsv.gz", mscaf = SCFS, stat = ["gerp", "fst"] ),
    busco_summaries = expand( "results/pinniped/{stat}/tsv/{stat}_busco_{mscaf}_summary.tsv.gz", mscaf = SCFS, stat = ["gerp", "fst"] )

def name_to_win_size(wildcards):
  pattern1 = re.compile(r'_w(.*?)k')
  pattern2 = re.compile(r'_s(.*?)k')
  out1 = re.findall(pattern1, wildcards.wsizes )
  out2 = re.findall(pattern2, wildcards.wsizes )
  return [ int(out1[0]) * 1000, int(out2[0]) * 1000 ]

rule create_window_bed:
    input:
      genome = "data/genomes/arcgaz_anc_h1.genome"
    output:
      bed = "data/genomes/arcgaz_anc_h1{wsizes}.bed.gz"
    params:
      w_sizes = lambda wc: name_to_win_size(wc)
    container: c_conda
    conda: "popgen_basics"
    shell:
      """
      bedtools makewindows \
        -s {params.w_sizes[1]} \
        -g {input.genome} \
        -w {params.w_sizes[0]} | \
        awk '{{print $0"\t"NR}}' | \
        gzip > {output.bed}
      """

rule create_busco_bed:
    input:
      dir_busco = "results/busco/arcgaz_anc_h1"
    output:
      bed = "results/pinniped/complete_buscos.bed.gz"
    params:
      tsv_busco = "run_carnivora_odb10/full_table.tsv"
    container: c_conda
    conda: "popgen_basics"
    shell:
      """
      tail -n +3  {input.dir_busco}/{params.tsv_busco} | \
        grep "Complete" | \
        awk -v OFS="\t" '{{ if ($4<$5) {{s = $4; e = $5}} else {{s =$5; e = $4}} {{print $3,s,e,$1}} }}' | \
        sort -k 1,1 -k2,2n | \
        grep "mscaf_a1" | \
        gzip > {output.bed}
      """

rule gerp_bed:
    input:
      tsv = "results/pinniped/maf/pinniped_set_{mscaf_nr}.maf.rates"
    output:
      bed = temp( "results/pinniped/gerp/beds/gerp_mscaf_a1_{mscaf_nr}.bed.gz" )
    params:
      outdir = "results/pinniped/gerp/beds/"
    log: "logs/gerp_to_bed_{mscaf_nr}.log"
    container: c_conda
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/gerp_to_bed.R {input.tsv} {output.bed} {params.outdir} &> {log}
      """

rule sliding_gerp:
    input:
      gerp = "results/pinniped/gerp/beds/gerp_{mscaf}.bed.gz",
      win = "data/genomes/arcgaz_anc_h1_w50k_s25k.bed.gz"
    output:
      tsv = "results/pinniped/gerp/beds/gerp_snps_{mscaf}.tsv.gz"
    container: c_conda
    conda: "popgen_basics"
    shell:
      """
      # chr, start(-1), end, win_idx, neutral_n, gerp_rs
      bedtools intersect \
        -a {input.win} \
        -b {input.gerp} \
        -wa -wb | \
        cut -f 1-3,4,8,9 | \
        gzip > {output.tsv}
      """

rule busco_gerp:
    input:
      gerp = "results/pinniped/gerp/beds/gerp_{mscaf}.bed.gz",
      busco = "results/pinniped/complete_buscos.bed.gz"
    output:
      tsv = "results/pinniped/gerp/beds/gerp_busco_{mscaf}.tsv.gz"
    container: c_conda
    conda: "popgen_basics"
    shell:
      """
      # chr, pos, neutral_n, gerp_rs, busco_id
      bedtools intersect \
        -a {input.gerp} \
        -b {input.busco} \
        -wa -wb | \
        cut -f 1,3,4,5,9 | \
        gzip > {output.tsv}
      """

rule summarize_gerp_busco:
    input:
      tsv = "results/pinniped/gerp/beds/gerp_busco_{mscaf}.tsv.gz"
    output:
      tsv = "results/pinniped/gerp/tsv/gerp_busco_{mscaf}_summary.tsv.gz"
    container: c_conda
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/summarize_gerp_busco.R {input.tsv} {output.tsv}
      """

rule summarize_gerp_win:
    input:
      tsv = "results/pinniped/gerp/beds/gerp_snps_{mscaf}.tsv.gz"
    output:
      tsv = "results/pinniped/gerp/tsv/gerp_snps_{mscaf}_summary.tsv.gz"
    container: c_conda
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/summarize_gerp_win.R {input.tsv} {output.tsv}
      """

rule fst_to_bed:
    input:
      tsv = "results/fst_bp.tsv.gz"
    output:
      bed = "results/pinniped/fst/bed/fst_bp_{mscaf}.bed.gz"
    shell:
      """
      zgrep {wildcards.mscaf} {input.tsv} | \
        grep -v "nan" | \
        awk -v OFS="\t" '{{print $1,$2-1,$2,$3}}' | \
        gzip > {output.bed}
      """

rule busco_fst:
    input:
      fst = "results/pinniped/fst/bed/fst_bp_{mscaf}.bed.gz",
      busco = "results/pinniped/complete_buscos.bed.gz"
    output:
      tsv = "results/pinniped/fst/beds/fst_busco_{mscaf}.tsv.gz"
    container: c_conda
    conda: "popgen_basics"
    shell:
      """
      # chr, pos, fst, busco_id
      bedtools intersect \
        -a {input.fst} \
        -b {input.busco} \
        -wa -wb | \
        cut -f 1,3,4,8 | \
        gzip > {output.tsv}
      """

rule summarize_fst_busco:
    input:
      tsv = "results/pinniped/fst/beds/fst_busco_{mscaf}.tsv.gz"
    output:
      tsv = "results/pinniped/fst/tsv/fst_busco_{mscaf}_summary.tsv.gz"
    container: c_conda
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/summarize_fst_busco.R {input.tsv} {output.tsv}
      """

rule sliding_fst_manual:
    input:
      fst = "results/pinniped/fst/bed/fst_bp_{mscaf}.bed.gz",
      win = "data/genomes/arcgaz_anc_h1_w50k_s25k.bed.gz"
    output:
      tsv = "results/pinniped/fst/beds/fst_snps_{mscaf}.tsv.gz"
    container: c_conda
    conda: "popgen_basics"
    shell:
      """
      # chr, start(-1), end, win_idx, fst
      bedtools intersect \
        -a {input.win} \
        -b {input.fst} \
        -wa -wb | \
        cut -f 1-3,4,8 | \
        gzip > {output.tsv}
      """

rule summarize_fst_win:
    input:
      tsv = "results/pinniped/fst/beds/fst_snps_{mscaf}.tsv.gz"
    output:
      tsv = "results/pinniped/fst/tsv/fst_snps_{mscaf}_summary.tsv.gz"
    container: c_conda
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/summarize_fst_win.R {input.tsv} {output.tsv}
      """