'''
-----------------------------------------------
NOTE ON THE FINAL FIGURES PART OF THE PIPELINE:
The scripts handled herein were in fact not run using snakmake but executed manually.
They were bundled here nontheless to conserve the dependency structure of the individual sctipts in relation to the rest of the analysis.
So, while this part was not run with snakemake when intially creating the figures, it should be possible to recreate theme with this smk file.

However, due to originally not being part of the snakemake pipeline, the computational environment was captured differently compared to the rest of the analysis.
In fact this part is simply a chain of individual R scripts that produce the final figures based on the output of the previous analysis.
The R environment in which they ran is captured in the file `renv.lock` and can be recreated by running `renv::restore()` within R.
-----------------------------------------------

snakemake -n --configfile workflow/config.yml -R final_plots
snakemake --rulegraph --configfile workflow/config.yml -R final_plots | dot -Tsvg > results/img/dag_final_plots.svg
'''

rule final_plots:
    input:
      tables = ["results/tab/aligned_genomes.tex"],
      plots = expand("results/img/{plt}", plt = ["arcgaz_a1_zalcal.pdf", "busco_go_term_2d_dens.pdf", "busco_go_term_2d_sub.pdf", "fig_s_anchoring.pdf", "fig_s_hirise_sizes.pdf", "hal_coverage.pdf", "neutral_tree.pdf", "win_coverage.pdf", "win_gerp_fst.pdf", "zoom_win_outlier_fst.pdf", "zoom_win_outlier_gerp.pdf", "zoom_win_outlier.pdf"])

rule figure_F1:
    input: "results/neutral_tree/rerooted.tree"
    output: "results/img/neutral_tree.pdf"
    shell:
      """
      Rscript R/plot_neutral_tree.R
      """

rule figure_F2:
    input:
      sizes = ["results/genome/arcgaz_anc_h1.size", "results/genome/zalcal_v1.size"],
      alignment = "results/psl/slim_zalcal_v1_on_arcgaz_anc_h1.psl.gz"
    output: "results/img/arcgaz_a1_zalcal.pdf"
    shell:
      """
      Rscript R/plot_genome_alignment.R
      """

rule input_for_F3:
    input:
      expand( "results/pinniped/gerp/tsv/gerp_snps_{scf}_summary.tsv.gz", scf = MSCAFS ),
      expand( "results/pinniped/gerp/tsv/gerp_busco_{scf}_summary.tsv.gz", scf = MSCAFS ),
      expand( "results/pinniped/fst/tsv/fst_snps_{scf}_summary.tsv.gz", scf = MSCAFS ),
      expand( "results/pinniped/fst/tsv/fst_busco_{scf}_summary.tsv.gz", scf = MSCAFS ),
      "data/genomes/arcgaz_anc_h1.genome",
      "results/busco/arcgaz_anc_h1/run_carnivora_odb10/full_table.tsv",
      "results/fst_win.tsv.gz"
    output:
      "results/pinniped/go_terms/thresholds.tsv",
      "results/pinniped/busco_gerp_fst.tsv",
      "results/pinniped/win_gerp_fst.tsv",
      "results/pinniped/win_outlier_summary.tsv"
    shell:
      """
      Rscript R/compile_window_and_busco_stats.R
      """

rule figure_F3:
    input:
      "data/genomes/arcgaz_anc_h1.genome",
      "results/pinniped/go_terms/thresholds.tsv",
      "results/pinniped/busco_gerp_fst.tsv",
      "results/pinniped/win_gerp_fst.tsv",
      "results/pinniped/win_outlier_summary.tsv"
    output: 
      plot = "results/img/win_gerp_fst.pdf",

    shell:
      """
      Rscript R/plot_gerp_and_fst.R
      """

rule link_busco_and_go_terms:
    input:
      "results/busco/arcgaz_anc_h1/run_carnivora_odb10/full_table.tsv"
    output:
      expand( "results/pinniped/go_terms/busco_id_to_go_term_{type}.tsv", type = [ "bp", "cc", "mf" ] )
    shell:
      """
      Rscript R/query_busco_go_terms.R
      """

rule figure_F4_and_SF7:
    input:
      "results/pinniped/go_terms/busco_id_to_go_term_bp.tsv",
      "results/pinniped/busco_gerp_fst.tsv",
      "results/pinniped/go_terms/thresholds.tsv"
    output:
      "results/pinniped/go_terms/go_term_info.tsv",
      "results/pinniped/go_terms/go_term_busco_stats.tsv",
      "results/img/go_sub_graph.pdf",
      "results/img/busco_go_term_2d_dens.pdf",
       "results/img/busco_go_term_2d_sub.pdf"
    shell:
      """
      Rscript R/go_term_enrichment.R
      """

rule figure_SF1:
    input:
      sizes = expand("results/genome/{genome}.size", genome = ["arcgaz_dt_h1_hardmasked", "arcgaz_dt_h2_hardmasked", "arcgaz_v3_hardmasked"])
    output: "results/img/fig_s_hirise_sizes.pdf"
    shell:
      """
      Rscript R/plot_pre_anchoring_sizes.R
      """

rule figure_SF2:
    input:
     bed = expand("results/anchoring/arcgaz_dt_h{ht}_hardmasked/anchored_arcgaz_dt_h{ht}.lifted.bed", ht = [1, 2]),
     faidx = expand("results/anchoring/arcgaz_dt_h{ht}_hardmasked/anchored_arcgaz_dt_h{ht}.fasta.gz.fai", ht = [1, 2])
    output: "results/img/fig_s_anchoring.pdf"
    shell:
      """
      Rscript R/plot_anchoring_assesment.R
      """

rule figure_SF3:
    input:
      expand("results/neutral_tree/cov/{scf}.collapsed.bed.gz", scf = MSCAFS),
      expand("results/neutral_tree/cov/fam/{fam}-{scf}.collapsed.bed.gz", fam = ["ota", "pho"], scf = MSCAFS)
    output: "results/img/hal_coverage.pdf"
    shell:
      """
      Rscript R/plot_hal_coverage.R
      """

rule figure_SF4:
    input: "results/pinniped/win_gerp_fst.tsv"
    output: "results/img/win_coverage.pdf"
    shell:
      """
      Rscript R/plot_coverage_windows.R
      """

rule figure_SF5_and_SF6:
    input:
      "data/genomes/arcgaz_anc_h1.genome",
      "results/pinniped/go_terms/thresholds.tsv",
      "results/pinniped/busco_gerp_fst.tsv",
      "results/pinniped/win_gerp_fst.tsv",
      "results/pinniped/win_outlier_summary.tsv"
    output:
      sf5 = "results/img/zoom_win_outlier_gerp.pdf",
      sf6 = "results/img/zoom_win_outlier_fst.pdf"
    shell:
      """
      Rscript R/plot_zoom.R
      """