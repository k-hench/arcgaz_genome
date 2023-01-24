// ---
// engine: knitr
// ---
//
// <!--- note that this nextflow script is converted into a quarto file
//   by running `bash sh/convert_nf_to_qmd.sh nf/00_template/template.nf`
//   from the `${params.base}` of the project --->
//
// # nf: Indexing the Reference Genomes
//
// :::{.callout-note}
//
// ## Minimal Summary
//
// This [nextflow](https://www.nextflow.io/) pipeline [@ditommaso17] ([`template.nf`](https://github.com/k-hench/ARCGAZ_GENOME/blob/master/nf/0_template/template.nf)) <span style:"color=red">does the followeing</span>
// :::
//
// ## About This Pipeline
//
// The pipeline is structured into a [configuration section](#config-section), a section specifying the individual
// [pipeline components](#workflow-components) and the workflow
// [execution section](#workflow-composition).
//
// As a shorthand to start the workflow, the run command is bundled in the
// accompanying `bash` sript `run_genotyping.sh`, so the workflow can be
// started by run_genotyping
//
// ```sh
// cd ${params.base}/nf/00_template
// bash run_nf.sh
// ```
//
// :::{.callout-warning}
//
//  ## Different *Base Directories*
//
// In all nextflow scripts, there might arise some confusion about
// the difference between the variables `${baseDir}` and `${params.base}`:
//
// - `${params.base}` points to the root of the project folder (which is idendical with the base of this `git` repository)
// - `${baseDir}` points to the path of the executed `nextflow` pipeline (typically `${params.base}/nf/<some_dir>`)
// :::
//
// ## Config Section
//
// The following pipeline parameters specify the base and
// output directory, which could alternatively be provided as
// command line options:
//
// ```groovy
// // ----- Config section -----
params.base = "${baseDir}/../.."
params.outdir = "${params.base}/results"
params.genomedir = "${params.base}/data/genomes"
// // color logging
( c0, c1, c2 ) = [ "\033[0m", "\033[0;32m", "\033[1;30m" ]
// ```
//
// The workflow parameters are logged at the start of the workflow execution.
//
// ```groovy
log.info"""
${c2}ARCGAZ_GENOME${c0}
===================================
${c1}Author:${c0} Kosmas Hench
-----------------------------------
${c1}base_dir${c0}        : ${params.base}
${c1}params.outdir${c0}   : ${params.outdir}
${c1}params.genomedir${c0}   : ${params.genomedir}
"""
// ```
//
// The the use of [nexflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html#)
// is enabled, which concludes the configuration.
//
// ```groovy
nextflow.enable.dsl = 2
// ```
//
// ## Workflow Components
//
//
// ```groovy
// // ----- workflow components -----
Channel
  .from( [ h1: "arcgaz_dt_h1",
          h2: "arcgaz_dt_h2",
          v1: "arcgaz_v1.fa.gz",
          zc: "zalcal_v1.fa.gz"] )
  .set{ all_genomes_ch }
// ```
//
// Before the alignment, we need to create a database for the reference genome using `lastdb`.
//
// ```groovy
process create_last_db {
  publishDir "${params.genomedir}/", mode: 'copy'
  label "Q_lastdb_c_map"
  memory '10. GB'

  input:
  val( ref )

  output:
  file( "${ref}_db*" )

  script:
  """
  lastdb -c ${ref}_db ${params.genomedir}/${ref}.fa.gz
  """
}
// ```
// ## Workflow Composition
//
//
//
// ```groovy
// // ----- run workflow -----
workflow {
  main:
  refrence_genome_ch = all_genomes_ch.map{ it["h1"] }
  mapped_genomes_ch = all_genomes_ch.map{ [ it["h2"], it["v1"], it["zc"] ] }
  refrence_genome_ch | create_last_db
}
// ```
//
// ------
