"""
snakemake --configfile workflow/config.yml --rerun-triggers mtime -n -R fsts

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
        -R fsts
"""

rule fsts:
  input:
    bp = "results/fst_bp.tsv.gz",
    win = "results/fst_win.tsv.gz"

rule vcf_header:
  input:
    vcf = "data/proto.vcf",
    genome = "data/genomes/arcgaz_anc_h1.genome"
  output:
    vcf = "results/vcf/header.vcf"
  params:
    date = datetime.datetime.today().strftime('%Y%m%d'),
    reference = "arcgaz_anc_h1.fa",
    species = '\\"dummy\\"'
  shell:
    """
    sed 's/XXgenomeXX/{params.reference}/; s/XXdateXX/{params.date}/' {input.vcf} > {output.vcf}
    awk '{{print "##contig=<ID="$1",length="$2",species={params.species}>"}}' {input.genome} >> {output.vcf}
    echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> {output.vcf}
    """

rule convert_scaf_vcf:
    input:
      tsv = lambda wc: "results/pinniped/snps/pinniped_set_" + scaf_to_nr(wc) + ".tsv"
    output:
      vcf = temp( "results/vcf/{mscaf}.vcf" )
    container: c_conda
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/snps_to_vcf.R {input.tsv} {output.vcf}
      """

rule combine_vcfs:
    input: 
      vcf_header = "results/vcf/header.vcf",
      vcf_body = expand( "results/vcf/{mscaf}.vcf", mscaf = SCFS )
    output:
      vcf = "results/vcf/genotpyes.vcf.gz"
    params:
      prefix = "results/vcf/genotpyes.vcf"
    shell:
      """
      cat {input.vcf_header} > {params.prefix}
      head -n 1 {input.vcf_body[0]} >> {params.prefix}
      tail -n +2 {input.vcf_body} | grep -v "==>" | grep -v "^$" >> {params.prefix}
      gzip {params.prefix}
      """

rule fst_win:
    input:
      vcf = "results/vcf/genotpyes.vcf.gz",
      pops = expand( "data/{p}.pop", p = [ "otariidae", "phocidae" ] )
    output:
      fst = "results/fst_win.tsv.gz"
    params:
      wsize = 50000,
      wstep = 25000
    container: c_vcfh
    log: "logs/fst_win.log"
    shell:
      """
      vcftools_haploid \
        --haploid \
        --gzvcf {input.vcf} \
        --weir-fst-pop {input.pops[0]} \
        --weir-fst-pop {input.pops[1]} \
        --fst-window-size {params.wsize} \
        --fst-window-step {params.wstep} \
        --stdout 2> {log} | gzip > {output.fst}
      """

rule fst_bp:
    input:
      vcf = "results/vcf/genotpyes.vcf.gz",
      pops = expand( "data/{p}.pop", p = [ "otariidae", "phocidae" ] )
    output:
      fst = "results/fst_bp.tsv.gz"
    container: c_vcfh
    log: "logs/fst_bp.log"
    shell:
      """
      vcftools_haploid \
        --haploid \
        --gzvcf {input.vcf} \
        --weir-fst-pop {input.pops[0]} \
        --weir-fst-pop {input.pops[1]} \
        --stdout 2> {log} | gzip > {output.fst}
      """