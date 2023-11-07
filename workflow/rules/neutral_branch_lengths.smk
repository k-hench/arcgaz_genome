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
        --export=ALL \
        -n {threads} \
        -e logs/{name}.{jobid}.err \
        -o logs/{name}.{jobid}.out \
        --mem={resources.mem_mb}' \
        --jn job_c.{name}.{jobid}.sh \
        -R create_neutral_tree
"""
WIN_SIZE = 1000
WIN_N = 5000
SCFS = expand( "mscaf_a1_{scf}", scf = MSCAFS)

rule create_neutral_tree:
    input:
      win_bed = expand( "results/neutral_tree/win/windows_{mscaf}.bed.gz", mscaf = SCFS ),
      tree = "results/neutral_tree/multifa/combined_windows.fa.treefile",
      gerp = expand( "results/pinniped/maf/pinniped_set_{mscaf_nr}.maf.rates", mscaf_nr = MSCAFS )

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
    container: c_conda
    conda: "r_tidy"
    shell:
      """
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
     min_cov = 9
   shell:
     """
     zcat {input.bed} | \
       awk '$4>{params.min_cov}{{print $0}}' | \
       gzip > {output.bed}
     """

rule negative_coverage_mask:
    input:
      genome = "data/genomes/arcgaz_anc_h1.genome",
      bed = expand( "results/neutral_tree/cov/filtered/{mscaf}.bed.gz", mscaf = SCFS )
    output:
      bed_cov = "results/neutral_tree/cov/filtered/whole_genome.bed.gz",
      bed_genome = "data/genomes/arcgaz_anc_h1.bed",
      bed_neg_cov = "results/neutral_tree/cov/filtered/whole_genome_exclude.bed.gz"
    container: c_conda
    conda: "popgen_basics"
    shell:
      """
      zcat {input.bed} | gzip > {output.bed_cov}
      awk '{{print $1"\t"0"\t"$2}}' {input.genome} > {output.bed_genome}
      bedtools subtract -a {output.bed_genome} -b {output.bed_cov} | gzip > {output.bed_neg_cov}
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

rule proto_windows:
    input:
      bed_cds = expand( "results/neutral_tree/masks/cds_{mscaf}.bed.gz", mscaf = SCFS )
    output:
      bed_cds = "results/neutral_tree/masks/cds.bed.gz",
      win_proto = "results/neutral_tree/win/proto.bed.gz"
    params:
      win_size = WIN_SIZE,
      win_n = WIN_N,
      win_proto_prefix = "results/neutral_tree/win/proto.bed",
      win_seed = 42
    log: "logs/win.log"
    shell:
      """
      zcat {input.bed_cds} | gzip > {output.bed_cds}

      mkdir -p results/neutral_tree/win/
      for k in $(seq 1 {params.win_n}); do echo -e "{SCFS[0]}\t0\t{params.win_size}" >> {params.win_proto_prefix}; done
      gzip {params.win_proto_prefix}
      """

rule shuffle_windows:
    input:
      genome = "data/genomes/arcgaz_anc_h1.genome",
      bed_cds = "results/neutral_tree/masks/cds.bed.gz",
      bed_neg_cov = "results/neutral_tree/cov/filtered/whole_genome_exclude.bed.gz",
      win_proto = "results/neutral_tree/win/proto.bed.gz"
    output:
      bed_tmp = temp( "results/neutral_tree/win/exclude_tmp.bed.gz" ),
      bed_exclude = "results/neutral_tree/win/exclude.bed.gz",
      bed_win = "results/neutral_tree/win/windows.bed.gz",
      win_n_scaf = "results/neutral_tree/win/win_n_scaf.txt"
    params:
      win_seed = 42
    container: c_conda
    conda: "popgen_basics"
    log: "logs/win.log"
    shell:
      """
      zcat {input.bed_cds} {input.bed_neg_cov} | \
        sort -k 1,1 -k2,2n | \
        gzip > {output.bed_tmp}

      bedtools merge -i {output.bed_tmp} | gzip > {output.bed_exclude}

      bedtools shuffle \
        -i {input.win_proto} \
        -g {input.genome} \
        -excl {output.bed_exclude} \
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

def scaf_to_nr(wildcards):
  pattern = re.compile(r'mscaf_a1_(.*?)$')
  out = re.findall(pattern, wildcards.mscaf)
  return out[0]

rule maf_to_fasta:
    input:
      maf = lambda wc: "results/pinniped/maf/pinniped_set_" + scaf_to_nr(wc) + ".maf",
      conf = "data/maffilter_templ.txt",
      windows = "results/neutral_tree/win/windows_{mscaf}.bed.gz",
      win_n_scaf = "results/neutral_tree/win/win_n_scaf.txt"
    output:
      fasta = "results/neutral_tree/multifa/{mscaf}.fa.gz"
    params:
      ref_spec = REF_SPEC,
      nr = lambda wc: scaf_to_nr(wc)
    log: "logs/maf2fa_{mscaf}.log"
    container: c_conda
    conda: "maffilter"
    shell:
      """
      N_WIN=$(grep {wildcards.mscaf} {input.win_n_scaf} | cut -f 2)
      BP_WIN=$((${{N_WIN}} * {WIN_SIZE}))
      maffilter param={input.conf} NR={params.nr} DATA={wildcards.mscaf} FASIZE=${{BP_WIN}} REF_SPEC={params.ref_spec} &> {log}
      """

rule single_multi_fasta:
    input:
      fas = expand( "results/neutral_tree/multifa/{mscaf}.fa.gz", mscaf = SCFS )
    output:
      fa = "results/neutral_tree/multifa/combined_windows.fa"
    container: c_conda
    conda: "seqkit"
    shell:
      """
      seqkit concat {input.fas} --full > {output.fa}
      """

rule estimate_branchlengths:
    input:
      fa = "results/neutral_tree/multifa/combined_windows.fa",
      tree = "data/topology.tree"
    output:
      tree = "results/neutral_tree/multifa/combined_windows.fa.treefile"
    container: c_conda
    conda: "phylo_ml"
    shell:
      """
      iqtree \
        -s {input.fa}\
        -m GTR+G \
        --threads-max 1 \
        --seed 42 \
        -g {input.tree}\
        --tree-fix 
      """

rule reroot_tree:
    input:
      tree = "results/neutral_tree/multifa/combined_windows.fa.treefile"
    output:
      tree = "results/neutral_tree/rerooted.tree"
    container: c_conda
    conda: "r_phytools"
    params:
      root_node = 15
    log: "logs/reroot_tree.log"
    shell:
      """
      Rscript --vanilla R/reroot_tree.R {input.tree} {params.root_node} {output.tree} &> {log}
      """

rule call_gerp:
    input:
      maf = "results/pinniped/maf/pinniped_set_{mscaf_nr}.maf",
      tree = "results/neutral_tree/rerooted.tree"
    output:
      rates = "results/pinniped/maf/pinniped_set_{mscaf_nr}.maf.rates"
    params:
      refname = "arcgaz"
    conda: "msa_phast"
    log: "logs/gerp_{mscaf_nr}.log"
    shell:
      """
      gerpcol -t {input.tree} -f {input.maf} -e {params.refname} -j -z -x ".rates" &> {log}
      """