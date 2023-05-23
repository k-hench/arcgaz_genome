"""
snakemake -n --configfile workflow/config.yml -R all_blast

snakemake -c 6 --configfile workflow/config.yml --use-conda -R all_blast
"""

B_REF = ["arcgaz_anc_h1", "arcgaz_anc_h2", "zalcal_v1", "arcgaz_v3", "arcgaz_v1_2", "arcgaz_v1_4"]
BLAST_SEQS = [ "arcgaz_mhc", "arcgaz_mt", "dog_mhc_transcripts", "DQB_Haps_Plus_new" ]

rule all_blast:
  input: 
    expand("results/blast/{que}_on_{ref}.csv", ref = B_REF, que = BLAST_SEQS),
    expand("results/busco/{ref}", ref = B_REF)

"""
# qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
# Query accesion.version
# Subject accession.version
# Percentage of identical matches
# Alignment length
# Number of mismatches
# Number of gap openings
# Start of alignment in query
# End of alignment in query
# Start of alignment in subject
# End of alignment in subject
# Expect value
# Bit score
"""
rule build_blast_db:
  input: "data/genomes/{ref}.fa.gz"
  output: "data/genomes/blastdb/{ref}.ndb"
  conda: "map_align"
  shell:
    """
    zcat {input} | \
      makeblastdb -in - \
      -title {wildcards.ref} \
      -dbtype nucl \
      -parse_seqids \
      -out data/genomes/blastdb/{wildcards.ref}
    """

rule blast_querry:
  input:
    blastndb = "data/genomes/blastdb/{ref}.ndb",
    querry = "data/ncbi_seqs/{que}.fa.gz"
  output: "results/blast/{que}_on_{ref}.csv"
  params: db_prefix = "data/genomes/blastdb/{ref}"
  conda: "map_align"
  shell:
    """
    zcat {input.querry} | \
      blastn -db {params.db_prefix} \
      -query /dev/stdin \
      -out - \
      -outfmt 10 > {output}
    """

rule busco_dl:
  output: directory("data/busco_downloads/lineages/carnivora_odb10")
  conda: "busco"
  shell:
    '''
    cd data/
    busco --download carnivora_odb10
    '''

rule busco:
  input: 
    ref = "data/genomes/{ref}.fa.gz",
    busco_db = "data/busco_downloads/lineages/carnivora_odb10"
  output: directory("results/busco/{ref}")
  log: "logs/busco/{ref}.log"
  conda: "busco"
  shell:
    """
    mkdir -p results/busco/
    zcat {input.ref} > tmp/{wildcards.ref}.fa

    busco -i tmp/{wildcards.ref}.fa \
        -l data/busco_downloads/lineages/carnivora_odb10 \
        --offline \
        -o {output} \
        -m genome \
        -c 3 &> {log}
    
    rm tmp/{wildcards.ref}.fa
    """
