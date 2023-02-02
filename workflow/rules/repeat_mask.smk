rule mask:
    input:
      'data/genomes/{species}.fa.gz'
    output:
      'data/genomes/{species}_masked.fa.gz'
    log:
      'logs/mask_plot_{species}.log'
    shell:
      """
      touch {output}
      """

rule kmer_cout:
    input:
      'data/genomes/{species}.fa.gz'
    output:
      'results/repeat_maksing/{species}.kat.m20.hist'
    log:
      'logs/kmer_hist_{species}.log'
    shell:
      """
      kat hist -t 12 {input} -m 20 -o {output}
      """