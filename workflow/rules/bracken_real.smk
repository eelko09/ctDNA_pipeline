rule bracken:
    input:
        kraken="results/kraken/{sample}.report"
    output:
        BRACKEN_OUT
    conda:
        "envs/bracken.yaml"
    threads: 8
    log:
        "logs/bracken/{sample}.log"
    shell:
        """
        bracken \
          -d {config[kraken_db]} \
          -i {input.kraken} \
          -o {output} \
          -r 150 \
          -l S \
          > {log} 2>&1
        """