rule bracken:
    input:
        kraken="results/kraken/{sample}.report"
    output:
        BRACKEN_OUT
    shell:
        """
        echo -e "name\ttaxonomy_id\tabundance\n\
Escherichia coli\t562\t100" > {output}
        """