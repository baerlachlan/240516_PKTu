# PK Tu analysis

## Pre-processing of SNP calls

The pre-processing of this dataset was executed via my Snakemake workflow:

1. Deploy the workflow

    - Clone smk-rnaseq-gatk-variants workflow from GitHub repository

        ```
        cd 240516_PKTu
        git clone git@github.com:baerlachlan/smk-rnaseq-gatk-variants.git
        mv smk-rnaseq-gatk-variants/ smk-variants/
        cd smk-variants/
        ## I used the following commit
        git reset --hard 3f4d2c460c22edba463dd2078603577f2d714f4b
        ```

    - Remove unneeded files and folders

        ```
        rm -rf .git .gitignore LICENSE README.md .snakemake-workflow-catalog.yml
        ```

1. Add the following snippet to the bottom of `workflow/Snakefile`

    ```python
    localrules:
        genome_get,
        annotation_get,
        known_variants_get,
    ```

    - These rules require internet access, but minimal compute, therefore run on the head node.

2. Create workflow-specific profile configuration for HPC execution in `workflow/profile/default`

   - This works in conjunction with a global profile configuration in my home directory. The repository is [here](https://github.com/baerlachlan/smk-cluster-generic-slurm).

   - The global profile contains the configuration for slurm execution in a HPC environment, while the workflow profile contains the configuration for workflow-specific resource requirements.

3. Create the directory `results/raw_data/fastq/` and copy raw FASTQ data into it.

   - `md5.txt` also added.

4. Modify the configurations as required in the following files:

    ```bash
    config/config.yaml
    config/samples.tsv
    config/units.tsv
    ```

5. Execute the workflow

    - Activate pre-made conda environment containing `snakemake v8.14.0`

    ```bash
    snakemake
    ```

    - This workflow executes using both the global and workflow-specific profile.
      - Global profile is defined using the `SNAKEMAKE_PROFILE` environmental variable in `~/.bashrc`, and Snakemake automatically locates the workflow-specific profile.
