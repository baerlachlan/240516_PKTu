# PK Tu analysis

## Pre-processing for gene-level DE analysis

The pre-processing of this dataset was executed via my Snakemake workflow:

1. Deploy the workflow

    - Download [smk-rnaseq-gatk-variants](https://github.com/baerlachlan/smk-rnaseq-star-featurecounts) `v1.1.2`

        ```bash
        wget https://github.com/baerlachlan/smk-rnaseq-star-featurecounts/archive/refs/tags/v1.1.2.zip
        unzip v1.1.2.zip && rm v1.1.2.zip
        cd smk-rnaseq-star-featurecounts-1.1.2
        ```

    - Remove uneeded files and directories

        ```bash
        rm -rf .gitignore .snakemake-workflow-catalog.yml .test LICENSE
        ```

1. Create workflow-specific profile configuration for HPC execution at `workflow/profiles/default/config.v8+.yaml`

   - This works in conjunction with a global profile configuration in my home directory. The repository is [here](https://github.com/baerlachlan/smk-cluster-generic-slurm).

   - The global profile contains the configuration for slurm execution in a HPC environment, while the workflow profile contains the configuration for workflow-specific resource requirements.

1. Create the directory `results/raw_data/fastq/` and copy raw FASTQ data into it.

   - `md5.txt` also added.

1. Modify the configurations as required in the following files:

    ```bash
    config/config.yaml
    config/samples.tsv
    config/units.tsv
    ```

1. Execute the workflow

    - Activate pre-made conda environment containing `snakemake v8.16.0`

    ```bash
    snakemake
    ```
    
    - This workflow executes using both the global and workflow-specific profile.
      - Global profile is added as an environmental variable in `~/.bashrc`, and Snakemake automatically locates the workflow-specific profile.
