# Vsearch Dereplication Script

The script is designed to dereplicate nucleotide sequences from multiple FASTQ files using the `vsearch` tool. It processes each sample, performs dereplication, and then aggregates the results into a single feature table and a FASTA file of unique sequences.

The core logic of this script is a direct adaptation of the dereplication functionality found in the `q2-vsearch` QIIME 2 plugin, refactored to operate as a standalone tool independent of the QIIME 2 framework. This allows it to produce standard, universally readable FASTA and TSV files.

## License

This project is licensed under the BSD 3-Clause License. See the accompanying `LICENSE` file for the full text.

## Citation

Since this script is a wrapper for `vsearch` and its logic is derived from a QIIME 2 plugin, please cite the original publications for these tools in your work:

* **VSEARCH**:
    Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahé, F. (2016). VSEARCH: a versatile open source tool for metagenomics. *PeerJ*, 4, e2584. <https://doi.org/10.7717/peerj.2584>

* **QIIME 2**:
    Bolyen, E., Rideout, J. R., Dillon, M. R., et al. (2019). Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. *Nature Biotechnology*, 37(8), 852–857. <https://doi.org/10.1038/s41587-019-0209-9>

## 1. Overview

The primary purpose of this script is to identify and count unique sequences from a collection of gzipped FASTQ files, where each file represents a single sample. It automates the process of:

1.  Converting each FASTQ file to a FASTA file.
2.  Running `vsearch` to dereplicate sequences within each sample.
3.  Parsing the `vsearch` output to determine the abundance of each unique sequence per sample.
4.  Aggregating the results from all samples into a master abundance table.
5.  Generating a single FASTA file containing all unique sequences across all samples, identified by their SHA1 hash.

This workflow is common in metabarcoding and metagenomics studies for creating a feature table (e.g., an ASV or OTU table) from raw sequencing data.

## 2. Dependencies

To run this script, you will need the following dependencies installed and available in your system's PATH.

### External Tools

* **vsearch**: A versatile open-source tool for sequence analysis. This script specifically uses its dereplication capabilities. It must be installed and accessible from the command line.

### Python Libraries

* **pandas**: Used for creating and managing the final abundance table.
* **Biopython**: Used for parsing FASTQ and FASTA files efficiently.

You can install the required Python libraries using pip:

```{bash, eval=FALSE}
pip install pandas biopython
```

## 3. Usage

The script is executed from the command line and accepts several arguments to control its behavior.

### Command-Line Arguments

| Argument            | Description                                                                                              | Required | Default |
|---------------------|----------------------------------------------------------------------------------------------------------|:--------:|:-------:|
| `--input_dir`       | Path to the directory containing your input `.fastq.gz` files.                                           |   Yes    |   N/A   |
| `--output_fasta`    | Path for the output FASTA file. This file will contain all unique dereplicated sequences.                |   Yes    |   N/A   |
| `--output_table`    | Path for the output TSV (tab-separated values) table, showing the abundance of each sequence per sample. |   Yes    |   N/A   |
| `--derep_prefix`    | If set, performs dereplication based on sequence prefixes. If not set, uses full-length dereplication.  |    No    | `False` |
| `--min_seq_length`  | The minimum length a sequence must have to be included in the analysis.                                  |    No    |    `1`    |
| `--min_unique_size` | The minimum number of times a unique sequence must appear to be considered a representative.             |    No    |    `1`    |
| `--threads`         | The number of parallel processes to use for processing samples.                                          |    No    |    `1`    |

### Example Command

```{bash, eval=FALSE}
python3 your_script_name.py \
    --input_dir ./my_fastq_files/ \
    --output_fasta unique_sequences.fasta \
    --output_table abundance_table.tsv \
    --min_seq_length 100 \
    --threads 8
```

## 4. Input and Output

### Input Files

The script expects a directory containing one or more FASTQ files compressed with `gzip` (i.e., ending in `.fastq.gz`). Each file is treated as a separate sample. The sample name is derived from the filename by removing the `.fastq.gz` extension. For example, `sample1.fastq.gz` will correspond to the sample ID `sample1`.

### Output Files

1.  **Output FASTA (`--output_fasta`)**: A standard FASTA file containing all unique sequences found across all input samples. The header for each sequence is its SHA1 hash, which ensures a unique and consistent identifier.

    *Example FASTA entry:*

    ```
    >0a1b2c3d4e5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b
    GATTACAGATTACAGATTACAGATTACAGATTACAGATTACA...
    ```

2.  **Output Table (`--output_table`)**: A tab-separated file representing the feature table.
    * The first column, `sequence_id`, contains the SHA1 hash corresponding to the sequences in the output FASTA file.
    * Subsequent columns correspond to each sample ID.
    * The values in the table are the counts (abundances) of each sequence in each sample.

    *Example Table:*

    ```
    sequence_id	sample1	sample2	sample3
    0a1b2c3d...	105	234	0
    1a2b3c4d...	50	0	112
    ...
    ```

## 5. Script Workflow

The script executes the following steps:

1.  **Argument Parsing**: Reads the command-line arguments provided by the user.
2.  **File Discovery**: Scans the `--input_dir` for all `.fastq.gz` files.
3.  **Parallel Processing**: It creates a temporary directory and processes each FASTQ file in parallel using a `ProcessPoolExecutor`. For each file:
    a.  **FASTQ to FASTA Conversion**: The `.fastq.gz` file is read, and its contents are converted into a FASTA file. The sample ID is prepended to each sequence header to track its origin (e.g., `>sample1_read1`).
    b.  **Vsearch Dereplication**: `vsearch` is called to perform dereplication on the generated FASTA file. This step produces a dereplicated FASTA file (containing only unique sequences for that sample) and a `.uc` file that maps reads to the unique sequences.
4.  **Result Aggregation**: After all samples are processed, the main process aggregates the results:
    a.  **Parse `.uc` Files**: For each sample, the `.uc` file is parsed to count how many original reads belong to each unique sequence.
    b.  **Build Global Data Structures**:
        * A dictionary `sequence_data` is created to store the actual sequence string for each unique sequence, keyed by its SHA1 hash.
        * A nested dictionary `abundance_data` is built to store the abundance of each sequence (`seq_id`) in each sample (`sample_id`).
5.  **Generate Outputs**:
    a.  **Write FASTA File**: The `sequence_data` dictionary is used to write the final dereplicated FASTA file.
    b.  **Write Abundance Table**: The `abundance_data` is converted into a pandas DataFrame and then written to the output TSV file. The table is structured with sequence IDs as rows and sample IDs as columns.
6.  **Cleanup**: The temporary directory and all intermediate files are automatically removed upon script completion.

## 6. Function Descriptions

* **`parse_args()`**: Defines and parses command-line arguments using `argparse`.
* **`convert_fastq_to_fasta(fastq_file, output_dir)`**: Reads a gzipped FASTQ file, converts records to FASTA format, and prepends the sample ID to headers.
* **`process_sample(fastq_file, temp_dir, args)`**: The main worker function for a single sample. It calls `convert_fastq_to_fasta` and then runs `vsearch`. It returns the sample ID and the paths to the dereplicated FASTA and `.uc` files.
* **`parse_uc_file(uc_file, sample_id)`**: Parses the `vsearch` `.uc` output file to count the number of reads clustered into each unique sequence. It handles both hit (`H`) and seed (`S`) lines.
* **`generate_sha1(sequence)`**: Computes the SHA1 hash of a given sequence string to create a unique identifier.
* **`main()`**: The main function that orchestrates the entire workflow, from setting up parallel processing to aggregating results and writing the final output files.
