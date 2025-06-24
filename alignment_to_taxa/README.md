# Alignment_to_taxa

### **1. General Description**

The script is a bioinformatics tool designed to resolve taxonomic classifications for sequences based on alignment data. It takes alignment results (BAM file), a table of sequence abundances, and a reference taxonomy map as input.

Its main purpose is to address the common challenge where a single sequence (or read) may align well to multiple different reference sequences (e.g., from closely related species). The script resolves these ambiguities using a Lowest Common Ancestor (LCA) algorithm.

Finally, it generates a new, clean abundance table where counts are aggregated by their final taxonomic assignment. To ensure stable and concise identifiers, each unique taxonomy string is given a permanent, 10-character MD5 hash ID.

### **2. Dependencies**

The script requires the following Python libraries:

* **pysam**: For reading and parsing the input BAM file.
* **pandas**: Used for handling and structuring data.

These can be installed using pip:

```{bash}
pip install pysam pandas
```

### **3. Core Functionality & Workflow**

The script follows a logical pipeline from raw alignments to a final, taxonomically-aggregated abundance table.

1.  **Argument Parsing**: The script begins by parsing required inputs (`--bam`, `--abundance`, `--taxonomy`) and optional parameters that control its behavior (`--score-threshold`, `--lca-threshold`, etc.).

2.  **Data Loading**:
    * **Taxonomy File**: It loads the reference taxonomy into a dictionary, mapping each reference sequence ID to its full taxonomic string.
    * **Abundance File**: It loads the per-sample counts for each unique sequence ID.

3.  **BAM File Processing**: The script iterates through the BAM file to build a dictionary where each `query_name` maps to a list of its alignment `(reference_name, alignment_score)`.

4.  **Taxonomy Assignment**: For each sequence, it filters alignments based on the `--score-threshold` relative to the best hit.
    * **One valid hit**: The sequence is directly assigned the taxonomy of that single reference.
    * **Multiple valid hits**: The script triggers the **Lowest Common Ancestor (LCA)** algorithm, which finds the deepest taxonomic rank shared by a consensus (`--lca-threshold`) of the hits.

5.  **Hash ID Generation**: To create a stable primary key, the script generates a 10-character MD5 hash (e.g., `a1b2c3d4e5`) for each unique taxonomic string produced.

6.  **Abundance Aggregation**: The script creates a new abundance table by summing the counts of all original sequences that resolve to the same final hash ID.

7.  **Output Generation**: The script writes the final results into four separate files.

### **4. Comprehensive Description of How Results are Handled**

The script produces a clean, verifiable, and easy-to-use set of outputs.

#### **Output Files**

Assuming no `--prefix` is used, the script generates the following files:

1.  **`abundance_table.tsv` - The Final Abundance Table**
    * **Purpose**: This is the main quantitative output. It provides the abundance of each final taxonomy across all samples.
    * **Format**: A tab-separated value (TSV) file.
    * **Columns**: `TaxID` (the 10-character hash), `Sample1`, `Sample2`, ...
    * **Rows**: Each row represents a unique taxonomic entity, and the values are the **sum of counts** from all original sequences assigned to that `TaxID`.

2.  **`taxonomy_table.tsv` - The Taxonomy Key**
    * **Purpose**: This file is the legend for `abundance_table.tsv`, mapping the hash ID back to a human-readable classification.
    * **Format**: A TSV file.
    * **Columns**: `TaxID`, `Kingdom`, `Phylum`, `Class`, `Order`, `Family`, `Genus`, `Species`.
    * **Rows**: Each row defines the full taxonomy for a given `TaxID`. Unassigned ranks are filled with "NA".

3.  **`sequence_hash_mapping.tsv` - The Traceability Map**
    * **Purpose**: This file provides full transparency, allowing a user to trace which original sequences were aggregated into which final `TaxID`.
    * **Format**: A two-column TSV file.
    * **Columns**: `Sequence`, `TaxID`.
    * **Rows**: A simple mapping of every successfully assigned sequence to its new taxonomic group ID.

4.  **`taxonomy_assignment_log.txt` - The Log File**
    * **Purpose**: Records the settings and summary statistics for the run, including parameters, assignment rates, and lookup failures.

### **5. Usage Example**

Below is an example of how to run the script from the command line.

```{bash}
./alignment_to_taxa.py \
    --bam processed_data/alignments.bam \
    --abundance counts/sequence_abundances.tsv \
    --taxonomy database/reference_taxonomy.tsv \
    --output-dir results/taxonomy_run_1 \
    --prefix my_project \
    --score-threshold 0.98 \
    --lca-threshold 0.80 \
    --verbose
```

This command would produce the following files inside the `results/taxonomy_run_1/` directory:
* `my_project_abundance_table.tsv`
* `my_project_taxonomy_table.tsv`
* `my_project_sequence_hash_mapping.tsv`
* `my_project_taxonomy_assignment_log.txt`

## License

This project is licensed under the GNU Lesser General Public License v3.0. See the `LICENSE` file for the full license text.
