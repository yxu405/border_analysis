# DiMeLo-Seq CENPA Analysis Pipeline

This repository contains Python scripts for analyzing DiMeLo-Seq data, with a focus on CENPA regions. The scripts process BAM and BED files, extract modified base information, calculate mod densities, and generate heatmaps for visualization.

## Table of Contents

- Requirements
- Installation
- Usage
  - BAM File Loading
  - Reference Genome Parsing
  - BED File Filtering
  - Mod Subset Calculation
  - Mod Density Calculation
  - Heatmap Plotting
- Output
- License

## Requirements

- Python 3.10+
- Required packages:
  - `pysam`
  - `numpy`
  - `biopython`
  - `pandas`
  - `matplotlib`
  - `scikit-learn`
  - `tabulate`

## Installation

1. Clone this repository:
    ```bash
    git clone <repository-url>
    cd <repository-directory>
    ```

2. Create and activate a virtual environment:
    ```bash
    python -m venv env
    source env/bin/activate  # For Linux/Mac
    env\Scripts\activate  # For Windows
    ```

3. Install the dependencies:
    ```bash
    pip install -r requirements.txt
    ```

## Usage

### BAM File Loading

The script loads a BAM file using the `pysam` library to extract aligned reads. Modified base information is fetched, and initial parsing is done to obtain insertion/deletion information for further processing.

```python
import pysam

bamfile = pysam.AlignmentFile("<path_to_bamfile>", "rb")
for read in bamfile.fetch():
    mod = read.modified_bases_forward
    # Analyze modified bases here
```

### Reference Genome Parsing

The reference genome is loaded from a FASTA file and parsed into a dictionary of chromosome sequences.

```python
from Bio import SeqIO

fasta_sequences = SeqIO.parse("<path_to_fasta>", "fasta")
assembly = {fasta.id: str(fasta.seq) for fasta in fasta_sequences}
```

### BED File Filtering

BED files are read into a Pandas DataFrame to extract and filter regions of interest, such as CDR and CDR_Transition regions.

```python
import pandas as pd

df = pd.read_csv("<path_to_bedfile>", sep="\t", names=columns)
selected_df = df[df['name'].isin(['CDR', 'CDR_Transition'])]
selected_df.to_csv("<output_bedfile>", sep="\t", index=False)
```

### Mod Subset Calculation

The function `mod_subset_producing_step()` isolates the desired regions (e.g., without dashes) within the modified bases numpy array.

```python
def mod_subset_producing_step(mod_no_dash, alignment_dash, target_start_no_dash, target_end_no_dash):
    # Implementation details here
    return mod_subset
```

### Mod Density Calculation

The `single_mod_density_calculator()` function calculates mod densities for specific genomic regions. It divides the region into intervals and computes both first and second derivatives for each read.

```python
def single_mod_density_calculator(bamfile, transition_coordinate, mod_tag, interval_num):
    # Implementation details here
    return read_first_derivatives, read_second_derivatives, read_id_list, forward_reverse_status, read_density_values
```

### Heatmap Plotting

The script generates heatmaps to visualize mod densities and their derivatives across genomic intervals.

- **Density Heatmap**: Plots densities of CpG and mA modifications.
- **First Derivative Heatmap**: Highlights changes in mod densities within specific intervals.

```python
def plot_heatmap(data1, data2, y_labels, target_region, derivative_mC, derivative_mA):
    # Implementation details here
    plt.show()
```

#### Example Usage

```python
plot_heatmap(read_density_values, read_density_values_mA, read_id_list, target_region, read_derivatives, read_derivatives_mA)
```

## Output

- **Filtered BED File**: Contains regions of interest (e.g., CDR and CDR_Transition).
- **Heatmaps**: Visual representations of mod densities and their derivatives, saved as visual outputs.

