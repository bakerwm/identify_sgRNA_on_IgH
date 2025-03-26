# Identify sgRNAs targeting chr21 of human genome

## Criteria
+ sgRNA is 23mer and ending with GG
+ sgRNA targets should be in tandom repeat (copy number > 10)
+ sgRNA targeting chr21 of human genome
+ sgRNA MUST not targeting **MOUSE** genome (off-target)

## Pipeline Overview

This pipeline performs the following steps:

1. Download and prepare the human reference genome (fasta, gtf) (GRCh38.chr21)
2. Detect tandom repeats using Tandem Repeats Finder (TRF)
3. Filter repeats based on criteria (length > 13bp, copy number > 10)
4. Extract candidate sgRNA loci, 23 bp in length and ending with GG (SpCas9)
5. Remove off-targets by genome-wide alignment to MOUSE genome
6. Generate summary reports and visualizations

## Requirements

The pipeline requires the following tools and packages:

- Python 3.6+
  - pybedtools
  - biopython
  - NumPy
  - matplotlib
- bowtie2
- samtools
- bedtools
- trf (Tandem Repeats Finder)

## Installation

### Conda Environment

The recommended way to install all dependencies is using conda:

```bash
conda env create -f find_sgrna.yml
conda activate find_sgrna
```

### Manual Installation

If you prefer to install dependencies manually:

1. Install Python packages:
   ```bash
   pip install biopython pybedtools numpy matplotlib
   ```

2. Install genome tools:
   ```bash
   # For Ubuntu/Debian
   sudo apt-get install bowtie2 samtools bedtools
   
   # For CentOS/RHEL
   sudo yum install bowtie2 samtools bedtools
   
   # For macOS
   brew install bowtie2 samtools bedtools
   ```

3. Install Tandem Repeats Finder (TRF):
   Download from https://tandem.bu.edu/trf/trf.html and add to your PATH.

## Usage

### Quick Start

To run the complete pipeline:

```bash
bash workflow.sh config.sh
```

### Configuration

The `config.sh` file contains all configurable parameters:

```bash
# System configuration
export N_CPU=12                   # Number of CPU cores to use

# Genome configuration
export MOUSE_BUILD="GRCm38"       # Mouse genome build
export HUMAN_BUILD="GRCh38"       # Human genome build 
export RELEASE=102                # Ensembl release number

# Analysis parameters
export MIN_COPY=10                # Minimum copy number
export MIN_LENGTH=13              # Minimum repeat length
```

### Individual Components

The pipeline is composed of multiple scripts that can be run independently:

- `01.download_genome.sh`: Download and prepare the human reference genome and mouse genome
- `03.find_repeats.sh`: Find repeat sequences in human chromosome 21
- `05.filter_repeats.py`: Filter repeats based on criteria
- `06.identify_sgrna.py`: Identify all possible sgRNAs within consensus sequence
- `07.remove_off_targets.sh`: Remove sgRNAs that map to mouse genome
- `08.generate_report.py`: Generate a summary report of the results

## Output

The pipeline generates the following outputs in the `results` directory:

- `repeats/`: Directory containing raw outputs from repeat-finding tools
- `repeats/all_repeats.filtered.bed`: BED6+1 file of tandom repeats
- `sgRNA/`: Directory containing sgRNA sequences and files
- `sgRNA/sgrna_raw.on_target.txt`: BED6+4 file of on target sgRNAs
- `report/`: Directory containing summary reports and visualizations
- `workflow_YYYYMMDD_HHMMSS.log`: Log file with timestamped entries

### Output Format

**sgRNA on target**

File: `results/sgRNA/sgrna_raw.on_target.txt`
Columns:
1. chromosome  
2. start  
3. end  
4. tandom_repeat_name  
5. copy_number  
6. strand  
7. consensus_sequence  
8. direction_of_sgRNA:fwd/rev  
9. sgRNA_id  
10. sgRNA_sequence

## Customization

You can modify the following parameters in the configuration file:

- `MOUSE_BUILD`: Mouse reference genome (e.g., GRCm38)
- `HUMAN_BUILD`: Human reference genome (e.g., GRCh38)
- `RELEASE`: Ensembl release number (e.g., 102)
- `MIN_LENGTH`: Minimum length of repeats
- `MIN_COPY`: Minimum copy number of repeats
- `N_CPU`: Number of CPUs to use

## Advanced Usage

### Resume Interrupted Runs

The pipeline is designed to automatically skip completed steps. If the workflow is interrupted, you can simply rerun the command and it will continue from the last successful step.

### Parallel Processing

The pipeline utilizes parallel processing for computationally intensive steps, particularly for genome alignment. You can adjust the number of CPU cores used by modifying the `N_CPU` parameter in the configuration file.

## Troubleshooting

### Common Issues

1. **Missing Dependencies**
   - **Symptom**: Error message about missing command
   - **Solution**: Ensure you've activated the conda environment with `conda activate find_sgrna`

2. **Memory Issues**
   - **Symptom**: Process killed or out of memory errors
   - **Solution**: Reduce `N_CPU` value or run on a machine with more RAM

3. **TRF Failure**
   - **Symptom**: "TRF Failed?" message in logs
   - **Solution**: Verify TRF is properly installed and executable

4. **No sgRNAs Found**
   - **Symptom**: Empty output file
   - **Solution**: Check if your target region is correct or try reducing the filtering criteria (`MIN_LENGTH`, `MIN_COPY`)

5. **Bowtie2 Index Error**
   - **Symptom**: "Could not locate a Bowtie index" error
   - **Solution**: Ensure the path to the Bowtie2 index is correct in your configuration

### Logging

The pipeline now generates detailed logs with timestamps for each step. If you encounter issues, check the log file in the `results` directory:

```bash
less results/workflow_YYYYMMDD_HHMMSS.log
```

## License

This software is provided under the MIT License.

## Citation

If you use this pipeline in your research, please cite:

`https://github.com/bakerwm/identify_sgRNA_on_IgH`

## Version History

- **v1.1** - Enhanced error handling, added resumable execution, improved logging
- **v1.0** - Initial release
