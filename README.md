# Seq2Feature

Seq2Feature is a comprehensive pipeline designed for annotating plasmid sequences and generating machine learning features using PLAnnotate. This tool automates feature extraction and annotation to support genomic research and machine learning model development.

This project is primarily inspired by [pLannotate](https://github.com/mmcguffi/pLannotate/) from [Barrick Lab](https://barricklab.org/twiki/bin/view/Lab), which provides the foundation for sequence annotation used in Seq2Feature. We extend and build upon this work to offer additional features and automation.

## Features

- **Automated Annotation**: Uses databases such as GenoLIB, FPbase, and Swiss PROT to annotate plasmid sequences.
- **Machine Learning Integration**: Extracts features from sequences to facilitate the development and training of ML models.
- **Scalable Processing**: Annotates large datasets efficiently with progress tracking.

## Installation

Seq2Feature can be installed using Conda or Docker. Choose the method that best fits your needs.

### Using Conda

1. **Clone the Repository**

   ```bash
   git clone https://github.com/yourusername/Seq2Feature.git
   cd Seq2Feature
   ```

2. **Create and Activate the Conda Environment**

   ```bash
   conda env create -f environment.yml
   conda activate seq2feature
   ```

### Using Docker

1. **Build the Docker Image**

   ```bash
   docker build -t seq2feature .
   ```

2. **Run the Docker Container**

   ```bash
   docker run -it --rm seq2feature
   ```

   By default, this command will run the `seqannotate.main` module. Adjust the Docker run command if you need to execute a different script or pass additional arguments.

## Usage

To use Seq2Feature, you need to run the annotation script which utilizes databases for feature extraction.

### Script Overview

- `Seq2Feature/seqannotate/main.py`:
  Runs the annotation script using databases such as GenoLIB, FPbase, and Swiss PROT. For more details, refer to [this paper](https://academic.oup.com/nar/article/49/W1/W516/6279845).

- `Seq2Feature/seqannotate/resources.py`:
  Contains annotated resources used for feature extraction and annotation.

- `Seq2Feature/gene_main.py`:
  Demonstrates how to use the `seqannotate` package to process and annotate gene data.

### Example

To process and annotate your gene data, follow these steps:

1. Update `read_loc` and `save_loc` in `Seq2Feature/gene_main.py` with the paths to your input CSV file and where you want to save the annotated files.

2. Run the script:

   ```bash
   python Seq2Feature/gene_main.py
   ```

   This script reads your input CSV, annotates each sequence, and saves the results in the specified location.

## Configuration

- **`read_loc`**: Path to the CSV file containing the sequences to be annotated.
- **`save_loc`**: Directory where annotated files will be saved.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for any enhancements or bug fixes.

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0). See the [LICENSE](LICENSE) file for details.

## Acknowledgements

This project is inspired by [pLannotate](https://github.com/mmcguffi/pLannotate/). Special thanks to the authors of pLannotate for their foundational work in sequence annotation.

## Contact

For questions or support, please contact [nikhilgnanavel@yahoo.com](mailto:nikhilgnanavel@yahoo.com).
