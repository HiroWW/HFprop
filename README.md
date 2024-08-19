# HFprop

# !!!! UNDER CONSTRUCTION !!!!
HFprop is a comprehensive software tool designed for both the analysis and design of propellers using Blade Element Momentum Theory (BEMT). This tool integrates two powerful modules: one for analyzing the performance of a propeller with a given geometry, and another for designing a propeller geometry that meets specified performance requirements. HFprop is aimed at engineers, researchers, and enthusiasts who need accurate and efficient propeller performance predictions and designs.

## Features

- **BEMT-Based Performance Analysis**: Analyze the aerodynamic performance of existing propeller geometries using the well-established Blade Element Momentum Theory.
- **Propeller Design Module**: Design propeller geometries to meet desired performance criteria, ensuring optimal efficiency and thrust.
- **User-Friendly Interface**: Intuitive and easy-to-use interface for both modules, making HFprop accessible to users of all levels.
- **Flexible and Extensible**: HFprop can be easily integrated into larger workflows and extended to suit specific use cases.

## Getting Started

### Prerequisites

- Python 3.8 or later
- Required Python libraries:
  - NumPy
  - SciPy
  - Matplotlib (for visualization)

You can install the required libraries using pip:

```bash
pip install numpy scipy matplotlib
```

### Installation

Clone the repository to your local machine:

```bash
git clone https://github.com/HiroWW/HFprop.git
cd HFprop
```

### Usage

#### 1. Performance Analysis

To analyze the performance of a propeller with a given geometry, use the `perform_analysis.py` script:

```bash
python perform_analysis.py --input geometry_file.txt --output results_file.txt
```

Options:
- `--input`: Path to the input file containing the propeller geometry data.
- `--output`: Path to the output file where the analysis results will be saved.

#### 2. Propeller Design

To design a propeller geometry that meets specific performance requirements, use the `design_propeller.py` script:

```bash
python design_propeller.py --requirements requirements_file.txt --output geometry_file.txt
```

Options:
- `--requirements`: Path to the file containing the performance requirements.
- `--output`: Path to the output file where the designed geometry will be saved.

## Contributing

We welcome contributions to HFprop! If you find a bug, have a feature request, or want to contribute code, please feel free to submit an issue or a pull request.

## License

HFprop is released under the MIT License. See the `LICENSE` file for more details.

## Contact

If you have any questions or need support, please contact me.
