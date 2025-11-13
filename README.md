# Variant Calling Analysis  

**Author:** Luis Aguilera

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

### Reproducing: 
**Circulating tumor DNA sequencing in colorectal cancer patients treated with first-line chemotherapy with anti-EGFR.**


# Installation

## Installation on a Local Computer

To install this repository and all its dependencies, we recommend using [Anaconda](https://www.anaconda.com).

* Clone the repository:
```sh
git clone --depth 1 https://github.com/luisub/circulatory_tumor_DNA_analysis.git
```

* Create a virtual environment and activate it:
```sh
conda create -n vca_env python=3.12 -y
conda activate vca_env
```

* Install the required dependencies:
```sh
pip install -r requirements.txt
```

To deactivate or remove the environment:

* Deactivate the environment:
```sh
conda deactivate
```

* Remove the environment:
```sh
conda env remove -n vca_env -y
```
