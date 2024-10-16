# Generate What You Can Make: Achieving in-house synthesizability with readily available resources in de novo drug design

## Link
[ChemRxiv (not online yet)](https://github.com/AlanHassen/led3score)

## Abstract
Molecules generated by Computer-Aided Drug Design often lack synthesizability to be valuable because Computer-Aided Synthesis Planning (CASP) and CASP-based approximated synthesizability scores have rarely been used as generation objectives, despite facilitating the in-silico generation of synthesizable molecules. Published scores approximate a general notion of CASP-based synthesizability with nearly unlimited building block resources. However, this approach is disconnected from the reality of small laboratory drug design, where building block resources are limited, making a notion of in-house synthesizability that uses already available resources highly desirable.

In this work, we show a successful *de novo* drug design workflow generating active and in-house synthesizable ligands of monoglyceride lipase (MGLL). We demonstrate the successful transfer of CASP from 17.4 million commercial building blocks to a small laboratory setting of roughly 6,000 building blocks with only a decrease of –12% in CASP success. Moreover, we present a rapidly retrainable in-house synthesizability score, successfully capturing our in-house synthesizability without relying on external building block resources. We show that including our in-house synthesizability score in a multi-objective *de novo* drug design workflow, alongside a simple QSAR model, provides thousands of potentially active and easily in-house synthesizable molecules. Further, we highlight differences between general and in-house synthesizability scores and demonstrate potential problems with the out-of-distribution predictive performance of synthesizability scores on generated molecules. Finally, we experimentally evaluate the synthesis and biochemical activity of three *de novo* candidates using their CASP-suggested synthesis routes using only in-house building blocks. We find one candidate with evident activity, suggesting potential new ligand ideas for MGLL inhibitors while showcasing the usefulness of our in-house synthesizability score.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Data](#data)
- [Citing](#citing)

## Installation

```bash
# Example installation command
make create_conda_env
conda activate led3_score
make install_packages
```

## Usage
Find the jupyter notebooks/scripts in the respective folders.

## Data

The models and data associated with this paper will be made available upon its publication or upon request.

## Citing
If you use our work in your research, please cite:

```bibtex
Will be updated...
}
```

## Funding

This study was partially funded by the European Union’s Horizon 2020 research and innovation program under the Marie Skłodowska-Curie Innovative Training Network European Industrial Doctorate grant agreement No. 956832 “Advanced machine learning for Innovative Drug Discovery”
