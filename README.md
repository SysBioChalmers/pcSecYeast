# pcSecYeast: The proteome constrained genome scale secretory model of _Saccharomyces cerevisiae_


* Brief Model Description:

This repository contains the current consensus consensus proteome constrained genome scale secretory model of _Saccharomyces cerevisiae_, which expands the metabolic model with protein translation, post-translational modifiaction and secretion process.

* Model KeyWords:

**GEM Category:** species; **Utilisation:** experimental data reconstruction, multi-omics integrative analysis, _in silico_ strain design; **Field:** metabolic-network reconstruction; **Type of Model:** reconstruction, curated; **Taxonomy:** _Saccharomyces cerevisiae_; 

* Last update: 2021-07-20

* Main Model Descriptors:

|Taxonomy | Template Model |
|:-------:|:--------------:|
|_Saccharomyces cerevisiae_|[Yeast 8.3.5](https://github.com/SysBioChalmers/yeast-GEM/blob/master/ModelFiles/xml/yeastGEM.xml)

This repository is administered by Feiran Li ([@feiranl](https://github.com/feiranl)), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology.


## Citation

- Currently, please cite the preprint paper "[Genome scale modeling of the protein secretory pathway reveals novel targets for improved recombinant protein production in yeast](https://doi.org/10.1101/2021.10.16.464630)"


## Installation

### Required Software - User:

* Matlab user:
  * A functional Matlab installation (MATLAB 7.3 or higher).
  * The [COBRA toolbox for MATLAB](https://github.com/opencobra/cobratoolbox).
  * solver [SoPlex](https://soplex.zib.de).



## Usage

For generate the figures in this paper, please download the intemediate file from the [Zenode](https://zenodo.org/record/5593654#.YXMoai8RrmF), and run the correponding function in [ComplementaryScripts/Simulation](https://github.com/SysBioChalmers/pcSecYeast/tree/main/ComplementaryScripts/Simulation)
For generatre the pcSecYeast model, run [buildModel](https://github.com/SysBioChalmers/pcSecYeast/tree/main/ComplementaryScripts).
