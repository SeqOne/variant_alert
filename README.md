# Variant Alert! framework

--------------------------------------------------------------------------------

![logo](img/variant-alert-logo.png)

## Overview

Collaborative variant interpretation databases such as ClinVar are updated weekly, resulting in changes that can impact diagnostic performance in every genetic centre worldwide. These frequent updates present at least two challenges: 
- (i) Promptly integrating the most up-to-date versions to ensure the patient receives the most relevant diagnostic
- (ii) Finding an effective way to re-interpret previous analyses when changes in databases alter interpretation results.

To do so, we provide Variant Alert!, a framework to monitor every significant alteration in variant classification and gene-disease association between two versions of ClinVar database.

**Variant Alert provide :**

- *compare-variant* : A tabulated or VCF file of variant with significant change in classification status
- *compare-gene* : A tabulated file of significant change in gene-disease association
- *clinvarome* : A list of genes with at least one pathogenic variant in ClinVar

If you want to cite this work :
> *Yauy et al.* Extension of the mendelian genome and continuous variant reinterpretation through ClinVar monitoring (2020).

## Installation

**Requirements**

```bash
python > 3.6
poetry
```

**Install**

```bash
poetry install
```

**Input recommended**

We recommend to use a VCF processed from Clinvcf, as VCF from ClinVar repository is not complete.

Clinvcf extract all variants from the monthly XML release of ClinVar into a VCF format (https://github.com/SeqOne/clinvcf).

## Run

```
poetry shell
```

```
Usage: variant-alert [OPTIONS] VCF_SOURCE_PATH VCF_TARGET_PATH COMMAND1
                     [ARGS]... [COMMAND2 [ARGS]...]...

  CLI for variant_alert package.

  VCF_SOURCE_PATH: Used as a reference VCF_TARGET_PATH: Compared to the
  source vcf

Options:
  -o, --output_directory PATH  Output directory of variant comparison file.
                               Default current working directory
  -f, -- output_format  tsv|vcf  Format of output file, either tsv or vcf for compare-variant function. Default = tsv
  --help                       Show this message and exit.

Commands:
  clinvarome       Output all pathogenic genes of target clinvar vcf
  compare-gene     Output genes with change of pathogenicity and/or review confidence status
  compare-variant  Output variants with change of pathogenicity
```

## Glossary

### Significant changes  

| Previous                                                               	| New                                                                    	| Clinical_impact 	|
|------------------------------------------------------------------------	|------------------------------------------------------------------------	|-----------------	|
| ..                                                                     	| ..                                                                     	| null            	|
| ..                                                                     	| Benign and/or likely_benign                                            	| null            	|
| ..                                                                     	| Uncertain_significance or Conflicting_interpretations_of_pathogenicity 	| unknown         	|
| ..                                                                     	| Likely_pathogenic                                                      	| major           	|
| ..                                                                     	| Pathogenic                                                             	| major           	|
| Benign and/or likely_benign                                            	| ..                                                                     	| warning         	|
| Benign and/or likely_benign                                            	| Uncertain_significance or Conflicting_interpretations_of_pathogenicity 	| unknown         	|
| Benign and/or likely_benign                                            	| Likely_pathogenic                                                      	| major           	|
| Benign and/or likely_benign                                            	| Pathogenic                                                             	| major           	|
| Uncertain_significance or Conflicting_interpretations_of_pathogenicity 	| ..                                                                     	| warning         	|
| Uncertain_significance or Conflicting_interpretations_of_pathogenicity 	| Benign and/or likely_benign                                            	| minor           	|
| Uncertain_significance or Conflicting_interpretations_of_pathogenicity 	| Likely_pathogenic                                                      	| major           	|
| Uncertain_significance or Conflicting_interpretations_of_pathogenicity 	| Pathogenic                                                             	| major           	|
| Likely_pathogenic or Pathogenic/Likely_pathogenic                      	| ..                                                                     	| warning         	|
| Likely_pathogenic or Pathogenic/Likely_pathogenic                      	| Benign and/or likely_benign                                            	| major           	|
| Likely_pathogenic or Pathogenic/Likely_pathogenic                      	| Uncertain_significance or Conflicting_interpretations_of_pathogenicity 	| major           	|
| Likely_pathogenic or Pathogenic/Likely_pathogenic                      	| Pathogenic                                                             	| minor           	|
| Pathogenic                                                             	| ..                                                                     	| warning         	|
| Pathogenic                                                             	| Benign and/or likely_benign                                            	| major           	|
| Pathogenic                                                             	| Uncertain_significance or Conflicting_interpretations_of_pathogenicity 	| major           	|
| Pathogenic                                                             	| Likely_pathogenic                                                      	| minor           	|


## Troubleshooting

- *Clinvarome*: you will also need to give 2 vcf in arguments for the clinvarome function, even if it will only use the lastest.
- *vcf output*: tsv output will take only few seconds, but print a VCF output could take 10 minutes. This option is only available for the compare-variant function.

--------------------------------------------------------------------------------
*Variant Alert! is a collaboration of :* 

[![SeqOne](img/logo-seqone.png)](https://seq.one/) 

[![Université Grenoble Alpes](img/logo-uga.png)](https://iab.univ-grenoble-alpes.fr/) 

[![CHU de Rouen](img/logo-CHU.png)](https://www.chu-rouen.fr/service/service-de-genetique/)
