# Optimization of a WIsH model 

[WIsH](https://github.com/soedinglab/wish) (Who Is the Host) is a software tool that builds a model for predicting the hosts of phages. 
The present scripts are to find the optimal model parameters such that the predictions are the most accurate for the particular dataset.

## Hypothetical experimental setup

In the present setting, there is a set of phages (e.g. phages from a metagenomic sample) for which we would like to predict the hosts. Using WIsH (Who Is the Host) is advisable for such analysis in case of short contigs (which is often the case in NGS shotgun viromics data). 

## Parameters to optimize

- Model order
- P-value cut off
- Evaluation method of prediction results (Best prediction, Majority vote of X best predictions, LCA of X best predictions)

For more information on the parameters, please see the [publication (and supplementary material)](https://academic.oup.com/bioinformatics/article/33/19/3113/3964377) on WIsH.

## Quick start

After installing R and RStudio (for version, see 'Software' bellow), clone the repository. <br />
From the master script, all other scripts can be run swiftly and easily. <br />
Keep in mind that some parts of the code may need to be taylored to the needs of your own dataset. 

## Software

The scripts were written using: <br />
R version 3.4.0 <br />
RStudio Version 1.0.143

## License

GPLv3 License - see the [LICENSE](https://github.com/eregenyi/WIsH-model-optimization/blob/master/LICENSE) file for details.

## To do
- construct dummy datasets for reproducibility

