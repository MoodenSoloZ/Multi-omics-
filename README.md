# Multi-view-clustering-for-multi-omics-data-by-latent-subspace-learning
## Table of contents
* [Descriptions](#descriptions)
* [Preparations](#preparations)
* [Debug the code files](#debug-the-code-files)
## Descriptions
### 1.Realworld datasets
All realworld datasets are saved in ``` "/Simulation_data/..." ``` including the following seven multi-omics datasets:
* Bladder
* Brainlower
* Breast
* Kidneyrenalclearcell
* Pcpg
* Skcm
* Thyroid
### 2.Sythetic dataset
The simulation dataset is save in ```/Pre_processed_multi_omics_data/...``` including three views named as ```v1_2.csv,v2_2.csv,v3_2.csv```, respectively.
### 3.Code files descriptions
#### 1.main.R
Complete multi-omics datasets clusteirng and incomplete multi-omics datasets clustering analysis. The default test data is set as the ```Brainlower```.
#### 2.run_data_col.R
Incomplete multi-omics data generation and integration.
#### 3.specClust.R
The spetcral clusteing method applied in the latent subspace to get the results.
#### 4.survive_analyze.R
Using clinical data to do survive analysis
#### 5.simulation_data_result.R
The intgration performance evalution on the sythetic dataset.
#### 6.lost_ratio.R
Genarate the incomplete multi-omics data with different loss ratio.


## Preparations


## Debug the code files
To run this project, install it locally using npm:

```
$ cd ../lorem
$ npm install
$ npm start
```
