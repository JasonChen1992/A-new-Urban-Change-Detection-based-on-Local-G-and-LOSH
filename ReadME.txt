# Urban Change Detection using Local G and LOSH

This repository contains the implementation of a new method for detecting urban changes based on Local G and Local Spatial Heteroscedasticity (LOSH) statistics. The method was proposed in the paper "A new urban change detection method based on the local G and local spatial heteroscedasticity statistics" by Yuzhou Chen and Ran Tao.

## Abstract
Accurate detection of urban changes is critical for guiding city planning that leads to smart and sustainable development. Few existing methods can detect urban changes in a timely manner while preserving great details of spatial heterogeneity. In this project, we propose a new method that combines two spatial statistics, namely the Local G and Local Spatial Heteroscedasticity (LOSH). By jointly analyzing the results of both statistics, we design the Urban Development Index (UDI) to assess the types of urban changes that each spatial unit has been experiencing. The experiments with both a synthetic dataset and a Rwanda population dataset demonstrate that our method can identify completed and ongoing phenomena of urban transition and unveil the heterogeneous nature of growth and/or shrinkage inside a city.

## Methodology

### Local G Statistic
The Local G statistic is used to identify clusters of high or low values in spatial data.


### Local Spatial Heteroscedasticity (LOSH)
LOSH is used to detect local heterogeneity in spatial data.


### Urban Development Index (UDI)
UDI combines the results of Local G and LOSH statistics to assess urban development. It categorizes spatial units into different stages of urban development based on the following criteria:
- If G is low and LOSH is low, UDI = 1
- If G is low and LOSH is not significant, UDI = 2
- If G is low and LOSH is high, UDI = 3
- If G is not significant and LOSH is high, UDI = 4
- If G is not significant and LOSH is not significant, UDI = 5
- If G is not significant and LOSH is low, UDI = 6
- If G is high and LOSH is high, UDI = 7
- If G is high and LOSH is not significant, UDI = 8
- If G is high and LOSH is low, UDI = 9

## Code Implementation
The code implementation is divided into two parts: calculating LOSH and calculating Local G. After calculating these values, we determine the UDI values for each spatial unit.

### Calculate LOSH
The LOSH calculation is implemented in the `calculateLOSH` function, which computes the LOSH value for a given key list and data dictionary. The significance of the LOSH values is determined using Monte Carlo permutation tests.

### Calculate Local G
The Local G calculation is implemented in the `calculateGetisG` function, which computes the Local G* value for a given key list and data dictionary. The significance of the Local G values is also determined using Monte Carlo permutation tests.
To reveal the absolute distributional dynamics instead of the relative ones, we adopt the resolution proposed by Tao and Chen (2022), which modifies the significance test by using
the first year of observation as the benchmark when performing the local G analysis on spatial panel data.
Tao, R., & Chen, Y. (2022). Applying local indicators of spatial association to analyze longitudinal data: The absolute perspective.
Geographical Analysis. https://doi.org/10.1111/gean.12323

### Determine UDI
The UDI values are determined based on the criteria mentioned above, using the calculated LOSH and Local G values and their significance.

## How to Use
1. Clone the repository to your local machine.
2. Ensure you have the necessary dependencies installed: `clusterpy`, `numpy`, `shapefile`, and `pandas`.
3. Place your input data files in the appropriate directory.
4. Run the scripts to calculate LOSH and Local G, and then determine the UDI values.
5. The results will be saved in a CSV file.

## Data
You can find the national data of Rwanda and city-level data of Kigali in the corresponding folders in Data

## Results
The results include the calculated LOSH and Local G values, their significance, and the UDI values for each spatial unit. The results are saved in a CSV file for further analysis.

## License
This project is licensed under the MIT License.

## Contact
For any questions or issues, please contact Yuzhou Chen at yuzhouchen@usf.edu or Ran Tao at rtao@usf.edu.
