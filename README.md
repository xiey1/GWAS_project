# GWAS_project
Final project for Quantitative Genomics and Genetics course 

# Aim:
Find the positions of causal polymorphisms for the five expressed genes.
### The data:
<br>`phenotypes.csv` contains the phenotype data for 344 samples and 5 genes.
<br>`genotypes.csv` contains the SNP data for 344 samples and 50000 genotypes.
<br>`gene_info.csv` contains information about each gene that was measured.
<br>`SNP info.csv` contains the additional information on the genotypes and has four columns.

# 1. Check the distribution of the phenotype file
## 1.1 check the distribution of relative gene expression for each gene
<br>Here the phenotype data includes the relative gene expression data for 5 genes: *ERAP2*, *PEX6*, *FAHD1*, *GFM1* and *MARCH7*. The density plot for the relative expression for each gene is shown as below:

<img src='https://github.com/xiey1/GWAS_project/blob/master/images/Figure1_density.png' width=900px>

<br>The gene expression has been normalized and centered around 0. It is possible that a Z score has been calculated and used to represent the relative gene expression in the phenotype file. Based on the density plot, all the gene expression conforms to normal distribution. Considering the phenotype data is continuous and follows normal distribution, linear regression model may be an option for the GWAS.

## 1.2 Check the distribution of genotype conditional on population structure
<br>The population structure including populations and genders are considered as the potential confounding factor in the genetic data. Covariates may need to be included in the modeling to account for the population structure. The following table summarizes the number of samples for each population/gender involved in this study. Approximately equal number of samples is included in each population/ gender.

<table>
  <thead>
    <tr>
      <th>Population</th>
      <th>CEU</th>
      <th>FIN</th>
      <th>GBR</th>
      <th>TSI</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>FEMALE</td>
      <td>37</td>
      <td>55</td>
      <td>46</td>
      <td>43</td>  
    </tr>
    <tr>
      <td>MALE</td>
      <td>41</td>
      <td>34</td>
      <td>39</td>
      <td>49</td>
    </tr>
  </tbody>
</table>

<br>To assess the impact of population structure on the distribution of genotype, PCA plots were also generated.

<img src='https://github.com/xiey1/GWAS_project/blob/master/images/Figure2_PCA_geno.png' width=900px>
<br>There is no obvious difference in genotype distribution between males and females, but different population groups form distinct clusters. CEU and GBR have similar genotype distribution and therefore are clustered together, while TSI and FIN each forms a separate cluster. However, PC1 accounts for 0.98% of the variance and PC2 accounts for 0.52% of the variance, which altogether account for only 1.5% of the total genotype variance. These suggest that although populations with different ancestry have distinct genotype structures, these don’t significantly affect genotype variation. Therefore, population structure may not be an important factor for this study.

## 1.3 Check the distribution of the phenotype conditional on population structure
<br>It is possible that different populations have large variance in the expression levels of genes measured in this study due to non-genetic factors such as environmental factors, which may confound the GWAS study. Besides, the phenotype (gene expression) data is obtained from RNA-seq, and technical issues such as batch effect during sequencing may also largely affect the expression data, further confounding our modeling process. To assess the potential impact of the other external or hidden factors on the phenotype data, it’s better to look at the distribution of the expression profiles.

<br>Density plot was generated for each population and gender for each gene.

<img src='https://github.com/xiey1/GWAS_project/blob/master/images/Figure3_density_sub.png' width=900px>
<br>Although several population/gender groups show deviation in gene expression, there is no obvious trend that individuals from different population/gender groups exhibit significant difference in gene expression. 

<br>PCA plots were also generated for the expression data.

<img src='https://github.com/xiey1/GWAS_project/blob/master/images/Figure4_PCA_pheno.png' width=900px>
<br>All data points group together with no distinct small clusters. There is no obvious separation of data points for either populations or gender groups.

<br>These altogether suggest that the phenotype data obtained in this study conforms to normal distribution with no significant variation conditional on population structure or gender groups. Thus, linear regression model can first be applied to this GWAS study.

## 1.4 Check minor allele frequency
<br>Since power is low when minor allele frequency (MAF) is low, MAF is calculated and genotypes with MAF < 0.05 should be eliminated. The histogram for the distribution of MAF is plotted and all genotypes have MAF > 0.05.

<img src='https://github.com/xiey1/GWAS_project/blob/master/images/Figure5_MAF.png' width=600px>

# 2. Linear regression model
<br>Since the phenotype data for each gene expression profile approximately conforms to the normal distribution and the population structure doesn’t have a significant impact on the genotype distribution, linear regression model without modeling covariates is first applied to the phenotype and genotype data. Xa and Xd matrices were generated, F-statistics and p-values were computed for each SNP per phenotype. Heatmap is generated for the p-values.

<img src='https://github.com/xiey1/GWAS_project/blob/master/images/Heatmap.png' width=900px>
<br>Population structure and hidden factors don’t have obvious effects on the p-value heatmap.

## 2.1 Manhattan plot
<br>Manhattan plot was generated for p-values of each phenotype. Controlling a study-wide type I error of 0.05, using the Bonferroni correction and the p-value cutoff is calculated as α=0.05/N=0.05/50000=1e-06. Using this cutoff, gene1_ ERAP2 has 73 significant SNPs, gene2_ PEX6 has 29 significant SNPs, gene3_ FAHD1 has 90 significant SNPs, and gene4_ GFM1 and gene5_MARCH7 don’t have significant SNPs. The p-value cutoff is indicated as the red line in each Manhattan plot. 

<img src='https://github.com/xiey1/GWAS_project/blob/master/images/Figure7_Manhattan_plot.png' width=900px>

## 2.2 QQ-plot
<br>QQ-plot is generated for each phenotype. For all 5 genes, the QQ-plot conforms to line y=x, suggesting that the null hypothesis is true for the majority of the SNPs. This indicates that the linear regression model works well for this study, consistent with the previous analysis that population structure may not largely affect the analysis.

<img src='https://github.com/xiey1/GWAS_project/blob/master/images/Figure8_QQ-plot.png' width=900px>

# 3. Linear regression model and include covariates to model the population structure and gender
<br>To exclude the potential impact of different populations, an Xz matrix was created to model the population structure. 6 columns are included in the Xz matrix, with 4 columns each representing a different population (1 indicates in the corresponding population group otherwise 0); 2 columns representing the gender.

<br>Both Manhattan plot and QQ-plot were generated for each phenotype. As expected, they look very similar to the plots generated without adding covariates. Choosing the same p-value cutoff, gene1_ ERAP2 has 74 significant SNPs, gene2_ PEX6 has 27 significant SNPs, gene3_FAHD1 has 90 significant SNPs, and gene4_GFM1 and gene5_ MARCH7 don’t have significant SNPs. The number of significant SNPs is comparable to the previous analysis without including covariates. However, adding covariates does slightly improve the QQ-plot. Thus, I choose the results obtained from the linear regression including covariates modeling as the model to perform further downstream analysis and focus on gene ERAP2, PEX6 and FAHD1 with potential causal SNPs.

# 4. Further analyze the potential causal SNPs obtained from the linear regression model
## 4.1 Hardy-Weinberg equilibrium test
<br>In the most circumstances, Hardy-Weinberg equilibrium applies to the human genome. SNPs that fail the Hardy-Weinberg equilibrium test need to be removed. With p-value=0.05 as cutoff, 3436 SNPs fail the Hardy-Weinberg equilibrium test. After removing these SNPs, gene1_ ERAP2 has 67 significant SNPs, gene2_ PEX6 has 26 significant SNPs, gene3_FAHD1 has 82 significant SNPs.

## 4.2 Check the localization of the potential causal SNPs for each phenotype – Local Manhattan plot
<br>To better study the genomic location of potential causal SNPs and the distribution of p-values for surrounding SNPs, the local Manhattan plot is generated for each gene covering the genomic area containing SNPs with significant p-values and 100 SNPs upstream and downstream of the potential causal genomic region as shown here:

<img src='https://github.com/xiey1/GWAS_project/blob/master/images/Figure9_Manhattan_Plot_Local.png' width=900px>

## 4.3 LD heatmap
<br>To better study the linkage disequilibrium for each potential causal SNP, the squared pearson correlation coefficient is calculated for each SNP in the causal SNP region as well as 50 SNPs upstream and downstream of the causal SNP region. A heatmap has been generated for the calculated r^2  value to visualize the LD among the SNPs in each specific genomic region as shown here:

<img src='https://github.com/xiey1/GWAS_project/blob/master/images/Figure9_LD_heatmap.png' width=900px>

# 5. Generate boxplot to view expression data conditional on each potential causal SNP
<br>To validate the causal SNPs identified in this GWAS study, boxplot is generated to visualize the relative gene expression data grouped by genotype ‘AA’, ‘AB’, ‘BB’ for each potential causal SNP.

# 6. Data Interpretation and Discussion
## 6.1 Gene1_ERAP2
<br>Based on calculation, the potential 67 causal SNPs for gene1_ERAP2 are located on chromosome 5, spanning from 96772432 to 97110808, covering 338376 bp. Based on the LD heatmap, there are two LD areas with the 1st LD region (chr5: 96772432- 96848898) exhibiting relatively weak association, while the 2nd LD region (chr5: 96868550-97110808) shows very strong correlation. It is also obvious from the Manhattan plot that the first 20 SNPs have p-values only slightly above the cutoff, and the majority of the SNPs in the second LD region are all well above the cutoff. Consistently, SNPs in the 1st LD region generally exhibits minor difference in ERAP2 expression (eg. rs26491) for different genotypes, while SNPs in the 2nd LD region (eg. snp rs1046395) show significant changes in gene expression affected by different genotypes as shown here:

<img src='https://github.com/xiey1/GWAS_project/blob/master/images/Figure11_Boxplot.png' width=900px>

<br>These altogether suggest that the potential causal SNP may be more likely localized in the 2nd LD region, which includes the 21th SNP (rs2549778) to the 67th SNP (rs72777610) (chr5: 96868550-97110808). ERAP2 is located on chr5: 96875939-96919702, and the identified potential SNP region covers ERAP2, suggesting the identified SNPs are “cis”-eQTL.

## 6.2 Gene2_PEX6
<br>The 26 potential causal SNPs identified for gene2_PEX6 are located on two chromosomes: ‘rs6854915’ on chr4 and the remaining 25 SNPs on chr6. Based on the local Manhattan plot, SNP ‘rs6854915’ has p-value only slightly above the cutoff and is the only SNP showing statistical significant on chr4 for PEX6. The squared pearson correlation coefficient calculated for ‘rs6854915’ and the other 25 SNPs on chr6 doesn’t show strong correlation, and boxplot of gene expression profile suggests that different genotypes for SNP ‘rs6854915’ only show mild difference in PEX6 expression. These altogether suggest that SNP ‘rs6854915’ might be a false positive target and only the 25 SNPs located on chr6 are further analyzed.

<br>Based on calculation, the potential 25 causal SNPs for gene2_PEX6 located on chromosome 6 span from 42889467 to 43108015, covering 218548 bp. There are two LD areas based on the LD heatmap with high correlation scores. The 1st LD region contains 13 SNPs (chr6: 42889467-42946612) with relatively low –log10(p-value). The 2nd LD region contains 11 SNPs (chr6: 42949278-42977844) and the 25th SNP is excluded from the 2nd LD region due to low –log10(p-value) and weak correlation with other SNPs. Detailed expression analysis also suggests that SNPs in the 1st LD region don’t have a significant impact on PEX6 expression (eg. rs2092554), while SNPs in the 2nd LD region has a larger impact on gene expression of PEX6 (eg. rs1129187).

<br>These altogether suggest that the potential causal SNP may be more likely localized in the 2nd LD region, which includes the 11 SNPs (rs4714638, rs9471976, rs58497441, rs13215983, rs6920547, rs7760250, rs1129187, rs3805953, rs10948061, rs9471985, rs199921136) (chr6: 42949278-42977844). PEX6 is located on chr6: 42963872- 42979242, and the identified potential SNP region partially overlaps with PEX6, suggesting the identified SNPs are “cis”-eQTL.

## 6.3 Gene3_FAHD1
<br>Based on calculation, the potential 82 causal SNPs for gene3_FAHD1 are located on chromosome 16, spanning from 1524250 to 1929366, covering 405116 bp. LD heatmap suggests all the SNPs are located in one LD region and the majority of the SNPs exhibit significant difference in gene expression of FAHD1 for different genotypes (eg. rs1657139). FAHD1 is located on chr16: 1827223-1840206. The identified potential SNP region (chr16: 1524250-1929366) covers FAHD1, suggesting the identified SNPs are “cis”-eQTL.
