> **DNA methylation signatures of left-handedness**<br>
> by Veronika V. Odintsova, Matthew Suderman, Fiona A. Hagenbeek, Doretta Caramaschi, Jouke-Jan Hottenga, René Pool, 
> BIOS consortium,
> Conor V. Dolan, Lannie Ligthart, Catharina E.M. van Beijsterveld, Gonneke Willemsen, Eco J.C. de Geus, Jeffrey J. Beck, Erik A.
> Ehli, Gabriel Cuellar-Partida, David M. Evans, Sara E. Medland, Caroline L. Relton, Dorret I. Boomsma, Jenny van Dongen

We report the first epigenome-wide association study (EWAS) of left-handedness, a trait with low heritability for which epigenetic mechanisms have been proposed as an underlying etiological mechanism. 
* The primary EWAS of handedness was performed in two cohorts with DNA methylation data in whole blood (Illumina, 450k): 
    1. [NTR](https://tweelingenregister.vu.nl/research) adults (N=2,682 individuals including twins, mean age at methylation 36.5, standard deviation (SD) 12.7)
        <br>`code`: [ntr/ewas_basic model_NTRadults.r](ntr/ewas_basic%20model_NTRadults.r) and [ntr/ewas_adjusted model_NTRadults.r](ntr/ewas_adjusted%20%model_NTRadults.r)
    2. [ALSPAC](http://www.bristol.ac.uk/alspac/) adults (N=1,232, mean age at methylation 48.98, SD 5.55)
        <br>`code`: [alspac/run-ewas.r](alspac/run-ewas.r)
* Summary statistics were combined in a meta-analysis (N=3,914) testing associations at 409,563 CpGs using METAL
    <br>`code`: [combined/metal.txt](combined/metal.txt) 
* We tested whether EWAS signal was enriched near loci detected in the previous GWAS on handedness [(Cuellar-Partida et al, 2021)](https://www.nature.com/articles/s41562-020-00956-y)
    <br>`code`: [combined/gwas_followup.r](combined/gwas_followup.r)
* We carried out within-pair twin analysis in MZ twins discordant for handedness (Nadults= 133 twin pairs, Nchildren= 86 twin pairs).
    <br>`code`: [ntr/ewas_discordant_MZ_adults.r](/ntr/ewas_discordant_MZ_adults.r), [/ntr/ewas_discordant_MZ_children.r](/ntr/ewas_discordant_MZ_children.r)
* Secondary EWAS were performed in different tissues collected from children:
    1. In cord blood and peripheral blood in [ALSPAC](http://www.bristol.ac.uk/alspac/) children (N=1,021 with DNA methylation data at birth, at 7, 17 years old, Illumina 450k chip, and/or at 24 years old, Illumina EPIC array)
        <br>`code`: [alspac/run-ewas.r](alspac/run-ewas.r) 
    2. In buccal cells in [NTR](https://tweelingenregister.vu.nl/research) children (N=946 twins, mean age 9.5, SD 1.85, Illumina EPIC array)
	    <br>`code`: [ntr/ewas_adjusted model_NTRchildren.r](ntr/ewas_adjusted%20model_NTRchildren.r)
* We then examined correlations between the effect sizes of top CpGs (defined as CpGs with the lowest p-value) in each analysis.
    <br>`code`: [combined/correlations.r](combined/correlations.r) 
* Finally, we created and tested polygenic and DNA methylation scores for left-handedness. 
    <br>`code`: [alspac/run-scores.r](alspac/run-scores.r), [ntr/methylation_scores_NTR.r](ntr/methylation_scores_NTR.r), [ntr/explained_variance_NTR.r](ntr/explained_variance_NTR.r)

*To cite:* Odintsova V., Suderman M., Hageenbeek F., Caramaschi D., Hottenga J.J., Pool R., BIOS consortium, Dolan C.V., Ligthart L., van Beijsterveld C.E.M., Willemsen G., de Geus E.J.C., Beck J.J., Ehli E.A., Cuellar-Partida G., Evans D.M., Medland S.E., Relton C.L., Boomsma D.I., van Dongen J. DNA-methylation signatures of left-handedness (under review)

