# handedness-ewas

> **DNA methylation signatures of left-handedness for different tissues and ages**
> by Veronika V. Odintsova, Matthew Suderman, Fiona A. Hagenbeek, Doretta Caramaschi, Jouke-Jan Hottenga, René Pool, 
> BIOS consortium,
> Conor V. Dolan, Lannie Ligthart, Catharina E.M. van Beijsterveld, Gonneke Willemsen, Eco J.C. de Geus, Jeffrey J. Beck, Erik A.
> Ehli, Gabriel Cuellar-Partida, David M. Evans, Sara E. Medland, Caroline L. Relton, Dorret I. Boomsma, Jenny van Dongen

We report the first epigenome-wide association study (EWAS) of left-handedness, a trait with low heritability for which epigenetic mechanisms have been proposed as an underlying etiological mechanism. 
* The **primary EWAS** of handedness was performed in two cohorts with DNA methylation data in whole blood (Illumina, 450k): 
    1. NTR adults (N=2,682 individuals including twins, mean age at methylation 36.5, standard deviation (SD) 12.7)
        *scripts*: [ntr/ewas_basic model_NTRadults.r] and [ntr/ewas_adjustrf model_NTRadults.r]
    2. ALSPAC adults (N=1,232, mean age at methylation 48.98, SD 5.55)
        *script*: [alspac/run-ewas.r]
* Summary statistics were combined in a meta-analysis (N=3,914) testing associations at 409,563 CpGs using METAL
    *script*: [combined/metal.txt] 
* We tested whether EWAS signal was enriched near loci detected in the previous GWAS on handedness (Cuellar-Partida et al, 2021). We carried out within-pair twin analysis in MZ twins discordant for handedness (Nadults= 133 twin pairs, Nchildren= 86 twin pairs). 
    *script*: [combined/gwas_followup.r]
* Secondary analyses were performed in different tissues
    1. In cord blood and peripheral blood in ALSPAC children (N=1,021 with DNA methylation data at birth, at 7, 17 years old, Illumina 450k chip, and/or at 24 years old, Illumina EPIC array)
        *script*: [alspac/run-ewas.r] 
    2. In buccal cells in NTR children (N=946 twins, mean age 9.5, SD 1.85, Illumina EPIC array)
* We then examined correlations between the effect sizes of top CpGs (defined as CpGs with the lowest p-value) in each analysis.
    *script*: [combined/correlations.r] 
* Finally, we created and tested polygenic and DNA methylation scores for left-handedness. 
    *scripts*: [alspac/run-scores.r], [ntr/methylation_scores_NTR.r], [ntr/explained_variance_NTR.r]


To cite: Odintsova V., Suderman M., Hageenbeek F., Caramaschi D., Hottenga J.J., Pool R., BIOS consortium, Dolan C.V., Ligthart L., van Beijsterveld C.E.M., Willemsen G., de Geus E.J.C., Beck J.J., Ehli E.A., Cuellar-Partida G., Evans D.M., Medland S.E., Relton C.L., Boomsma D.I., van Dongen J. DNA-methylation signatures of left-handedness for different tissues and ages (under revision)

