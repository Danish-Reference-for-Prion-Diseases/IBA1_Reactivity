# IBA1_Protein_Morph_Reactivity

Code and data for article: "Combined neocortical protein and morphological profiling of reactive microglia across Alzheimer’s and Creutzfeldt-Jakob disease".

Raw data can be found here:
GeoMx protein source data ->  /Data/GeoMx (protein)/1. IBA1 Raw Geomx Dataset.xlsx
Morphometrical Data extracted from CellProfiler -> /Data/Cellprofiler (morphology)/

If not done beforehand, please extract ALL contents of the .zip file before running R-scripts or R-markdown files through Rstudio using the R-project.

Please install missing packages/dependencies before running each R-script (using default install.packages() function in R):
Note, some packages are from Bioconductor that will need to be installed through it (https://www.bioconductor.org/install/). The author could at least identify these packages from Bioconductor: 
- ComplexHeatmap (https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
- variancePartition (https://www.bioconductor.org/packages/release/bioc/html/variancePartition.html)
- MOFA2 (https://www.bioconductor.org/packages/release/bioc/html/MOFA2.html)


Figure_2 employs a package ("JLutils") that needs to fetched remotely from GitHub. Please make sure Git is installed on your system (https://git-scm.com/downloads). Then install "JLutils" using the "remotes" package (see https://rdrr.io/github/larmarange/JLutils/ or try install.packages("remotes") and then remotes::install_github("larmarange/JLutils")).


The R-scripts utilize relative pathfinding with the "here" package and should be able to run with the correct corresponding data. 

The recommended order of executed R-scripts to see the complete data analysis and results:

0. QC_Segments_Targets.R (utilizing Nanostring's own Geoscript for Segment and Target QC Evaluation - optional)
1. GeoMx Evaluation and Clean-up.R
2. Figure_1+2.R
3. Cellprofiler Clean-up, normalization and feature selection.R
4. Figure_3.R
(5. Add_ID for MoBIE.R)
6. SupplementaryFigure_6
7. Figure_7
   
Baseline Tables.R is standalone and can be executed whenever.

To inspect the CellProfiler pipeline (found at /Scripts/Cellprofiler – IBA1 Morphology Extraction.cpproj), CellProfiler needs to be installed (https://cellprofiler.org/). Version 4.2.6 was used to create this particular pipeline.


**Session infos:**

  **QC_Segments_Targets**

     R version 4.5.0 (2025-04-11 ucrt)
     Platform: x86_64-w64-mingw32/x64
     Running under: Windows 11 x64 (build 22631)

    Matrix products: default
    LAPACK version 3.12.1

    locale:
    [1] LC_COLLATE=Danish_Denmark.utf8  LC_CTYPE=Danish_Denmark.utf8    LC_MONETARY=Danish_Denmark.utf8
    [4] LC_NUMERIC=C                    LC_TIME=Danish_Denmark.utf8    

    time zone: Europe/Copenhagen
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
    [1] here_1.0.1

    loaded via a namespace (and not attached):
    [1] compiler_4.5.0     rprojroot_2.1.1    tools_4.5.0        rstudioapi_0.17.1  performance_0.15.1
    [6] insight_1.4.2     
  

  **GeoMx evaluation and Clean-up**
  
     R version 4.5.0 (2025-04-11 ucrt)
     Platform: x86_64-w64-mingw32/x64
     Running under: Windows 11 x64 (build 22631)

    Matrix products: default
    LAPACK version 3.12.1

    locale:
    [1] LC_COLLATE=Danish_Denmark.utf8  LC_CTYPE=Danish_Denmark.utf8    LC_MONETARY=Danish_Denmark.utf8
    [4] LC_NUMERIC=C                    LC_TIME=Danish_Denmark.utf8    

    time zone: Europe/Copenhagen
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
    [1] ggh4x_0.3.1       emmeans_1.11.2-8  lme4_1.1-37       Matrix_1.7-4      TidyDensity_1.5.1 gridExtra_2.3    
    [7] readxl_1.4.5      lubridate_1.9.4   forcats_1.0.0     stringr_1.5.1     dplyr_1.1.4       purrr_1.1.0      
    [13] readr_2.1.5       tidyr_1.3.1       tibble_3.3.0      ggplot2_3.5.2     tidyverse_2.0.0   here_1.0.1       

    loaded via a namespace (and not attached):
    [1] generics_0.1.4     stringi_1.8.7      lattice_0.22-7     hms_1.1.3          magrittr_2.0.3    
    [6] estimability_1.5.1 grid_4.5.0         timechange_0.3.0   RColorBrewer_1.1-3 mvtnorm_1.3-3     
    [11] cellranger_1.1.0   rprojroot_2.1.1    scales_1.4.0       Rdpack_2.6.4       reformulas_0.4.1  
    [16] cli_3.6.5          rlang_1.1.6        rbibutils_2.3      performance_0.15.1 splines_4.5.0     
    [21] withr_3.0.2        tools_4.5.0        tzdb_0.5.0         nloptr_2.2.1       minqa_1.2.8       
    [26] boot_1.3-32        vctrs_0.6.5        R6_2.6.1           lifecycle_1.0.4    MASS_7.3-65       
    [31] insight_1.4.2      pkgconfig_2.0.3    pillar_1.11.0      gtable_0.3.6       glue_1.8.0        
    [36] data.table_1.17.8  Rcpp_1.1.0         tidyselect_1.2.1   rstudioapi_0.17.1  xtable_1.8-4      
    [41] farver_2.1.2       nlme_3.1-168       compiler_4.5.0    


  **Figure_1+2:**

     R version 4.5.0 (2025-04-11 ucrt)
     Platform: x86_64-w64-mingw32/x64
     Running under: Windows 11 x64 (build 22631)

    Matrix products: default
      LAPACK version 3.12.1

    locale:
    [1] LC_COLLATE=Danish_Denmark.utf8  LC_CTYPE=Danish_Denmark.utf8    LC_MONETARY=Danish_Denmark.utf8
    [4] LC_NUMERIC=C                    LC_TIME=Danish_Denmark.utf8    

    time zone: Europe/Copenhagen
    tzcode source: internal

    attached base packages:
    [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] openxlsx_4.2.8           kableExtra_1.4.0         glmmSeq_0.5.5            ggh4x_0.3.1             
     [5] ggrepel_0.9.6            performance_0.15.1       emmeans_1.11.2-8         lmerTest_3.1-3          
     [9] lme4_1.1-37              Matrix_1.7-4             variancePartition_1.39.1 BiocParallel_1.42.1     
     [13] limma_3.64.1             gridExtra_2.3            GGally_2.4.0             circlize_0.4.16         
     [17] RColorBrewer_1.1-3       ComplexHeatmap_2.25.2    plotly_4.11.0            factoextra_1.0.7        
     [21] lubridate_1.9.4          forcats_1.0.0            stringr_1.5.1            dplyr_1.1.4             
     [25] purrr_1.1.0              readr_2.1.5              tidyr_1.3.1              tibble_3.3.0            
     [29] ggplot2_3.5.2            tidyverse_2.0.0          here_1.0.1              

    loaded via a namespace (and not attached):
      [1] rstudioapi_0.17.1   jsonlite_2.0.0      shape_1.4.6.1       magrittr_2.0.3      estimability_1.5.1 
      [6] rmarkdown_2.29      farver_2.1.2        nloptr_2.2.1        GlobalOptions_0.1.2 vctrs_0.6.5        
      [11] minqa_1.2.8         rstatix_0.7.2       progress_1.2.3      htmltools_0.5.8.1   broom_1.0.9        
      [16] Formula_1.2-5       sass_0.4.10         bslib_0.9.0         KernSmooth_2.23-26  htmlwidgets_1.6.4  
      [21] pbkrtest_0.5.5      plyr_1.8.9          sandwich_3.1-1      cachem_1.1.0        zoo_1.8-14         
      [26] TMB_1.9.17          lifecycle_1.0.4     iterators_1.0.14    pkgconfig_2.0.3     R6_2.6.1           
      [31] fastmap_1.2.0       rbibutils_2.3       clue_0.3-66         digest_0.6.37       numDeriv_2016.8-1.1
      [36] colorspace_2.1-1    S4Vectors_0.46.0    rprojroot_2.1.1     crosstalk_1.2.2     textshaping_1.0.2  
      [41] ggpubr_0.6.1        labeling_0.4.3      timechange_0.3.0    httr_1.4.7          abind_1.4-8        
      [46] mgcv_1.9-3          compiler_4.5.0      aod_1.3.3           withr_3.0.2         doParallel_1.0.17  
      [51] S7_0.2.0            backports_1.5.0     carData_3.0-5       ggstats_0.10.0      gplots_3.2.0       
      [56] ggsignif_0.6.4      MASS_7.3-65         rjson_0.2.23        corpcor_1.6.10      gtools_3.9.5       
      [61] caTools_1.18.3      tools_4.5.0         zip_2.3.3           remaCor_0.0.20      glue_1.8.0         
      [66] nlme_3.1-168        cluster_2.1.8.1     reshape2_1.4.4      generics_0.1.4      glmmTMB_1.1.12     
      [71] gtable_0.3.6        tzdb_0.5.0          data.table_1.17.8   hms_1.1.3           xml2_1.4.0         
      [76] car_3.1-3           BiocGenerics_0.55.1 foreach_1.5.2       pillar_1.11.0       splines_4.5.0      
      [81] lattice_0.22-7      tidyselect_1.2.1    pbapply_1.7-4       knitr_1.50          reformulas_0.4.1   
      [86] IRanges_2.42.0      svglite_2.2.1       xfun_0.52           RhpcBLASctl_0.23-42 stats4_4.5.0       
      [91] Biobase_2.68.0      statmod_1.5.0       matrixStats_1.5.0   stringi_1.8.7       yaml_2.3.10        
      [96] lazyeval_0.2.2      boot_1.3-32         evaluate_1.0.5      codetools_0.2-20    qvalue_2.41.0      
      [101] cli_3.6.5           systemfonts_1.2.3   pbmcapply_1.5.1     xtable_1.8-4        Rdpack_2.6.4       
      [106] jquerylib_0.1.4     Rcpp_1.1.0          EnvStats_3.1.0      png_0.1-8           parallel_4.5.0     
      [111] prettyunits_1.2.0   bitops_1.0-9        viridisLite_0.4.2   mvtnorm_1.3-3       scales_1.4.0       
      [116] insight_1.4.2       crayon_1.5.3        fANCOVA_0.6-1       GetoptLong_1.0.5    rlang_1.1.6



  **Cellprofiler Clean-up, normalization and feature selection**
  
      R version 4.5.0 (2025-04-11 ucrt)
     Platform: x86_64-w64-mingw32/x64
     Running under: Windows 11 x64 (build 22631)

    Matrix products: default
      LAPACK version 3.12.1

    locale:
    [1] LC_COLLATE=Danish_Denmark.utf8  LC_CTYPE=Danish_Denmark.utf8    LC_MONETARY=Danish_Denmark.utf8
    [4] LC_NUMERIC=C                    LC_TIME=Danish_Denmark.utf8    

    time zone: Europe/Copenhagen
    tzcode source: internal

    attached base packages:
    [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] circlize_0.4.16       ComplexHeatmap_2.25.2 collapse_2.1.3        caret_7.0-1           lattice_0.22-7       
     [6] sjmisc_2.8.11         VIM_6.2.2             colorspace_2.1-1      lubridate_1.9.4       forcats_1.0.0        
     [11] stringr_1.5.1         dplyr_1.1.4           purrr_1.1.0           readr_2.1.5           tidyr_1.3.1          
     [16] tibble_3.3.0          ggplot2_3.5.2         tidyverse_2.0.0       here_1.0.1           

    loaded via a namespace (and not attached):
     [1] pROC_1.19.0.1        rlang_1.1.6          magrittr_2.0.3       clue_0.3-66          GetoptLong_1.0.5    
     [6] snakecase_0.11.1     matrixStats_1.5.0    e1071_1.7-16         compiler_4.5.0       png_0.1-8           
     [11] vctrs_0.6.5          reshape2_1.4.4       pkgconfig_2.0.3      shape_1.4.6.1        crayon_1.5.3        
     [16] prodlim_2025.04.28   tzdb_0.5.0           recipes_1.3.1        cluster_2.1.8.1      parallel_4.5.0      
     [21] R6_2.6.1             stringi_1.8.7        vcd_1.4-13           RColorBrewer_1.1-3   ranger_0.17.0       
     [26] parallelly_1.45.1    car_3.1-3            boot_1.3-32          rpart_4.1.24         lmtest_0.9-40       
     [31] Rcpp_1.1.0           iterators_1.0.14     future.apply_1.20.0  zoo_1.8-14           IRanges_2.42.0      
     [36] Matrix_1.7-4         splines_4.5.0        nnet_7.3-20          timechange_0.3.0     tidyselect_1.2.1    
     [41] rstudioapi_0.17.1    abind_1.4-8          timeDate_4041.110    doParallel_1.0.17    codetools_0.2-20    
     [46] sjlabelled_1.2.0     listenv_0.9.1        plyr_1.8.9           withr_3.0.2          future_1.67.0       
     [51] survival_3.8-3       proxy_0.4-27         pillar_1.11.0        carData_3.0-5        foreach_1.5.2       
     [56] stats4_4.5.0         insight_1.4.2        generics_0.1.4       rprojroot_2.1.1      sp_2.2-0            
     [61] hms_1.1.3            S4Vectors_0.46.0     scales_1.4.0         laeken_0.5.3         globals_0.18.0      
     [66] class_7.3-23         glue_1.8.0           tools_4.5.0          robustbase_0.99-4-1  data.table_1.17.8   
     [71] ModelMetrics_1.2.2.2 gower_1.0.2          ipred_0.9-15         nlme_3.1-168         Formula_1.2-5       
     [76] cli_3.6.5            lava_1.8.1           gtable_0.3.6         DEoptimR_1.1-4       digest_0.6.37       
     [81] BiocGenerics_0.55.1  rjson_0.2.23         farver_2.1.2         lifecycle_1.0.4      hardhat_1.4.2       
     [86] GlobalOptions_0.1.2  MASS_7.3-65          


  **Figure_3**
  
     R version 4.5.0 (2025-04-11 ucrt)
     Platform: x86_64-w64-mingw32/x64
     Running under: Windows 11 x64 (build 22631)

    Matrix products: default
      LAPACK version 3.12.1

    locale:
    [1] LC_COLLATE=Danish_Denmark.utf8  LC_CTYPE=Danish_Denmark.utf8    LC_MONETARY=Danish_Denmark.utf8 LC_NUMERIC=C                   
    [5] LC_TIME=Danish_Denmark.utf8    

    time zone: Europe/Copenhagen
    tzcode source: internal

    attached base packages:
    [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] scales_1.4.0          viridis_0.6.5         viridisLite_0.4.2     plotly_4.11.0         umap_0.2.10.0         caroline_0.9.9       
     [7] ggcorrplot_0.1.4.1    rmcorr_0.7.0          psych_2.5.6           emmeans_1.11.2-8      performance_0.15.1    glmmTMB_1.1.12       
     [13] lme4_1.1-37           Matrix_1.7-4          fmsb_0.7.6            JLutils_1.24.0        rstatix_0.7.2         circlize_0.4.16      
     [19] ComplexHeatmap_2.25.2 Boruta_9.0.0          caret_7.0-1           lattice_0.22-7        MASS_7.3-65           factoextra_1.0.7     
     [25] lubridate_1.9.4       forcats_1.0.0         stringr_1.5.1         dplyr_1.1.4           purrr_1.1.0           readr_2.1.5          
     [31] tidyr_1.3.1           tibble_3.3.0          ggplot2_3.5.2         tidyverse_2.0.0       here_1.0.1           

    loaded via a namespace (and not attached):
      [1] splines_4.5.0        later_1.4.4          datawizard_1.2.0     hardhat_1.4.2        pROC_1.19.0.1        rpart_4.1.24        
      [7] lifecycle_1.0.4      Rdpack_2.6.4         doParallel_1.0.17    rprojroot_2.1.1      globals_0.18.0       insight_1.4.2       
      [13] backports_1.5.0      magrittr_2.0.3       httpuv_1.6.16        askpass_1.2.1        reticulate_1.43.0    minqa_1.2.8         
      [19] RColorBrewer_1.1-3   abind_1.4-8          DHARMa_0.4.7         BiocGenerics_0.55.1  nnet_7.3-20          sandwich_3.1-1      
      [25] ipred_0.9-15         labelled_2.14.1      lava_1.8.1           IRanges_2.42.0       S4Vectors_0.46.0     ggrepel_0.9.6       
      [31] listenv_0.9.1        RSpectra_0.16-2      parallelly_1.45.1    codetools_0.2-20     tidyselect_1.2.1     shape_1.4.6.1       
      [37] farver_2.1.2         matrixStats_1.5.0    stats4_4.5.0         jsonlite_2.0.0       GetoptLong_1.0.5     e1071_1.7-16        
      [43] Formula_1.2-5        survival_3.8-3       iterators_1.0.14     foreach_1.5.2        tools_4.5.0          Rcpp_1.1.0          
      [49] glue_1.8.0           mnormt_2.1.1         prodlim_2025.04.28   gridExtra_2.3        mgcv_1.9-3           withr_3.0.2         
      [55] numDeriv_2016.8-1.1  fastmap_1.2.0        ggh4x_0.3.1          GGally_2.4.0         boot_1.3-32          openssl_2.3.3       
      [61] digest_0.6.37        timechange_0.3.0     R6_2.6.1             mime_0.13            estimability_1.5.1   colorspace_2.1-1    
      [67] see_0.11.0           utf8_1.2.6           generics_0.1.4       data.table_1.17.8    recipes_1.3.1        class_7.3-23        
      [73] httr_1.4.7           htmlwidgets_1.6.4    ggstats_0.10.0       ModelMetrics_1.2.2.2 pkgconfig_2.0.3      gtable_0.3.6        
      [79] timeDate_4041.110    S7_0.2.0             htmltools_0.5.8.1    carData_3.0-5        TMB_1.9.17           clue_0.3-66         
      [85] png_0.1-8            gower_1.0.2          reformulas_0.4.1     rstudioapi_0.17.1    tzdb_0.5.0           reshape2_1.4.4      
      [91] rjson_0.2.23         nlme_3.1-168         nloptr_2.2.1         proxy_0.4-27         zoo_1.8-14           GlobalOptions_0.1.2 
      [97] parallel_4.5.0       miniUI_0.1.2         pillar_1.11.0        vctrs_0.6.5          promises_1.3.3       car_3.1-3           
      [103] xtable_1.8-4         cluster_2.1.8.1      mvtnorm_1.3-3        cli_3.6.5            compiler_4.5.0       rlang_1.1.6         
      [109] crayon_1.5.3         future.apply_1.20.0  labeling_0.4.3       plyr_1.8.9           stringi_1.8.7        questionr_0.8.1     
      [115] lazyeval_0.2.2       bayestestR_0.17.0    patchwork_1.3.2      hms_1.1.3            future_1.67.0        shiny_1.11.1        
      [121] highr_0.11           haven_2.5.5          rbibutils_2.3        broom_1.0.9         
    

  **Add_ID for MoBIE**
  
     R version 4.5.0 (2025-04-11 ucrt)
     Platform: x86_64-w64-mingw32/x64
     Running under: Windows 11 x64 (build 22631)

    Matrix products: default
      LAPACK version 3.12.1

    locale:
    [1] LC_COLLATE=Danish_Denmark.utf8  LC_CTYPE=Danish_Denmark.utf8    LC_MONETARY=Danish_Denmark.utf8 LC_NUMERIC=C                   
    [5] LC_TIME=Danish_Denmark.utf8    

    time zone: Europe/Copenhagen
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] caroline_0.9.9    viridis_0.6.5     viridisLite_0.4.2 lubridate_1.9.4   forcats_1.0.0     stringr_1.5.1     dplyr_1.1.4       purrr_1.1.0      
     [9] readr_2.1.5       tidyr_1.3.1       tibble_3.3.0      ggplot2_3.5.2     tidyverse_2.0.0   here_1.0.1       

    loaded via a namespace (and not attached):
     [1] gtable_0.3.6       compiler_4.5.0     tidyselect_1.2.1   gridExtra_2.3      scales_1.4.0       R6_2.6.1           generics_0.1.4    
     [8] rprojroot_2.1.1    pillar_1.11.0      RColorBrewer_1.1-3 tzdb_0.5.0         rlang_1.1.6        stringi_1.8.7      timechange_0.3.0  
     [15] cli_3.6.5          withr_3.0.2        magrittr_2.0.3     grid_4.5.0         rstudioapi_0.17.1  hms_1.1.3          lifecycle_1.0.4   
     [22] vctrs_0.6.5        glue_1.8.0         farver_2.1.2       tools_4.5.0        pkgconfig_2.0.3    

              

  **SupplementaryFigure_6**
  
     R version 4.5.0 (2025-04-11 ucrt)
      Platform: x86_64-w64-mingw32/x64
     Running under: Windows 11 x64 (build 22631)

    Matrix products: default
      LAPACK version 3.12.1

    locale:
    [1] LC_COLLATE=Danish_Denmark.utf8  LC_CTYPE=Danish_Denmark.utf8    LC_MONETARY=Danish_Denmark.utf8 LC_NUMERIC=C                   
    [5] LC_TIME=Danish_Denmark.utf8    

    time zone: Europe/Copenhagen
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] cowplot_1.2.0            ggh4x_0.3.1              ggrepel_0.9.6            kableExtra_1.4.0         performance_0.15.1      
     [6] glmmSeq_0.5.5            emmeans_1.11.2-8         lme4_1.1-37              Matrix_1.7-4             variancePartition_1.39.1
     [11] BiocParallel_1.42.1      limma_3.64.1             gridExtra_2.3            GGally_2.4.0             circlize_0.4.16         
     [16] RColorBrewer_1.1-3       plotly_4.11.0            factoextra_1.0.7         caroline_0.9.9           viridis_0.6.5           
     [21] viridisLite_0.4.2        lubridate_1.9.4          forcats_1.0.0            stringr_1.5.1            dplyr_1.1.4              
     [26] purrr_1.1.0              readr_2.1.5              tidyr_1.3.1              tibble_3.3.0             ggplot2_3.5.2            
     [31] tidyverse_2.0.0          here_1.0.1              

    loaded via a namespace (and not attached):
      [1] rstudioapi_0.17.1   jsonlite_2.0.0      shape_1.4.6.1       datawizard_1.2.0    magrittr_2.0.3      estimability_1.5.1  farver_2.1.2       
      [8] nloptr_2.2.1        rmarkdown_2.29      GlobalOptions_0.1.2 vctrs_0.6.5         minqa_1.2.8         rstatix_0.7.2       progress_1.2.3     
      [15] htmltools_0.5.8.1   broom_1.0.9         Formula_1.2-5       sass_0.4.10         bslib_0.9.0         KernSmooth_2.23-26  htmlwidgets_1.6.4  
      [22] pbkrtest_0.5.5      plyr_1.8.9          sandwich_3.1-1      cachem_1.1.0        zoo_1.8-14          TMB_1.9.17          lifecycle_1.0.4    
      [29] iterators_1.0.14    pkgconfig_2.0.3     R6_2.6.1            fastmap_1.2.0       rbibutils_2.3       digest_0.6.37       numDeriv_2016.8-1.1
      [36] colorspace_2.1-1    patchwork_1.3.2     rprojroot_2.1.1     crosstalk_1.2.2     textshaping_1.0.2   ggpubr_0.6.1        labeling_0.4.3     
      [43] timechange_0.3.0    httr_1.4.7          abind_1.4-8         mgcv_1.9-3          compiler_4.5.0      aod_1.3.3           withr_3.0.2        
      [50] S7_0.2.0            backports_1.5.0     carData_3.0-5       ggstats_0.10.0      gplots_3.2.0        ggsignif_0.6.4      MASS_7.3-65        
      [57] corpcor_1.6.10      gtools_3.9.5        caTools_1.18.3      tools_4.5.0         remaCor_0.0.20      glue_1.8.0          nlme_3.1-168       
      [64] grid_4.5.0          reshape2_1.4.4      see_0.11.0          generics_0.1.4      glmmTMB_1.1.12      gtable_0.3.6        tzdb_0.5.0         
      [71] data.table_1.17.8   hms_1.1.3           xml2_1.4.0          car_3.1-3           BiocGenerics_0.55.1 pillar_1.11.0       splines_4.5.0      
      [78] lattice_0.22-7      tidyselect_1.2.1    pbapply_1.7-4       knitr_1.50          reformulas_0.4.1    svglite_2.2.1       RhpcBLASctl_0.23-42
      [85] xfun_0.52           Biobase_2.68.0      statmod_1.5.0       matrixStats_1.5.0   stringi_1.8.7       yaml_2.3.10         lazyeval_0.2.2     
      [92] boot_1.3-32         evaluate_1.0.5      codetools_0.2-20    qvalue_2.41.0       cli_3.6.5           xtable_1.8-4        pbmcapply_1.5.1    
      [99] systemfonts_1.2.3   Rdpack_2.6.4        jquerylib_0.1.4     Rcpp_1.1.0          EnvStats_3.1.0      parallel_4.5.0      prettyunits_1.2.0  
      [106] bayestestR_0.17.0   bitops_1.0-9        mvtnorm_1.3-3       lmerTest_3.1-3      scales_1.4.0        crayon_1.5.3        insight_1.4.2      
      [113] fANCOVA_0.6-1       rlang_1.1.6    


  **Figure_7**
  
     R version 4.5.0 (2025-04-11 ucrt)
     Platform: x86_64-w64-mingw32/x64
     Running under: Windows 11 x64 (build 22631)

    Matrix products: default
      LAPACK version 3.12.1

    locale:
    [1] LC_COLLATE=Danish_Denmark.utf8  LC_CTYPE=Danish_Denmark.utf8    LC_MONETARY=Danish_Denmark.utf8 LC_NUMERIC=C                   
    [5] LC_TIME=Danish_Denmark.utf8    

    time zone: Europe/Copenhagen
    tzcode source: internal

    attached base packages:
    [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] MOFA2_1.18.0          viridis_0.6.5         viridisLite_0.4.2     ComplexHeatmap_2.25.2 SCORPIUS_1.0.9        cowplot_1.2.0      
     [7] ggh4x_0.3.1           ggrepel_0.9.6         gridExtra_2.3         GGally_2.4.0          circlize_0.4.16       RColorBrewer_1.1-3   
     [13] plotly_4.11.0         factoextra_1.0.7      lubridate_1.9.4       forcats_1.0.0         stringr_1.5.1         dplyr_1.1.4          
     [19] purrr_1.1.0           readr_2.1.5           tidyr_1.3.1           tibble_3.3.0          ggplot2_3.5.2         tidyverse_2.0.0      
     [25] here_1.0.1           

    loaded via a namespace (and not attached):
      [1] rstudioapi_0.17.1     jsonlite_2.0.0        shape_1.4.6.1         magrittr_2.0.3        farver_2.1.2          corrplot_0.95        
      [7] GlobalOptions_0.1.2   vctrs_0.6.5           rstatix_0.7.2         htmltools_0.5.8.1     S4Arrays_1.8.1        broom_1.0.9          
      [13] Rhdf5lib_1.31.0       SparseArray_1.8.0     dynparam_1.0.2        Formula_1.2-5         rhdf5_2.53.4          htmlwidgets_1.6.4    
      [19] basilisk_1.21.5       princurve_2.1.6       desc_1.4.3            plyr_1.8.9            igraph_2.1.4          lifecycle_1.0.4      
      [25] iterators_1.0.14      pkgconfig_2.0.3       Matrix_1.7-4          R6_2.6.1              fastmap_1.2.0         babelwhale_1.2.0     
      [31] MatrixGenerics_1.21.0 clue_0.3-66           digest_0.6.37         colorspace_2.1-1      S4Vectors_0.46.0      ps_1.9.1              
      [37] rprojroot_2.1.1       irlba_2.3.5.1         ggpubr_0.6.1          filelock_1.0.3        labeling_0.4.3        timechange_0.3.0     
      [43] lmds_0.1.0            httr_1.4.7            abind_1.4-8           mgcv_1.9-3            compiler_4.5.0        remotes_2.5.0        
      [49] withr_3.0.2           doParallel_1.0.17     S7_0.2.0              backports_1.5.0       carData_3.0-5         ggstats_0.10.0       
      [55] HDF5Array_1.37.0      ggsignif_0.6.4        MASS_7.3-65           proxyC_0.5.2          DelayedArray_0.34.1   rjson_0.2.23         
      [61] tools_4.5.0           ranger_0.17.0         dynutils_1.0.11       glue_1.8.0            h5mread_1.1.1         nlme_3.1-168         
      [67] rhdf5filters_1.21.0   Rtsne_0.17            cluster_2.1.8.1       reshape2_1.4.4        generics_0.1.4        isoband_0.2.7        
      [73] dynwrap_1.2.4         gtable_0.3.6          tzdb_0.5.0            data.table_1.17.8     hms_1.1.3             utf8_1.2.6           
      [79] carrier_0.2.0         car_3.1-3             XVector_0.48.0        BiocGenerics_0.55.1   RANN_2.6.2            foreach_1.5.2        
      [85] pillar_1.11.0         splines_4.5.0         lattice_0.22-7        tidyselect_1.2.1      pbapply_1.7-4         IRanges_2.42.0       
      [91] stats4_4.5.0          matrixStats_1.5.0     pheatmap_1.0.13       stringi_1.8.7         lazyeval_0.2.2        yaml_2.3.10          
      [97] codetools_0.2-20      cli_3.6.5             uwot_0.2.3            reticulate_1.43.0     processx_3.8.6        Rcpp_1.1.0           
      [103] dir.expiry_1.17.0     png_0.1-8             parallel_4.5.0        assertthat_0.2.1      mclust_6.1.1          scales_1.4.0         
      [109] crayon_1.5.3          GetoptLong_1.0.5      rlang_1.1.6           TSP_1.2-5     

        
  **Baseline Tables**
  
     R version 4.5.0 (2025-04-11 ucrt)
     Platform: x86_64-w64-mingw32/x64
     Running under: Windows 11 x64 (build 22631)

    Matrix products: default
      LAPACK version 3.12.1

    locale:
    [1] LC_COLLATE=Danish_Denmark.utf8  LC_CTYPE=Danish_Denmark.utf8    LC_MONETARY=Danish_Denmark.utf8 LC_NUMERIC=C                   
    [5] LC_TIME=Danish_Denmark.utf8    

    time zone: Europe/Copenhagen
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] ggalluvial_0.12.5 table1_1.5.0      kableExtra_1.4.0  lubridate_1.9.4   forcats_1.0.0     stringr_1.5.1     dplyr_1.1.4       purrr_1.1.0      
     [9] readr_2.1.5       tidyr_1.3.1       tibble_3.3.0      ggplot2_3.5.2     tidyverse_2.0.0   here_1.0.1       

    loaded via a namespace (and not attached):
     [1] gtable_0.3.6       compiler_4.5.0     tidyselect_1.2.1   xml2_1.4.0         textshaping_1.0.2  systemfonts_1.2.3  scales_1.4.0      
     [8] fastmap_1.2.0      R6_2.6.1           labeling_0.4.3     generics_0.1.4     Formula_1.2-5      knitr_1.50         rprojroot_2.1.1   
     [15] svglite_2.2.1      pillar_1.11.0      RColorBrewer_1.1-3 tzdb_0.5.0         rlang_1.1.6        stringi_1.8.7      xfun_0.52         
     [22] viridisLite_0.4.2  timechange_0.3.0   cli_3.6.5          withr_3.0.2        magrittr_2.0.3     digest_0.6.37      grid_4.5.0        
     [29] rstudioapi_0.17.1  hms_1.1.3          lifecycle_1.0.4    vctrs_0.6.5        evaluate_1.0.5     glue_1.8.0         farver_2.1.2      
     [36] rmarkdown_2.29     tools_4.5.0        pkgconfig_2.0.3    htmltools_0.5.8.1  
