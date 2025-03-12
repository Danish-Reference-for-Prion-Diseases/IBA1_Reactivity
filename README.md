# IBA1_Protein_Morph_Reactivity

Code and Data for article: 


The order of executed R-scripts to get the complete data and results:

0. QC_Segments_Targets (Nanostring's own Geoscript for Segment and Target QC Evaluation)
1. GeoMx Evaluation and Clean-up
2. Figure_1
3. Cellprofiler Clean-up, normalization and feature selection
4. Figure_2
5. Add_ID for MoBIE
6. Figure_4

(Baseline Tables is standalone and can be executed whenever)



 



**Session infos:**

  **QC_Segments_Targets**

    R version 4.3.0 (2023-04-21 ucrt)
    Platform: x86_64-w64-mingw32/x64 (64-bit)
    Running under: Windows 10 x64 (build 19045)
    
    Matrix products: default
    
    
    locale:
    [1] LC_COLLATE=Danish_Denmark.utf8  LC_CTYPE=Danish_Denmark.utf8    LC_MONETARY=Danish_Denmark.utf8 LC_NUMERIC=C                   
    [5] LC_TIME=Danish_Denmark.utf8    
    
    time zone: Europe/Copenhagen
    tzcode source: internal
    
    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     
    
    other attached packages:
    [1] here_1.0.1
    
    loaded via a namespace (and not attached):
    [1] compiler_4.3.0    rprojroot_2.0.4   tools_4.3.0       rstudioapi_0.17.1
  

  **GeoMx evaluation and Clean-up**
  
    R version 4.3.0 (2023-04-21 ucrt)
    Platform: x86_64-w64-mingw32/x64 (64-bit)
    Running under: Windows 10 x64 (build 19045)
    
    Matrix products: default
    
    
    locale:
    [1] LC_COLLATE=Danish_Denmark.utf8  LC_CTYPE=Danish_Denmark.utf8    LC_MONETARY=Danish_Denmark.utf8 LC_NUMERIC=C                   
    [5] LC_TIME=Danish_Denmark.utf8    
    
    time zone: Europe/Copenhagen
    tzcode source: internal
    
    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     
    
    other attached packages:
     [1] TidyDensity_1.5.0 gridExtra_2.3     readxl_1.4.3      lubridate_1.9.3   forcats_1.0.0     stringr_1.5.1     dplyr_1.1.4       purrr_1.0.2      
     [9] readr_2.1.4       tidyr_1.3.0       tibble_3.2.1      ggplot2_3.5.0     tidyverse_2.0.0   here_1.0.1       
    
    loaded via a namespace (and not attached):
     [1] gtable_0.3.6       compiler_4.3.0     tidyselect_1.2.1   Rcpp_1.0.11        xml2_1.3.6         scales_1.3.0       fastmap_1.2.0      R6_2.5.1          
     [9] plyr_1.8.9         generics_0.1.3     insight_0.20.5     munsell_0.5.1      rprojroot_2.0.4    pillar_1.9.0       tzdb_0.4.0         rlang_1.1.4       
    [17] utf8_1.2.4         stringi_1.8.3      performance_0.12.3 timechange_0.2.0   cli_3.6.2          withr_3.0.2        magrittr_2.0.3     digest_0.6.33     
    [25] grid_4.3.0         rstudioapi_0.17.1  hms_1.1.3          lifecycle_1.0.4    vctrs_0.6.5        data.table_1.14.10 glue_1.8.0         cellranger_1.1.0  
    [33] gt_0.11.1          fansi_1.0.6        colorspace_2.1-0   tools_4.3.0        pkgconfig_2.0.3    htmltools_0.5.8.1 


  **Figure_1:**

    R version 4.3.0 (2023-04-21 ucrt)
    Platform: x86_64-w64-mingw32/x64 (64-bit)
    Running under: Windows 11 x64 (build 22631)

    Matrix products: default


    locale:
    [1] LC_COLLATE=Danish_Denmark.utf8  LC_CTYPE=Danish_Denmark.utf8    LC_MONETARY=Danish_Denmark.utf8 LC_NUMERIC=C                   
    [5] LC_TIME=Danish_Denmark.utf8    

    time zone: Europe/Copenhagen
    tzcode source: internal

    attached base packages:
     [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] ggh4x_0.2.8              ggrepel_0.9.6            performance_0.12.3       emmeans_1.10.0           lmerTest_3.1-3          
     [6] lme4_1.1-35.1            Matrix_1.6-5             variancePartition_1.32.5 BiocParallel_1.36.0      limma_3.58.1            
     [11] gridExtra_2.3            GGally_2.2.1             circlize_0.4.16          RColorBrewer_1.1-3       ComplexHeatmap_2.18.0   
     [16] plotly_4.10.4            factoextra_1.0.7         lubridate_1.9.3          forcats_1.0.0            stringr_1.5.1           
     [21] dplyr_1.1.4              purrr_1.0.2              readr_2.1.4              tidyr_1.3.0              tibble_3.2.1            
     [26] ggplot2_3.5.0            tidyverse_2.0.0          here_1.0.1               MOFA2_1.12.1            

    loaded via a namespace (and not attached):
     [1] rstudioapi_0.17.1     jsonlite_1.8.9        shape_1.4.6.1         magrittr_2.0.3        TH.data_1.1-2        
     [6] estimability_1.5.1    corrplot_0.95         nloptr_2.0.3          GlobalOptions_0.1.2   zlibbioc_1.48.2      
     [11] vctrs_0.6.5           minqa_1.2.6           htmltools_0.5.8.1     S4Arrays_1.2.1        broom_1.0.7          
     [16] Rhdf5lib_1.24.2       SparseArray_1.2.4     rhdf5_2.46.1          KernSmooth_2.23-20    htmlwidgets_1.6.4    
     [21] basilisk_1.14.3       sandwich_3.1-1        pbkrtest_0.5.3        plyr_1.8.9            zoo_1.8-12           
     [26] lifecycle_1.0.4       iterators_1.0.14      pkgconfig_2.0.3       R6_2.5.1              fastmap_1.2.0        
     [31] rbibutils_2.2.16      MatrixGenerics_1.14.0 clue_0.3-65           digest_0.6.33         numDeriv_2016.8-1.1  
     [36] colorspace_2.1-0      S4Vectors_0.40.2      rprojroot_2.0.4       filelock_1.0.3        fansi_1.0.6          
     [41] timechange_0.2.0      httr_1.4.7            abind_1.4-8           compiler_4.3.0        aod_1.3.3            
     [46] withr_3.0.2           doParallel_1.0.17     backports_1.4.1       ggstats_0.7.0         HDF5Array_1.30.1     
     [51] gplots_3.2.0          MASS_7.3-60.0.1       DelayedArray_0.28.0   rjson_0.2.21          corpcor_1.6.10       
     [56] gtools_3.9.5          caTools_1.18.3        tools_4.3.0           remaCor_0.0.18        glue_1.8.0           
     [61] nlme_3.1-164          rhdf5filters_1.14.1   Rtsne_0.17            cluster_2.1.6         reshape2_1.4.4       
     [66] generics_0.1.3        gtable_0.3.6          tzdb_0.4.0            data.table_1.14.10    hms_1.1.3            
     [71] utf8_1.2.4            XVector_0.42.0        BiocGenerics_0.48.1   foreach_1.5.2         pillar_1.9.0         
     [76] splines_4.3.0         lattice_0.21-8        survival_3.5-5        tidyselect_1.2.1      IRanges_2.36.0       
     [81] RhpcBLASctl_0.23-42   stats4_4.3.0          Biobase_2.62.0        statmod_1.5.0         matrixStats_1.2.0    
     [86] pheatmap_1.0.12       stringi_1.8.3         lazyeval_0.2.2        boot_1.3-31           codetools_0.2-19     
     [91] cli_3.6.2             uwot_0.2.2            xtable_1.8-4          reticulate_1.39.0     Rdpack_2.6.1         
     [96] munsell_0.5.1         Rcpp_1.0.11           EnvStats_3.0.0        dir.expiry_1.10.0     coda_0.19-4.1        
     [101] png_0.1-8             parallel_4.3.0        basilisk.utils_1.14.1 bitops_1.0-9          mvtnorm_1.2-4        
     [106] viridisLite_0.4.2     scales_1.3.0          insight_0.20.5        crayon_1.5.3          fANCOVA_0.6-1        
     [111] GetoptLong_1.0.5      rlang_1.1.4           multcomp_1.4-26       cowplot_1.1.3


  **Cellprofiler Clean-up, normalization and feature selection**
  
    R version 4.3.0 (2023-04-21 ucrt)
    Platform: x86_64-w64-mingw32/x64 (64-bit)
    Running under: Windows 10 x64 (build 19045)
    
    Matrix products: default
    
    
    locale:
    [1] LC_COLLATE=Danish_Denmark.utf8  LC_CTYPE=Danish_Denmark.utf8    LC_MONETARY=Danish_Denmark.utf8 LC_NUMERIC=C                   
    [5] LC_TIME=Danish_Denmark.utf8    
    
    time zone: Europe/Copenhagen
    tzcode source: internal
    
    attached base packages:
    [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
    
    other attached packages:
     [1] circlize_0.4.16       ComplexHeatmap_2.18.0 collapse_2.0.16       caret_6.0-94          lattice_0.21-8        sjmisc_2.8.10        
     [7] VIM_6.2.2             colorspace_2.1-0      lubridate_1.9.3       forcats_1.0.0         stringr_1.5.1         dplyr_1.1.4          
    [13] purrr_1.0.2           readr_2.1.4           tidyr_1.3.0           tibble_3.2.1          ggplot2_3.5.0         tidyverse_2.0.0      
    [19] here_1.0.1           
    
    loaded via a namespace (and not attached):
     [1] pROC_1.18.5          rlang_1.1.4          magrittr_2.0.3       clue_0.3-65          GetoptLong_1.0.5     matrixStats_1.2.0    e1071_1.7-16        
     [8] compiler_4.3.0       png_0.1-8            vctrs_0.6.5          reshape2_1.4.4       pkgconfig_2.0.3      shape_1.4.6.1        crayon_1.5.3        
    [15] utf8_1.2.4           prodlim_2023.08.28   tzdb_0.4.0           recipes_1.1.0        cluster_2.1.6        parallel_4.3.0       R6_2.5.1            
    [22] stringi_1.8.3        vcd_1.4-13           RColorBrewer_1.1-3   ranger_0.16.0        parallelly_1.38.0    car_3.1-2            boot_1.3-31         
    [29] rpart_4.1.19         lmtest_0.9-40        Rcpp_1.0.11          iterators_1.0.14     future.apply_1.11.3  zoo_1.8-12           IRanges_2.36.0      
    [36] Matrix_1.6-5         splines_4.3.0        nnet_7.3-18          timechange_0.2.0     tidyselect_1.2.1     rstudioapi_0.17.1    abind_1.4-8         
    [43] timeDate_4041.110    doParallel_1.0.17    codetools_0.2-19     sjlabelled_1.2.0     listenv_0.9.1        plyr_1.8.9           withr_3.0.2         
    [50] future_1.34.0        survival_3.5-5       proxy_0.4-27         pillar_1.9.0         carData_3.0-5        foreach_1.5.2        stats4_4.3.0        
    [57] insight_0.20.5       generics_0.1.3       rprojroot_2.0.4      sp_2.1-4             hms_1.1.3            S4Vectors_0.40.2     munsell_0.5.1       
    [64] scales_1.3.0         laeken_0.5.3         globals_0.16.3       class_7.3-21         glue_1.8.0           tools_4.3.0          robustbase_0.99-4-1 
    [71] data.table_1.14.10   ModelMetrics_1.2.2.2 gower_1.0.1          ipred_0.9-14         nlme_3.1-164         cli_3.6.2            fansi_1.0.6         
    [78] lava_1.8.0           gtable_0.3.6         DEoptimR_1.1-3       digest_0.6.33        BiocGenerics_0.48.1  rjson_0.2.21         lifecycle_1.0.4     
    [85] hardhat_1.4.0        GlobalOptions_0.1.2  MASS_7.3-60.0.1


  **Figure_2**
  
    R version 4.3.0 (2023-04-21 ucrt)
    Platform: x86_64-w64-mingw32/x64 (64-bit)
    Running under: Windows 10 x64 (build 19045)
    
    Matrix products: default
    
    
    locale:
    [1] LC_COLLATE=Danish_Denmark.utf8  LC_CTYPE=Danish_Denmark.utf8    LC_MONETARY=Danish_Denmark.utf8 LC_NUMERIC=C                   
    [5] LC_TIME=Danish_Denmark.utf8    
    
    time zone: Europe/Copenhagen
    tzcode source: internal
    
    attached base packages:
    [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
    
    other attached packages:
     [1] caroline_0.9.9        ggcorrplot_0.1.4.1    psych_2.4.6.26        emmeans_1.10.0        performance_0.12.3    glmmTMB_1.1.10       
     [7] lme4_1.1-35.1         Matrix_1.6-5          fmsb_0.7.6            rstatix_0.7.2         circlize_0.4.16       ComplexHeatmap_2.18.0
    [13] Boruta_8.0.0          caret_6.0-94          lattice_0.21-8        MASS_7.3-60.0.1       JLutils_1.24.0        factoextra_1.0.7     
    [19] lubridate_1.9.3       forcats_1.0.0         stringr_1.5.1         dplyr_1.1.4           purrr_1.0.2           readr_2.1.4          
    [25] tidyr_1.3.0           tibble_3.2.1          ggplot2_3.5.0         tidyverse_2.0.0       here_1.0.1           
    
    loaded via a namespace (and not attached):
      [1] RColorBrewer_1.1-3   rstudioapi_0.17.1    shape_1.4.6.1        magrittr_2.0.3       TH.data_1.1-2        estimability_1.5.1   nloptr_2.0.3        
      [8] GlobalOptions_0.1.2  vctrs_0.6.5          minqa_1.2.6          htmltools_0.5.8.1    haven_2.5.4          broom_1.0.7          pROC_1.18.5         
     [15] parallelly_1.38.0    sandwich_3.1-1       plyr_1.8.9           zoo_1.8-12           TMB_1.9.15           mime_0.12            lifecycle_1.0.4     
     [22] iterators_1.0.14     pkgconfig_2.0.3      R6_2.5.1             fastmap_1.2.0        rbibutils_2.2.16     future_1.34.0        shiny_1.9.1         
     [29] clue_0.3-65          digest_0.6.33        numDeriv_2016.8-1.1  colorspace_2.1-0     GGally_2.2.1         S4Vectors_0.40.2     rprojroot_2.0.4     
     [36] fansi_1.0.6          timechange_0.2.0     abind_1.4-8          mgcv_1.8-42          compiler_4.3.0       withr_3.0.2          doParallel_1.0.17   
     [43] backports_1.4.1      carData_3.0-5        ggstats_0.7.0        highr_0.11           lava_1.8.0           rjson_0.2.21         ModelMetrics_1.2.2.2
     [50] tools_4.3.0          httpuv_1.6.15        future.apply_1.11.3  nnet_7.3-18          glue_1.8.0           questionr_0.7.8      nlme_3.1-164        
     [57] promises_1.3.0       cluster_2.1.6        reshape2_1.4.4       generics_0.1.3       recipes_1.1.0        gtable_0.3.6         labelled_2.13.0     
     [64] tzdb_0.4.0           class_7.3-21         data.table_1.14.10   hms_1.1.3            car_3.1-2            utf8_1.2.4           BiocGenerics_0.48.1 
     [71] ggrepel_0.9.6        foreach_1.5.2        pillar_1.9.0         later_1.3.2          splines_4.3.0        survival_3.5-5       tidyselect_1.2.1    
     [78] miniUI_0.1.1.1       reformulas_0.3.0     IRanges_2.36.0       stats4_4.3.0         hardhat_1.4.0        timeDate_4041.110    matrixStats_1.2.0   
     [85] stringi_1.8.3        boot_1.3-31          codetools_0.2-19     cli_3.6.2            rpart_4.1.19         xtable_1.8-4         Rdpack_2.6.1        
     [92] munsell_0.5.1        Rcpp_1.0.11          globals_0.16.3       coda_0.19-4.1        png_0.1-8            parallel_4.3.0       gower_1.0.1         
     [99] listenv_0.9.1        mvtnorm_1.2-4        ipred_0.9-14         scales_1.3.0         prodlim_2023.08.28   insight_0.20.5       crayon_1.5.3        
    [106] GetoptLong_1.0.5     rlang_1.1.4          mnormt_2.1.1         multcomp_1.4-26


  **Add_ID for MoBIE**
  
    R version 4.3.0 (2023-04-21 ucrt)
    Platform: x86_64-w64-mingw32/x64 (64-bit)
    Running under: Windows 10 x64 (build 19045)
    
    Matrix products: default
    
    
    locale:
    [1] LC_COLLATE=Danish_Denmark.utf8  LC_CTYPE=Danish_Denmark.utf8    LC_MONETARY=Danish_Denmark.utf8 LC_NUMERIC=C                   
    [5] LC_TIME=Danish_Denmark.utf8    
    
    time zone: Europe/Copenhagen
    tzcode source: internal
    
    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     
    
    other attached packages:
     [1] lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.4     tidyr_1.3.0     tibble_3.2.1    ggplot2_3.5.0  
    [10] tidyverse_2.0.0 caroline_0.9.9  here_1.0.1     
    
    loaded via a namespace (and not attached):
     [1] vctrs_0.6.5       cli_3.6.2         rlang_1.1.4       stringi_1.8.3     generics_0.1.3    glue_1.8.0        colorspace_2.1-0  rprojroot_2.0.4  
     [9] hms_1.1.3         scales_1.3.0      fansi_1.0.6       grid_4.3.0        munsell_0.5.1     tzdb_0.4.0        lifecycle_1.0.4   compiler_4.3.0   
    [17] timechange_0.2.0  pkgconfig_2.0.3   rstudioapi_0.17.1 R6_2.5.1          tidyselect_1.2.1  utf8_1.2.4        pillar_1.9.0      magrittr_2.0.3   
    [25] tools_4.3.0       withr_3.0.2       gtable_0.3.6     

  **Figure_4**
  
    R version 4.3.1 (2023-06-16 ucrt)
    Platform: x86_64-w64-mingw32/x64 (64-bit)
    Running under: Windows 11 x64 (build 22631)
    
    Matrix products: default
    
    
    locale:
    [1] LC_COLLATE=English_Denmark.utf8  LC_CTYPE=English_Denmark.utf8   
    [3] LC_MONETARY=English_Denmark.utf8 LC_NUMERIC=C                    
    [5] LC_TIME=English_Denmark.utf8    
    
    time zone: Europe/Berlin
    tzcode source: internal
    
    attached base packages:
    [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
    
    other attached packages:
     [1] MOFA2_1.12.1             viridis_0.6.5            viridisLite_0.4.2       
     [4] ComplexHeatmap_2.18.0    SCORPIUS_1.0.9           cowplot_1.1.3           
     [7] ggh4x_0.2.8              ggrepel_0.9.5            performance_0.12.4      
    [10] emmeans_1.10.3           lme4_1.1-35.3            Matrix_1.6-5            
    [13] variancePartition_1.32.5 BiocParallel_1.36.0      limma_3.58.1            
    [16] gridExtra_2.3            GGally_2.2.1             circlize_0.4.16         
    [19] RColorBrewer_1.1-3       plotly_4.10.4            factoextra_1.0.7        
    [22] lubridate_1.9.3          forcats_1.0.0            stringr_1.5.1           
    [25] dplyr_1.1.4              purrr_1.0.2              readr_2.1.5             
    [28] tidyr_1.3.0              tibble_3.2.1             ggplot2_3.5.1           
    [31] tidyverse_2.0.0          here_1.0.1              
    
    loaded via a namespace (and not attached):
      [1] splines_4.3.1         bitops_1.0-8          filelock_1.0.3        basilisk.utils_1.14.1
      [5] lifecycle_1.0.4       Rdpack_2.6.1          doParallel_1.0.17     rprojroot_2.0.4      
      [9] processx_3.8.4        lattice_0.22-5        MASS_7.3-60.0.1       insight_0.20.5       
     [13] backports_1.5.0       magrittr_2.0.3        yaml_2.3.10           remotes_2.5.0        
     [17] reticulate_1.35.0     pbapply_1.7-2         minqa_1.2.8           zlibbioc_1.48.2      
     [21] multcomp_1.4-26       abind_1.4-5           Rtsne_0.17            EnvStats_2.8.1       
     [25] BiocGenerics_0.48.1   TH.data_1.1-2         sandwich_3.1-0        IRanges_2.36.0       
     [29] S4Vectors_0.40.2      pbkrtest_0.5.3        irlba_2.3.5.1         pheatmap_1.0.12      
     [33] codetools_0.2-19      DelayedArray_0.28.0   tidyselect_1.2.1      babelwhale_1.2.0     
     [37] shape_1.4.6.1         dynwrap_1.2.4         TSP_1.2-4             matrixStats_1.2.0    
     [41] stats4_4.3.1          jsonlite_1.8.8        GetoptLong_1.0.5      survival_3.5-7       
     [45] iterators_1.0.14      foreach_1.5.2         tools_4.3.1           Rcpp_1.0.12          
     [49] glue_1.6.2            SparseArray_1.2.4     ranger_0.16.0         MatrixGenerics_1.14.0
     [53] HDF5Array_1.30.1      withr_3.0.1           numDeriv_2016.8-1.1   fastmap_1.2.0        
     [57] basilisk_1.14.3       rhdf5filters_1.14.1   boot_1.3-28.1         fansi_1.0.6          
     [61] caTools_1.18.2        digest_0.6.37         timechange_0.3.0      R6_2.5.1             
     [65] estimability_1.5.1    colorspace_2.1-0      gtools_3.9.5          RhpcBLASctl_0.23-42  
     [69] utf8_1.2.4            generics_0.1.3        data.table_1.15.4     corpcor_1.6.10       
     [73] httr_1.4.7            htmlwidgets_1.6.4     S4Arrays_1.2.1        ggstats_0.6.0        
     [77] uwot_0.1.16           pkgconfig_2.0.3       gtable_0.3.5          XVector_0.42.0       
     [81] lmds_0.1.0            remaCor_0.0.18        htmltools_0.5.8.1     clue_0.3-65          
     [85] scales_1.3.0          Biobase_2.62.0        png_0.1-8             corrplot_0.94        
     [89] fANCOVA_0.6-1         rstudioapi_0.16.0     tzdb_0.4.0            reshape2_1.4.4       
     [93] rjson_0.2.21          coda_0.19-4.1         nlme_3.1-164          nloptr_2.1.1         
     [97] rhdf5_2.46.1          zoo_1.8-12            GlobalOptions_0.1.2   KernSmooth_2.23-22   
    [101] parallel_4.3.1        desc_1.4.3            pillar_1.9.0          proxyC_0.4.1         
    [105] vctrs_0.6.5           RANN_2.6.1            gplots_3.1.3.1        dynparam_1.0.2       
    [109] xtable_1.8-4          cluster_2.1.6         princurve_2.1.6       mvtnorm_1.2-6        
    [113] cli_3.6.1             compiler_4.3.1        rlang_1.1.2           crayon_1.5.3         
    [117] mclust_6.1.1          carrier_0.1.1         ps_1.7.7              plyr_1.8.9           
    [121] stringi_1.8.3         assertthat_0.2.1      lmerTest_3.1-3        munsell_0.5.1        
    [125] lazyeval_0.2.2        aod_1.3.3             dir.expiry_1.10.0     hms_1.1.3            
    [129] Rhdf5lib_1.24.2       statmod_1.5.0         rbibutils_2.2.16      igraph_2.0.3         
    [133] broom_1.0.6           dynutils_1.0.11


**Baseline Tables**

    R version 4.3.0 (2023-04-21 ucrt)
    Platform: x86_64-w64-mingw32/x64 (64-bit)
    Running under: Windows 10 x64 (build 19045)
    
    Matrix products: default
    
    
    locale:
    [1] LC_COLLATE=Danish_Denmark.utf8  LC_CTYPE=Danish_Denmark.utf8    LC_MONETARY=Danish_Denmark.utf8 LC_NUMERIC=C                   
    [5] LC_TIME=Danish_Denmark.utf8    
    
    time zone: Europe/Copenhagen
    tzcode source: internal
    
    attached base packages:
    [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
    
    other attached packages:
     [1] PRISMAstatement_1.1.1 gtsummary_2.0.3       gt_0.11.1             readxl_1.4.3          ggalluvial_0.12.5     table1_1.4.3         
     [7] kableExtra_1.4.0      lubridate_1.9.3       forcats_1.0.0         stringr_1.5.1         dplyr_1.1.4           purrr_1.0.2          
    [13] readr_2.1.4           tidyr_1.3.0           tibble_3.2.1          ggplot2_3.5.0         tidyverse_2.0.0       here_1.0.1           
    
    loaded via a namespace (and not attached):
     [1] utf8_1.2.4        generics_0.1.3    xml2_1.3.6        stringi_1.8.3     hms_1.1.3         digest_0.6.33     magrittr_2.0.3    evaluate_1.0.1   
     [9] timechange_0.2.0  fastmap_1.2.0     cellranger_1.1.0  rprojroot_2.0.4   Formula_1.2-5     fansi_1.0.6       viridisLite_0.4.2 scales_1.3.0     
    [17] cli_3.6.2         rlang_1.1.4       munsell_0.5.1     withr_3.0.2       tools_4.3.0       tzdb_0.4.0        colorspace_2.1-0  vctrs_0.6.5      
    [25] R6_2.5.1          lifecycle_1.0.4   pkgconfig_2.0.3   pillar_1.9.0      gtable_0.3.6      glue_1.8.0        systemfonts_1.1.0 xfun_0.48        
    [33] tidyselect_1.2.1  rstudioapi_0.17.1 knitr_1.48        htmltools_0.5.8.1 rmarkdown_2.28    svglite_2.1.3     compiler_4.3.0 
