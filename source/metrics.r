box::use(
  dplyr[d_filter = filter, d_select = select, ...],
  # for unite
  tidyr[...],
  # for enframe  
  tibble[...],
  stringr[...],
  readr[...],
  mltools[...], 
  epiR[...], 
 )

compute_metrics_for_discrete = function (tool,
                                         label_p, label_b, 
                                         predictions, 
                                         pathogenic,
                                         benign){
   
   TP = pathogenic %>%
      # filter column of tool in current call by the label it uses for pathogenic
      # variants
      d_filter(str_detect(!!as.name(tool), label_p)) %>% 
      summarise(TP = n()) %>%
      mutate(TP = ifelse(TP == 0, 1, TP)) %>%
      .[['TP']]
   
   TN = benign %>%
      # filter column of tool in current call by the label it uses for benign
      # variants
      d_filter(str_detect(!!as.name(tool), label_b)) %>% 
      summarise(TN = n()) %>%
      mutate(TN = ifelse(TN == 0, 1, TN)) %>%
      .[['TN']]
   
   FP = benign %>%
      # filter column of tool in current call by the label it uses for pathogenic
      # variants
      d_filter(str_detect(!!as.name(tool), label_p)) %>% 
      summarise(FP = n()) %>%
      mutate(FP = ifelse(FP == 0, 1, FP)) %>%
      .[['FP']]
   
   FN = pathogenic %>%
      # filter column of tool in current call by the label it uses for benign
      # variants
      d_filter(str_detect(!!as.name(tool), label_b)) %>% 
      summarise(FN = n()) %>%
      mutate(FN = ifelse(FN == 0, 1, FN)) %>%
      .[['FN']]
   
   # TP, FP etc as matrix, was the input of the epi.test below
   # df = as.table(matrix(c(TP,FN,FP,TN)), nrow = 2, byrow = TRUE)
   # and it generated a one column matrix. The order of the stats
   # is wrong according to the documentation but somehow the function freturned
   # correct results. For my sanity using a vector of the measures in the 
   # correct order according to documentation which is: TP, FP, FN, TN
   stats = c(TP, FP, FN, TN)
   rval = epi.tests(stats, conf.level = 0.95)

   YR = predictions %>%
      d_filter(!is.na(!!as.name(tool))) %>%
      summarise("YR" = n() / (predictions %>% nrow())) %>%
      .[['YR']]

   #OPP = sqrt(((rval$detail$pv.pos$est^2) + (rval$detail$pv.neg$est^2) + (YR^2))/3)
   OPP = sqrt((((rval$detail %>% d_filter(statistic == "pv.pos") %>% .[['est']]) ^ 2) + ((rval$detail %>% d_filter(statistic == "pv.neg") %>% .[['est']]) ^ 2) + (YR^2))/3)

   MCC = mcc(TP = TP, FP = FP, TN = TN, FN = FN)
   ACC = (TP + TN)/(TP + FN + FP + TN)
   
   results = tibble(
      Tool = tool, 
      All_SNVs = predictions %>% nrow(), 
      # number of variants for which annotated pathogenicity is not VUS
      # annotated pathogenic variants
      APV = predictions %>% 
         d_filter(str_detect(Observed_pathogenicity, "VUS", negate = TRUE)) %>% 
         nrow(), 
      # number of non VUS variants for which tool gives predictions
      # predicted variants, non VUS annotation
      PV_NV = predictions %>% 
         d_filter(str_detect(Observed_pathogenicity, "VUS", negate = TRUE), !is.na(!!as.name(tool))) %>% 
         summarise(tmp = n()) %>% .[['tmp']], 
      # number of variants (including VUS) for which tool gives predictions
      # predicted variants
      PV = predictions %>% d_filter(!is.na(!!as.name(tool))) %>%
         summarise(tmp = n()) %>% .[['tmp']],
      TP = TP,
      FN = FN,
      FP = FP,
      TN = TN,
      YR = YR,
      Accuracy = ACC,
      Sensitivity = rval$detail %>% d_filter(statistic == "se") %>% .[['est']],
      Specificity = rval$detail %>% d_filter(statistic == "sp") %>% .[['est']],
      PPV = rval$detail %>% d_filter(statistic == "pv.pos") %>% .[['est']],
      NPV = rval$detail %>% d_filter(statistic == "pv.neg") %>% .[['est']],
      OPP = OPP,
      MCC = MCC,
      `LR+` = rval$detail %>% d_filter(statistic == "lr.pos") %>% .[['est']],
      `LR+_CI_lower` = rval$detail %>% d_filter(statistic == "lr.pos") %>% .[['lower']],
      `LR+_CI_upper` = rval$detail %>% d_filter(statistic == "lr.pos") %>% .[['upper']],
      `LR-` = rval$detail %>% d_filter(statistic == "lr.neg") %>% .[['est']],
      `LR-_CI_lower` = rval$detail %>% d_filter(statistic == "lr.neg") %>% .[['lower']],
      `LR-_CI_upper` = rval$detail %>% d_filter(statistic == "lr.neg") %>% .[['upper']],
   )
   
   return(results)
}


compute_metrics_for_continuous = function(tool,
                                          threshold, th_class, 
                                          predictions, 
                                          pathogenic,
                                          benign){
   
   if(th_class == "below"){
      TP = pathogenic %>%
         # filter by threshold
         d_filter(!is.na(!!as.name(tool)) & !!as.name(tool) < threshold) %>% 
         summarise(TP = n()) %>%
         mutate(TP = ifelse(TP == 0, 1, TP)) %>%
         .[['TP']]
      
      TN = benign %>%
         d_filter(!is.na(!!as.name(tool)) & !!as.name(tool) >= threshold) %>% 
         summarise(TN = n()) %>%
         mutate(TN = ifelse(TN == 0, 1, TN)) %>%
         .[['TN']]
      
      FP = benign %>%
         d_filter(!is.na(!!as.name(tool)) & !!as.name(tool) < threshold) %>% 
         summarise(FP = n()) %>%
         mutate(FP = ifelse(FP == 0, 1, FP)) %>%
         .[['FP']]
      
      FN = pathogenic %>%
         d_filter(!is.na(!!as.name(tool)) & !!as.name(tool) >= threshold) %>% 
         summarise(FN = n()) %>%
         mutate(FN = ifelse(FN == 0, 1, FN)) %>%
         .[['FN']]
   } else if(th_class == "below_closed") {
      TP = pathogenic %>%
         # filter by threshold
         d_filter(!is.na(!!as.name(tool)) & !!as.name(tool) <= threshold) %>% 
         summarise(TP = n()) %>%
         mutate(TP = ifelse(TP == 0, 1, TP)) %>%
         .[['TP']]
      
      TN = benign %>%
         d_filter(!is.na(!!as.name(tool)) & !!as.name(tool) > threshold) %>% 
         summarise(TN = n()) %>%
         mutate(TN = ifelse(TN == 0, 1, TN)) %>%
         .[['TN']]
      
      FP = benign %>%
         d_filter(!is.na(!!as.name(tool)) & !!as.name(tool) <= threshold) %>% 
         summarise(FP = n()) %>%
         mutate(FP = ifelse(FP == 0, 1, FP)) %>%
         .[['FP']]
      
      FN = pathogenic %>%
         d_filter(!is.na(!!as.name(tool)) & !!as.name(tool) > threshold) %>% 
         summarise(FN = n()) %>%
         mutate(FN = ifelse(FN == 0, 1, FN)) %>%
         .[['FN']]
   } else if(th_class == "above"){
      TP = pathogenic %>%
         # filter by threshold
         d_filter(!is.na(!!as.name(tool)) & !!as.name(tool) > threshold) %>% 
         summarise(TP = n()) %>%
         mutate(TP = ifelse(TP == 0, 1, TP)) %>%
         .[['TP']]
      
      TN = benign %>%
         d_filter(!is.na(!!as.name(tool)) & !!as.name(tool) <= threshold) %>% 
         summarise(TN = n()) %>%
         mutate(TN = ifelse(TN == 0, 1, TN)) %>%
         .[['TN']]
      
      FP = benign %>%
         d_filter(!is.na(!!as.name(tool)) & !!as.name(tool) > threshold) %>% 
         summarise(FP = n()) %>%
         mutate(FP = ifelse(FP == 0, 1, FP)) %>%
         .[['FP']]
      
      FN = pathogenic %>%
         d_filter(!is.na(!!as.name(tool)) & !!as.name(tool) <= threshold) %>% 
         summarise(FN = n()) %>%
         mutate(FN = ifelse(FN == 0, 1, FN)) %>%
         .[['FN']]
   } else if(th_class == "above_closed"){
      TP = pathogenic %>%
         # filter by threshold
         d_filter(!is.na(!!as.name(tool)) & !!as.name(tool) >= threshold) %>% 
         summarise(TP = n()) %>%
         mutate(TP = ifelse(TP == 0, 1, TP)) %>%
         .[['TP']]
      
      TN = benign %>%
         d_filter(!is.na(!!as.name(tool)) & !!as.name(tool) < threshold) %>% 
         summarise(TN = n()) %>%
         mutate(TN = ifelse(TN == 0, 1, TN)) %>%
         .[['TN']]
      
      FP = benign %>%
         d_filter(!is.na(!!as.name(tool)) & !!as.name(tool) >= threshold) %>% 
         summarise(FP = n()) %>%
         mutate(FP = ifelse(FP == 0, 1, FP)) %>%
         .[['FP']]
      
      FN = pathogenic %>%
         d_filter(!is.na(!!as.name(tool)) & !!as.name(tool) < threshold) %>% 
         summarise(FN = n()) %>%
         mutate(FN = ifelse(FN == 0, 1, FN)) %>%
         .[['FN']]
   }
      
   
   # TP, FP etc as matrix, was the input of the epi.test below
   # df = as.table(matrix(c(TP,FN,FP,TN)), nrow = 2, byrow = TRUE)
   # and it generated a one column matrix. The order of the stats
   # is wrong according to the documentation but somehow the function freturned
   # correct results. For my sanity using a vector of the measures in the 
   # correct order according to documentation which is: TP, FP, FN, TN
   stats = c(TP, FP, FN, TN)
   rval = epi.tests(stats, conf.level = 0.95)

   YR = predictions %>%
      d_filter(!is.na(!!as.name(tool))) %>%
      summarise("YR" = n() / (predictions %>% nrow())) %>%
      .[['YR']]

   OPP = sqrt((((rval$detail %>% d_filter(statistic == "pv.pos") %>% .[['est']]) ^2) + ((rval$detail %>% d_filter(statistic == "pv.neg") %>% .[['est']]) ^2) + (YR^2))/3)

   MCC = mcc(TP = TP, FP = FP, TN = TN, FN = FN)
   ACC = (TP + TN)/(TP + FN + FP + TN)
   
   results = tibble(
      Tool = tool, 
      Threshold = threshold,
      All_SNVs = predictions %>% nrow(), 
      # number of variants for which annotated pathogenicity is not VUS
      # annotated pathogenic variants
      APV = predictions %>% 
         d_filter(str_detect(Observed_pathogenicity, "VUS", negate = TRUE)) %>% 
         nrow(), 
      # number of non VUS variants for which tool gives predictions
      # predicted variants, non VUS annotation
      PV_NV = predictions %>% 
         d_filter(str_detect(Observed_pathogenicity, "VUS", negate = TRUE) & !is.na(!!as.name(tool))) %>% 
         nrow(), 
      # number of variants (including VUS) for which tool gives predictions
      # predicted variants
      PV = predictions %>% d_filter(!is.na(!!as.name(tool))) %>%
         nrow(),
      TP = TP,
      FN = FN,
      FP = FP,
      TN = TN,
      YR = YR,
      Accuracy = ACC,
      Sensitivity = rval$detail %>% d_filter(statistic == "se") %>% .[['est']],
      Specificity = rval$detail %>% d_filter(statistic == "sp") %>% .[['est']],
      PPV = rval$detail %>% d_filter(statistic == "pv.pos") %>% .[['est']],
      NPV = rval$detail %>% d_filter(statistic == "pv.neg") %>% .[['est']],
      OPP = OPP,
      MCC = MCC,
      `LR+` = rval$detail %>% d_filter(statistic == "lr.pos") %>% .[['est']],
      `LR+_CI_lower` = rval$detail %>% d_filter(statistic == "lr.pos") %>% .[['lower']],
      `LR+_CI_upper` = rval$detail %>% d_filter(statistic == "lr.pos") %>% .[['upper']],
      `LR-` = rval$detail %>% d_filter(statistic == "lr.neg") %>% .[['est']],
      `LR-_CI_lower` = rval$detail %>% d_filter(statistic == "lr.neg") %>% .[['lower']],
      `LR-_CI_upper` = rval$detail %>% d_filter(statistic == "lr.neg") %>% .[['upper']],
   )
   
   return(results)
}

# ST: MaxentScan score = diff*100/ref && if diff > 3 && score > 30% => set as pathogenic
compute_metrics_for_maxentscan = function(tool_1, tool_2,
                                th1, th2, 
                                predictions, 
                                pathogenic,
                                benign){
   
   TP = pathogenic %>%
      # filter by threshold
      d_filter(!is.na(!!as.name(tool_1)) & abs(!!as.name(tool_1)) > th1 &
                  !is.na(!!as.name(tool_2)) & abs(!!as.name(tool_2)) > th2 ) %>% 
      summarise(TP = n()) %>%
      mutate(TP = ifelse(TP == 0, 1, TP)) %>%
      .[['TP']]
   
   TN = benign %>%
      d_filter( ! (!is.na(!!as.name(tool_1)) & abs(!!as.name(tool_1)) > th1 &
                  !is.na(!!as.name(tool_2)) & abs(!!as.name(tool_2)) > th2 ) &
                   (!is.na(!!as.name(tool_1)) & !is.na(!!as.name(tool_2)))
                ) %>% 
      summarise(TN = n()) %>%
      mutate(TN = ifelse(TN == 0, 1, TN)) %>%
      .[['TN']]
   
   # if it passes the filter (> both thresholds) and is predicted pathogenic
   # but the observed phenotype is benign
   # then it is a FP
   FP = benign %>%
      d_filter(!is.na(!!as.name(tool_1)) & abs(!!as.name(tool_1)) > th1 &
                  !is.na(!!as.name(tool_2)) & abs(!!as.name(tool_2)) > th2 ) %>% 
      summarise(FP = n()) %>%
      mutate(FP = ifelse(FP == 0, 1, FP)) %>%
      .[['FP']]
   
   # if it passes the filter (< both filters) the filter and is predicted 
   # pathogenic but the observed phenotype is benign
   # then it is a FP
   FN = pathogenic %>%
      d_filter( ! (!is.na(!!as.name(tool_1)) & abs(!!as.name(tool_1)) > th1 &
                  !is.na(!!as.name(tool_2)) & abs(!!as.name(tool_2)) > th2 ) &
                   (!is.na(!!as.name(tool_1)) & !is.na(!!as.name(tool_2)))
                ) %>% 
      summarise(FN = n()) %>%
      mutate(FN = ifelse(FN == 0, 1, FN)) %>%
      .[['FN']]

   # TP, FP etc as matrix, was the input of the epi.test below
   # df = as.table(matrix(c(TP,FN,FP,TN)), nrow = 2, byrow = TRUE)
   # and it generated a one column matrix. The order of the stats
   # is wrong according to the documentation but somehow the function freturned
   # correct results. For my sanity using a vector of the measures in the 
   # correct order according to documentation which is: TP, FP, FN, TN
   stats = c(TP, FP, FN, TN)
   rval = epi.tests(stats, conf.level = 0.95)

   YR = predictions %>%
      d_filter(!is.na(!!as.name(tool_1)) & !is.na(!!as.name(tool_2))) %>%
      summarise("YR" = n() / (predictions %>% nrow())) %>%
      .[['YR']]

   OPP = sqrt((((rval$detail %>% d_filter(statistic == "pv.pos") %>% .[['est']]) ^2) + ((rval$detail %>% d_filter(statistic == "pv.neg") %>% .[['est']]) ^2) + (YR^2))/3)

   MCC = mcc(TP = TP, FP = FP, TN = TN, FN = FN)
   ACC = (TP + TN)/(TP + FN + FP + TN)
   
   results = tibble(
      Tool_1 = tool_1, 
      Tool_2 = tool_2, 
      Threshold_1 = th1,
      Threshold_2 = th2,
      All_SNVs = predictions %>% nrow(), 
      # number of variants for which annotated pathogenicity is not VUS
      # annotated pathogenic variants
      APV = predictions %>%
         d_filter(str_detect(Observed_pathogenicity, "VUS", negate = TRUE)) %>% 
         nrow(), 
      # number of non VUS variants for which tool gives predictions
      # predicted variants, non VUS annotation
      PV_NV = predictions %>% 
         d_filter(str_detect(Observed_pathogenicity, "VUS", negate = TRUE) &
                  !is.na(!!as.name(tool_1)) & !is.na(!!as.name(tool_2))) %>% 
         nrow(), 
      # number of variants (including VUS) for which tool gives predictions
      # predicted variants
      PV = predictions %>% 
         d_filter(!is.na(!!as.name(tool_1)) & !is.na(!!as.name(tool_2)))%>%
         nrow(),
      TP = TP,
      FN = FN,
      FP = FP,
      TN = TN,
      YR = YR,
      Accuracy = ACC,
      Sensitivity = rval$detail %>% d_filter(statistic == "se") %>% .[['est']],
      Specificity = rval$detail %>% d_filter(statistic == "sp") %>% .[['est']],
      PPV = rval$detail %>% d_filter(statistic == "pv.pos") %>% .[['est']],
      NPV = rval$detail %>% d_filter(statistic == "pv.neg") %>% .[['est']],
      OPP = OPP,
      MCC = MCC,
      `LR+` = rval$detail %>% d_filter(statistic == "lr.pos") %>% .[['est']],
      `LR+_CI_lower` = rval$detail %>% d_filter(statistic == "lr.pos") %>% .[['lower']],
      `LR+_CI_upper` = rval$detail %>% d_filter(statistic == "lr.pos") %>% .[['upper']],
      `LR-` =  rval$detail %>% d_filter(statistic == "lr.neg") %>% .[['est']],
      `LR-_CI_lower` = rval$detail %>% d_filter(statistic == "lr.neg") %>% .[['lower']],
      `LR-_CI_upper` = rval$detail %>% d_filter(statistic == "lr.neg") %>% .[['upper']],
   )
}

