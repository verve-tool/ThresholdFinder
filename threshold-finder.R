
directories <- c("fullpage", "horizontal", "horizontal-vertical", "limited-horizontal-vertical", "maxAOC");
mdirectories <- c("fullpage-mutated", "horizontal-mutated", "horizontal-vertical-mutated", "limited-horizontal-vertical-mutated", "maxAOC-mutated");
edirectories <- c("fullpage-ext", "horizontal-ext", "horizontal-vertical-ext", "limited-horizontal-vertical-ext", "maxAOC-ext");
ext_site_names <- c("EatThisMuch","Forvo","GMapStreetViewPlayer","HoursOf","RetailMeNot","SimilarSites","Tiiime","startbootstrap-agency","startbootstrap-business-casual","startbootstrap-clean-blog","startbootstrap-coming-soon","startbootstrap-creative","startbootstrap-freelancer","startbootstrap-grayscale","startbootstrap-landing-page","startbootstrap-new-age","startbootstrap-one-page-wonder","startbootstrap-resume","startbootstrap-sb-admin-2","startbootstrap-stylish-portfolio")

mTypes <-c("alt"="Alternative Chi-Square", "bha"="Bhattacharyya Distance", "chi"="Chi-Square", "cor"="Correlation", "int"="Intersection", "kld"="Kullback-Leibler Divergence")

step <- 0.01;

input_dir <- file.path(getwd(),"input", "second-submission")
output_dir <- file.path(getwd(),"output", "second-submission", "threshold-finder")
if(!file.exists(output_dir)){
  dir.create(output_dir)
}

absolute_columns <- function(df){
  for(c in seq(from=5, to=ncol(df)-1, by=3)){
    df[,c+1] <-  abs(df[,c+1]);
    df[,c+2] <-  abs(df[,c+2]);
  }
  return(df)
}
get_possible_thresholds <- function(df, step, type){
  thresholds <- c();
  for(c in seq(from=5, to=ncol(df)-1, by=3)){
    if(type == "int"){
      base <- df[,c];
      min_oracle <- df[,c+1];
      max_oracle <- df[,c+2];
      int_value <- (base-min_oracle)/base;
      thresholds <- append(thresholds, int_value);
      int_value <- (base-max_oracle)/base;
      thresholds <- append(thresholds, int_value);
    }else{
      thresholds <- append(thresholds, df[,c]);
      thresholds <- append(thresholds, df[,c+1]);
      thresholds <- append(thresholds, df[,c+2]);
    }
  }
  thresholds <- append(thresholds, 0);
  thresholds <- append(thresholds, 1);
  
  uniq_df <- unique(thresholds);
  u_plus_df<- uniq_df+step;
  u_minus_df<- uniq_df-step;
  uniq_df <- append(uniq_df, u_plus_df);
  uniq_df <- append(uniq_df, u_minus_df);
  uniq_df <- unique(round(uniq_df,digits = 2));
  if(type == "cor" || type == "bha" || type == "int"){
    uniq_df <- uniq_df[!uniq_df > 1];
  }
  uniq_df <- uniq_df[!uniq_df < 0];
  
  return(sort(uniq_df));
}
min_oracle <- function(df){
  min <- 9999999999;
  for(c in seq(from=5, to=ncol(df)-1, by=3)){
    min <- min(df[,c+1], df[,c+2]);
  }
  return(min)
}
max_oracle <- function(df){
  max <- -9999999999;
  for(c in seq(from=5, to=ncol(df)-1, by=3)){
    max <- max(df[,c+1], df[,c+2]);
  }
  return(max)
}
alt_bha_chi_kld_classify <- function(df,threshold){
  for(r in 1:nrow(df)){
    df[r,ncol(df)] = "FP";
    for(c in seq(from=5, to=ncol(df)-1, by=3)){
      base <- df[r,c];
      min_oracle <- df[r,c+1];
      max_oracle <- df[r,c+2];
      oracle <-min(min_oracle,max_oracle);
      if(oracle >= threshold){
        df[r,ncol(df)] = "TP";
        break;
      }
    }
  }
  return(df);
}
cor_classify <- function(df,threshold){
  for(r in 1:nrow(df)){
    df[r,ncol(df)] = "FP";
    for(c in seq(from=5, to=ncol(df)-1, by=3)){
      base <- df[r,c];
      min_oracle <- df[r,c+1];
      max_oracle <- df[r,c+2];
      oracle <-max(min_oracle,max_oracle);
      if(oracle < threshold){
        df[r,ncol(df)] = "TP";
        break;
      }
    }
  }
  return(df);
}
int_classify <- function(df,threshold){
  for(r in 1:nrow(df)){
    df[r,ncol(df)] = "FP";
    for(c in seq(from=5, to=ncol(df)-1, by=3)){
      base <- df[r,c];
      min_oracle <- df[r,c+1];
      max_oracle <- df[r,c+2];
      oracle <-max(min_oracle,max_oracle);
      int_value <- (base-oracle)/base;
      if( int_value >= threshold){
        df[r,ncol(df)] = "TP";
        break;
      }
    }
  }
  return(df);
}

get_accuracy <- function(auto, man){
  #Last column is the classification
  # ac <- ncol(auto)
  # mc <- ncol(man)
  # matched <- 0
  # for(r in 1:nrow(auto)){
  #   a <- auto[r, ac]
  #   m <- man[r, mc]
  #   if(a == m){
  #     matched <- matched + 1
  #   }
  # }
  auto[,ncol(auto)] <- factor(auto[,ncol(auto)], levels=c("TP", "FP", "NOI"));
  man[,ncol(man)] <- factor(man[,ncol(man)], levels=c("TP", "FP", "NOI"));
  matched <- nrow(auto[auto[,ncol(auto)] == man[,ncol(man)],]);
  return(round((matched/nrow(auto)) * 100, digits = 2));
}

get_best_thresholds <- function(auto, man, thresholds, step, type){
  minimum <- min(thresholds);
  maximum <- max(thresholds);
  #print(paste("Test-Thresholds:", length(thresholds),"Min:", min(thresholds), "Max:", max(thresholds)));
  results_df <- data.frame(Threshold=numeric(),Accuracy=numeric());
  for(threshold in thresholds){
    auto <- classify(auto, threshold, type);
    results_df[nrow(results_df) + 1,] = c(threshold , get_accuracy(auto,man));
  }
  best_df <- results_df[results_df$Accuracy == max(results_df$Accuracy) ,];
  return(best_df);
}
run_for_file <- function(auto, man, step, type){


  start_time <- Sys.time();
  #print(paste("--------------", type,"--------------"));

  thresholds <- get_possible_thresholds(auto, step, type);
  result <- get_best_thresholds(auto, man, thresholds, step, type);
  threshold <- result[ceiling(nrow(result)/2),1];
  accuracy <- max(result[,2]);
  #accuracy <- paste(accuracy, "%", sep = "");
  end_time <- Sys.time();
  #print(paste("Threshold:", threshold));
  #print(paste("Accuracy:", accuracy));
  #print(paste("Threshold:", threshold, "Accuracy:", accuracy,"%", sep = ""), "Threshold-Score:", paste(round((nrow(result)/length(thresholds))*100, digits = 2),"%", sep = ""),"Count:", nrow(result), "Time:", paste(round(end_time-start_time,digits = 0), "s", sep = "")));
  return(c(threshold,accuracy));
}
run_for_file_using_threshold <- function(auto, man, threshold, type){
  auto <- classify(auto, threshold, type);
  accuracy <- get_accuracy(auto,man);
  #accuracy <- paste(accuracy, "%", sep = "");
  #print(paste("Accuracy:", accuracy));
  #print(paste("Mutation-Threshold:", threshold, " Accuracy:", accuracy,"%", sep = ""));
  return(accuracy);
}
classify <- function(auto, threshold, type){
  if(type == "alt" || type == "bha" || type == "chi" || type == "kld"){
    classify_method <- alt_bha_chi_kld_classify;
  }
  else if(type == "cor"){
    classify_method <- cor_classify;
    
  }else if(type == "int"){
    classify_method <- int_classify;
  }
  auto <- classify_method(auto, threshold);
  return(auto);
}
run_for_directory_by_parsing <- function(directory, mdirectory, approachName){
  alt_full <- read.csv(file.path(input_dir, directory, "Alternative-Chi-Square.csv"), stringsAsFactors = FALSE);
  bha_full <- read.csv(file.path(input_dir, directory, "Bhattacharyya.csv"), stringsAsFactors = FALSE);
  chi_full <- read.csv(file.path(input_dir, directory, "Chi-Square.csv"), stringsAsFactors = FALSE);
  cor_full <- read.csv(file.path(input_dir, directory, "Correlation.csv"), stringsAsFactors = FALSE);
  int_full <- read.csv(file.path(input_dir, directory, "Intersection.csv"), stringsAsFactors = FALSE);
  kld_full <- read.csv(file.path(input_dir, directory, "Kullback-Leibler-Divergence.csv"), stringsAsFactors = FALSE);
  man_full <- read.csv(file.path(input_dir, "manual-report-all-failures.csv"), stringsAsFactors = FALSE);
  man_full <- man_full[man_full$Failure == "Small-Range" | man_full$Failure == "small-range",]
  kld_full <- absolute_columns(kld_full);
  
  alt_m <- read.csv(file.path(input_dir, mdirectory, "Alternative-Chi-Square.csv"), stringsAsFactors = FALSE);
  bha_m <- read.csv(file.path(input_dir, mdirectory, "Bhattacharyya.csv"), stringsAsFactors = FALSE);
  chi_m <- read.csv(file.path(input_dir, mdirectory, "Chi-Square.csv"), stringsAsFactors = FALSE);
  cor_m <- read.csv(file.path(input_dir, mdirectory, "Correlation.csv"), stringsAsFactors = FALSE);
  int_m <- read.csv(file.path(input_dir, mdirectory, "Intersection.csv"), stringsAsFactors = FALSE);
  kld_m <- read.csv(file.path(input_dir, mdirectory, "Kullback-Leibler-Divergence.csv"), stringsAsFactors = FALSE);
  man_m <- read.csv(file.path(input_dir, "manual-report-mutated-pages.csv"), stringsAsFactors = FALSE);
  kld_m <- absolute_columns(kld_m);
  
  alt_m <- alt_m[order(alt_m$UID),]
  bha_m <- bha_m[order(bha_m$UID),]
  chi_m <- chi_m[order(chi_m$UID),]
  cor_m <- cor_m[order(cor_m$UID),]
  int_m <- int_m[order(int_m$UID),]
  kld_m <- kld_m[order(kld_m$UID),]
  man_m <- man_m[order(man_m$UID),]
  
  alt_e <- subset(alt_full, subset = (Webpage %in% ext_site_names))
  bha_e <- subset(bha_full, subset = (Webpage %in% ext_site_names))
  chi_e <- subset(chi_full, subset = (Webpage %in% ext_site_names))
  cor_e <- subset(cor_full, subset = (Webpage %in% ext_site_names))
  int_e <- subset(int_full, subset = (Webpage %in% ext_site_names))
  kld_e <- subset(kld_full, subset = (Webpage %in% ext_site_names))
  man_e <- subset(man_full, subset = (Webpage %in% ext_site_names))
  
  alt_e <- alt_e[order(alt_e$UID),]
  bha_e <- bha_e[order(bha_e$UID),]
  chi_e <- chi_e[order(chi_e$UID),]
  cor_e <- cor_e[order(cor_e$UID),]
  int_e <- int_e[order(int_e$UID),]
  kld_e <- kld_e[order(kld_e$UID),]
  man_e <- man_e[order(man_e$UID),]
  
  
  alt_r <- subset(alt_full, subset = !(Webpage %in% ext_site_names))
  bha_r <- subset(bha_full, subset = !(Webpage %in% ext_site_names))
  chi_r <- subset(chi_full, subset = !(Webpage %in% ext_site_names))
  cor_r <- subset(cor_full, subset = !(Webpage %in% ext_site_names))
  int_r <- subset(int_full, subset = !(Webpage %in% ext_site_names))
  kld_r <- subset(kld_full, subset = !(Webpage %in% ext_site_names))
  man_r <- subset(man_full, subset = !(Webpage %in% ext_site_names))
  
  alt_r <- alt_r[order(alt_r$UID),]
  bha_r <- bha_r[order(bha_r$UID),]
  chi_r <- chi_r[order(chi_r$UID),]
  cor_r <- cor_r[order(cor_r$UID),]
  int_r <- int_r[order(int_r$UID),]
  kld_r <- kld_r[order(kld_r$UID),]
  man_r <- man_r[order(man_r$UID),]

  #Combine all sites and mutated
  alt_a <- rbind(setNames(alt_r, names(alt_e)),setNames(alt_m, names(alt_e)),alt_e);
  bha_a <- rbind(setNames(bha_r, names(bha_e)),setNames(bha_m, names(bha_e)),bha_e);
  chi_a <- rbind(setNames(chi_r, names(chi_e)),setNames(chi_m, names(chi_e)),chi_e);
  cor_a <- rbind(setNames(cor_r, names(cor_e)),setNames(cor_m, names(cor_e)),cor_e);
  int_a <- rbind(setNames(int_r, names(int_e)),setNames(int_m, names(int_e)),int_e);
  kld_a <- rbind(setNames(kld_r, names(kld_e)),setNames(kld_m, names(kld_e)),kld_e);
  man_a <- rbind(setNames(man_r, names(man_e)),setNames(man_m, names(man_e)),man_e);



  #Threshold for all sites and mutated
  # alt <- run_for_file(alt_a,man_a,step,"alt");
  # bha <- run_for_file(bha_a,man_a,step,"bha");
  # chi <- run_for_file(chi_a,man_a,step,"chi");
  # cor <- run_for_file(cor_a,man_a,step,"cor");
  # int <- run_for_file(int_a,man_a,step,"int");
  # kld <- run_for_file(kld_a,man_a,step,"kld");
  
  
  results_df <- data.frame("Approach"=character(),"Measure"=character(),"Threshold"=numeric(),"Optimal Threshold"=numeric(),"Initial Webpages Accuracy"=numeric(),"Additional Webpages Accuracy"=numeric(), "Mutated Webpages Accuracy"=numeric(), "All Webpages Accuracy"=numeric(), "All Webpages Optimal Accuracy"=numeric(), stringsAsFactors = FALSE);
  print("*********************************");
  print(paste("Approach:",approachName));
  print(paste("Original--Directory:",directory));
  print(paste("Mutatuted-Directory:",mdirectory));
  print("*********************************");
  
  
  result <- process("alt", alt_r, alt_a, alt_e, alt_m, man_r, man_a, man_e, man_m);
  results_df[nrow(results_df) + 1,] = c(approachName, mTypes["alt"], result);
  
  result <- process("bha", bha_r, bha_a, bha_e, bha_m, man_r, man_a, man_e, man_m);
  results_df[nrow(results_df) + 1,] = c(approachName, mTypes["bha"], result);
  
  result <- process("chi", chi_r, chi_a, chi_e, chi_m, man_r, man_a, man_e, man_m);
  results_df[nrow(results_df) + 1,] = c(approachName, mTypes["chi"], result);
  
  result <- process("cor", cor_r, cor_a, cor_e, cor_m, man_r, man_a, man_e, man_m);
  results_df[nrow(results_df) + 1,] = c(approachName, mTypes["cor"], result);
  
  result <- process("int", int_r, int_a, int_e, int_m, man_r, man_a, man_e, man_m);
  results_df[nrow(results_df) + 1,] = c(approachName, mTypes["int"], result);
  
  result <- process("kld", kld_r, kld_a, kld_e, kld_m, man_r, man_a, man_e, man_m);
  results_df[nrow(results_df) + 1,] = c(approachName, mTypes["kld"], result);
  
  
  
  return(results_df);
  
}
run_for_directory_by_parsing_using_threshold <- function(directory, mdirectory){
  alt_t <- 0.84
  bha_t <- 0.24
  chi_t <- 0.46
  cor_t <- 1;
  int_t <- 0.15;
  kld_t <- 1.85;
  
  alt_full <- read.csv(file.path(input_dir, directory, "Alternative-Chi-Square.csv"), stringsAsFactors = FALSE);
  bha_full <- read.csv(file.path(input_dir, directory, "Bhattacharyya.csv"), stringsAsFactors = FALSE);
  chi_full <- read.csv(file.path(input_dir, directory, "Chi-Square.csv"), stringsAsFactors = FALSE);
  cor_full <- read.csv(file.path(input_dir, directory, "Correlation.csv"), stringsAsFactors = FALSE);
  int_full <- read.csv(file.path(input_dir, directory, "Intersection.csv"), stringsAsFactors = FALSE);
  kld_full <- read.csv(file.path(input_dir, directory, "Kullback-Leibler-Divergence.csv"), stringsAsFactors = FALSE);
  man_full <- read.csv(file.path(input_dir, "manual-report-all-failures.csv"), stringsAsFactors = FALSE);
  man_full <- man_full[man_full$Failure == "Small-Range" | man_full$Failure == "small-range",]
  kld_full <- absolute_columns(kld_full);

  alt_m <- read.csv(file.path(input_dir, mdirectory, "Alternative-Chi-Square.csv"), stringsAsFactors = FALSE);
  bha_m <- read.csv(file.path(input_dir, mdirectory, "Bhattacharyya.csv"), stringsAsFactors = FALSE);
  chi_m <- read.csv(file.path(input_dir, mdirectory, "Chi-Square.csv"), stringsAsFactors = FALSE);
  cor_m <- read.csv(file.path(input_dir, mdirectory, "Correlation.csv"), stringsAsFactors = FALSE);
  int_m <- read.csv(file.path(input_dir, mdirectory, "Intersection.csv"), stringsAsFactors = FALSE);
  kld_m <- read.csv(file.path(input_dir, mdirectory, "Kullback-Leibler-Divergence.csv"), stringsAsFactors = FALSE);
  man_m <- read.csv(file.path(input_dir, "manual-report-mutated-pages.csv"), stringsAsFactors = FALSE);
  kld_m <- absolute_columns(kld_m);
  
  alt_m <- alt_m[order(alt_m$UID),]
  bha_m <- bha_m[order(bha_m$UID),]
  chi_m <- chi_m[order(chi_m$UID),]
  cor_m <- cor_m[order(cor_m$UID),]
  int_m <- int_m[order(int_m$UID),]
  kld_m <- kld_m[order(kld_m$UID),]
  man_m <- man_m[order(man_m$UID),]
  
  alt_e <- subset(alt_full, subset = (Webpage %in% ext_site_names))
  bha_e <- subset(bha_full, subset = (Webpage %in% ext_site_names))
  chi_e <- subset(chi_full, subset = (Webpage %in% ext_site_names))
  cor_e <- subset(cor_full, subset = (Webpage %in% ext_site_names))
  int_e <- subset(int_full, subset = (Webpage %in% ext_site_names))
  kld_e <- subset(kld_full, subset = (Webpage %in% ext_site_names))
  man_e <- subset(man_full, subset = (Webpage %in% ext_site_names))
  
  alt_e <- alt_e[order(alt_e$UID),]
  bha_e <- bha_e[order(bha_e$UID),]
  chi_e <- chi_e[order(chi_e$UID),]
  cor_e <- cor_e[order(cor_e$UID),]
  int_e <- int_e[order(int_e$UID),]
  kld_e <- kld_e[order(kld_e$UID),]
  man_e <- man_e[order(man_e$UID),]
  
  
  alt_r <- subset(alt_full, subset = !(Webpage %in% ext_site_names))
  bha_r <- subset(bha_full, subset = !(Webpage %in% ext_site_names))
  chi_r <- subset(chi_full, subset = !(Webpage %in% ext_site_names))
  cor_r <- subset(cor_full, subset = !(Webpage %in% ext_site_names))
  int_r <- subset(int_full, subset = !(Webpage %in% ext_site_names))
  kld_r <- subset(kld_full, subset = !(Webpage %in% ext_site_names))
  man_r <- subset(man_full, subset = !(Webpage %in% ext_site_names))
  
  alt_r <- alt_r[order(alt_r$UID),]
  bha_r <- bha_r[order(bha_r$UID),]
  chi_r <- chi_r[order(chi_r$UID),]
  cor_r <- cor_r[order(cor_r$UID),]
  int_r <- int_r[order(int_r$UID),]
  kld_r <- kld_r[order(kld_r$UID),]
  man_r <- man_r[order(man_r$UID),]
  
  
  #Combine all sites and mutated
  alt_a <- rbind(setNames(alt_r, names(alt_e)),setNames(alt_m, names(alt_e)),alt_e);
  bha_a <- rbind(setNames(bha_r, names(bha_e)),setNames(bha_m, names(bha_e)),bha_e);
  chi_a <- rbind(setNames(chi_r, names(chi_e)),setNames(chi_m, names(chi_e)),chi_e);
  cor_a <- rbind(setNames(cor_r, names(cor_e)),setNames(cor_m, names(cor_e)),cor_e);
  int_a <- rbind(setNames(int_r, names(int_e)),setNames(int_m, names(int_e)),int_e);
  kld_a <- rbind(setNames(kld_r, names(kld_e)),setNames(kld_m, names(kld_e)),kld_e);
  man_a <- rbind(setNames(man_r, names(man_e)),setNames(man_m, names(man_e)),man_e);
  
  
  #Threshold for all sites and mutated
  # alt <- run_for_file(alt_a,man_a,step,"alt");
  # bha <- run_for_file(bha_a,man_a,step,"bha");
  # chi <- run_for_file(chi_a,man_a,step,"chi");
  # cor <- run_for_file(cor_a,man_a,step,"cor");
  # int <- run_for_file(int_a,man_a,step,"int");
  # kld <- run_for_file(kld_a,man_a,step,"kld");
  
  
  results_df <- data.frame(Approach=character(),Measure=character(),Threshold=numeric(),Threshold_Optimal=numeric(),Accuracy=numeric(), Accuracy_Ext=numeric(), Accuracy_Mut=numeric(), Accuracy_All=numeric(), Accuracy_All_Optimal=numeric(), stringsAsFactors = FALSE);
  print("*********************************");
  print(paste("Original--Directory:",directory));
  print(paste("Mutatuted-Directory:",mdirectory));
  print("*********************************");
  
  
  result <- process_using_threshold("alt", alt_r, alt_a, alt_e, alt_m, man_r, man_a, man_e, man_m, alt_t);
  results_df[nrow(results_df) + 1,] = c(directory, "alt", result);
  
  result <- process_using_threshold("bha", bha_r, bha_a, bha_e, bha_m, man_r, man_a, man_e, man_m, bha_t);
  results_df[nrow(results_df) + 1,] = c(directory, "bha", result);
  
  result <- process_using_threshold("chi", chi_r, chi_a, chi_e, chi_m, man_r, man_a, man_e, man_m, chi_t);
  results_df[nrow(results_df) + 1,] = c(directory, "chi", result);
  
  result <- process_using_threshold("cor", cor_r, cor_a, cor_e, cor_m, man_r, man_a, man_e, man_m, cor_t);
  results_df[nrow(results_df) + 1,] = c(directory, "cor", result);
  
  result <- process_using_threshold("int", int_r, int_a, int_e, int_m, man_r, man_a, man_e, man_m, int_t);
  results_df[nrow(results_df) + 1,] = c(directory, "int", result);
  
  result <- process_using_threshold("kld", kld_r, kld_a, kld_e, kld_m, man_r, man_a, man_e, man_m, kld_t);
  results_df[nrow(results_df) + 1,] = c(directory, "kld", result);
  
  
  
  return(results_df);
  
}
run_for_directory <- function(directory, mdirectory, edirectory){
  alt_r <- read.csv(file.path(input_dir, directory, "Alternative-Chi-Square.csv"), stringsAsFactors = FALSE);
  bha_r <- read.csv(file.path(input_dir, directory, "Bhattacharyya.csv"), stringsAsFactors = FALSE);
  chi_r <- read.csv(file.path(input_dir, directory, "Chi-Square.csv"), stringsAsFactors = FALSE);
  cor_r <- read.csv(file.path(input_dir, directory, "Correlation.csv"), stringsAsFactors = FALSE);
  int_r <- read.csv(file.path(input_dir, directory, "Intersection.csv"), stringsAsFactors = FALSE);
  kld_r <- read.csv(file.path(input_dir, directory, "Kullback-Leibler-Divergence.csv"), stringsAsFactors = FALSE);
  man_r <- read.csv(file.path(input_dir, directory, "Manual-Classifications.csv"), stringsAsFactors = FALSE);
  kld_r <- absolute_columns(kld_r);
  
  alt_m <- read.csv(file.path(input_dir, mdirectory, "Alternative-Chi-Square.csv"), stringsAsFactors = FALSE);
  bha_m <- read.csv(file.path(input_dir, mdirectory, "Bhattacharyya.csv"), stringsAsFactors = FALSE);
  chi_m <- read.csv(file.path(input_dir, mdirectory, "Chi-Square.csv"), stringsAsFactors = FALSE);
  cor_m <- read.csv(file.path(input_dir, mdirectory, "Correlation.csv"), stringsAsFactors = FALSE);
  int_m <- read.csv(file.path(input_dir, mdirectory, "Intersection.csv"), stringsAsFactors = FALSE);
  kld_m <- read.csv(file.path(input_dir, mdirectory, "Kullback-Leibler-Divergence.csv"), stringsAsFactors = FALSE);
  man_m <- read.csv(file.path(input_dir, mdirectory, "Manual-Classifications.csv"), stringsAsFactors = FALSE);
  kld_m <- absolute_columns(kld_m);
  
  alt_e <- read.csv(file.path(input_dir, edirectory, "Alternative-Chi-Square.csv"), stringsAsFactors = FALSE);
  bha_e <- read.csv(file.path(input_dir, edirectory, "Bhattacharyya.csv"), stringsAsFactors = FALSE);
  chi_e <- read.csv(file.path(input_dir, edirectory, "Chi-Square.csv"), stringsAsFactors = FALSE);
  cor_e <- read.csv(file.path(input_dir, edirectory, "Correlation.csv"), stringsAsFactors = FALSE);
  int_e <- read.csv(file.path(input_dir, edirectory, "Intersection.csv"), stringsAsFactors = FALSE);
  kld_e <- read.csv(file.path(input_dir, edirectory, "Kullback-Leibler-Divergence.csv"), stringsAsFactors = FALSE);
  man_e <- read.csv(file.path(input_dir, edirectory, "Manual-Classifications.csv"), stringsAsFactors = FALSE);
  kld_e <- absolute_columns(kld_e);

  #Combine all sites and mutated

  alt_a <- rbind(setNames(alt_r, names(alt_e)),setNames(alt_m, names(alt_e)),alt_e);
  bha_a <- rbind(setNames(bha_r, names(bha_e)),setNames(bha_m, names(bha_e)),bha_e);
  chi_a <- rbind(setNames(chi_r, names(chi_e)),setNames(chi_m, names(chi_e)),chi_e);
  cor_a <- rbind(setNames(cor_r, names(cor_e)),setNames(cor_m, names(cor_e)),cor_e);
  int_a <- rbind(setNames(int_r, names(int_e)),setNames(int_m, names(int_e)),int_e);
  kld_a <- rbind(setNames(kld_r, names(kld_e)),setNames(kld_m, names(kld_e)),kld_e);
  man_a <- rbind(setNames(man_r, names(man_e)),setNames(man_m, names(man_e)),man_e);

  #Threshold for all sites and mutated
  # alt <- run_for_file(alt_a,man_a,step,"alt");
  # bha <- run_for_file(bha_a,man_a,step,"bha");
  # chi <- run_for_file(chi_a,man_a,step,"chi");
  # cor <- run_for_file(cor_a,man_a,step,"cor");
  # int <- run_for_file(int_a,man_a,step,"int");
  # kld <- run_for_file(kld_a,man_a,step,"kld");


  results_df <- data.frame(Approach=character(),Measure=character(),Threshold=numeric(),Threshold_Optimal=numeric(),Accuracy=numeric(), Accuracy_Ext=numeric(), Accuracy_Mut=numeric(), Accuracy_All=numeric(), Accuracy_All_Optimal=numeric(), stringsAsFactors = FALSE);
  print("*********************************");
  print(paste("Original--Directory:",directory));
  print(paste("Extension-Directory:",edirectory));
  print(paste("Mutatuted-Directory:",mdirectory));
  print("*********************************");


  result <- process("alt", alt_r, alt_a, alt_e, alt_m, man_r, man_a, man_e, man_m);
  results_df[nrow(results_df) + 1,] = c(directory, "alt", result);

  result <- process("bha", bha_r, bha_a, bha_e, bha_m, man_r, man_a, man_e, man_m);
  results_df[nrow(results_df) + 1,] = c(directory, "bha", result);

  result <- process("chi", chi_r, chi_a, chi_e, chi_m, man_r, man_a, man_e, man_m);
  results_df[nrow(results_df) + 1,] = c(directory, "chi", result);

  result <- process("cor", cor_r, cor_a, cor_e, cor_m, man_r, man_a, man_e, man_m);
  results_df[nrow(results_df) + 1,] = c(directory, "cor", result);

  result <- process("int", int_r, int_a, int_e, int_m, man_r, man_a, man_e, man_m);
  results_df[nrow(results_df) + 1,] = c(directory, "int", result);

  result <- process("kld", kld_r, kld_a, kld_e, kld_m, man_r, man_a, man_e, man_m);
  results_df[nrow(results_df) + 1,] = c(directory, "kld", result);



  return(results_df);

}
process_using_threshold <- function(type, df_r, df_a, df_e, df_m, man_r, man_a, man_e, man_m, threshold){
  print(paste("--------------", type,"--------------"));
  tmp <- run_for_file(df_a,man_a,step,type);
  
  opt_threshold <- tmp[1]
  opt_a_accuracy <- tmp[2]
  
  a_accuracy <- paste(format(round(run_for_file_using_threshold(df_a, man_a, threshold, type), digits = 1), nsmall = 1), "%");
  e_accuracy <- paste(format(round(run_for_file_using_threshold(df_e, man_e, threshold, type), digits = 1), nsmall = 1), "%");
  m_accuracy <- paste(format(round(run_for_file_using_threshold(df_m, man_m, threshold, type), digits = 1), nsmall = 1), "%");
  r_accuracy <- paste(format(round(run_for_file_using_threshold(df_r, man_r, threshold, type), digits = 1), nsmall = 1), "%");
  
  print(paste("Threshold-orig     :", threshold));
  print(paste("Threshold-all-optim:", opt_threshold));
  print(paste("Accuracy-orig      :", r_accuracy));
  print(paste("Accuracy-ext       :", e_accuracy));
  print(paste("Accuracy-mut       :", m_accuracy));
  print(paste("Accuracy-all       :", a_accuracy));
  print(paste("Accuracy-all-optim :", opt_a_accuracy));
  return(c(threshold,opt_threshold,r_accuracy, e_accuracy, m_accuracy, a_accuracy, opt_a_accuracy));
}
process <- function(type, df_r, df_a, df_e, df_m, man_r, man_a, man_e, man_m){
  print(paste("--------------", type,"--------------"));
  vec_a <- run_for_file(df_a,man_a,step,type);
  a_accuracy <- paste(format(round(vec_a[2], digits = 1), nsmall = 1), "%");
  vec <- run_for_file(df_r,man_r,step,type);
  accuracy <- paste(format(round(vec[2], digits = 1), nsmall = 1), "%");
  e_accuracy <- paste(format(round(run_for_file_using_threshold(df_e, man_e, vec[1], type), digits = 1), nsmall = 1), "%");
  m_accuracy <- paste(format(round(run_for_file_using_threshold(df_m, man_m, vec[1], type), digits = 1), nsmall = 1), "%");
  r_a_accuracy <- paste(format(round(run_for_file_using_threshold(df_a, man_a, vec[1], type), digits = 1), nsmall = 1), "%");

  #Invert Correlation Threshold
  if(type == "cor"){
    vec[1] <- 1 - vec[1];
    vec_a[1] <- 1- vec_a[1];
  }
  print(paste("Threshold-orig     :", vec[1]));
  print(paste("Threshold-all-optim:", vec_a[1]));
  print(paste("Accuracy-orig      :", accuracy));
  print(paste("Accuracy-ext       :", e_accuracy));
  print(paste("Accuracy-mut       :", m_accuracy));
  print(paste("Accuracy-all       :", r_a_accuracy));
  print(paste("Accuracy-all-optim :", a_accuracy));
  return(c(format(vec[1], nsmall = 2),format(vec_a[1], nsmall = 2),accuracy, e_accuracy, m_accuracy, r_a_accuracy, a_accuracy));
}
test_stvr_first_thresholds_on_mutations <- function(mdirectory){
  alt_m <- read.csv(file.path(input_dir, mdirectory, "Alternative-Chi-Square.csv"), stringsAsFactors = FALSE);
  bha_m <- read.csv(file.path(input_dir, mdirectory, "Bhattacharyya.csv"), stringsAsFactors = FALSE);
  chi_m <- read.csv(file.path(input_dir, mdirectory, "Chi-Square.csv"), stringsAsFactors = FALSE);
  cor_m <- read.csv(file.path(input_dir, mdirectory, "Correlation.csv"), stringsAsFactors = FALSE);
  int_m <- read.csv(file.path(input_dir, mdirectory, "Intersection.csv"), stringsAsFactors = FALSE);
  kld_m <- read.csv(file.path(input_dir, mdirectory, "Kullback-Leibler-Divergence.csv"), stringsAsFactors = FALSE);
  man_m <- read.csv(file.path(input_dir, mdirectory, "Manual-Classifications.csv"), stringsAsFactors = FALSE);
  kld_m <- absolute_columns(kld_m);

  bha_t <- 0.23;
  chi_t <- 1.85;
  alt_t <- 1.2;
  cor_t <- 1;
  int_t <- 0.22;
  kld_t <- 2;
  
  print("** alt **");
  run_for_file_using_threshold(alt_m, man_m, alt_t, "alt");
  print("** bha **")
  run_for_file_using_threshold(bha_m, man_m, bha_t, "bha");
  print("** chi **");
  run_for_file_using_threshold(chi_m, man_m, chi_t, "chi");
  print("** cor **");
  run_for_file_using_threshold(cor_m, man_m, cor_t, "cor");
  print("** int **");
  run_for_file_using_threshold(int_m, man_m, int_t, "int");
  print("** kld **");
  run_for_file_using_threshold(kld_m, man_m, kld_t, "kld");
}


latexTable <- function(df, caption, fileName){
  table <- kable(df, "latex",  booktabs = TRUE, caption = caption, row.names = FALSE, linesep = "", align=c("l", "l", rep("r", 7)),col.names  = c("Approach", "Measure", "Experimental", "Optimal", "Initial Pages", "Additional Pages", "Mutated Pages", "All Pages", "Optimal")) %>% 
    row_spec(0, bold = TRUE) %>%
    add_header_above(escape=F, c(rep(" ", 2),"\\\\textbf{Threshold}" = 2,  "\\\\textbf{Accuracy}" = 5)) %>%
    kable_styling(latex_options = c("striped", "scale_down"), position = "center", stripe_index = c(seq(1,nrow(df),2)))
  table <- gsub("\\{tab:\\}",paste("{table:",fileName,"}", sep = ""),table)
  table <- gsub("\\\\begin\\{tabular\\}",paste("\\\\fontsize\\{",firstFontSize,"\\}\\{",secondFontSize,"\\}\\\\selectfont\n\\\\begin\\{tabular\\}", sep = ""),table)
  table <- addAfterLine(table, "begin{table}", "\\captionsetup{justification=normal,singlelinecheck=false}")
  write_file(table, file.path(output_dir, paste(fileName, ".tex", sep = "")))
}

#Main program
print(paste("==========","PROGRAM START","=========="));
  # print("--- STVR Real ---");
  # test_stvr_first_thresholds_on_mutations("horizontal");
  # print("--- STVR Mutated ---");
  # test_stvr_first_thresholds_on_mutations("horizontal-mutated");

# for(dirIndex in 1:length(directories_30)){
#     start_time <- Sys.time()
#     if(dirIndex == 1){
#       timing = "..."
#     }else{
#       timing <- paste(round((round((end_time-start_time),digits = 3) * (length(directories_30)-dirIndex)),digits = 3), "s", sep = "")
#     }
#     print(paste("Progress:", paste(dirIndex,length(directories_30),sep = "/"),"Completes-in:", timing))
#     dir <- file.path("horizontal-vertical-30", directories_30[dirIndex])
#     if(dirIndex == 1){
#       df <- run_for_directory_by_parsing(dir,"stvr-horizontal-vertical-mutation-pages");
#     }else{
#       temp <- run_for_directory_by_parsing(dir, "stvr-horizontal-vertical-mutation-pages");
#       df <- rbind(df,temp);
#     }
#     end_time <- Sys.time()
# }
# write.csv(df, file.path(output_dir, "best-threshold-30.csv"));

#df <-run_for_directory_by_parsing_using_threshold("stvr-horizontal-vertical-all-pages-min","stvr-horizontal-vertical-mutation-pages")

#Horizontal Vertical threshold finder
dfHV <-run_for_directory_by_parsing(file.path("horizontal-vertical","to-calculate-threshold"),file.path("horizontal-vertical","mutated"), "Horizontal & Vertical")
write.csv(df, file.path(output_dir, "horizontal-vertical-best-threshold.csv"));
#Horizontal ONLY threshold finder
dfH <-run_for_directory_by_parsing(file.path("horizontal","to-calculate-threshold"),file.path("horizontal","mutated"), "Horizontal")
write.csv(df, file.path(output_dir, "horizontal-best-threshold.csv"));
df <- rbind(dfH,dfHV)

df <- df[order(df$All.Webpages.Optimal.Accuracy, df$All.Webpages.Accuracy, decreasing = T),]

print(df);
latexTable(df, "Comparing \\smallrangehorizontal vs \\smallrangehorizontalvertical approaches using automatically determined thresholds","thresholds")
# for(i in 1:length(directories)){
#   directory <- directories[i];
#   mdirectory <- mdirectories[i];
#   edirectory <- edirectories[i];
#   if(i == 1){
#     df <- run_for_directory(directory, mdirectory, edirectory);
#   }else{
#     temp <- run_for_directory(directory, mdirectory, edirectory);
#     df <- rbind(df,temp);
#   }
# }



print(paste("==========","PROGRAM END","=========="))

