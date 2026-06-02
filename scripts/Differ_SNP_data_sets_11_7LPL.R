
########################################################################################################################################################################
# 日期：2024/11/02
# 内容：最新修改的在训练群体进行GWAS分析的8个策略且是P.level构建SNP矩阵
# LPL: 修改
#******************************************************************************#
########################################################################################################################################################################
source("./gapit_functionsw.txt")

library(openxlsx)
library(tidyverse)
library(data.table)
library(ggplot2)
library(caret)
library(dplyr)
library(caret)
library(Matrix)

GD = read.table(file = "Wheat_geno.txt", header = TRUE)
myY = read.table(file = "Wheat_pheno.txt", header = TRUE)[ ,c(1,2)]   #每次只输入一个性状
GM = read.table(file = "SNPinfor.txt", header = TRUE)

# 最佳模型的确定

myGAPIT = GAPIT(
  Y=myY[ ,c(1,2)], # 第一列是 ID
  GD = GD,
  GM = GM,
  PCA.total = 3,
  model = c("GLM", "MLM", "CMLM", "MLMM", "FarmCPU", "BLINK", "SUPER"),#
  Multiple_analysis = F)


files <- list.files(pattern = "^GAPIT\\.Association\\.Filter_GWAS_results\\.csv$")# 查找以 GAPIT.Association.Filter_GWAS_results.csv 开头的文件，在确定有显著位点情况下才能进行

if (length(files) > 0) {    #检查是否找到文件
  file_to_read <- files[1]  #如果有多个匹配文件，选择第1个
  p.sig_sum <- read.csv(file_to_read)
}

summary_counts <- p.sig_sum %>%  #统计第7列的类别数量
  group_by(across(7)) %>%
  summarise(count = n())

most_common_category <- summary_counts %>% slice(1) %>% pull(1)  #找到数量最多的类别

best_model <- sub("\\..*$", "", most_common_category[[1]])       #从第一列中提取.号之前的字符

dir.create("differ-stra")  #建立新的文件夹，设置不同策略的工作路径，查看路径，执行不同策略
setwd("differ-stra")
getwd()

############################################################
###############编写基于GWAS辅助的GS预测函数#################
############################################################
GWAS_assisted_GS <- function(GD, GM, myY, models, p.levels, n_repeats, n_folds) {
  results <- data.frame(models = character(),repeats = integer(), folds = integer(), PandSNPnumber = numeric(), strategy = character(), r = numeric())  #创建一个数据框来存储结果
  set.seed(123)  # 保证结果可重复
  folds_indices <- createMultiFolds(myY[, 2], k = n_folds, times = n_repeats)
  names(folds_indices) <- gsub("Fold0+([1-9])", "Fold\\1", names(folds_indices)) #正则表达式操作来格式化结果名称。例如，始终移除Fold后面多余的零。
  
 for (m in 1:length(models)) {  #外部循环：重复次数
   #m=2
    model <- models[m] 
  for (repeats in 1:n_repeats) {
    #repeats=1
    for (fold in 1:n_folds) {
      #fold=1
      Fold.Rep <- as.name(paste0("Fold", fold, ".Rep", repeats))  #获取当前折的索引作为测试集
      train_indices <- folds_indices[[Fold.Rep]]
      test_indices <- setdiff(1:nrow(myY),  train_indices)        #除去测试集的索引，剩余的样本作为训练集
      source("./gapit_functions王佳博.txt") 
      myGAPIT = GAPIT(                                            #在训练群体中进行GWAS分析，7个模型均进行，其中有一个是最佳模型
        Y = myY[train_indices, c(1,2)], # 第一列是 ID
        GD = GD,
        GM = GM,
        PCA.total = 3,
        model = c("MLMM", "FarmCPU", "BLINK", "GLM", "MLM", "CMLM", "SUPER"),
        Multiple_analysis = F)
      
      #策略1，随机选取同等数量的SNP，作为整个对照（可以考虑不放在整个策略的循环中）
      #整个群体进行GWAS分析，根据GAPIT默认阈值下的显著位点的数量确定最佳模型（楸树的为BLINK模型）
      #为了简化模型运行的时间，每次循环下都进行7个模型，然后从中筛选最佳模型
      
      #策略2，以上一步最佳模型作为训练群体的GWAS分析模型，根据P值构建不同P.level的SNP子集
      #策略3，以上一步最佳模型作为训练群体的GWAS分析模型，将P最小的作为固定效应，随机选取同等数量的SNP子集
      #策略4，以上一步最佳模型作为训练群体的GWAS分析模型，将P最小的作为固定效应，根据P值构建不同P.level的SNP子集
      #策略5，在训练群体中进行7个模型的GWAS分析，将每个模型下P值最小的位点的并集作为固定效应，随机选取同等数量的SNP子集
      #策略6，在训练群体中进行7个模型的GWAS分析，将每个模型下P值最小的位点的并集作为固定效应，根据P值构建不同P.level的SNP子集
      #策略7，在训练群体中进行7个模型的GWAS分析，将每个模型下P值最小的前5个位点的并集作为固定效应，随机选取同等数量的SNP子集
      #策略8，在训练群体中进行7个模型的GWAS分析，将每个模型下P值最小的前5个位点的并集作为固定效应，根据P值构建不同P.level的SNP子
      #@*特别说明这里同等数量是指：根据不同P值下的SNP数量，在所有标记中随机选取与P值下相同数量的标记，作为对照*
      source("./gapit_functions.txt")
      for (p.level in p.levels) {    #针对不同的 p.level 进行策略选择
        #p.level=0.01
        all_files <- list.files()    #从当前文件夹读取最佳模型的结果
        pattern <- paste0("GAPIT.Association.GWAS_Results.", best_model, ".*")
        target_files <- grep(pattern, all_files, value = TRUE)
        
        for (file in target_files) {  #读取每个符合条件的文件
          file_path <- file.path(getwd(), file)
          best_model_GWAS <- read.csv(file_path)
        }
       
##策略1,随机选取同等数量的SNP，作为整个对照（可以考虑不放在整个策略的循环中）
        p01 = best_model_GWAS[, c(1, 4)]
        selected_indices01 <- which(p01$P.value < p.level)    #根据不同P值筛选标记的索引
        num_selected01 <- length(selected_indices01)               #从第2列开始随机抽取，计算 selected_indices 中的元素数量
        selected_snps01 <- sample(names(GD)[-1], num_selected01)   #提取标记的名称
        #subset_matrix_GD01 <- GD[ , c("Taxa", selected_snps01)]     #构建新的数据框，包括第1列和随机抽取的列
        subset_matrix_GD01 <- GD %>% dplyr::select(Taxa, all_of(selected_snps01))
        subset_matrix_GM01 <- GM[GM[ , 1] %in% colnames(subset_matrix_GD01), ]# 根据筛选的 SNP 构建子集GM矩阵

        stra1 <- GAPIT(
          Y=myY[train_indices,c(1,2)],
          GD=subset_matrix_GD01, 
          GM=subset_matrix_GM01,
          PCA.total=3,
          model="model", 
          file.out=F)
        m1=merge(myY[test_indices,c(1,2)],stra1$Pred[,c(1,3,5,8)],by.x="Taxa",by.y="Taxa")
        R1=cor(m1[,5],m1[,2])^2
        r1=sqrt(R1)
        results <- rbind(results, data.frame(models = model, repeats = repeats, folds = fold, PandSNPnumber =  num_selected01, strategy = "str1", r = r1))  #将策略1的结果添加到数据框中
        
##策略2，根据GWAS的结果选择不同P水平的SNP子集
        p02 = best_model_GWAS[, c(1, 4)]
        selected_indices02 <- which(p02$P.value < p.level)    #根据不同P值筛选标记的索引
        selected_snps_p02 <- p02[selected_indices02, 1]         #根据索引提取标记的名称
        #subset_matrix_p_GD02 <- GD[, c("Taxa", selected_snps_p02)]                  #提取GD矩阵
        subset_matrix_p_GD02 <- GD %>% dplyr::select(Taxa, all_of(selected_snps_p02))
        subset_matrix_p_GM02 <- GM[GM[, 1] %in% colnames(subset_matrix_p_GD02), ]   #提取GM矩阵
       
        stra2 <- GAPIT(
          Y=myY[train_indices,c(1,2)],
          GD=subset_matrix_p_GD02, 
          GM=subset_matrix_p_GM02,
          PCA.total=3,
          model="model", 
          file.out=F)
        
        m2=merge(myY[test_indices,c(1,2)],stra2$Pred[,c(1,3,5,8)],by.x="Taxa",by.y="Taxa")
        R2=cor(m2[,5],m2[,2])^2
        r2=sqrt(R2)
        results <- rbind(results, data.frame(models = model, repeats = repeats, folds = fold, PandSNPnumber = p.level, strategy = "str2", r = r2)) 
        
##策略3,以上一步最佳模型作为训练群体的GWAS分析模型，将P最小的作为固定效应，随机选取同等数量的SNP子集
        #为策略3和策略4准备，挑选出P值最小的位点，作为固定效应
        best_sig_snp_1 <- best_model_GWAS[which.min(best_model_GWAS$P.value), ] #筛选最显著的SNP--筛选最小的 P.value 所对应的第一列的SNP位点编号
        snp_value <- best_sig_snp_1$SNP
        best_sig_matrix_GD_1 <- as.data.frame(GD %>% dplyr::select(Taxa, all_of(snp_value)))   #根据筛选的 SNP 构建子集GD矩阵作为固定效应
        
        p03 = best_model_GWAS[, c(1, 4)]
        selected_indices03 <- which(p03$P.value < p.level)    #根据不同P值筛选标记的索引
        num_selected03 <- length(selected_indices03)               #从第2列开始随机抽取，计算 selected_indices 中的元素数量
        selected_snps03 <- sample(names(GD)[-1], num_selected03)   #提取标记的名称
        #subset_matrix_GD03 <- GD[, c("Taxa", selected_snps03)]     #构建新的数据框，包括第1列和随机抽取的20列
        subset_matrix_GD03 <- GD %>% dplyr::select(Taxa, all_of(selected_snps03))
        subset_matrix_GM03 <- GM[GM[, 1] %in% colnames(subset_matrix_GD03), ]# 根据筛选的 SNP 构建子集GM矩阵
        
        stra3 <- GAPIT(
          Y=myY[train_indices,c(1,2)],
          GD=subset_matrix_GD03, 
          GM=subset_matrix_GM03,
          PCA.total=3,
          CV=best_sig_matrix_GD_1,    #将PCA的结果和显著的SNP位点构成的矩阵作为协变量加入到模型的固定效应中
          model="model",
          SNP.test=FALSE,
          memo="MAS+model",file.out=F)
        
        m3=merge(myY[test_indices,c(1,2)],stra3$Pred[,c(1,3,5,8)],by.x="Taxa",by.y="Taxa")
        R3=cor(m3[,5],m3[,2])^2
        r3=sqrt(R3)
        results <- rbind(results, data.frame(models = model, repeats = repeats, folds = fold, PandSNPnumber = num_selected03, strategy = "str3", r = r3))
        
        
##策略4，以上一步最佳模型作为训练群体的GWAS分析模型，将P最小的作为固定效应，根据P值构建不同P.level的SNP子集
        best_sig_snp_1 <- best_model_GWAS[which.min(best_model_GWAS$P.value), ] #筛选最显著的SNP--筛选最小的 P.value 所对应的第一列的SNP位点编号
        snp_value <- best_sig_snp_1$SNP
        best_sig_matrix_GD_1 <- as.data.frame(GD %>% dplyr::select(Taxa, all_of(snp_value)))   #根据筛选的 SNP 构建子集GD矩阵作为固定效应
        
        p04 = best_model_GWAS[, c(1, 4)]
        selected_indices04 <- which(p04$P.value < p.level)    #根据不同P值筛选标记的索引
        selected_snps_p04 <- p04[selected_indices04, 1]         #根据索引提取标记的名称
        subset_matrix_p_GD04 <- GD %>% dplyr::select(Taxa, all_of(selected_snps_p04))
        subset_matrix_p_GM04 <- GM[GM[, 1] %in% colnames(subset_matrix_p_GD04), ]   #提取GM矩阵
        
        stra4 <- GAPIT(
          Y=myY[train_indices,c(1,2)],
          GD=subset_matrix_p_GD04, 
          GM=subset_matrix_p_GM04,
          PCA.total=3,
          CV=best_sig_matrix_GD_1,    #将PCA的结果和显著的SNP位点构成的矩阵作为协变量加入到模型的固定效应中
          model="model",
          SNP.test=FALSE,
          memo="MAS+model",file.out=F)
        
        m4=merge(myY[test_indices,c(1,2)],stra4$Pred[,c(1,3,5,8)],by.x="Taxa",by.y="Taxa")
        R4=cor(m4[,5],m4[,2])^2
        r4=sqrt(R4)
        results <- rbind(results, data.frame(models = model, repeats = repeats, folds = fold, PandSNPnumber = p.level, strategy = "str4", r = r4))
        
##策略5,在训练群体中进行7个模型的GWAS分析，将每个模型下P值最小的位点的并集作为固定效应，随机选取同等数量的SNP子集
        #为策略5和策略6准备，选取7个模型下P最小的SNP的并集
        p05 = best_model_GWAS[, c(1, 4)]
        selected_indices05 <- which(p05$P.value < p.level)    
        num_selected05 <- length(selected_indices05)               
        selected_snps05 <- sample(names(GD)[-1], num_selected05) 
        subset_matrix_GD05 <- GD %>% dplyr::select(Taxa, all_of(selected_snps05))
        subset_matrix_GM05 <- GM[GM[, 1] %in% colnames(subset_matrix_GD05), ]
        
        all_files <- list.files()                           #获取当前工作目录中的所有文件名
        pattern <- "^GAPIT.Association.GWAS_Results\\..*"   #创建匹配模式，匹配所有以 GAPIT.Association.GWAS_Results. 开头的文件
        U_target_files <- grep(pattern, all_files, value = TRUE) #过滤出符合命名规则的文件
        U_snp_names <- c()                                  #存储所有文件中 P.value 最小的位点名称
        
        for (file in U_target_files) {             #该for循环是读取每个符合条件的文件并提取P.value最小的位点名称
          file_path <- file.path(getwd(), file)
          U_target_data <- read.csv(file_path)                    #假设文件是CSV格式，可以根据实际情况修改读取函数
          min_p_value_index <- which.min(U_target_data$P.value)   #找到P.value最小值对应的 SNP 位点名称
          min_p_value_snp <- U_target_data[min_p_value_index, 1]
          U_snp_names <- c(U_snp_names, min_p_value_snp)          #将SNP位点名称加入集合
        }
        snp_names_unique <- unique(U_snp_names)                  #取所有文件中P.value最小的位点名称的并集
        U_subset_matrix_GD <- as.data.frame(GD %>% dplyr::select(Taxa, all_of(snp_names_unique)))  #构建新的数据框，包括第1列和并集列
        stra5 <- GAPIT(
          Y=myY[train_indices,c(1,2)],
          GD=subset_matrix_GD05, 
          GM=subset_matrix_GM05,
          PCA.total=3,
          CV=U_subset_matrix_GD,    #将PCA的结果和显著的SNP位点构成的矩阵作为协变量加入到模型的固定效应中
          model="model",
          SNP.test=FALSE,
          memo="MAS+model",file.out=F)
        
        m5=merge(myY[test_indices,c(1,2)],stra5$Pred[,c(1,3,5,8)],by.x="Taxa",by.y="Taxa")
        R5=cor(m5[,5],m5[,2])^2
        r5=sqrt(R5)
        
        # 将策略5的结果添加到数据框中
        results <- rbind(results, data.frame(models = model, repeats = repeats, folds = fold, PandSNPnumber =  num_selected05 , strategy = "str5", r = r5))
        
##策略6,在训练群体中进行7个模型的GWAS分析，将每个模型下P值最小的位点的并集作为固定效应，根据P值构建不同P.level的SNP子集
        p06 = best_model_GWAS[, c(1, 4)]
        selected_indices06 <- which(p06$P.value < p.level)    #根据不同P值筛选标记的索引
        selected_snps_p06 <- p06[selected_indices06, 1]         #根据索引提取标记的名称
        subset_matrix_p_GD06 <- GD %>% dplyr::select(Taxa, all_of(selected_snps_p06)) #提取GD矩阵
        subset_matrix_p_GM06 <- GM[GM[, 1] %in% colnames(subset_matrix_p_GD06), ]   #提取GM矩阵
        
        all_files <- list.files()                           #获取当前工作目录中的所有文件名
        pattern <- "^GAPIT.Association.GWAS_Results\\..*"   #创建匹配模式，匹配所有以 GAPIT.Association.GWAS_Results. 开头的文件
        U_target_files <- grep(pattern, all_files, value = TRUE) #过滤出符合命名规则的文件
        U_snp_names <- c()                                  #存储所有文件中 P.value 最小的位点名称
        
        for (file in U_target_files) {             #该for循环是读取每个符合条件的文件并提取P.value最小的位点名称
          file_path <- file.path(getwd(), file)
          U_target_data <- read.csv(file_path)                    #假设文件是CSV格式，可以根据实际情况修改读取函数
          min_p_value_index <- which.min(U_target_data$P.value)   #找到P.value最小值对应的 SNP 位点名称
          min_p_value_snp <- U_target_data[min_p_value_index, 1]
          U_snp_names <- c(U_snp_names, min_p_value_snp)          #将SNP位点名称加入集合
        }
        snp_names_unique <- unique(U_snp_names)                  #取所有文件中P.value最小的位点名称的并集
        U_subset_matrix_GD <- as.data.frame(GD %>% dplyr::select(Taxa, all_of(snp_names_unique)))  #构建新的数据框，包括第1列和并集列
        
       stra6 <- GAPIT(
          Y=myY[train_indices,c(1,2)],
          GD=subset_matrix_p_GD06, 
          GM=subset_matrix_p_GM06,
          PCA.total=3,
          CV=U_subset_matrix_GD,    #将PCA的结果和显著的SNP位点构成的矩阵作为协变量加入到模型的固定效应中
          model="model",
          SNP.test=FALSE,
          memo="MAS+model",file.out=F)
        
        m6=merge(myY[test_indices,c(1,2)],stra6$Pred[,c(1,3,5,8)],by.x="Taxa",by.y="Taxa")
        R6=cor(m6[,5],m6[,2])^2
        r6=sqrt(R6)
        
        # 将策略6的结果添加到数据框中
        results <- rbind(results, data.frame(models = model, repeats = repeats, folds = fold, PandSNPnumber = p.level, strategy = "str6", r = r6))
        
##策略7,在训练群体中进行7个模型的GWAS分析，将每个模型下P值最小的前5个位点的并集作为固定效应，随机选取同等数量的SNP子集
        p07 = best_model_GWAS[, c(1, 4)]
        selected_indices07 <- which(p07$P.value < p.level)    
        num_selected07 <- length(selected_indices07)               
        selected_snps07 <- sample(names(GD)[-1], num_selected07)   
        subset_matrix_GD07 <- GD %>% dplyr::select(Taxa, all_of(selected_snps07))
        subset_matrix_GM07 <- GM[GM[, 1] %in% colnames(subset_matrix_GD07), ]
        
        selected_snps <- c()                           #存储所有文件中选取的 SNP 位点名称
        for (file in U_target_files) {                 #该for循环是读取每个符合条件的文件并选取SNP位点名称
          file_path <- file.path(getwd(), file)
          data <- read.csv(file_path)                  #假设文件是CSV格式，可以根据实际情况修改读取函数
          sorted_data <- data[order(data$P.value), ]   #根据 P.value 从小到大排序
          top_5_snps <- head(sorted_data$SNP, 5)       #选取前5个 SNP 位点名称
          selected_snps <- c(selected_snps, top_5_snps)#将选取的 SNP 位点名称加入集合
        }
        
        snp_names_unique_5 <- unique(selected_snps)    #取所有文件中选取的 SNP 位点名称的并集
        U_5_subset_matrix_GD <- as.data.frame(GD %>% dplyr::select(Taxa, all_of(snp_names_unique_5)))  #构建新的数据框，包括第1列和并集列
        
        stra7 <- GAPIT(
          Y=myY[train_indices,c(1,2)],
          GD=subset_matrix_GD07, 
          GM=subset_matrix_GM07,
          PCA.total=3,
          CV=U_5_subset_matrix_GD,    #将PCA的结果和显著的SNP位点构成的矩阵作为协变量加入到模型的固定效应中
          model="model",
          SNP.test=FALSE,
          memo="MAS+model",file.out=F)
        
        m7=merge(myY[test_indices,c(1,2)],stra7$Pred[,c(1,3,5,8)],by.x="Taxa",by.y="Taxa")
        R7=cor(m7[,5],m7[,2])^2
        r7=sqrt(R7)
        results <- rbind(results, data.frame(models = model, repeats = repeats, folds = fold, PandSNPnumber = num_selected07, strategy = "str7", r = r7))
        
##策略8,在训练群体中进行7个模型的GWAS分析，将每个模型下P值最小的前5个位点的并集作为固定效应，根据P值构建不同P.level的SNP子
        p08 = best_model_GWAS[, c(1, 4)]
        selected_indices08 <- which(p08$P.value < p.level)    #根据不同P值筛选标记的索引
        selected_snps_p08 <- p08[selected_indices08, 1]         #根据索引提取标记的名称
        subset_matrix_p_GD08 <- GD %>% dplyr::select(Taxa, all_of(selected_snps_p08))
        subset_matrix_p_GM08 <- GM[GM[, 1] %in% colnames(subset_matrix_p_GD08), ]   #提取GM矩阵
        
        
        stra8 <- GAPIT(
          Y=myY[train_indices,c(1,2)],
          GD=subset_matrix_p_GD08, 
          GM=subset_matrix_p_GM08,
          PCA.total=3,
          CV=U_5_subset_matrix_GD,    #将PCA的结果和显著的SNP位点构成的矩阵作为协变量加入到模型的固定效应中
          model="model",
          SNP.test=FALSE,
          memo="MAS+model",file.out=F)
        
        m8=merge(myY[test_indices,c(1,2)],stra8$Pred[,c(1,3,5,8)],by.x="Taxa",by.y="Taxa")
        R8=cor(m8[,5],m8[,2])^2
        r8=sqrt(R8)
        results <- rbind(results, data.frame(models = model, repeats = repeats, folds = fold, PandSNPnumber = p.level, strategy = "str8", r = r8))
       
       }
     }
   }
  
  }
  return(results)
}

#执行函数
result <- GWAS_assisted_GS(GD, GM, myY,
                           models =c("gBLUP"), #, "sBLUP", "cBLUP"
                           p.levels <- c(0.01, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5),
                           n_repeats = 5, 
                           n_folds =10)
 write.csv(result, "trait1_stra8_fold10-20250327.csv")




















