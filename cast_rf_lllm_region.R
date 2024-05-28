## create conda env
# conda create --name R_4.1.1_cast
## activate conda env
# conda activate R_4.1.1_cast
## install R 4.1.1
# conda install -c conda-forge r-base=4.1.1
## start R kernel
# R 
# install.packages("CAST")
# install.packages("caret")
# install.packages("randomForest")

library(CAST); packageVersion("CAST")
library(caret); packageVersion("caret")
library(randomForest); packageVersion("randomForest")

train_data <- read.csv("output/stats/Random_forest_train_LLLM_region.csv", header = TRUE)
columns <- colnames(train_data)
predictors <- columns[! columns %in% c('longest_blade_length_mean', 'region')]

# create indices for targeted cross validation (leave location out)
indices <- CreateSpacetimeFolds(train_data, spacevar = "region", k=5)

# evaluate mtry with tune length of 11; ~30% of predictors in each bag gives lowest RMSE
model_LLM_RMSE_mtry <- train(train_data[,predictors],train_data$longest_blade_length_mean,
      method="rf", ntree = 500, importance=TRUE, tuneLength = 11,
      trControl=trainControl(method="cv", index = indices$index,number=5))

set.seed(10)
model_LLM <- train(train_data[,predictors],train_data$longest_blade_length_mean,
                   method="rf",tuneGrid=data.frame("mtry"=35), ntree = 500, importance=TRUE,
                   trControl=trainControl(method="cv", index = indices$index, number=5))

plot(varImp(model_LLM))

test_data <- read.csv("output/stats/Random_forest_test_LLLM_region.csv", header = TRUE)
test_data <- subset(test_data, select = -region)

pred_LLLM_Region <- predict(model_LLM, test_data)

test_data$pred_LLLM_Region = pred_LLLM_Region
write.csv(test_data, "output/stats/Random_forest_test_LLLM_region_plus_CAST_rf_prediction.csv", row.names = FALSE)
plot(test_data$longest_blade_length_mean, test_data$pred_LLLM_Region)
cor.test(test_data$longest_blade_length_mean, test_data$pred_LLLM_Region) # 87% t = 3.4518, df = 4, p-value = 0.02601
cor_test_spearman <- cor.test(test_data$longest_blade_length_mean, test_data$pred_LLLM_Region, method = "spearman") 

cor_test_df <- data.frame(correlation=cor_test_spearman$estimate,p=cor_test_spearman$p.value)
write.csv(cor_test_df, "output/stats/Random_forest_LLLM_predict_test_BB_cortest_spearman_cast.csv", row.names = FALSE)

Imp_rf_LLLM <- varImp(model_LLM)
Imp_rf_LLLM_df <- as.data.frame(Imp_rf_LLLM$importance)
Imp_rf_LLLM_df$Species_percent_sim <- row.names(Imp_rf_LLLM_df)

# sort by overall importance scores
Imp_rf_LLLM_df <- Imp_rf_LLLM_df[order(Imp_rf_LLLM_df$Overall, decreasing = TRUE), ]
write.csv(Imp_rf_LLLM_df, "output/stats/Random_forest_LLLM_importance_cast.csv", row.names = FALSE)
# subset by 10 most important
Imp_rf_LLLM_top_10 <- Imp_rf_LLLM_df[1:10, ]
write.csv(Imp_rf_LLLM_top_10, "output/stats/Random_forest_LLLM_top_10_importance_cast.csv", row.names = FALSE)
