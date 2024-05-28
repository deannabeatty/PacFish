# activate conda env
# conda activate R_4.1.1_cast

library(CAST); packageVersion("CAST")
library(caret); packageVersion("caret")
library(randomForest); packageVersion("randomForest")

train_data <- read.csv("output/stats/Random_forest_train_SST_region.csv", header = TRUE)
columns <- colnames(train_data)
predictors <- columns[! columns %in% c('sst', 'region')]

# create indices for targeted cross validation (leave location out)
indices <- CreateSpacetimeFolds(train_data, spacevar = "region", k=5)

# evaluate mtry with tune length of 11
model_SST_RMSE_mtry <- train(train_data[,predictors],train_data$sst,
      method="rf", ntree = 500, importance=TRUE, tuneLength = 11,
      trControl=trainControl(method="cv", index = indices$index, number=5))

set.seed(10)
model_sst <- train(train_data[,predictors],train_data$sst,
                   method="rf",tuneGrid=data.frame("mtry"=35), ntree = 500, importance=TRUE,
                   trControl=trainControl(method="cv", index = indices$index, number=5))
plot(varImp(model_sst))

test_data <- read.csv("output/stats/Random_forest_test_sst_region.csv", header = TRUE)
pred_SST_Region <- predict(model_sst, test_data)

# performs poorly on new data; not significant
test_data$pred_SST_Region = pred_SST_Region
write.csv(test_data, "output/stats/Random_forest_test_sst_region_plus_CAST_rf_prediction.csv", row.names = FALSE)
plot(test_data$sst, test_data$pred_SST_Region)
cor_test_spearman <- cor.test(test_data$sst, test_data$pred_SST_Region, method = "spearman") 
cor_test_df <- data.frame(correlation=cor_test_spearman$estimate,p=cor_test_spearman$p.value)
write.csv(cor_test_df, "output/stats/Random_forest_sst_predict_test_BB_cortest_spearman_cast.csv", row.names = FALSE)

Imp_rf_SST <- varImp(model_sst)
Imp_rf_SST_df <- as.data.frame(Imp_rf_SST$importance)
Imp_rf_SST_df$Species_percent_sim <- row.names(Imp_rf_SST_df)

# sort by overall importance scores
Imp_rf_SST_df <- Imp_rf_SST_df[order(Imp_rf_SST_df$Overall, decreasing = TRUE), ]
write.csv(Imp_rf_SST_df, "output/stats/Random_forest_SST_importance_cast.csv", row.names = FALSE)
# subset by 10 most important
Imp_rf_SST_top_10 <- Imp_rf_SST_df[1:10, ]
write.csv(Imp_rf_SST_top_10, "output/stats/Random_forest_SST_top_10_importance_cast.csv", row.names = FALSE)
