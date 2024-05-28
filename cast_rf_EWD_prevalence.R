# activate conda env
# conda activate R_4.1.1_cast

library(CAST); packageVersion("CAST")
library(caret); packageVersion("caret")
library(randomForest); packageVersion("randomForest")

train_data <- read.csv("output/stats/Random_forest_train_EWD_prevalence_region.csv", header = TRUE)
columns <- colnames(train_data)
predictors <- columns[! columns %in% c('PrevalenceMean', 'Region')]

# create indices for targeted cross validation (leave location out)
indices <- CreateSpacetimeFolds(train_data, spacevar = "Region", k=5)

# evaluate mtry with tune length of 11
model_Prev_RMSE_mtry <- train(train_data[,predictors],train_data$PrevalenceMean,
      method="rf", ntree = 500, importance=TRUE, tuneLength = 11,
      trControl=trainControl(method="cv", index = indices$index, number=5))

set.seed(10)
model_EWD_Prev <- train(train_data[,predictors],train_data$PrevalenceMean,
                   method="rf",tuneGrid=data.frame("mtry"=35), ntree = 500, importance=TRUE,
                   trControl=trainControl(method="cv", index = indices$index, number=5))

plot(varImp(model_EWD_Prev))

test_data <- read.csv("output/stats/Random_forest_test_EWD_prevalence_region.csv", header = TRUE)
test_data <- subset(test_data, select = -Region)

pred_EWD_Prev_Region <- predict(model_EWD_Prev, test_data)

# negative correlation; performs poorly on new data
test_data$pred_EWD_Prev_Region = pred_EWD_Prev_Region
write.csv(test_data, "output/stats/Random_forest_test_EWD_prev_region_plus_CAST_rf_prediction.csv", row.names = FALSE)
plot(test_data$PrevalenceMean, test_data$pred_EWD_Prev_Region)
cor_test_spearman <- cor.test(test_data$PrevalenceMean, test_data$pred_EWD_Prev_Region, method = "spearman") 
cor_test_df <- data.frame(correlation=cor_test_spearman$estimate,p=cor_test_spearman$p.value)
write.csv(cor_test_df, "output/stats/Random_forest_EWD_prev_predict_test_BB_cortest_spearman_cast.csv", row.names = FALSE)
