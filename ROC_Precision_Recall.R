## This program is based on a blog based on Anshul Kundaje tweet
## http://davemcg.github.io/post/are-you-in-genomics-stop-using-roc-use-pr/
## Better to look ROC and PR call to make conclusions.
setwd("C:/Users/nitish.mishra/Desktop")
library(tidyverse)
clinvar <- read_tsv('clinvar_one_hot_CS_toy_set.tsv')
clinvar $Status <- factor(clinvar$Status, levels=c('Pathogenic','NotPathogenic'))
clinvar$Status %>% table()
clinvar %>% select(-Status) %>% colnames()

clinvar$index <- seq(1:nrow(clinvar))
set.seed(9253)
train_set <- clinvar %>% group_by(Status) %>% sample_frac(0.5)
test_set <- clinvar %>% filter(!index %in% train_set$index)
# remove index so the models don't use them to classify
train_set <- train_set %>% select(-index)
test_set <- test_set %>% select(-index)

library(caret)
fitControl_naive <- trainControl(
  classProbs=T, # we want probabilites returned for each prediction
  method = "cv",
  number = 5,
  summaryFunction = twoClassSummary
)


bglmFit <- train(Status ~ ., data=train_set, 
                 method = 'bayesglm',
                 trControl = fitControl_naive)
rfFit <- train(Status ~ ., data=train_set, 
               method = 'rf',
               trControl = fitControl_naive)
# let's see how a popular pathogenicity score does alone
caddFit <- train(Status ~ ., data=train_set %>% select(Status, cadd_phred), 
                 method = 'glm',
                 trControl = fitControl_naive)


my_models <- list()
my_models$bglm <- bglmFit
my_models$rfFit <- rfFit
my_models$cadd <- caddFit

library(PRROC)
# AUROC
roc_maker <- function(model, data) {
  # new predictions on test set
  # don't use the training set - if you are overfitting you will not get accurate idea of your models merit
  new_predictions <- predict(model, data, type = 'prob') %>%
    mutate(Answers = data$Status, 
           Prediction = case_when(Pathogenic > 0.5 ~ 'Pathogenic', 
                                  TRUE ~ 'NotPathogenic'))
  roc.curve(scores.class0 = new_predictions %>% filter(Answers=='Pathogenic') %>% pull(Pathogenic),
            scores.class1 = new_predictions %>% filter(Answers=='NotPathogenic') %>% pull(Pathogenic),
            curve = T)
}

#bglm
plot(roc_maker(bglmFit, test_set))

roc_data <- data.frame()
for (i in names(my_models)){
  print(my_models[[i]]$method)
  roc <- roc_maker(my_models[[i]], test_set)
  out <- roc$curve[,1:2] %>% data.frame()
  colnames(out) <- c('FPR','Sensitivity')
  out$model <- i
  out$AUC <- roc$auc
  out$'Model (AUC)' <- paste0(i, ' (',round(roc$auc, 2),')' )
  roc_data <- rbind(roc_data, out)
}


roc <- roc_data %>% 
  ggplot(aes(x=FPR, y=Sensitivity, colour=`Model (AUC)`)) + 
  geom_line() + 
  theme_minimal()  + 
  ggtitle('AUROC') +
  ggsci::scale_color_startrek()

roc


cm_maker <- function(model, data, cutoff=0.5, mode = 'sens_spec') {
  new_predictions <- predict(model, data, type='prob') %>%
    mutate(Answers = as.factor(data$Status), Prediction = as.factor(case_when(Pathogenic > cutoff ~ 'Pathogenic', TRUE ~ 'NotPathogenic')))
  confusionMatrix(data = new_predictions$Prediction, reference = as.factor(new_predictions$Answers), mode= mode)
}

for (i in names(my_models)){
  print(i)
  print(cm_maker(my_models[[i]], test_set)$table)
}


# precision recall AUC
pr_maker <- function(model, data) {
  new_predictions <- predict(model, data, type = 'prob') %>%
    mutate(Answers = data$Status, Prediction = case_when(Pathogenic > 0.5 ~ 'Pathogenic', TRUE ~ 'NotPathogenic'))
  pr.curve(scores.class0 = new_predictions %>% filter(Answers=='Pathogenic') %>% pull(Pathogenic),
           scores.class1 = new_predictions %>% filter(Answers=='NotPathogenic') %>% pull(Pathogenic),
           curve = T)
}

pr_data <- data.frame()
for (i in names(my_models)){
  print(my_models[[i]]$method)
  pr <- pr_maker(my_models[[i]], test_set)
  out <- pr$curve[,1:2] %>% data.frame()
  colnames(out) <- c('Recall','Precision')
  out$AUC <- pr$auc.integral
  out$model <- i
  out$'Model (AUC)' <- paste0(i, ' (',round(pr$auc.integral,2),')' )
  pr_data <- rbind(pr_data, out)
}


pr <- pr_data %>% 
  ggplot(aes(x=Recall, y=Precision, colour=`Model (AUC)`)) + 
  geom_line() + 
  theme_minimal() + 
  ggtitle('Precision Recall Curve') +
  ggsci::scale_color_startrek()

pr
