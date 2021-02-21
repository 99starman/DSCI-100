# supress warning messages
options(warn=-1)

library(tidyverse)
library(repr)
library(caret)
library(GGally)
install.packages('e1071', dependencies=TRUE)


# read dataset
yeast_url <- 'http://archive.ics.uci.edu/ml/machine-learning-databases/yeast/yeast.data'
yeast_data <- read_delim(yeast_url, delim = '  ', col_name = F)


yeast_data = as.data.frame(yeast_data)

# preprocess the dataset
colnames(yeast_data) <- c('name', 'mcg', 'gvh', 'alm', 'mit', 'erl', 'pox', 'vac', 'nuc', 'class')
cols <- 2:9
yeast_data[cols] <- as.numeric(unlist(yeast_data[cols]))

clean_yeast_data <- yeast_data %>%
    select(-name) %>%
    mutate(class = as.factor(class))

# split into training and testing set
set.seed(1234)
training_rows <- clean_yeast_data %>% 
  select(class) %>% 
  unlist() %>%
  createDataPartition(p = 0.70, list = FALSE)

training_set <- clean_yeast_data %>% slice(training_rows)
testing_set <- clean_yeast_data %>% slice(-training_rows)


# Scale Data (Using Only Training Data)
scale_transformer <- preProcess(training_set, method = c("range"))
training_set <- predict(scale_transformer, training_set)
testing_set <- predict(scale_transformer, testing_set)


# exploratory data analysis 1

# Make Summary Table
summary_table <- training_set %>%
    group_by(class) %>%
    summarize(n = n(), mcg_avg = mean(mcg), gvh_avg = mean(gvh), alm_avg = mean(alm), mit_avg = mean(mit), 
              erl_avg = mean(erl), pox_avg = mean(pox), vac_avg = mean(vac), nuc_avg = mean(nuc))
    
head(training_set)
head(summary_table)


# exploratory data analysis 2

# Visualise variable relationships with ggpairs
options(repr.plot.width = 15, repr.plot.height = 15) 
ggp <- ggpairs(training_set, ggplot2::aes(colour=class), progress = F)
ggp


training_set <- training_set[,c(1:9)]
training_set <- as.data.frame(training_set)


# Number of each class
table(training_set$class )

# Up-sampling to balance class
balanced_training_set <- upSample(x=select(training_set, -class),
                                  y= select(training_set, class) %>% unlist())  
balanced_training_set %>% 
    group_by(Class) %>%
    summarize(n = n())    
      

balanced_X_train <- balanced_training_set %>% select(-Class) %>% data.frame()
balanced_Y_train <- balanced_training_set %>% select(Class) %>% unlist()
X_train <- training_set %>% select(-class) %>% data.frame()
Y_train <- training_set %>% select(class) %>% unlist()

# Parameter value selection
ks <- data.frame(k = seq(from = 1, to = 99, by = 2))

# Cross-validation
train_control <- trainControl(method="cv", number = 5)

# training using balanced data
set.seed(1234)
classifier_ks <- train(x= X_train, y= Y_train, method="knn", tuneGrid= ks, trControl=train_control)
classifier_ks

accuracy_upsampling <- classifier_ks$results
ks_plot_upsampling <- ggplot(accuracy_upsampling, aes(x = k, y = Accuracy)) +
           geom_point() +
           geom_line()
ks_plot_upsampling

classifier_k2 <- train(x= balanced_X_train, y= balanced_Y_train, method="knn", tuneGrid= data.frame(k=2), trControl=train_control)
classifier_k2

X_test <- testing_set %>% select(-class) %>% data.frame()
Y_test <- testing_set %>% select(class) %>% unlist()

# prediction based on k=2
Y_test_predicted <- predict(object = classifier_k2, X_test)

model_quality_k2 <- confusionMatrix(data = Y_test_predicted, reference = Y_test)
model_quality_k2

# prediction accuracy
model_quality_k2$overall[1]

# Visualization of the result--- a visualized confusion matrix
library(ggplot2)
library(scales)

options(repr.plot.width=16, repr.plot.height=12)

ggplotConfusionMatrix <- function(m){
  mytitle <- paste("Accuracy", percent_format()(m$overall[1]))
                   
  p <-
    ggplot(data = as.data.frame(m$table) ,
           aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = Freq), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    geom_text(aes(x = Reference, y = Prediction, label = Freq)) +
    theme(legend.position = "none") +
    ggtitle(mytitle)
  return(p)
}

ggplotConfusionMatrix(model_quality_k2)

# Visualization of the result--- a filled bar chart

# get balanced accuracy(b_acc)
b_acc=list()
for (i in (1:10)){
    b_acc[i]= model_quality_k2$byClass[i,11]
}

# get predicted instances of each class
count=list()
for (i in (1:10)){
     count[i]= sum(model_quality_k2$table[(10*i-9) : (10*i)])
}

# get number of correct predictions in each class
cls=list("CYT","ERL", "EXC","ME1", "ME2", "ME3", "MIT", "NUC", "POX", "VAC")
stat <- do.call(rbind, Map(data.frame, Class=cls, Balanced_Accuracy= b_acc, Count= count))
stat <- mutate(stat, correct= round(Balanced_Accuracy * Count))

df1 <- stat %>% select(Class, Count)
df2 <- stat %>% select(Class, correct)

# plot showing proportion of correct prediction in each class

ggplot()+
  geom_bar(aes(x= Class, y= Count, fill="a"),stat="identity",data= df1, position ="identity",alpha=.4 ) +
  geom_bar(aes(x= Class, y= correct , fill= "b"),stat="identity", data= df2, position ="identity",alpha=.8  )+
  scale_fill_manual(name = 'data source', 
                     values =c('a'='lightblue', 'b'='red'),
                     labels = c('wrong prediction','correct prediction')) + 
  ggtitle("Correct prediction count in each class based on balanced accuracy")+
  theme(plot.title = element_text(size = 24, face = "bold"),
       axis.title=element_text(size=18,face="bold"),
       legend.text=element_text(size=18))
 

# comparison(baseline1)
# An alternative approach to balancing data: SMOTE

install.packages("DMwR")
library(DMwR)

smo_df  <- training_set %>% mutate_if(is.factor, as.character) 

group1 <- filter(smo_df, class == " CYT" | class == " VAC") %>% mutate(class= as.factor(class))
group2 <- filter(smo_df, class == " NUC" | class == " POX") %>% mutate(class= as.factor(class))
group3 <- filter(smo_df, class == " MIT" | class == " ERL") %>% mutate(class= as.factor(class))
group4 <- filter(smo_df, class == " EXC" | class == " ME1") %>% mutate(class= as.factor(class))
group5 <- filter(smo_df,  class == " ME2" | class == " ME3") %>% mutate(class= as.factor(class))

smo_1 <- SMOTE(class ~ ., group1, perc.over = 600, perc.under = 200) %>% na.omit()
smo_2 <- SMOTE(class ~ ., group2, perc.over = 600, perc.under = 200) %>% na.omit()
smo_3 <- SMOTE(class ~ ., group3, perc.over = 1000, perc.under = 300) %>% na.omit()
smo_4 <- SMOTE(class ~ ., group4, perc.over = 200) %>% na.omit()
smo_5 <- SMOTE(class ~ ., group5, perc.over = 100) %>% na.omit()

# new balanced training set
smo_training_set <- do.call("rbind", list(smo_1, smo_2, smo_3, smo_4, smo_5))
table(smo_training_set$class)

smo_X_train <- smo_training_set %>% select(-class) %>% data.frame() 
smo_Y_train <- smo_training_set %>% select(class) %>% unlist() 


set.seed(1234)
classifier_smo_ks <- train(x= smo_X_train, y= smo_Y_train, method="knn", tuneGrid= ks, trControl=train_control)


set.seed(14)
classifier_smo <- train(x= smo_X_train, y= smo_Y_train, method="knn", tuneGrid= data.frame(k=2), trControl=train_control)


# prediction based on k=2, preprocessed using SMOTE
Y_test_predicted_smo <- predict(object = classifier_smo, X_test)

model_quality_k2_smo <- confusionMatrix(data = Y_test_predicted_smo, reference = Y_test)
model_quality_k2_smo

# comparison(baseline2)

# prediction based on k=15, preprocessed using upSampling from caret
set.seed(1234)
classifier_k15 <- train(x= balanced_X_train, y= balanced_Y_train, method="knn", tuneGrid= data.frame(k=15), trControl=train_control)

Y_test_predicted_15 <- predict(object = classifier_k15, X_test)

model_quality_k15 <- confusionMatrix(data = Y_test_predicted_15, reference = Y_test)
model_quality_k15$overall[1]
