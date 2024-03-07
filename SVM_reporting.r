#Import required packages (install first)
library(caret)
library(corrplot)
library(kernlab)

#Find and set working directory
getwd()
setwd("C:\\Users\\file_path\\file_directory")

#Ensure csv files are in the working directory
degron <- read.csv("training_dataset.csv")
query <- read.csv("degron_data_to_predict.csv")

#Generate correlation plot between column values
cor_degron = cor(degron[c(21:27, 28)])
corrplot(cor_degron, type="upper",method="circle", tl.col = 'black', tl.cex=0.75,
         number.cex=0.4, addCoef.col='blue')
corrplot.mixed(cor_degron, lower = 'number', upper = 'circle', tl.col = 'black',
               tl.cex=0.75, number.cex=0.4, tl.pos ="lt")

#Split data into training and testing dataset
#p = fraction to use for training dataset
cutoff <- createDataPartition(degron$Protein.Stability.Index, p=0.010, list=FALSE) 
training <- degron[cutoff,]
testing <- degron[-cutoff,]

#Set training controls
control = trainControl(method="cv", number=10, classProbs = TRUE)

#Train model, may be time intensive
model = train(Protein.Stability.Index ~., data = training, method="svmRadial", 
              tuneLength = 8, preProc = c("center","scale"), trControl=control)
model
plot(model)
print("Training Done")

#Predict testing data with trained model
test=predict(model, degron)

#Plot correlation plot and calculate Pearson correlation coefficient
#for prediction vs observed PSI
plot(testing$Protein.Stability.Index, test)
cor(testing$Protein.Stability.Index, test)

#Export predictions
write.csv(test, "output.csv")

