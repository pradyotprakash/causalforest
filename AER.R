
library(causalTree)

# set global parameters
minsize.temp = 40
split.Bucket.temp = T
bucketNum.temp = 10
bucketMax.temp = 100

# Do causal tree estimation
split.Rule.temp = "CT" #CT
cv.option.temp = "CT" #CT
split.Honest.temp = T
cv.Honest.temp = T
split.alpha.temp = .5
cv.alpha.temp = .5
xvalsize = 10
cp.temp = 0

get_causal_tree <- function(forml, dataTrain){
    if(.Platform$OS.type == "unix"){
        sink("/dev/null");
    }
    else{
        sink("nul");
    }

    tree <- causalTree(as.formula(forml), data = dataTrain, treatment = dataTrain$w,
                     split.Rule=split.Rule.temp, split.Honest=T, split.Bucket=split.Bucket.temp, bucketNum = bucketNum.temp, 
                     bucketMax = bucketMax.temp, cv.option=cv.option.temp, cv.Honest=cv.Honest.temp, minsize = minsize.temp, 
                     split.alpha = split.alpha.temp, cv.alpha = cv.alpha.temp, xval=xvalsize, cp = cp.temp)
    sink();
    tree
}

remove_na <- function(data, vec){
    for(v in vec){
        tmp <- data[[v]]
        data <- data[!is.na(tmp),]
    }
    data
}

get_treatment_effect <- function(data_, indices, outcome_var){
    relevant_data <- data_[indices,]

    # the block probabilities can be introduced here to normalize things
    # will need to supply them based on some measure
    treatment_values <- relevant_data[relevant_data$w == 1, outcome_var]
    control_values <- relevant_data[relevant_data$w == 0, outcome_var]
    mean(treatment_values) - mean(control_values)
}

get_model <- function(data, fml, outcome_var, frac_tr, frac_e, B, precision=4, adaptive_sampling=FALSE){

    N <- nrow(data)
    ntr <- frac_tr*N
    ne <- frac_e*N
    nb <- ntr + ne

    trees <- list()
    treatments <- list()
    adaptive_samples <- matrix(1, nrow=B, ncol=nrow(data))

    for(jj in 1:B){
        if(jj %% 50 == 0){
            cat(sprintf("Constructing tree %d\n", jj))
        }

        sub_sample <- sample(N, nb)

        if(adaptive_sampling){
            adaptive_samples[jj, sub_sample] <- 0    
        }

        train_indices <- sample(sub_sample, ntr)
        estimation_indices <- setdiff(sub_sample, train_indices)
        sampled_tr <- data[train_indices, ]
        sampled_e <- data[estimation_indices, ]
        cf <- get_causal_tree(fml, sampled_tr)

        predictions <- predict(cf, newdata=sampled_e)
        predictions <- sapply(predictions, function(t){as.character(round(t, precision))})
        all_unique <- unique(predictions)

        l <- list()
        for(i in all_unique){
            l[[i]] <- vector()
        }

        count = 1
        for(i in predictions){
            l[[i]] <- c(l[[i]], count)
            count <- count + 1
        }

        treatment_effects <- list()
        for(i in names(l)){
            treatment_effects[[i]] <- get_treatment_effect(sampled_e, l[[i]], outcome_var)
        }

        trees[[jj]] <- cf
        treatments[[jj]] <- treatment_effects
    }
    list(trees=trees, treatments=treatments, adaptive_samples=adaptive_samples)
}

make_prediction <- function(model, newdata, precision=4){
    predictions <- vector(length=nrow(newdata))
    trees <- model$trees
    treatments <- model$treatments
    adaptive_samples <- model$adaptive_samples
    B <- length(trees)

    for(i in 1:B){
        cf <- trees[[i]]
        treatment <- treatments[[i]]
        all_predictions <- predict(cf, new_data=newdata)
        all_predictions <- sapply(all_predictions, function(t){v <- as.character(round(t, precision)); as.numeric(treatment[[v]])})
        predictions <- predictions + adaptive_samples[i]*all_predictions
    }

    predictions/colSums(adaptive_samples)
}

do_regression <- function(data, treatment_effects, outcome_var, features){
    threshold <- mean(treatment_effects)
    positive_indices <- treatment_effects > threshold
    negative_indices <- !positive_indices
    positive_data <- data[positive_indices, ]
    negative_data <- data[negative_indices, ]

    modified_data <- cbind(data, positive=as.numeric(positive_indices), negative=as.numeric(negative_indices))
    fml1 <- paste(outcome_var, " ~ 0 + w + w:positive + ", paste(features, collapse=" + "))
    fml2 <- paste(outcome_var, " ~ 0 + w + w:negative + ", paste(features, collapse=" + "))
    h1 <- lm(fml1, data=modified_data)
    h2 <- lm(fml2, data=modified_data)
    list(positive=h1, negative=h2)
}

library(foreign)

features <- c("gender", "race", "g1surban", "g1freelunch", "g1tgen", "g1trace", "w")

data_source <- "kindergarten.dta"
data <- read.dta(data_source)
data <- data[data$flagsg1 == "yes",]
data$y <- (data$g1tlistss + data$g1treadss + data$g1tmathss + data$g1wordskillss)
data$w <- sapply(as.character(data$g1classtype), function(x){if(x == 'SMALL CLASS') return (1) else return (0)})
                                                                     
data <- remove_na(data, c("y", "w", features))

scores <- list()
interaction <- list()
beta <- 5
feat <- head(features, -1)

for(f in feat){
    l1 <- list()
    l2 <- list()
    c <- sample(5:10, 1)

    for(level in levels(data[[f]])){
        l1[[level]] <- sample(seq(1, 10), 1)
        l2[[level]] <- l1[[level]]/c
    }
    scores[[f]] <- l1
    interaction[[f]] <- l2
}

get_scores <- function(row){
#     feat <- head(features, -1)
    feat <- c("gender", "race", "g1surban", "g1freelunch", "g1tgen", "g1trace")

    x <- 0
    for(f in feat){
        x <- x + scores[[f]][[row[[f]]]]
    }
    
    x <- x + as.numeric(row[["w"]])*(beta + rnorm(1, sd=1))
    
#     # interaction terms
#     pairs <- list(c("gender", "race"), c("race", "g1freelunch"), c("gender", "g1freelunch"))
    
#     for(pair in pairs){
#         x <- x + scores[[pair[1]]][[row[[pair[1]]]]] * scores[[pair[2]]][[row[[pair[2]]]]]/(interaction[[pair[1]]][[row[[pair[1]]]]]*interaction[[pair[2]]][[row[[pair[2]]]]])
#     }

    x
}

data$y_tilde <- apply(data, 1, get_scores)
data$noise <- sample(20:30, nrow(data), replace=TRUE)

features <- c("gender", "race", "g1surban", "g1freelunch", "g1tgen", "g1trace")
outcome_var <- "y_tilde"

out <- paste(features, collapse="+")
causal_fml <- paste(outcome_var, " ~ ", out, sep="")
model <- get_model(data, causal_fml, outcome_var, 0.1, 0.1, 100)

treatment_effects <- make_prediction(model, data)
h <- do_regression(data, treatment_effects, outcome_var, features)

h

library(randomForest)

tau <- treatment_effects - mean(treatment_effects)
labels <- factor(2*as.numeric(tau > 0) - 1)

rhs <- paste(features, collapse="+")
forml <- paste("labels ~", rhs)

forest <- randomForest(as.formula(forml), data=data.frame(data, labels=labels), importance=TRUE)

sum(predict(forest, newdata=data) == labels)/length(labels)
# sqrt(sum((predict(forest, newdata=data) - data$y_tilde)^2)/length(labels))

importance(forest)
options(repr.plot.width=4, repr.plot.height=3)
varImpPlot(forest, type=1)

table(data[,"race"])
table(data[,"g1freelunch"])

forml <- paste("y_tilde ~ 0+labels+w+labels:w+", paste(features, collapse="+"))
reg <- lm(as.formula(forml), data=data.frame(data, labels=as.numeric(labels)))
reg
forml

library(caret)
varImp(reg)

class1 <- data[tau > 0, ]
class2 <- data[tau <= 0, ]

# imp_features <- rownames(importance(forest))
l_ <- list()

for(f in features){
    t1 <- table(class1[, f])
    t2 <- table(class2[, f])
    l_[[f]] <- cbind(t1, t1/(t1+t2), t2, t2/(t1+t2))
}

l_
