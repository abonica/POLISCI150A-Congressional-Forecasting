library(Matrix)
library(ggplot2)
library(boot)
library(caret)
library(GGally)
library(doMC)
library(gbm)
library(randomForest)
library(RCurl)

##MULTI-CORE
registerDoMC() ##for MAC OS and Linux

##Set Team Name (Required to output model predictions)
team.name <- 'team_name' ##team name with underscores in place of spaces.

###############################################################################
##PLOTTING/EVALUATION FUNCTIONS
###############################################################################
pred.plot <- function(y,pv1,mtitle=''){
    if(max(y,na.rm=T) < 10){
        y <- y * 100
    }
    if(max(pv1,na.rm=T) < 10){
        pv1 <- pv1 * 100
    }
    if(min(y,na.rm=T) > -1){
        pv1 <- pv1 - 50
        y <- y - 50
    }

    require('ggplot2')
    tt <- table(ifelse(pv1 >= 0,1,0),
                ifelse(y >= 0,1,0)
                )
    cc <- sum(diag(tt))/sum(tt)
    if(sum(tt)<435){
        cc <- 1-(sum(tt)-sum(diag(tt)))/435
    }

    gg <- cbind(pv1,
                y,
                ifelse(y >= 0,
                        ifelse(pv1 >= 0,'blue','cyan'),
                       ifelse(pv1 < 0,'red','pink')))
    gg <- as.data.frame(gg)
    colnames(gg) <- c('pv','tv','col')
    gg$pv <- as.numeric(as.character(gg$pv))
    gg$tv <- as.numeric(as.character(gg$tv))
    mode(gg$tv) <- 'numeric'
    p <- qplot(x = pv, y=tv,data= gg,alpha=.8,
               main=mtitle,size=1.5,
               colour=as.character(gg$col))
    p <- p + ylab('Dem. Two Party Vote Share (Actual)')
    p <- p + xlab(paste('Dem. Two Party Vote Share (Predicted) \n',
                        '% Correctly Predicted:',round(cc*100,1)))
    p <- p + scale_colour_manual(values = c("red"="darkred",
                                     "pink" = "red",
                                     "blue"= "blue",
                                     "cyan" = "skyblue"))
    p <- p + scale_size_identity()
    p <- p + theme_bw()
    p <- p + theme(legend.position="none")

    mmin <- min(c(pv1,y))
    mmax <- max(c(pv1,y))
    mmin <- -50
    mmax <- 50
    p <- p + xlim(mmin,mmax) + ylim(mmin,mmax)

    print(p)
    return(p)
}


get.cv.pred.plots <- function(mod.in=NA,
                              y.train.in=NA,x.train.in=NA,
                              y.test.in=NA,x.test.in=NA,
                              ttl='',
                              add.top.features=FALSE,
                              pdfout=FALSE
                              ){


    if(typeof(y.train.in)!='double'){
        print('Requires regression model to generate plots.')
        return(NULL)
    }

    ##CROSS-VALIDATED
    pv.train <- predict(mod.in,
                        as.matrix(x.train.in))
    tt <- table(y.train.in>=0, pv.train>=0 )

    print(paste0('Pct. Correctly Predicted (Training Set): ',round(sum(diag(tt))/sum(tt),3)))
    q1.rf <- pred.plot(y.train.in,pv.train,paste0('Crossvalidated Training Set: ',paste(train.years,collapse=', ')))
    if(pdfout){
        pdf(file=paste0('figures/insample_',ttl,
                        '_num_features_',ncol(x.train.in),'.pdf'),
            width=8,height=8)
        print(q1.rf)
        dev.off()
    }

    ##TEST SET
    if(sum(!is.na(y.test.in))>10){
        pv.test <- predict(mod.in,x.test.in)
        tt <- table(y.test.in>=0 ,pv.test>=0 )
        print(paste0('Pct. Correctly Predicted (Test Set): ',round(sum(diag(tt))/sum(tt),3)))
        print(cat(c('Pct. Correctly Predicted (Out-of-Sample): ',sum(diag(tt))/sum(tt),'\n')))
        q2.rf <- pred.plot(y.test.in,pv.test,paste0('Test Set: ',test.years))

        if(pdfout){
            pdf(file=paste0('figures/oob_',ttl,
                    '_num_features_',ncol(x.train.in),'.pdf'),
                width=8,height=8)
            print(q2.rf)
            dev.off()
        }

        require(gridExtra)
        print(grid.arrange(q1.rf, q2.rf, ncol=2))

        X.in <- X
        X.in$pv.test <- pv.test[match(X.in$ICPSR,rownames(x.test.in))]

        oo <- cbind(X.in$pv.test,X.in$dem.vote.share,X.in)
        o1 = oo[which(X$ICPSR %in% rownames(x.test.in[which(((y.test.in>= 0)  !=  (pv.test >=0 ))),])),]

        m1 <- match(rownames(o1),X.in$dem.icpsr)
        cn.use <- c('pv.test','dem.vote.share','dem.gen.pct','cycle','district','dem.icpsr','rep.icpsr','dem.gen.winner')
        o2 <- X.in[m1,cn.use]
        o2$dem.name <- cands$Name[match(o2$dem.icpsr,cands$ICPSR)]
        o2$rep.name <- cands$Name[match(o2$rep.icpsr,cands$ICPSR)]

        if(add.top.features){
            ii <- importance(mod.in$finalModel)
            cn.use2 <- names(ii[rev(order(ii))[1:50],])
        cn.use2 <- cn.use2[cn.use2 %in% colnames(o1)][1:10]
            oo <- cbind(o2,o1[,cn.use2])
        }else{
            oo <- o2
        }
        od <- abs(oo[,1]-oo[,2])
        oo <- oo[rev(order(od)),]
        print(oo)
        return(oo)
    }
}


##This function will return the top n PACs/donors ranked by the number of contributions
get.n.top.donors <- function(X.in=X,diff.cm.in=diff.cm,n=100,subset.on.bipart.donors=F){
    if(n==0){
        return(X.in)
    }

    if(subset.on.bipart.donors){
        tmp.d <- diff.cm.in
        tmp.d[tmp.d <0] <- 0
        tmp.r <- diff.cm.in
        tmp.r[tmp.r > 0] <- 0
        gave.to.both.parties <- colSums(tmp.d)>0 & colSums(abs(tmp.r))>0
        diff.cm.in <- diff.cm.in[,gave.to.both.parties]
    }


    diff.cm.out <- data.matrix(diff.cm.in[match(rownames(X.in),rownames(diff.cm.in)),1:n])
    return(cbind(X.in,diff.cm.out))
}

##This function generates predictions for the 2018 midterm elections
pred.contests <- function(mod.in=NA){

    test.set <- x.train[test.index,]
    test.set <- test.set[complete.cases(test.set),]
    preds <- predict(mod.in,test.set)

    if(typeof(preds)=='double'){
        wins <- ifelse(as.numeric(preds>0),'D','R')
        rval <- data.frame(pred.dem.vote.share=round(preds,3),winner=wins)
    }else{
        pred.prob <- predict(mod.in, test.set, type="prob")[,1]
        wins <- preds
        rval <- data.frame(pred.prob.dem.wins=round(pred.prob,3),pred.winner=wins)
    }
    print(table(wins))
    rownames(rval) <- rownames(test.set)
    return(rval)
}


##This function generates the output in a consistent format
output.predictions <- function(preds.in=NA, ct =cong.trainset){
    x <- cong.trainset$X
    x <- x[match(rownames(preds.in),rownames(x)),]
    pwin <- preds.in[,2]

    pred.women <- ifelse(pwin =='D',x$dem.female,x$rep.female)
    pred.dime <- ifelse(pwin =='D',x$dem.ip,x$rep.ip)

    print(paste0('Predicted Number of Women in 116th Congress: ', sum(pred.women,na.rm=T) ))

    print(aggregate(pred.dime,list(pwin),mean,na.rm=T))

    preds.in[,1] <- round(preds.in[,1],3)

    rval <- data.frame(cbind(x[,c('dem.name','dem.rid','rep.name','rep.rid')],preds.in))
    rval <- rval[rev(order(rval[,5])),]
    return(rval)

}

##Print the cross-tab of Cook expert ratings against the model predictions
compare.predictions.with.cook.expert.ratings <- function(preds.in=NA, ct =cong.trainset){
    x <- cong.trainset$expert.ratings
    x <- x[match(rownames(preds.in),x$distcyc),]
    pwin <- preds.in[,2]
    print(table(x$Cook.Political.Report,pwin))
}


##Drop any races where a candidate is running unopposed
drop.unopposed <- function(X.in){
    X.in <- X.in[which(is.na(X.in$dem.vote.share) | abs(X.in$dem.vote.share) <= .49),]
    if(!is.null(X.in$uncontested.dem)){
        X.in <- X.in[(substr(rownames(X.in),1,4) %in% test.years |
                      (X.in$uncontested.dem==0 & X.in$uncontested.rep==0)),]
    }
    return(X.in)
}


##Check for new version of dataset and load into R
get.cong.trainset <- function(){
    if(file.exists('cong_training_set.rdata')){##load from disk if file exists
        url.in <- "https://www.dropbox.com/s/jqrpzlbw7i7e0q2/cong_training_set.rdata?dl=1"
        res <- url.exists(url.in, .header=TRUE)
        ##Check for newer file
        if(abs(as.numeric(res['Content-Length']) - as.numeric(file.size('cong_training_set.rdata')))>5){
            print('Newer version of data set available. Download disk? [y/n]')
            ch <- scan(what = character(),n=1)
            if(length(ch)==0){ch='n'}
            if(ch=='y'){
                print('Downloading newest dataset.')
                if(file.exists('cong_training_set.rdata')){
                    load('cong_training_set.rdata')
                    cn1 <- colnames(cong.trainset$house)
                }
                load(url(url.in))
                cn2 <- colnames(cong.trainset$house)
                print('Names of new variables:')
                print(t(cn2[!(cn2 %in% cn1)]))
                save(cong.trainset,file='cong_training_set.rdata')
            }else{
                print('Using local version of dataset.')
                load('cong_training_set.rdata')
            }
        }else{
            print('Using local version of dataset.')
            load('cong_training_set.rdata')
        }
    }else{
        load(url(url.in))
    }
    return(cong.trainset)
}



###############################################################################
##Load data set
###############################################################################
setwd('~/midterm/') #Set working directory

##Check for newer version of data set online and load
cong.trainset <- get.cong.trainset()

X1 <- cong.trainset$house  #House data set
diff.cm <- cong.trainset$diff.cm #Individual donor columns

###############################################################################
##Construct/tranform predictor variables
###############################################################################


##Put any code used to transform, create, or add variables in this function.
create.add.features <- function(x=X1){
    ##add/create new features to x which is returned at end of function.

    ##The difference in log total spending
    x$log.spending.diff <- log1p(x$total.disbursements.dem)-log1p(x$total.disbursements.rep)

    ##The difference between log number of donors giving to the D and R candidates
    x$log.donor.diff <- log1p(x$num.distinct.donors.dem)-log1p(x$num.distinct.donors.rep)

    ##Difference between the number of D and R candidates running for House seats
    x$n.cands.diff <- x$num.cands.in.cycle.dem-x$num.cands.in.cycle.rep

    ##Proportion of votes in primaries going to Democrats
    x$dem.tot.prim.prop <- x$total.prim.votes.cast.dem/(x$total.prim.votes.cast.dem +x$total.prim.votes.cast.rep)

    ##ADD CODE HERE

    return(x)
}

##e.g. constructing new variables
X1 <- create.add.features(X1)

###############################################################################
##Select predictor variables
###############################################################################

##Which election cycles to include in training set?
train.years <- c(2004,2006,2008,2010,2012,2014,2016)
test.years <- 2018
train.years <- train.years[!(train.years %in% test.years)]

##Target variable to predict
cols.main <- c("dem.vote.share") #Target vector--outcome to predict

##Specify with high-level features to use
cols <- c(cols.main,
          'generic.ballot',
          'dist.pres.vs',
          'dem.inc',
          'rep.inc',
          'openseat',
          'log.tot.diff',
          'log.pac.diff',
          'log.donor.diff',
          'log.spending.diff',
          'dem.tot.prim.prop',
          'dem.maj.party.h',
          'dem.pres',
          "uncontested.dem",
          "uncontested.rep",
          'midterm'
          )
cols.use <- cols[cols %in% colnames(X1)]

##Subset data set
X <- as.data.frame(cbind(X1[X1$cycle %in% c(train.years,test.years), cols.use]) )

##Specify the number of top donors to include (Set n=0 to exclude individual donors)
X <- get.n.top.donors(n=0, #how many top donors to include
                      subset.on.bipart.donors=T #only include donors that have given to both parties
                      )

##Drop contests that are unopposed. (Usually counterproductive when predicticing binary win/loss outcomes).
##X <- drop.unopposed(X)

###############################################################################
##CREATE TRAINING SET
###############################################################################

##FEATURES
x.train <- data.matrix(X[,!c(colnames(X) %in% c('dem.vote.share')) ])
mode(x.train) <- 'numeric'

##Target vector
y.train <- X[,'dem.vote.share']
names(y.train) <- rownames(x.train)

##Alternative binary target vector
y.train.win <- factor(ifelse(y.train >= 0,'D','R'))

#Create indices for year to include in training/test sets
train.index <- complete.cases(cbind(y.train,x.train)) & substr(rownames(X),1,4) %in% train.years
test.index <- substr(rownames(X),1,4) %in% test.years

##Removing features with little variation
nz <- nearZeroVar(x.train,freqCut=100,saveMetrics=T)
x.train <- x.train[,!nz$nzv]

###############################################################################
##Train model
###############################################################################


set.seed(1234)##for replicability
caret.out <- train(x=x.train[train.index,],
                   y=y.train[train.index],
                   method='gbm',
                   preProcess=c('scale','center'),##preprocessing
                   tuneGrid=expand.grid(interaction.depth = c(1,3, 5, 9),##Tuning parameters
                                        n.trees = seq(1,30,3)*50,
                                        shrinkage = 0.1,
                                        n.minobsinnode = 20),
                   trControl=trainControl(method="repeatedcv", ##Control arguments
                                          number=10,
                                          repeats = 3,
                                          verbose=TRUE))
print(caret.out)

##Note: Pay close attention to variable importance when deciding which features to include.
var.importance <- varImp(caret.out, scale=FALSE)
print(var.importance)
plot(var.importance)

##Predict 2018 contests
pred <- pred.contests(caret.out)

if(test.years==2018){

    ##compare prediction with Cook Report expert ratings
    compare.predictions.with.cook.expert.ratings(pred)

    ##generate output matrix
    output <- output.predictions(preds.in=pred)

    ##Write model predictions to .csv
    write.csv(output,file=paste0('sml_',test.years,'_pred_winners_',team.name,'.csv'))

    ##Save caret output to disk
    save(caret.out,file=paste0('sml_',test.years,'_caret_',team.name,'.rdata'))

}else{
    ##Plot scatterplot comparison. (Only works when predicticting vote shares.)
    cv.pp.out <- get.cv.pred.plots(mod.in=caret.out,
                                   y.train.in = y.train[train.index],
                                   x.train.in = x.train[train.index,],
                                   y.test.in = y.train[test.index],
                                   x.test.in = x.train[test.index,],
                                   ttl=paste(''),
                                   add.top.features=TRUE)
}


