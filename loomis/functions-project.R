
NAToZero = function(x)
{
  for(row in 1:nrow(x)){
    
     for (col in 1:length(x)){
       if(is.na(x[row,col])){
         x[row,col]= 0
       }
     }
  }
return(x)
}
  

#Convert units

ConvertCmToInches = function(x)
{
  for (row in 1:nrow(x)){
    if(x[row, 27] == "cm"){
      for(i in 4:26){
      x[row,i] = x[row,i]*2.5
      }
    }
  }
  return(x)
}

normalDiagnosticPlot = function(x,  normalityTest=TRUE,
                                    showDensity=TRUE,
                                    showNormal=TRUE,
                                    showSDs=FALSE,
                                    showAxis=TRUE
                                )
  {
  xx = na.omit(x);
  x.stats = doStatsSummary(x);
  # x.table = table(x);
  
  # library(KernSmooth); # install.packages("KernSmooth", dependencies=TRUE);
  # bin.count = dpih(xx);
  # mybreaks = 100 * bin.count;
  
  mxlim = c(x.stats$mean - 3.5 * x.stats$sd , 
            x.stats$mean + 3.5 * x.stats$sd );
  h = hist(xx, breaks="Sturges", plot=F);
  mylim = c(0, max(h$counts));
  
  myMain = paste0(  "Histogram (mean: ",
                  round(x.stats$mean,digits=3), 
                  ", sd: ",
                  round(x.stats$sd,digits=3),
                  ")"
                  );
  
  
  
mxlab = "";  
  if(normalityTest)
    {
    isNormal = NULL;
    if(x.stats$shapiro.is.normal$`0.10`) { isNormal = 0.10; }
    if(x.stats$shapiro.is.normal$`0.05`) { isNormal = 0.05; }
    if(x.stats$shapiro.is.normal$`0.01`) { isNormal = 0.01; }
    
    isNormalResult = FALSE;
    if(!is.null(isNormal)) { isNormalResult = TRUE;}
    if(is.null(isNormal)) { isNormal = 0.05;}
    
    mxlab = paste0("Shapiro Normality test at (alpha = ",
                isNormal, ") is ... ",isNormalResult);
    }
  
  
### Histogram  
  hist(xx, breaks="Sturges",  xlim=mxlim, ylim=mylim,
      xlab=mxlab, xaxt='n', main=myMain);
  
  if(showDensity)
    {
    par(new=T); # overlay
  ### Density Plot (remember first reading?)
    plot( density(xx, kernel="epanechnikov") ,
            xlim=mxlim, 
            main="", 
            xlab="", 
            ylab="", 
            xaxt='n', 
            yaxt='n'  
        );
    }
    
  
  if(showNormal)
    {    
    par(new=T); # overlay  
  ### Normal Curve
    xt = seq(-3.5,3.5, length=100);
  			yt = dnorm(xt);
  
  	plot( xt, yt, 
  	      type="l", 
  	      lwd=2, 
  	      col = "red",
  	      axes=F, 
  	      xlab="",
  	      ylab=""
  	    );	
    }
  	
    
  if(showSDs)
    {
  ### vertical lines at sd's of data ...	
    abline(v=x.stats$mean,lwd=4,col="blue");
      abline(v=x.stats$mean - 1 * x.stats$sd , col="green",lwd=3);
      abline(v=x.stats$mean + 1 * x.stats$sd , col="green",lwd=3);
      abline(v=x.stats$mean - 2 * x.stats$sd , col="green",lwd=2);
      abline(v=x.stats$mean + 2 * x.stats$sd , col="green",lwd=2);
      abline(v=x.stats$mean - 3 * x.stats$sd , col="green",lwd=1);
      abline(v=x.stats$mean + 3 * x.stats$sd , col="green",lwd=1);
    }
  
    
  if(showAxis)
    {
  ### axis labels showing the ability to use expression			
  	axis(1, at = -3:3, labels = c( expression("-3"~hat(sigma) ), expression("-2"~sigma ), expression("-1"~hat(s) ), expression(bar(x)), expression("1"~hat(s) ), "2s", c( expression("3"~hat(sigma) ))) );
  		#axis(1, at = -3:3, labels = c("-3s", "-2s", "-1s", "hat(mu)", "1s", "2s", "3s"))		
    }
  			
  
  }

cleanUpData = function(measure) {

measure = measure[ ,1:27]

#convert some columns to factors and standardize case
measure$writing = (sapply(measure$writing, function(x) gsub("\"", "", x)))
measure$writing = factor(tolower(measure$writing));

measure$eye = (sapply(measure$eye, function(x) gsub("\"", "", x)))
measure$eye = factor(tolower(measure$eye));

measure$eye.color = (sapply(measure$eye.color, function(x) gsub("\"", "", x)))
measure$eye.color = factor(tolower(measure$eye.color));
measure$eye[measure$eye=="equal"] = "both";

measure$gender = (sapply(measure$gender, function(x) gsub("\"", "", x)))
measure$gender = factor(tolower(measure$gender));
measure$gender[measure$gender=="f"] = "female";
measure$gender[measure$gender=="m"] = "male";

measure$swinging = (sapply(measure$swinging, function(x) gsub("\"", "", x)))
measure$swinging [startsWith(measure$swinging , "l")] <- "left"
measure$swinging [startsWith(measure$swinging , "r")] <- "right"
measure$swinging = factor(tolower(measure$swinging));

measure$ethnicity = (sapply(measure$ethnicity, function(x) gsub("\"", "", x)))
measure$ethnicity = (tolower(measure$ethnicity));
measure$ethnicity[measure$ethnicity=="asain"] = "asian";
measure$ethnicity[measure$ethnicity=="caucasain"] = "caucasian";
measure$ethnicity[measure$ethnicity=="white italian"] = "caucasian";
measure$ethnicity[measure$ethnicity=='white non-hispanic'] = "caucasian";
measure$ethnicity[measure$ethnicity=='white'] = "caucasian";
measure$ethnicity[measure$ethnicity=='korean'] = "asian";
measure$ethnicity[measure$ethnicity=='filipino'] = "asian";
measure$ethnicity[measure$ethnicity=='japanese'] = "asian";
measure$ethnicity[measure$ethnicity=='chinese'] = "asian";
measure$ethnicity[measure$ethnicity=='laotian'] = "asian";
measure$ethnicity[measure$ethnicity=='black'] = "african american";
measure$ethnicity[measure$ethnicity=='anglo'] = "caucasian";
measure$ethnicity[measure$ethnicity=='caucasian/asian'] = "mixed race";
measure$ethnicity[measure$ethnicity=='asian/latino'] = "mixed race";
measure$ethnicity[measure$ethnicity=='white-filipino'] = "mixed race";
measure$ethnicity[measure$ethnicity=='japanese italian'] = "mixed race";
measure$ethnicity[measure$ethnicity=='latino'] = "hispanic";
measure$ethnicity[measure$ethnicity=='latin american'] = "hispanic";
measure$ethnicity[measure$ethnicity=='pacific islander'] = "other";
measure$ethnicity[measure$ethnicity=='indian'] = "other";
measure$ethnicity[measure$ethnicity=='native american'] = "other";
measure$ethnicity = factor(measure$ethnicity);

return(measure)
}