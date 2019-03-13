#Temporal analysis of contaminats function, taCont()

###The function tacont() has the variables:

* **date** is given as date or year
* **y** is either the contaminant, stable isotope or biometric value
* **censored** is a logical vector (TRUE/FALSE) stating if y is censored or not. TRUE means that the y is censored i.e. <LOQ
* **plot** TRUE/FALSE, TRUE by default.
* **pub** TRUE/FALSE, if TRUE the key statistical values are printed on the plot
* **palmost** TRUE/FALSE, models with p<0.05 are plotted as solid lines, p (p almost) between 0.05 and 0.1 are plotted as dashed lines. palmost is TRUE by default
* **cenPerc** is the percentage of allowed cencored data per year, the default is 80
* **onlyRes** TRUE/FALSE, FALSE by default, if true then the list results only contains the dataframe and not the plot

###The result from the function is **a list with up to four elements**, depending on dataset and the chosen variables in the function *plot* and *onlyRes*

1. ggplot
2. data.frame with key results
3. log-linear regression model for the whole dataset
4. log-linear regression model for the last 10 years of the dataset
