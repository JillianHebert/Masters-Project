# Masters Project

For my Masters Project at University of Wisconsin - La Crosse, I investigated the impacts of variable specification and model selection. In many fields, predominantly the sciences, experiments include both fixed and random variables. Typically, within a nested design, many researchers have found their papers returned from peer reviewers due to their variable specification or model selection. This project aims to investigate the differences in results based on variable specification and model selection.

To investigate the differences in results, four model types were compared. The four models differed in the way they handled fixed and random variables. The first model was a linear model that pooled over the random variable. This meant that instead of including the random variable in the model, an average of the response variable was taken by pooling the observations over the random effect. The second model was another linear model that used all variables, but they were all specified as fixed variables in the model, even if the variable was actually a random effect. The third model was a generalized linear model that specified variables as fixed or random and used a maximum likelihood algorithm to calculate results. The fourth and final model was also a generalized linear model that specified variables as fixed or random, but this model used a restricted maximum likelihood algorithm to produce results. 

To further dive into model differences, two simulated data sets were used. These data sets emulate a nested design structure with one significant difference, the treatment effect. One simulated data set has no treatment effect, while the other has 0.5 times the standard deviation of an effect. These two simulated data sets were used to compare the power of the four models. Comparisons were also made between the model's confidence intervals to identify if certain models had smaller interval widths.

The two empirical data sets used for this project both come from U.S. Geological Survey. Due to privacy restrictions, the data sets used could not be posted in this repository. The data file that has been included was used to change the inputs for the simulated data to determine how models change under conditions like sample size, within subject size, and between subject size. 

