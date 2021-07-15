
fit <- lm(mpg ~ cyl + hp, data = mtcars)
summary(fit)
##Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept) 36.90833    2.19080  16.847  < 2e-16 ***
## cyl         -2.26469    0.57589  -3.933  0.00048 ***
## hp          -0.01912    0.01500  -1.275  0.21253 


plot(mpg ~ cyl, data = mtcars, xlab = "Cylinders", ylab = "Miles per gallon")
abline(coef(fit)[1:2])

## rounded coefficients for better output
cf <- round(coef(fit), 2) 
eq <- paste0("mpg = ", cf[1], sign(unname(coef(fit)[2])), abs(cf[2]), " cyl ", sign(unname(coef(fit)[3])), abs(cf[3]), " hp")

# printing of the equation
mtext(eq, 3, line=-2)
