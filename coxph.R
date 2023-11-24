library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

########################################################
data(ovarian , package="survival")
View(heart)

fit <- coxph(Surv(futime, fustat) ~ age + ecog.ps,  
             data=ovarian) 
temp <- cox.zph(fit) 
print(temp)                   
plot(temp)                   
#######################################################

data1 <- heart_transplant
View(data1)
names(data1)
summary(data1)
data1$new_status <- ifelse(data1$survived == "dead", 0, 1)

data1 <- data1 %>% 
  mutate(new_status = recode(survived, "dead" = 0, "alive" = 1))

data1$new_status <- as.numeric(factor(data1$survived, levels = c("dead", "alive")))
View(data1)
data1[is.na(data1)] <- 0
View(data1)
frequency(transplant$treatment)
count(transplant)
table(data1$transplant)

#s(t) from cox model
cox.mod <- coxph(Surv(survtime,new_status) ~ age + transplant, data = data1)
summary(cox.mod)

coefficients <- coef(cox.mod)
View(coefficients)
cal_lin_pred <- function(coefs, covariate_values) {
  linear_predictor <- sum(coefs * covariate_values)
  return(linear_predictor)
}
View(cal_lin_pred)

# Assuming LP is the calculated linear predictor for the specific individual/group
LP <- c(-0.02324734,-0.3650623)
baseline_survival <- survfit(cox.mod, data1 = data.frame())$surv
survival_probability <- exp(-baseline_survival * exp(LP))
print(survival_probability)

###########################################################################
cox.mod2 <- coxph(Surv(survtime,new_status) ~ age, data=data1)
anova(cox.mod, cox.mod2, test= "LRT" )
cox.mod3 <- coxph(Surv(survtime,new_status) ~ transplant, data=data1)
anova(cox.mod, cox.mod3, test= "LRT" )
###########################################################################

temp <- cox.zph(cox.mod) 
print(temp)                  
plot(temp)
plot(predict(cox.mod), residual(cox.mod, type= "martingale"),
     xlab = "fitted value", ylab= "Martingale residuals",
     main = "Residual plots", las = 1)
martingale_res <- residuals(cox.mod, type = "martingale")
print(martingale_res)

plot(martingale_res, ylab = "Martingale Residuals", xlab = "Observation Order", main = "Martingale Residual Plot")
abline(h = 0, col = "red", lty = 2)


fitted_values <- predict(cox.mod, type = "expected")
martingale_data <- data.frame(Fitted = fitted_values, Martingale = martingale_res)
ggplot(martingale_data, aes(x = Fitted, y = Martingale)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(
    x = "Fitted Values",
    y = "Martingale Residuals",
    title = "Martingale Residual Plot"
  )
ggplot(martingale_data, aes(x = Fitted, y = "deviance")) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(
    x = "Fitted Values",
    y = "Martingale Residuals",
    title = "Martingale Residual Plot"
  )


covariate_name <- "age"

# Compute the scaled Schoenfeld residuals for the selected covariate
sch_res <- residuals(cox.mod, type = "scaledsch")
View(sch_res)

# Get the values of the selected covariate from the dataset
selected_covariate <- data1[[covariate_name]]
View(selected_covariate)

###############################################################################################
# Create a data frame for plotting
plot_data <- data.frame(
  Covariate = selected_covariate,
  ScaledSchoenfeldResiduals = sch_res
)

plot_data <- data.frame(Covariate = selected_covariate, ScaledSchoenfeldResiduals = sch_res)
#############################################################################################

# Create the scaled Schoenfeld residuals vs. covariate plot
ggplot(plot_data, aes(x = Covariate, y = ScaledSchoenfeldResiduals)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(
    x = covariate_name,
    y = "Scaled Schoenfeld Residuals",
    title = paste("Scaled Schoenfeld Residuals vs.", covariate_name),
    subtitle = "Check for linearity"
  )

cox.zph(cox.mod)
par(mfrow=c(2,1))
plot(cox.zph(cox.mod))

par(mfrow=c(1,1))
plot(cox.zph(cox.mod)[1])
abline(h=0, col=2)
plot(cox.zph(cox.mod)[2])
abline(h=1, col=2)

------------------------------------------------------------------
#convert data to a format used by the survival package
heartSurv <- Surv(time=data1$survtime,
                  event=data1$new_status)

#calculate survival curves by fitting survival data to treatments
ht <- survfit(heartSurv~transplant, data= data1,
              type="kaplan-meier")

#survival plots
ggsurvplot(ht, conf.int = TRUE, pval = FALSE,
           risk.table = FALSE,
           legend.labs=c("control", "treatment"),
           legend= c(0.10,0.10), break.time.by=200,
           censor.shape= "1", censor.size= 4,
           palette = c("firebrick","goldenrod1"),
           xlab="Time since diagonsis (in days)",
           ylab="Proportion surviving")
print(ht)

# hazard ratio quantities using survdiff function
htDiff <- survdiff(heartSurv ~ transplant, data= data1)
View(htDiff)
D1 <- htDiff$obs[1] # 4
D2 <- htDiff$obs[2] # 24
E1 <- htDiff$exp[1] # 2.33
E2 <- htDiff$exp[2] # 25.63

#hazard ratio, std error, c.i
HR <- (D1/D2)/(E1/E2) ; HR
SE_lnHR <- sqrt(1/E1 + 1/E2) ; SE_lnHR
L <- log(HR) ; L
Zcrit <- qnorm(0.975) ; Zcrit
lower <- exp(L - Zcrit*SE_lnHR) ; lower
upper <- exp(L + Zcrit*SE_lnHR) ; upper
ci95 <- c(lower=lower, upper=upper) ; ci95

#log rank test
X1 <- (D1-E1)^2/E1 + (D2-E2)^2/E2 ; X1
1 - pchisq(X1,df=1)

#log rank test - using survdiff function (chi-sq calc different)
print(htDiff, digits = 4)

#assess pH graphically
plot(ht,fun="cloglog",col = c("firebrick","goldenrod1"),
     las=1, lwd = 2, ylab= "ln(-ln(Proportional surviving)",
     xlab= "ln(Time since diagnosis)")
  
