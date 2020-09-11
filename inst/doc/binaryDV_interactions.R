## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
knitr::opts_chunk$set(echo = TRUE, echo=TRUE, message=FALSE, warning=FALSE)
library(tibble)
library(dplyr)
library(DAMisc)
library(ggplot2)

## -----------------------------------------------------------------------------
set.seed(1234)
df1 <- tibble(
  x1 = as.factor(rbinom(1000, 1, .5)), 
  x2 = as.factor(rbinom(1000, 1, .4)), 
  z = rnorm(1000),
  ystar = 0 + as.numeric(x1 == "1") - as.numeric(x2 == "1") + 
    2*as.numeric(x1=="1")*as.numeric(x2=="1") + z, 
  p = plogis(ystar), 
  y = rbinom(1000, 1, p)
)

mod1 <- glm(y ~ x1*x2 + z, data=df1, family=binomial)

## -----------------------------------------------------------------------------
## make the model matrix for all conditions
X11 <- X10 <- X01 <- X00 <- model.matrix(mod1)
## set the conditions for each of the four different
## scenarios above

## x1 = 1, x2=1
X11[,"x11"] <- X11[,"x21"] <- X11[,"x11:x21"] <- 1 

## x1=1, x2=0
X10[,"x11"] <- 1
X10[,"x21"] <- X10[,"x11:x21"] <- 0

## x1=0, x2=1
X01[,"x21"] <- 1
X01[,"x11"] <- X01[,"x11:x21"] <- 0

## x1=0, x2=0
X00[,"x11"] <- X00[,"x21"] <-  X00[,"x11:x21"] <- 0
 
## calculate the probabilities
p11 <- plogis(X11 %*% coef(mod1))
p10 <- plogis(X10 %*% coef(mod1))
p01 <- plogis(X01 %*% coef(mod1))
p00 <- plogis(X00 %*% coef(mod1))

eff1 <- p11 - p10 - p01 + p00

## -----------------------------------------------------------------------------
i1 <- intEff(mod1, c("x1", "x2"), df1)

## ---- fig.height=6, fig.width=6, out.width="65%", fig.align="center"----------
tibble(e1 = eff1, i1 = i1$byobs$int$int_eff) %>% 
  ggplot(mapping= aes(x=e1, y=i1)) + 
  geom_point(pch=1) + 
  theme_bw() + 
  labs(x="Calculated by Hand", y="intEff Function Output")

## ---- fig.height=6, fig.width=6, out.width="65%", fig.align="center"----------
i1$byobs$int %>% 
  ggplot(mapping=aes(x=int_eff)) + 
  geom_histogram() + 
  theme_bw() + 
  labs(x="Interaction Effect")

## ---- fig.height=6, fig.width=6, out.width="65%", fig.align="center"----------
i1$byobs$int %>% 
  mutate(sig = ifelse(abs(i1$byobs$int$zstat) > 1.96, 1, 0), 
         sig = factor(sig, levels=c(0,1), labels=c("No", "Yes"))) %>% 
  ggplot(mapping=aes(x=int_eff)) + 
  geom_histogram() + 
  theme_bw() + 
  facet_wrap(~sig) + 
  labs(x="Interaction Effect")

## -----------------------------------------------------------------------------
sd1 <- secondDiff(mod1, c("x1", "x2"), df1)
summary(sd1)

## -----------------------------------------------------------------------------
i1$byobs$int[1,]

## -----------------------------------------------------------------------------
i1$atmean

## -----------------------------------------------------------------------------
set.seed(1234)
df2 <- tibble(
  x2 = as.factor(rbinom(1000, 1, .5)), 
  x1 = runif(1000, -2,2), 
  z = rnorm(1000),
  ystar = 0 + as.numeric(x2 == "1") - x1 + 
    .75*as.numeric(x2=="1")*x1 + z, 
  p = plogis(ystar), 
  y = rbinom(1000, 1, p)
)

mod2 <- glm(y ~ x1*x2 + z, data=df2, family=binomial)

## -----------------------------------------------------------------------------
X0 <- X1 <- model.matrix(mod2)
## set the conditions for each of the four different
## scenarios above

## x1 = 1, x2=1
X1[,"x21"] <- 1
X1[,"x1:x21"] <- X1[,"x1"]

## x1=1, x2=0
X0[,"x21"] <- 0
X0[,"x1:x21"] <-  0


b <- coef(mod2)

## print the coefficients to show that the two coefficients
## we want are the second and fifth ones. 
b

## calculate the first effect
e1 <- (b[2] + b[5])*dlogis(X1 %*% b)

## calculate the second effect
e2 <- (b[2] )*dlogis(X0 %*% b)


## calculate the probabilities
eff2 <- e1 - e2

## ---- fig.height=6, fig.width=6, out.width="65%", fig.align="center"----------
i2 <- intEff(mod2, c("x1", "x2"), df2)
ggplot(mapping=aes(y = i2$byobs$int[,1], x=eff2)) + 
  geom_point() + 
  theme_bw() + 
  labs(x="Calculated by hand", y= "intEff Output")

## -----------------------------------------------------------------------------
i2$byobs$int[1,]

## ---- fig.height=6, fig.width=6, out.width="65%", fig.align="center"----------
tmpX1 <- X1[rep(1, 51), ]
tmpX0 <- X0[rep(1, 51), ]
tmpX1[,"x1"] <- tmpX1[,"x1:x21"] <- c(seq(-2, 2, length=50), 1.350536)
tmpX0[,"x1"] <- c(seq(-2, 2, length=50), 1.350536)
tmpX0[,"x1:x21"] <- 0

phat1 <- plogis(tmpX1 %*% b)
phat0 <- plogis(tmpX0 %*% b)
plot.df <- tibble(phat = c(phat0[1:50], phat1[1:50]), 
                  x = rep(seq(-2,2,length=50), 2), 
                  x2 = factor(rep(c(0,1), each=50), levels=c(0,1), labels=c("x2 = 0", "x2 = 1")))

yint1 <- phat1[51] - e1[1]*tmpX1[51, "x1"]
yint0 <- phat0[51] - e2[1]*tmpX0[51, "x1"]


plot.df %>% 
  ggplot(aes(x=x, y=phat, colour = x2)) + 
  geom_line() + 
  scale_colour_manual(values=c("#0072B2", "#D55E00")) + 
  geom_abline(slope=e1[1], intercept=yint1, colour="#D55E00", lty=3) + 
  geom_abline(slope=e2[1], intercept=yint0, colour="#0072B2", lty=3) + 
  theme_bw() + 
  labs(x="x1", y="Predicted Probability of y=1", colour="x2")
  



## ---- fig.height=6, fig.width=12, out.width="100%", fig.align="center"--------
i2$byobs$int %>% 
  mutate(sig = ifelse(abs(zstat) > 1.96, 1, 0), 
         sig = factor(sig, levels=c(0,1), labels=c("No", "Yes"))) %>% 
  ggplot(mapping=aes(x=int_eff)) + 
  geom_histogram() + 
  theme_bw() + 
  facet_wrap(~sig) + 
  labs(x="Interaction Effect") + 
  theme(aspect.ratio=1)

## -----------------------------------------------------------------------------
s2 <- secondDiff(mod2, c("x1", "x2"), df2, vals=list(x1=c(-2,-1), x2=c("0", "1")))
summary(s2)

## -----------------------------------------------------------------------------
set.seed(1234)
df3 <- tibble(
  x2 = runif(1000, -2,2), 
  x1 = runif(1000, -2,2), 
  z = rnorm(1000),
  ystar = 0 + as.numeric(x2 == "1") - x1 + 
    .75*as.numeric(x2=="1")*x1 + z, 
  p = plogis(ystar), 
  y = rbinom(1000, 1, p)
)

mod3 <- glm(y ~ x1*x2 + z, data=df3, family=binomial)

## -----------------------------------------------------------------------------
X <- model.matrix(mod3)
b <- coef(mod3)
e3 <- b[5]*dlogis(X%*% b) + (b[2] + b[5]*X[,"x2"])*(b[3] + b[5]*X[,"x1"])*dlogis(X%*%b)*(1-(2*plogis(X %*% b)))

## -----------------------------------------------------------------------------
i3 <- intEff(mod3, c("x1", "x2"), data=df3)

## -----------------------------------------------------------------------------
e3[1]
i3$byobs$int[1,]

## -----------------------------------------------------------------------------
s3 <- secondDiff(mod3, c("x1", "x2"), data=df3, vals =list(x1=c(-1,0), x2=c(-2,2)))
summary(s3)

