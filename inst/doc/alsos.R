## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  message=FALSE, 
  warning=FALSE, 
  fig.retina=3
)

## -----------------------------------------------------------------------------
data(elec92, package="optiscale")

## -----------------------------------------------------------------------------
library(DAMisc)
mod <- alsos(choice ~ party + ideol + econ4yr, data=elec92, scale_dv=TRUE)

## ---- fig.height=8, fig.width=8, fig.align="center", out.width="90%"----------
plot(mod)

## -----------------------------------------------------------------------------
library(boot)
boot.mod <- boot.alsos(choice ~ party + ideol + econ4yr, data=elec92, level=2, R=50)

## ---- fig.height=8, fig.width=8, fig.align="center", out.width="90%"----------
library(ggplot2)
ggplot(boot.mod$data, aes(x=raw_vals, y=os_vals)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=.15) + 
  geom_line() + 
  facet_wrap(~variable, scales="free") + 
  theme_bw() + 
  labs(x="Raw values", y="Optimally scaled values")

