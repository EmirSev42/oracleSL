
### SimKFData.R --- 
#----------------------------------------------------------------------
## Author: Emir S
## Created: Dec 8, 2023
## Version: 2
## Last-Updated: Dec 12, 2023
##           By: 
##     Update #: Dec 12, 2023
# Comments: This function generates competing risk survival data
# for KF (primary, event 1) and mortality (secondary, event 2).
# By default, male, age, lACRc, eGFR, dm and cvd are generated
# with their coefficients taken from empirical data.

# "god" represents variability that we do not capture that causes 
# the primary event, created as N(0,1).
# It's used to artificially increase the number
# of events, which is otherwise far too few.

# the numbers on the input formula correspond to log(effect ratio)
# feel free to adjust in function call as needed.

# if you wish to increase/decrease the number of primary events,
# simply adjust the input for variable "god".

# the included vairables are: 
# male ( if male, 0 else)
# age (numeric)
# dm (diabetes, 1 yes 0 no)
# eGFR (blood GFR measurement, cont.)
# lACRc (log ACR, cont.)
# cvd (cardiovascular disease, 1 eys 0 no)
#----------------------------------------------------------------------


simKFData <- function (n = 500, outcome = "competing.risks", 
  formula = ~f(male, 0.2700271) + f(age,-0.03045921)
  + f(lACRc,0.5364934) + f(eGFR,-0.08338161) + f(dm,0.1823216) + f(cvd,0.1823216) 
  + f(god,5.5), intercept = 0)
{
  require(riskRegression)
  require(survival)
  require(prodlim)
  require(lava)
  
  if(!is.numeric(n)){stop("number of data points n must be numeric")}
  if(!(outcome %in% c("survival", "competing.risks", "binary"))){
    stop("outcome must be one of: survival, competing.risks, binary")}
  
  
  
  male = age = lACRc = eGFR = dm = NULL
  outcome <- match.arg(outcome, c("survival", "competing.risks", 
                                  "binary"))
  m <- lava::lvm()
  # GFR
  lava::distribution(m, ~eGFR) <- lava::normal.lvm(mean = 35, 
                                                 sd = 5)
  # lACRc
  lava::distribution(m, ~lACRc) <- lava::normal.lvm(mean = 1.5, 
                                                 sd = 0.5)
  # AGE
  lava::distribution(m, ~age) <- lava::normal.lvm(mean = 65, 
                                                 sd = 5)
  # PROP MALE
  lava::distribution(m, ~male) <- lava::binomial.lvm(p = c(0.45))
  
  # dmETES
  lava::distribution(m, ~dm) <- lava::binomial.lvm(p = c(0.4))
  
  # CVD
  lava::distribution(m, ~cvd) <- lava::binomial.lvm(p = c(0.4))
  
  # GOD
  # this is a variable to artificially create more event 1's
  # it represents whatever causes of mortality we're not capturing
  lava::distribution(m, ~god) <- lava::normal.lvm(mean = 0, 
                                                   sd = 1)
  
  if ("binary" %in% outcome) {
    lava::distribution(m, ~Y) <- lava::binomial.lvm()
    lava::regression(m) <- stats::update(formula, "Y~.")
    lava::intercept(m, ~Y) <- intercept
  }
  if ("survival" %in% outcome) {
    lava::distribution(m, "eventtime") <- lava::coxWeibull.lvm(scale = 1/100)
    lava::distribution(m, "censtime") <- lava::coxWeibull.lvm(scale = 1/100)
    m <- lava::eventTime(m, time ~ min(eventtime = 1, censtime = 0), 
                         "event")
    lava::regression(m) <- stats::update(formula, "eventtime~.")
  }
  if ("competing.risks" %in% outcome) {
    lava::distribution(m, "eventtime1") <- lava::coxWeibull.lvm(scale = 1/100)
    lava::distribution(m, "eventtime2") <- lava::coxWeibull.lvm(scale = 1/100)
    lava::distribution(m, "censtime") <- lava::coxWeibull.lvm(scale = 1/100)
    m <- lava::eventTime(m, time ~ min(eventtime1 = 1, eventtime2 = 2, 
                                       censtime = 0), "event")
    lava::regression(m) <- stats::update(formula, "eventtime1~.")
  }
  out <- data.table::as.data.table(lava::sim(m, n))
  out[, `:=`(male, factor(male, levels = c("0", "1"), 
                        labels = c("0", "1")))]
  out[, `:=`(dm, factor(dm, levels = c("0", "1"), 
                        labels = c("0", "1")))]
  out[, `:=`(cvd, factor(cvd, levels = c("0", "1"), 
                        labels = c("0", "1")))]
  out[]
}


# # example:
# library(tidyverse)
# df <- simKFData() %>% select(-c("eventtime1","eventtime2","censtime"))
# df %>% head()
# # see number of events
# table(df$event)
# # prop male
# table(df$male)
# # dm
# table(df$dm)
# # cvd
# table(df$cvd)
# # see distributions of cont. variables
# ggplot(df) + geom_histogram(aes(x = eGFR), col = "white", fill = "black")
# ggplot(df) + geom_histogram(aes(x = lACRc), col = "white", fill = "black")
# ggplot(df) + geom_histogram(aes(x = age), col = "white", fill = "black")
# # fit CSC and see
# riskRegression::CSC(Hist(time,event)~eGFR + lACRc + age + male + dm + cvd, data = df)


