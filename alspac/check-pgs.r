## expect full 
coef(summary(lm(handed.mom=="left"~pgs.mom, alspac.table)))[2,]
coef(summary(lm(handed=="left"~pgs.child, alspac.table)))[2,]

## expect half
coef(summary(lm(handed=="left"~pgs.mom, alspac.table)))[2,]
coef(summary(lm(handed.mom=="left"~pgs.child, alspac.table)))[2,]
coef(summary(lm(handed.partner=="left"~pgs.child, alspac.table)))[2,]

## expect none
coef(summary(lm(handed.partner=="left"~pgs.mom, alspac.table)))[2,]

## > ## expect full 
## > coef(summary(lm(handed.mom=="left"~pgs.mom, alspac.table)))[2,]
##    Estimate  Std. Error     t value    Pr(>|t|) 
## 0.010956810 0.003365612 3.255518297 0.001138108 
## > coef(summary(lm(handed=="left"~pgs.child, alspac.table)))[2,]
##    Estimate  Std. Error     t value    Pr(>|t|) 
## 0.009986924 0.004019210 2.484797507 0.012988935 
## > ## expect half
## > coef(summary(lm(handed=="left"~pgs.mom, alspac.table)))[2,]
##    Estimate  Std. Error     t value    Pr(>|t|) 
## 0.007112932 0.004114507 1.728744635 0.083905784 
## > coef(summary(lm(handed.mom=="left"~pgs.child, alspac.table)))[2,]
##    Estimate  Std. Error     t value    Pr(>|t|) 
## 0.005276617 0.003416938 1.544253144 0.122578262 
## > coef(summary(lm(handed.partner=="left"~pgs.child, alspac.table)))[2,]
##    Estimate  Std. Error     t value    Pr(>|t|) 
## 0.007452921 0.004287903 1.738127207 0.082271549 
## > ## expect none
## > coef(summary(lm(handed.partner=="left"~pgs.mom, alspac.table)))[2,]
##     Estimate   Std. Error      t value     Pr(>|t|) 
## -0.003808645  0.004456123 -0.854699271  0.392775118 

