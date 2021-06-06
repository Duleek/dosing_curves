#non-parametric test for coefficient of cocordance
#use this to see if the experiment's ranks of drug sensitivity for patients/cell lines
#matches actual drug sensitivity rankings 
#(reported rankings by cancerxgene in the case of cell lines)

install.packages("synchrony")
library(synchrony)


#input data where each sample/cell line is a row
#each column is a measurement (either dose or reported IC50)
#each cell can be a cell viability or dose
#each column should be consistent (ex. all cells in column 2 are cell viability %)
#different columns can use different units (ex. column 2 uses cv%, column 3 uses dose)


#doxorubicin; day 1 exposure and 2 reported IC50s from cancerxgene

dox1 <- read.csv("doxorubicin cell viability day 1.csv", head = TRUE)
head(dox1)

dox_concordance1 = kendall.w (dox1, nrands = 0, type = 1, quiet = FALSE)

#Poor agreement = Less than 0.20
#Fair agreement = 0.21 to 0.40
#Moderate agreement = 0.41 to 0.60
#Good agreement = 0.61 to 0.80
#Very good agreement = 0.81 to 1.00
#Spearman's ranked correlation values can also go to -1

#doxorubicin; day 2 exposure and 2 reported IC50s from cancerxgene

dox2 <- read.csv("doxorubicin cell viability day 2.csv", head = TRUE)
head(dox2)

dox_concordance2 = kendall.w (dox2, nrands = 0, type = 1, quiet = FALSE)



#paclitaxel; day 1 exposure and 1 reported IC50 from cancerxgene

pac1 <- read.csv("paclitaxel cell viability day 1.csv", head = TRUE)
head(pac1)

pac_concordance1 = kendall.w (pac1, nrands = 0, type = 1, quiet = FALSE)


#paclitaxel; day 2 exposure and 1 reported IC50 from cancerxgene

pac2 <- read.csv("paclitaxel cell viability day 2.csv", head = TRUE)
head(pac2)

pac_concordance2 = kendall.w (pac2, nrands = 0, type = 1, quiet = FALSE)


print(dox_concordance1)
print(dox_concordance2)
print(pac_concordance1)
print(pac_concordance2)

#to also get p value, print (c(dox_concordance1))
#to also get p value, print (c(dox_concordance2))
#to also get p value, print (c(pac_concordance1))
#to also get p value, print (c(pac_concordance2))

#compare the 2 IC50s that are reported for dox

ic50 = dox1[,c('IC50_1', 'IC50_2')]
ic50_concordance = kendall.w (ic50, nrands = 0, type = 1, quiet = FALSE)
print(ic50_concordance)
#to also get p value, print (c(ic50_concordance))

#compare individual dox day 1 doses to IC50s

dox1_40 = dox1[,c('X40', 'IC50_1', 'IC50_2')]
dox1_40_concordance = kendall.w (dox1_40, nrands = 0, type = 1, quiet = FALSE)

dox1_80 = dox1[,c('X80', 'IC50_1', 'IC50_2')]
dox1_80_concordance = kendall.w (dox1_80, nrands = 0, type = 1, quiet = FALSE)

dox1_100 = dox1[,c('X100', 'IC50_1', 'IC50_2')]
dox1_100_concordance = kendall.w (dox1_100, nrands = 0, type = 1, quiet = FALSE)

dox1_200 = dox1[,c('X200', 'IC50_1', 'IC50_2')]
dox1_200_concordance = kendall.w (dox1_200, nrands = 0, type = 1, quiet = FALSE)

dox1_500 = dox1[,c('X500', 'IC50_1', 'IC50_2')]
dox1_500_concordance = kendall.w (dox1_500, nrands = 0, type = 1, quiet = FALSE)

dox1_1000 = dox1[,c('X1000', 'IC50_1', 'IC50_2')]
dox1_1000_concordance = kendall.w (dox1_1000, nrands = 0, type = 1, quiet = FALSE)

dox1_5000 = dox1[,c('X5000', 'IC50_1', 'IC50_2')]
dox1_5000_concordance = kendall.w (dox1_5000, nrands = 0, type = 1, quiet = FALSE)

print (c(dox1_40_concordance))
print (c(dox1_80_concordance))
print (c(dox1_100_concordance))
print (c(dox1_200_concordance))
print (c(dox1_500_concordance))
print (c(dox1_1000_concordance))
print (c(dox1_5000_concordance))

#compare individual dox day 2 doses to IC50s

dox2_40 = dox2[,c('X40', 'IC50_1', 'IC50_2')]
dox2_40_concordance = kendall.w (dox2_40, nrands = 0, type = 1, quiet = FALSE)

dox2_80 = dox2[,c('X80', 'IC50_1', 'IC50_2')]
dox2_80_concordance = kendall.w (dox2_80, nrands = 0, type = 1, quiet = FALSE)

dox2_100 = dox2[,c('X100', 'IC50_1', 'IC50_2')]
dox2_100_concordance = kendall.w (dox2_100, nrands = 0, type = 1, quiet = FALSE)

dox2_200 = dox2[,c('X200', 'IC50_1', 'IC50_2')]
dox2_200_concordance = kendall.w (dox2_200, nrands = 0, type = 1, quiet = FALSE)

dox2_500 = dox2[,c('X500', 'IC50_1', 'IC50_2')]
dox2_500_concordance = kendall.w (dox2_500, nrands = 0, type = 1, quiet = FALSE)

dox2_1000 = dox2[,c('X1000', 'IC50_1', 'IC50_2')]
dox2_1000_concordance = kendall.w (dox2_1000, nrands = 0, type = 1, quiet = FALSE)

dox2_5000 = dox2[,c('X5000', 'IC50_1', 'IC50_2')]
dox2_5000_concordance = kendall.w (dox2_5000, nrands = 0, type = 1, quiet = FALSE)

print (c(dox2_40_concordance))
print (c(dox2_80_concordance))
print (c(dox2_100_concordance))
print (c(dox2_200_concordance))
print (c(dox2_500_concordance))
print (c(dox2_1000_concordance))
print (c(dox2_5000_concordance))


#compare individual pac day 1 doses to IC50

pac1_20 = pac1[,c('X20', 'IC50_1')]
pac1_20_concordance = kendall.w (pac1_20, nrands = 0, type = 1, quiet = FALSE)

pac1_50 = pac1[,c('X50', 'IC50_1')]
pac1_50_concordance = kendall.w (pac1_50, nrands = 0, type = 1, quiet = FALSE)

pac1_80 = pac1[,c('X80', 'IC50_1')]
pac1_80_concordance = kendall.w (pac1_80, nrands = 0, type = 1, quiet = FALSE)

pac1_120 = pac1[,c('X120', 'IC50_1')]
pac1_120_concordance = kendall.w (pac1_120, nrands = 0, type = 1, quiet = FALSE)

pac1_200 = pac1[,c('X200', 'IC50_1')]
pac1_200_concordance = kendall.w (pac1_200, nrands = 0, type = 1, quiet = FALSE)

pac1_500 = pac1[,c('X500', 'IC50_1')]
pac1_500_concordance = kendall.w (pac1_500, nrands = 0, type = 1, quiet = FALSE)

pac1_1000 = pac1[,c('X1000', 'IC50_1')]
pac1_1000_concordance = kendall.w (pac1_1000, nrands = 0, type = 1, quiet = FALSE)

print (c(pac1_20_concordance))
print (c(pac1_50_concordance))
print (c(pac1_80_concordance))
print (c(pac1_120_concordance))
print (c(pac1_200_concordance))
print (c(pac1_500_concordance))
print (c(pac1_1000_concordance))

#compare individual pac day 2 doses to IC50

pac2_20 = pac2[,c('X20', 'IC50_1')]
pac2_20_concordance = kendall.w (pac2_20, nrands = 0, type = 1, quiet = FALSE)

pac2_50 = pac2[,c('X50', 'IC50_1')]
pac2_50_concordance = kendall.w (pac2_50, nrands = 0, type = 1, quiet = FALSE)

pac2_80 = pac2[,c('X80', 'IC50_1')]
pac2_80_concordance = kendall.w (pac2_80, nrands = 0, type = 1, quiet = FALSE)

pac2_120 = pac2[,c('X120', 'IC50_1')]
pac2_120_concordance = kendall.w (pac2_120, nrands = 0, type = 1, quiet = FALSE)

pac2_200 = pac2[,c('X200', 'IC50_1')]
pac2_200_concordance = kendall.w (pac2_200, nrands = 0, type = 1, quiet = FALSE)

pac2_500 = pac2[,c('X500', 'IC50_1')]
pac2_500_concordance = kendall.w (pac2_500, nrands = 0, type = 1, quiet = FALSE)

pac2_1000 = pac2[,c('X1000', 'IC50_1')]
pac2_1000_concordance = kendall.w (pac2_1000, nrands = 0, type = 1, quiet = FALSE)

print (c(pac2_20_concordance))
print (c(pac2_50_concordance))
print (c(pac2_80_concordance))
print (c(pac2_120_concordance))
print (c(pac2_200_concordance))
print (c(pac2_500_concordance))
print (c(pac2_1000_concordance))


#cyclophosphamide; day 1 exposure and 1 reported IC50 from cancerxgene

cyclo1 <- read.csv("cyclophosphamide cell viability day 1.csv", head = TRUE)
head(cyclo1)

cyclo1_concordance = kendall.w (cyclo1, nrands = 0, type = 1, quiet = FALSE)
print (c(cyclo1_concordance))

#cyclophosphamide; day 2 exposure and 1 reported IC50 from cancerxgene

cyclo2 <- read.csv("cyclophosphamide cell viability day 2.csv", head = TRUE)
head(cyclo2)

cyclo2_concordance = kendall.w (cyclo2, nrands = 0, type = 1, quiet = FALSE)
print (c(cyclo2_concordance))

#compare individual cyclo day 1 doses to IC50

cyclo1_600 = cyclo1[,c('X600', 'IC50_1')]
cyclo1_600_concordance = kendall.w (cyclo1_600, nrands = 0, type = 1, quiet = FALSE)

cyclo1_10000 = cyclo1[,c('X10000', 'IC50_1')]
cyclo1_10000_concordance = kendall.w (cyclo1_10000, nrands = 0, type = 1, quiet = FALSE)

cyclo1_50000 = cyclo1[,c('X50000', 'IC50_1')]
cyclo1_50000_concordance = kendall.w (cyclo1_50000, nrands = 0, type = 1, quiet = FALSE)

cyclo1_100000 = cyclo1[,c('X100000', 'IC50_1')]
cyclo1_100000_concordance = kendall.w (cyclo1_100000, nrands = 0, type = 1, quiet = FALSE)

cyclo1_200000 = cyclo1[,c('X200000', 'IC50_1')]
cyclo1_200000_concordance = kendall.w (cyclo1_200000, nrands = 0, type = 1, quiet = FALSE)

cyclo1_500000 = cyclo1[,c('X500000', 'IC50_1')]
cyclo1_500000_concordance = kendall.w (cyclo1_500000, nrands = 0, type = 1, quiet = FALSE)

cyclo1_1000000 = cyclo1[,c('X1000000', 'IC50_1')]
cyclo1_1000000_concordance = kendall.w (cyclo1_1000000, nrands = 0, type = 1, quiet = FALSE)

print (c(cyclo1_600_concordance))
print (c(cyclo1_10000_concordance))
print (c(cyclo1_50000_concordance))
print (c(cyclo1_100000_concordance))
print (c(cyclo1_200000_concordance))
print (c(cyclo1_500000_concordance))
print (c(cyclo1_1000000_concordance))

#compare individual cyclo day 2 doses to IC50

cyclo2_600 = cyclo2[,c('X600', 'IC50_1')]
cyclo2_600_concordance = kendall.w (cyclo2_600, nrands = 0, type = 1, quiet = FALSE)

cyclo2_10000 = cyclo2[,c('X10000', 'IC50_1')]
cyclo2_10000_concordance = kendall.w (cyclo2_10000, nrands = 0, type = 1, quiet = FALSE)

cyclo2_50000 = cyclo2[,c('X50000', 'IC50_1')]
cyclo2_50000_concordance = kendall.w (cyclo2_50000, nrands = 0, type = 1, quiet = FALSE)

cyclo2_100000 = cyclo2[,c('X100000', 'IC50_1')]
cyclo2_100000_concordance = kendall.w (cyclo2_100000, nrands = 0, type = 1, quiet = FALSE)

cyclo2_200000 = cyclo2[,c('X200000', 'IC50_1')]
cyclo2_200000_concordance = kendall.w (cyclo2_200000, nrands = 0, type = 1, quiet = FALSE)

cyclo2_500000 = cyclo2[,c('X500000', 'IC50_1')]
cyclo2_500000_concordance = kendall.w (cyclo2_500000, nrands = 0, type = 1, quiet = FALSE)

cyclo2_1000000 = cyclo2[,c('X1000000', 'IC50_1')]
cyclo2_1000000_concordance = kendall.w (cyclo2_1000000, nrands = 0, type = 1, quiet = FALSE)

print (c(cyclo2_600_concordance))
print (c(cyclo2_10000_concordance))
print (c(cyclo2_50000_concordance))
print (c(cyclo2_100000_concordance))
print (c(cyclo2_200000_concordance))
print (c(cyclo2_500000_concordance))
print (c(cyclo2_1000000_concordance))


#carboplatin; day 1 exposure and 2 reported IC50s for cisplatin from cancerxgene

carbo1 <- read.csv("carboplatin cell viability day 1.csv", head = TRUE)
head(carbo1)

carbo1_concordance = kendall.w (carbo1, nrands = 0, type = 1, quiet = FALSE)
print (c(carbo1_concordance))

#carboplatin; day 2 exposure and 2 reported IC50s for cisplatin from cancerxgene

carbo2 <- read.csv("carboplatin cell viability day 2.csv", head = TRUE)
head(carbo2)

carbo2_concordance = kendall.w (carbo2, nrands = 0, type = 1, quiet = FALSE)
print (c(carbo2_concordance))

#compare the 2 IC50s that are reported for cisplatin

ic50_cisplatin = carbo1[,c('IC50_1', 'IC50_2')]
ic50_cisplatin_concordance = kendall.w (ic50_cisplatin, nrands = 0, type = 1, quiet = FALSE)
print (c(ic50_cisplatin_concordance))

#compare individual carbo day 1 doses to the 2 cisplatin IC50s

carbo1_1000 = carbo1[,c('X1000', 'IC50_1', 'IC50_2')]
carbo1_1000_concordance = kendall.w (carbo1_1000, nrands = 0, type = 1, quiet = FALSE)

carbo1_5000 = carbo1[,c('X5000', 'IC50_1', 'IC50_2')]
carbo1_5000_concordance = kendall.w (carbo1_5000, nrands = 0, type = 1, quiet = FALSE)

carbo1_10000 = carbo1[,c('X10000', 'IC50_1', 'IC50_2')]
carbo1_10000_concordance = kendall.w (carbo1_10000, nrands = 0, type = 1, quiet = FALSE)

carbo1_20000 = carbo1[,c('X20000', 'IC50_1', 'IC50_2')]
carbo1_20000_concordance = kendall.w (carbo1_20000, nrands = 0, type = 1, quiet = FALSE)

carbo1_80000 = carbo1[,c('X80000', 'IC50_1', 'IC50_2')]
carbo1_80000_concordance = kendall.w (carbo1_80000, nrands = 0, type = 1, quiet = FALSE)

carbo1_100000 = carbo1[,c('X100000', 'IC50_1', 'IC50_2')]
carbo1_100000_concordance = kendall.w (carbo1_100000, nrands = 0, type = 1, quiet = FALSE)

carbo1_500000 = carbo1[,c('X500000', 'IC50_1', 'IC50_2')]
carbo1_500000_concordance = kendall.w (carbo1_500000, nrands = 0, type = 1, quiet = FALSE)

print (c(carbo1_1000_concordance))
print (c(carbo1_5000_concordance))
print (c(carbo1_10000_concordance))
print (c(carbo1_20000_concordance))
print (c(carbo1_80000_concordance))
print (c(carbo1_100000_concordance))
print (c(carbo1_500000_concordance))

#compare individual carbo day 2 doses to the 2 cisplatin IC50s

carbo2_1000 = carbo2[,c('X1000', 'IC50_1', 'IC50_2')]
carbo2_1000_concordance = kendall.w (carbo2_1000, nrands = 0, type = 1, quiet = FALSE)

carbo2_5000 = carbo2[,c('X5000', 'IC50_1', 'IC50_2')]
carbo2_5000_concordance = kendall.w (carbo2_5000, nrands = 0, type = 1, quiet = FALSE)

carbo2_10000 = carbo2[,c('X10000', 'IC50_1', 'IC50_2')]
carbo2_10000_concordance = kendall.w (carbo2_10000, nrands = 0, type = 1, quiet = FALSE)

carbo2_20000 = carbo2[,c('X20000', 'IC50_1', 'IC50_2')]
carbo2_20000_concordance = kendall.w (carbo2_20000, nrands = 0, type = 1, quiet = FALSE)

carbo2_80000 = carbo2[,c('X80000', 'IC50_1', 'IC50_2')]
carbo2_80000_concordance = kendall.w (carbo2_80000, nrands = 0, type = 1, quiet = FALSE)

carbo2_100000 = carbo2[,c('X100000', 'IC50_1', 'IC50_2')]
carbo2_100000_concordance = kendall.w (carbo2_100000, nrands = 0, type = 1, quiet = FALSE)

carbo2_500000 = carbo2[,c('X500000', 'IC50_1', 'IC50_2')]
carbo2_500000_concordance = kendall.w (carbo2_500000, nrands = 0, type = 1, quiet = FALSE)

print (c(carbo2_1000_concordance))
print (c(carbo2_5000_concordance))
print (c(carbo2_10000_concordance))
print (c(carbo2_20000_concordance))
print (c(carbo2_80000_concordance))
print (c(carbo2_100000_concordance))
print (c(carbo2_500000_concordance))


