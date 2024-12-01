rm(list = ls())

datanul <- read.csv(file='C:/Users/damie/desktop/MASTER/ADM/tp/Datagenus.csv', sep=';')
data <- datanul[1:1000,]

espece <- paste0("gen", 1:27)
densite <- data[espece] / data$surface
densite <- as.matrix(densite)

n <- nrow(densite) 
p <- ncol(densite) 

moyennes_especes <- colMeans(densite)
sd_especes <- sqrt(colSums((densite - matrix(moyennes_especes, n, p, byrow = TRUE))^2) / (n))

tableau <- (densite - matrix(moyennes_especes, n, p, byrow = TRUE)) / matrix(sd_especes, n, p, byrow = TRUE)

dp = dist(tableau, method="euclidean")
CAHDP = hclust(d=dp, method = "ward.D")
plot(CAHDP)

PDP2 = cutree(tree = CAHDP, k=4)
rect.hclust(CAHDP, 4, border="blue")

R2_PDP2 = cbind(rep(0 , ncol(tableau)))
for (i in 1:ncol(tableau)) {
  R2_PDP2[i] = summary(lm(tableau[,i]~as.factor(PDP2)))$r.squared
}
row.names(R2_PDP2) = colnames(tableau)                        
R2G_PDP2 = mean(R2_PDP2)

# Partie K-means
IC2DP = data.frame(model.matrix(~as.factor(PDP2)-1))
mIC2DP = as.matrix(IC2DP)
mDP = as.matrix(tableau)
CentresC2 = solve(t(mIC2DP) %*% mIC2DP) %*% t(mIC2DP)%*% mDP 
KMDP2 = kmeans(tableau, CentresC2)
KMDP2$cluster

boxplot(tableau[, 1]~as.factor(KMDP2$cluster), main="Boxplot pour une variable spécifique par Cluster K-means")

R2_KMDP2 = cbind(rep(0, ncol(tableau)))
for (i in 1:ncol(tableau)) {
  R2_KMDP2[i] = summary(lm(tableau[, i] ~ as.factor(KMDP2$cluster)))$r.squared
}
row.names(R2_KMDP2) = colnames(tableau)
R2G_KMDP2 = mean(R2_KMDP2)


print(paste("Le R^2 global après K-means est de :", R2G_KMDP2))



data$cluster <- KMDP2$cluster
forest_by_cluster <- table(data$forest, data$cluster)
print(addmargins(forest_by_cluster))
proportions <- prop.table(forest_by_cluster, 2)  
print(proportions)


#Tshuprow

chisq_test <- chisq.test(table(type_forestier, KMDP2$cluster))
Tschuprow_T <- sqrt(chisq_test$statistic / (n * sqrt((length(unique(type_forestier))-1)*(4-1))))
print(Tschuprow_T)

chisq_test <- chisq.test(table(type_geo, KMDP2$cluster))
Tschuprow_T <- sqrt(chisq_test$statistic / (n * sqrt((length(unique(type_geo))-1)*(4-1))))
print(Tschuprow_T)




