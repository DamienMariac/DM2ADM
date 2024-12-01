rm(list = ls())
library(factoextra)

# On charge le dataframe.
tab <- read.csv("./Datagenus.csv", sep = ";")
# On enleve la ligne qui pose probleme.
data <- tab[1 :1000,]

# Sélectionner les colonnes des especes
espece <- paste0("gen", 1:27)   # Calcule de la densité de peuplement pour chaque espece
densite <- data[espece] / data$surface
densite <- as.matrix(densite)   # Convertir la dataframe densité en matrice densité

#Partie 1 preparation des donnees

# Statistiques descriptives pour chaque espece
summary(densite)
# On s'apercoit que les moyennes, le min et le max differe entre chaque espece
# donc la distance euclidienne entre chaque variable peut etre consequente. 

# Calcul direct de la matrice des distances
# Initialiser une matrice pour stocker les contributions par variable
n = 1000 # Nombre de parcelle
p = 27   # Nombre de caractéristique
contributions <- array(0, dim = c(n, n, p)) 
# Calcul des distances et contributions
for (i in 1:n) {
  for (j in i:n) {
    # Différence au carré pour chaque variable
    diff_carre <- (densite[i, ] - densite[j, ])^2
    # Distance totale (somme des différences au carré)
    dist_carre <- sum(diff_carre)
    # Contribution de chaque variable
    if (dist_carre != 0) {
      contributions[i, j, ] <- round((diff_carre / dist_carre)*100,2)
      
    }
  }
}
print(contributions[1, 2,])
# On standardise !
# Centrage et réduction 
n <- nrow(densite) # Nombre de parcelles
p <- ncol(densite) # Nombre d'espèces
# Calcul des moyennes et écarts-types par espèce
moyennes_especes <- colMeans(densite)
sd_especes <- sqrt(colSums((densite - matrix(moyennes_especes, n, p, byrow = TRUE))^2) / (n))
# Centrage et réduction des densités
densite_centree_reduite <- (densite - matrix(moyennes_especes, n, p, byrow = TRUE)) / matrix(sd_especes, n, p, byrow = TRUE)





# Partie 2 CAH des parcelles sur les densites de peuplement.

#-----------------------------------Pour l'indice de Ward----------------------------------
#a) création de la matrice des distances euclidiennes:
dp=dist(densite_centree_reduite, method="euclidean")
#b) CAH avec Ward 
CAHDP = hclust(d=dp, method = "ward.D")
plot(CAHDP , main="",xlab ="Parcelles", ylab="Indice" ,sub="")
# Nombre de classe observe
k=4
rect.hclust(CAHDP, k, border="blue")  # Pour k classes
#Coupure de l'arbre et fabrication de la variable de classe correspondant à la partition obtenue.
PDP2 = cutree(tree = CAHDP, k)
WARD = cutree(tree=CAHDP,k)
#Calcul du R2 des variables avec la variable de classe. On va stocker tous les R2 dans un seul vecteur
R2_PDP2 = cbind(rep(0 , ncol(densite)))
#Puis, on calcule les R2 de toutes les variables avec la variable de classe et on met les résultats dans R2:
for (i in cbind(1:ncol(densite))) {
  R2_PDP2[i] = summary(lm(densite[,i]~as.factor(PDP2)))$r.squared
}
#On peut réassigner les noms des variables aux éléments de ce vecteur:
row.names(R2_PDP2) = colnames(densite)
#f) Calcul du R2 de la partition:
R2G_PV2 = mean(R2_PDP2)

# Calcul de la matrice des distances
dp = dist(densite_centree_reduite, method = "euclidean")

# CAH avec la méthode de Ward
CAHDP = hclust(d = dp, method = "ward.D")


# Création d'une liste pour stocker les parcelles dans chaque classe
classes_parcelles_ward <- list()
# Remplir la liste avec les indices des parcelles dans chaque classe
for (i in 1:k) {
  classes_parcelles_ward[[i]] <- which(PDP2 == i)
}
# Création d'une liste pour stocker les résultats par classe
resultats_forest_par_classe_ward <- list()

# Parcourir chaque classe et compter le nombre de parcelles par type de géologie
for (i in 1:k) {
  # Obtenir les indices des parcelles dans la classe i
  parcelles_classe_i <- classes_parcelles_ward[[i]]
  
  # Extraire les types de géologie correspondants pour ces parcelles
  forest_parcelles_i <- data$forest[parcelles_classe_i]
  
  # Calculer la fréquence des types de géologie dans cette classe
  resultats_forest_par_classe_ward[[i]] <- table(forest_parcelles_i)
}

fviz_dend(CAHDP, 
          k = 4, 
          show_labels = FALSE, 
          rect = TRUE, 
          xlab = "Parcelles",  # Nom de l'axe X
          ylab = "Indice",      # Nom de l'axe Y
          main =""
)

#-------------------------------------Pour l'indice du saut maximal
#CAH avec saut maximal
CAHDP2 = hclust(d=dp, method = "complete")
plot(CAHDP2, main="Dendogramme (indice du saut maximal)",xlab ="Parcelles", ylab="Indice" ,sub="")
rect.hclust(CAHDP2, k, border="blue")  # Pour k classes
#Coupure de l'arbre et fabrication de la variable de classe correspondant à la partition obtenue.
PDP2_2 = cutree(tree = CAHDP2, k)
MAX = cutree(tree = CAHDP2, k)
#Calcul du R2 des variables avec la variable de classe. On va stocker tous les R2 dans un seul vecteur
R2_PDP2_2 = cbind(rep(0 , ncol(densite)))
#Puis, on calcule les R2 de toutes les variables avec la variable de classe et on met les résultats dans R2:
for (i in cbind(1:ncol(densite))) {
  R2_PDP2_2[i] = summary(lm(densite[,i]~as.factor(PDP2_2)))$r.squared
}
#On peut réassigner les noms des variables aux éléments de ce vecteur:
row.names(R2_PDP2_2) = colnames(densite)
#f) Calcul du R2 de la partition:
R2G_PV2_2 = mean(R2_PDP2_2)

# Création d'une liste pour stocker les parcelles dans chaque classe
classes_parcelles_max <- list()
# Remplir la liste avec les indices des parcelles dans chaque classe
for (i in 1:k) {
  classes_parcelles_max[[i]] <- which(PDP2_2 == i)
}
# Création d'une liste pour stocker les résultats par classe
resultats_forest_par_classe_max <- list()

# Parcourir chaque classe et compter le nombre de parcelles par type de géologie
for (i in 1:k) {
  # Obtenir les indices des parcelles dans la classe i
  parcelles_classe_i <- classes_parcelles_max[[i]]
  
  # Extraire les types de géologie correspondants pour ces parcelles
  forest_parcelles_i <- data$forest[parcelles_classe_i]
  
  # Calculer la fréquence des types de géologie dans cette classe
  resultats_forest_par_classe_max[[i]] <- table(forest_parcelles_i)
}
fviz_dend(CAHDP2, 
          k = 4, 
          show_labels = FALSE, 
          rect = TRUE, 
          xlab = "Parcelles",  # Nom de l'axe X
          ylab = "Indice",      # Nom de l'axe Y
          main=""
)
#-------------------------------------Pour l'indice du saut minimal
CAHDP3 = hclust(d=dp, method = "single")
plot(CAHDP3, main="Dendogramme (indice du saut minimal)",xlab ="Parcelles", ylab="Indice" ,sub="")
#Coupure de l'arbre et fabrication de la variable de classe correspondant à la partition obtenue.
PDP2_3= cutree(tree = CAHDP3, k)
MIN= cutree(tree = CAHDP3, k)
#Calcul du R2 des variables avec la variable de classe. On va stocker tous les R2 dans un seul vecteur
R2_PDP2_3 = cbind(rep(0 , ncol(densite)))
#Puis, on calcule les R2 de toutes les variables avec la variable de classe et on met les résultats dans R2:
for (i in cbind(1:ncol(densite))) {
  R2_PDP2_3[i] = summary(lm(densite[,i]~as.factor(PDP2_3)))$r.squared
}
#On peut réassigner les noms des variables aux éléments de ce vecteur:
row.names(R2_PDP2_3) = colnames(densite)
#f) Calcul du R2 de la partition:
R2G_PV2_3 = mean(R2_PDP2_3)
hist(CAHDP3$height)

# Création d'une liste pour stocker les parcelles dans chaque classe
classes_parcelles_min <- list()
# Remplir la liste avec les indices des parcelles dans chaque classe
for (i in 1:k) {
  classes_parcelles_min[[i]] <- which(PDP2_3 == i)
}
# Création d'une liste pour stocker les parcelles dans chaque classe
classes_parcelles_min <- list()
# Remplir la liste avec les indices des parcelles dans chaque classe
for (i in 1:k) {
  classes_parcelles_min[[i]] <- which(PDP2_3 == i)
}
# Création d'une liste pour stocker les résultats par classe
resultats_forest_par_classe_min <- list()

# Parcourir chaque classe et compter le nombre de parcelles par type de géologie
for (i in 1:k) {
  # Obtenir les indices des parcelles dans la classe i
  parcelles_classe_i <- classes_parcelles_min[[i]]
  
  # Extraire les types de géologie correspondants pour ces parcelles
  forest_parcelles_i <- data$forest[parcelles_classe_i]
  
  # Calculer la fréquence des types de géologie dans cette classe
  resultats_forest_par_classe_min[[i]] <- table(forest_parcelles_i)
}

#-------------------------------------Pour l'indice du saut moyen
CAHDP4 = hclust(d=dp, method = "average")
plot(CAHDP4, main="Dendogramme (indice du saut moyen)",xlab ="Parcelles", ylab="Indice" ,sub="")
#Coupure de l'arbre et fabrication de la variable de classe correspondant à la partition obtenue.
PDP2_4= cutree(tree = CAHDP4, k)
MOYEN= cutree(tree = CAHDP4, k)
#Calcul du R2 des variables avec la variable de classe. On va stocker tous les R2 dans un seul vecteur
R2_PDP2_4 = cbind(rep(0 , ncol(densite)))
#Puis, on calcule les R2 de toutes les variables avec la variable de classe et on met les résultats dans R2:
for (i in cbind(1:ncol(densite))) {
  R2_PDP2_4[i] = summary(lm(densite[,i]~as.factor(PDP2_4)))$r.squared
}
#On peut réassigner les noms des variables aux éléments de ce vecteur:
row.names(R2_PDP2_4) = colnames(densite)
#f) Calcul du R2 de la partition:
R2G_PV2_4 = mean(R2_PDP2_4)
hist(CAHDP4$height)
# Création d'une liste pour stocker les parcelles dans chaque classe
classes_parcelles_moy <- list()
# Remplir la liste avec les indices des parcelles dans chaque classe
for (i in 1:k) {
  classes_parcelles_moy[[i]] <- which(PDP2_4 == i)
}
# Création d'une liste pour stocker les parcelles dans chaque classe
classes_parcelles_moy <- list()
# Remplir la liste avec les indices des parcelles dans chaque classe
for (i in 1:k) {
  classes_parcelles_moy[[i]] <- which(PDP2_4 == i)
}
# Création d'une liste pour stocker les résultats par classe
resultats_forest_par_classe_moy <- list()

# Parcourir chaque classe et compter le nombre de parcelles par type de géologie
for (i in 1:k) {
  # Obtenir les indices des parcelles dans la classe i
  parcelles_classe_i <- classes_parcelles_moy[[i]]
  
  # Extraire les types de géologie correspondants pour ces parcelles
  forest_parcelles_i <- data$forest[parcelles_classe_i]
  
  # Calculer la fréquence des types de géologie dans cette classe
  resultats_forest_par_classe_moy[[i]] <- table(forest_parcelles_i)
}
fviz_dend(CAHDP4, 
          k = 4, 
          show_labels = FALSE, 
          rect = TRUE, 
          xlab = "Parcelles",  # Nom de l'axe X
          ylab = "Indice",      # Nom de l'axe Y
          main=""
)
#-------------------INDICE DE RAND---------------------------
# Fonction pour calculer l'indice de Rand
rand_index <- function(partition1, partition2) {
  n <- length(partition1)
  
  # Initialisation des compteurs
  C1 <- 0 
  C2 <- 0 
  D1 <- 0 
  D2 <- 0 
  
  # Boucle sur toutes les paires possibles
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      same_in_P <- (partition1[i] == partition1[j])
      same_in_P_prime <- (partition2[i] == partition2[j])
      
      if (same_in_P && same_in_P_prime) {
        C1 <- C1 + 1
      } else if (!same_in_P && !same_in_P_prime) {
        C2 <- C2 + 1
      } else if (same_in_P && !same_in_P_prime) {
        D1 <- D1 + 1
      } else if (!same_in_P && same_in_P_prime) {
        D2 <- D2 + 1
      }
    }
  }
  
  # Calcul de l'indice de Rand
  rand_index_value <- (C1 + C2) / (C1 + C2 + D1 + D2)
  
  return(rand_index_value)
}

indice_rand_1 <- rand_index(WARD, MAX )
cat(indice_rand_1)
#indice_rand_2 <- rand_index(WARD, MIN)
indice_rand_3 <- rand_index(WARD, MOYEN)
cat(indice_rand_3)
#indice_rand_4 <- rand_index(MAX, MIN)
indice_rand_5 <- rand_index(MAX, MOYEN)
cat(indice_rand_5)
indice_rand_6 <- rand_index(MIN, MOYEN)
cat(indice_rand_6)

# Nous allons maintenant optimiser avec la methode du Kmeans

#- Transformation d'une variable qualitative en matrice d'indicatrices:
IC2DP = data.frame(model.matrix(~as.factor(PDP2)-1))
#- Calcul matriciel des centres de gravité de classes de la CAH:
mIC2DP = as.matrix(IC2DP)
mDP = as.matrix(densite)
CentresC2 = solve(t(mIC2DP) %*% mIC2DP) %*% t(mIC2DP)%*% mDP
#- K-means à partir de ces centres initiaux:
KMDP2 = kmeans(densite_centree_reduite, CentresC2)
#- La variable de classe ainsi produite est dans:
KMDP2$cluster

# Extraire les indices des parcelles appartenant au cluster 
liste_cluster_indices <- list()
for(i in 1:k){
  liste_cluster_indices[[i]] <- which(KMDP2$cluster == i) # indice des parcelles appartenant au i eme cluster
}

#- Boxplot d'une variable xj conditionnellement à la variable de classe:
boxplot(densite_centree_reduite[,8]~as.factor(PDP2), main="",xlab ="classes")
boxplot(densite_centree_reduite[,8]~as.factor(KMDP2$cluster), xlab ="classes")



# Voir l'opti des classes et ce qu'il y a dedans.
# Création d'une liste pour stocker les parcelles dans chaque classe
classes_parcelles_ward_kmeans <- list()
# Remplir la liste avec les indices des parcelles dans chaque classe
for (i in 1:k) {
  classes_parcelles_ward_kmeans[[i]] <- which(KMDP2$cluster == i)
}
# Création d'une liste pour stocker les résultats par classe
resultats_forest_par_classe_ward_kmeans <- list()

# Parcourir chaque classe et compter le nombre de parcelles par type de géologie
for (i in 1:k) {
  # Obtenir les indices des parcelles dans la classe i
  parcelles_classe_i <- classes_parcelles_ward_kmeans[[i]]
  
  # Extraire les types de géologie correspondants pour ces parcelles
  forest_parcelles_i <- data$forest[parcelles_classe_i]
  
  # Calculer la fréquence des types de géologie dans cette classe
  resultats_forest_par_classe_ward_kmeans[[i]] <- table(forest_parcelles_i)
}