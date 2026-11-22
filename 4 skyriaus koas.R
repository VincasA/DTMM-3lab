# =====================================================================
# 3 užduotis. 4 skyrius: klasterizavimas originalioje erdvėje
# Algoritmai: k-means, hierarchinis, DBSCAN
# Duomenys: A25_normuotas.csv (RobustScaler, 1050 objektų, 25 požymiai)
# Reikia paketų: cluster, dbscan   ->  install.packages(c("cluster", "dbscan"))
# =====================================================================
library(cluster)
library(dbscan)

# ---------------- Parametrai ----------------
POZ_AIBE    <- "A"          # "A" (25 požymiai) arba "B" (16 požymių)
K_MAIN      <- 5            # pagrindinis klasterių skaičius (3 skyrius)
K_ALT       <- 3            # papildomas k (atitinka 3 žinomas klases)
K_LIST      <- c(K_ALT, K_MAIN)
NSTART      <- 25           # k-means pradžių skaičius
SEED        <- 0            # užfiksuota sėkla
N_SEEDS     <- 10           # sėklų skaičius k-means jautrumo tyrimui (nstart = 1)
MINPTS_GRID <- c(5, 10, 25) # DBSCAN minPts reikšmės
MINPTS_MAIN <- 10           # DBSCAN minPts galutiniam variantui
EPS_GRID    <- seq(1, 3.5, by = 0.25)  # DBSCAN eps tinklelis
EPS_FINAL   <- NULL         # NULL = eps pagal kelio (knee) metodą; arba įrašykite skaičių

isv_dir <- paste0("klasterizavimo_rezultatai_", POZ_AIBE)
dir.create(isv_dir, showWarnings = FALSE)
isv <- function(f) file.path(isv_dir, f)

# ---------------- Duomenys ----------------
n_df     <- read.csv("A25_normuotas.csv")
pozymiai <- setdiff(names(n_df), c("class", "id"))
X_A      <- as.matrix(n_df[, pozymiai])
izmesti_B <- c("LogOfAreas", "SigmoidOfAreas", "Sum_of_Luminosity",
               "X_Perimeter", "Y_Perimeter", "Log_Y_Index", "Log_X_Index",
               "Y_Maximum", "X_Maximum")
X <- if (POZ_AIBE == "A") X_A else X_A[, !colnames(X_A) %in% izmesti_B]
stopifnot(nrow(X) == 1050, !anyNA(X))
klase <- factor(n_df$class)                 # tik interpretacijai
ids   <- n_df$id
cat("Požymių aibė:", POZ_AIBE, "-", ncol(X), "požymių\n")

D_euk <- dist(X)                            # euklidinis atstumas

# PCA 2D tik klasterių VIZUALIZAVIMUI (klasterizavimas vyksta originalioje erdvėje)
pca <- prcomp(X, center = TRUE, scale. = FALSE)
pc  <- pca$x[, 1:2]
pc_proc <- round(100 * pca$sdev[1:2]^2 / sum(pca$sdev^2), 1)

# ---------------- Pagalbinės funkcijos ----------------
# Pakoreguotas Rand indeksas (ARI)
ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  sc <- function(x) sum(choose(x, 2))
  sij <- sc(tab); sa <- sc(rowSums(tab)); sb <- sc(colSums(tab))
  ex <- sa * sb / choose(n, 2); mx <- (sa + sb) / 2
  if (mx == ex) return(1)
  (sij - ex) / (mx - ex)
}

# Vidutinis siluetas (triukšmo taškai, žymimi 0, neįtraukiami)
sil_mean <- function(lab, D) {
  nz <- lab > 0
  if (length(unique(lab[nz])) < 2) return(NA_real_)
  if (all(nz)) return(mean(silhouette(lab, D)[, "sil_width"]))
  Dm <- as.matrix(D)[nz, nz]
  mean(silhouette(lab[nz], as.dist(Dm))[, "sil_width"])
}

# WCSS euklidinėje erdvėje (palyginamumui tarp visų metodų)
wcss_of <- function(Xm, lab) {
  sum(sapply(split(seq_len(nrow(Xm)), lab), function(i)
    sum(scale(Xm[i, , drop = FALSE], scale = FALSE)^2)))
}

# Viena rezultatų lentelės eilutė
eilute <- function(metodas, param, lab, D, bazine = NULL) {
  nz <- lab > 0
  sz <- as.integer(table(lab[nz]))
  data.frame(metodas = metodas, parametrai = param,
             n_klasteriu = length(sz), triuksas = sum(!nz),
             WCSS = round(wcss_of(X[nz, , drop = FALSE], lab[nz]), 1),
             silhouette = round(sil_mean(lab, D), 3),
             min_dydis = min(sz), max_dydis = max(sz),
             ARI_vs_baze = if (is.null(bazine)) NA else round(ari(lab, bazine), 3),
             stringsAsFactors = FALSE)
}

# Spalvos ir PCA grafikas (0 = triukšmas)
pal <- c("steelblue", "firebrick", "goldenrod", "darkgreen", "purple",
         "orange3", "cyan4", "brown", "deeppink", "olivedrab")
spalvos <- function(lab) ifelse(lab == 0, "gray60", pal[(pmax(lab, 1) - 1) %% length(pal) + 1])
braizyti_pca <- function(lab, pav) {
  plot(pc[, 1], pc[, 2], col = adjustcolor(spalvos(lab), 0.7),
       pch = ifelse(lab == 0, 4, 16), cex = 0.6, asp = 1,
       xlab = paste0("PC1 (", pc_proc[1], " %)"),
       ylab = paste0("PC2 (", pc_proc[2], " %)"), main = pav)
  ks <- sort(unique(lab[lab > 0]))
  legend("topright", bty = "n", cex = 0.75,
         legend = c(paste("Klasteris", ks), if (any(lab == 0)) "Triukšmas"),
         col = c(spalvos(ks), if (any(lab == 0)) "gray60"),
         pch = c(rep(16, length(ks)), if (any(lab == 0)) 4))
}
png_pca <- function(failas, lab_list, pavadinimai) {
  png(isv(failas), width = 1500 * length(lab_list), height = 1400, res = 200)
  par(mfrow = c(1, length(lab_list)), mar = c(4, 4, 3, 1))
  for (i in seq_along(lab_list)) braizyti_pca(lab_list[[i]], pavadinimai[i])
  dev.off()
}
# Kelio (knee) taškas didėjančiai išgaubtai kreivei
knee_eps <- function(y) {
  y <- sort(y); x <- seq_along(y)
  xn <- (x - min(x)) / (max(x) - min(x)); yn <- (y - min(y)) / (max(y) - min(y))
  y[which.max(xn - yn)]
}

# =====================================================================
# 4.1. K-MEANS
# =====================================================================
cat("\n===== K-MEANS =====\n")
km_baze <- list(); km_eil <- list()
for (k in K_LIST) {
  set.seed(SEED)
  km <- kmeans(X, centers = k, nstart = NSTART, iter.max = 300,
               algorithm = "Hartigan-Wong")
  km_baze[[as.character(k)]] <- km$cluster
}

for (k in K_LIST) {
  baze <- km_baze[[as.character(k)]]
  # a) nstart poveikis (algoritmas visur Hartigan-Wong)
  for (ns in c(1, 25, 100)) {
    set.seed(SEED)
    l <- kmeans(X, centers = k, nstart = ns, iter.max = 300,
                algorithm = "Hartigan-Wong")$cluster
    km_eil[[length(km_eil) + 1]] <- eilute("k-means",
                                           sprintf("k=%d, euklidinis, Hartigan-Wong, nstart=%d", k, ns), l, D_euk, baze)
  }
  # b) panašumo matas: k-medoids (PAM) su skirtingais atstumais
  for (m in c("euclidean", "manhattan")) {
    Dm <- if (m == "euclidean") D_euk else dist(X, method = "manhattan")
    l <- pam(Dm, k = k, diss = TRUE)$clustering
    km_eil[[length(km_eil) + 1]] <- eilute("k-medoids (PAM)",
                                           sprintf("k=%d, atstumas=%s", k, m), l, Dm, baze)
  }
}
km_tab <- do.call(rbind, km_eil)
write.csv(km_tab, isv("4_1_kmeans_parametrai.csv"), row.names = FALSE)
print(km_tab, row.names = FALSE)

# d) sėklų poveikis (nstart = 1, N_SEEDS sėklų)
cat("\nK-means jautrumas sėklai (nstart = 1):\n")
seed_eil <- list()
for (k in K_LIST) {
  baze <- km_baze[[as.character(k)]]
  ls <- lapply(0:(N_SEEDS - 1), function(s) { set.seed(s); kmeans(X, k, nstart = 1, iter.max = 300, algorithm = "Hartigan-Wong")$cluster })
  aris <- sapply(ls, function(l) ari(l, baze))
  wc   <- sapply(ls, function(l) wcss_of(X, l))
  seed_eil[[length(seed_eil) + 1]] <- data.frame(k = k,
                                                 ARI_vid = round(mean(aris), 3), ARI_sd = round(sd(aris), 3), ARI_min = round(min(aris), 3),
                                                 WCSS_min = round(min(wc), 1), WCSS_max = round(max(wc), 1))
}
seed_tab <- do.call(rbind, seed_eil)
write.csv(seed_tab, isv("4_1_kmeans_seklos.csv"), row.names = FALSE)
print(seed_tab, row.names = FALSE)

# Klasterių dydžiai ir sankirta su tikromis klasėmis (tik interpretacijai)
for (k in K_LIST) {
  cat("\nk-means, k =", k, "- klasterių sankirta su klasėmis:\n")
  print(table(klasteris = km_baze[[as.character(k)]], klase = klase))
}
png_pca("4_1_kmeans_pca.png", lapply(K_LIST, function(k) km_baze[[as.character(k)]]),
        paste0("k-means, k = ", K_LIST))

# =====================================================================
# 4.2. HIERARCHINIS KLASTERIZAVIMAS
# =====================================================================
cat("\n===== HIERARCHINIS =====\n")
metrikos <- c("euclidean", "manhattan", "maximum")
D_met <- lapply(metrikos, function(m) dist(X, method = m)); names(D_met) <- metrikos
junginiai <- list(euclidean = c("ward.D2", "complete", "average", "single"),
                  manhattan = c("complete", "average", "single"),
                  maximum   = c("complete", "average", "single"))

h_baze <- hclust(D_met$euclidean, method = "ward.D2")      # atskaitos variantas
hc_baze <- lapply(K_LIST, function(k) cutree(h_baze, k)); names(hc_baze) <- K_LIST

hc_eil <- list()
for (m in metrikos) for (lk in junginiai[[m]]) {
  h <- hclust(D_met[[m]], method = lk)
  kof <- cor(D_met[[m]], cophenetic(h))
  for (k in K_LIST) {
    r <- eilute("hierarchinis", sprintf("k=%d, atstumas=%s, jungimas=%s", k, m, lk),
                cutree(h, k), D_met[[m]], hc_baze[[as.character(k)]])
    r$kofenetine_r <- round(kof, 3)
    hc_eil[[length(hc_eil) + 1]] <- r
  }
}
hc_tab <- do.call(rbind, hc_eil)
write.csv(hc_tab, isv("4_2_hierarchinis_parametrai.csv"), row.names = FALSE)
print(hc_tab, row.names = FALSE)

# (Dendrograma klasterių skaičiui nustatyti braižoma 3 skyriaus skripte.)
png_pca("4_2_hierarchinis_pca.png", unname(hc_baze), paste0("Hierarchinis (Ward), k = ", K_LIST))

for (k in K_LIST) {
  cat("\nHierarchinis (Ward), k =", k, "- sankirta su klasėmis:\n")
  print(table(klasteris = hc_baze[[as.character(k)]], klase = klase))
}

# =====================================================================
# 4.3. DBSCAN
# =====================================================================
cat("\n===== DBSCAN =====\n")
# a) k-NN atstumų grafikai ir eps pagal kelio (knee) metodą
eps_knee <- list()
png(isv("4_3_dbscan_kdist.png"), width = 600 * length(MINPTS_GRID), height = 800, res = 200)
par(mfrow = c(1, length(MINPTS_GRID)), mar = c(4, 4, 3, 1))
for (mp in MINPTS_GRID) {
  kd <- sort(kNNdist(X, k = mp - 1))
  ek <- knee_eps(kd); eps_knee[[as.character(mp)]] <- ek
  plot(kd, type = "l", xlab = "Taškai (surūšiuoti)", ylab = paste0(mp - 1, "-NN atstumas"),
       main = paste0("minPts = ", mp))
  abline(h = ek, lty = 2, col = "firebrick")
  legend("topleft", bty = "n", cex = 0.8, legend = paste("eps (knee) =", round(ek, 2)),
         lty = 2, col = "firebrick")
}
dev.off()
cat("eps pagal knee metodą:\n"); print(round(unlist(eps_knee), 3))

# b) parametrų tinklelis (minPts x eps)
db_eil <- list()
for (mp in MINPTS_GRID) for (eps in EPS_GRID) {
  l <- dbscan(X, eps = eps, minPts = mp)$cluster
  nz <- l > 0
  db_eil[[length(db_eil) + 1]] <- data.frame(
    minPts = mp, eps = eps, n_klasteriu = length(unique(l[nz])),
    triuksas = sum(!nz), triuksas_proc = round(100 * mean(!nz), 1),
    didziausias_klasteris = if (any(nz)) max(table(l[nz])) else 0,
    silhouette = round(sil_mean(l, D_euk), 3))
}
db_tab <- do.call(rbind, db_eil)
write.csv(db_tab, isv("4_3_dbscan_parametrai.csv"), row.names = FALSE)
print(db_tab, row.names = FALSE)

png(isv("4_3_dbscan_parametrai.png"), width = 2400, height = 1000, res = 200)
par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))
sp <- c("steelblue", "firebrick", "darkgreen")
plot(NA, xlim = range(EPS_GRID), ylim = range(db_tab$n_klasteriu),
     xlab = "eps", ylab = "Klasterių skaičius", main = "Klasterių skaičius pagal eps")
for (i in seq_along(MINPTS_GRID)) {
  d <- db_tab[db_tab$minPts == MINPTS_GRID[i], ]
  lines(d$eps, d$n_klasteriu, type = "b", pch = 16, col = sp[i])
}
legend("topright", bty = "n", legend = paste("minPts =", MINPTS_GRID), col = sp, pch = 16, lty = 1)
plot(NA, xlim = range(EPS_GRID), ylim = c(0, 100),
     xlab = "eps", ylab = "Triukšmo taškų dalis (%)", main = "Triukšmo dalis pagal eps")
for (i in seq_along(MINPTS_GRID)) {
  d <- db_tab[db_tab$minPts == MINPTS_GRID[i], ]
  lines(d$eps, d$triuksas_proc, type = "b", pch = 16, col = sp[i])
}
dev.off()

# c) galutinis variantas: eps pagal knee (arba rankiniu būdu)
eps_use <- if (is.null(EPS_FINAL)) round(eps_knee[[as.character(MINPTS_MAIN)]], 2) else EPS_FINAL
cat(sprintf("\nGalutinis DBSCAN: minPts = %d, eps = %.2f\n", MINPTS_MAIN, eps_use))
db <- dbscan(X, eps = eps_use, minPts = MINPTS_MAIN)
lab_db <- db$cluster                                      # 0 = triukšmas
cat("Klasterių dydžiai (0 = triukšmas):\n"); print(table(lab_db))
cat("\nDBSCAN - sankirta su klasėmis:\n")
print(table(klasteris = lab_db, klase = klase))
png_pca("4_3_dbscan_pca.png", list(lab_db),
        sprintf("DBSCAN, minPts = %d, eps = %.2f", MINPTS_MAIN, eps_use))

# =====================================================================
# 4.4. ALGORITMŲ PALYGINIMAS (originali erdvė) - atskaitos rezultatai
# =====================================================================
cat("\n===== PALYGINIMAS =====\n")
lab_km <- km_baze[[as.character(K_MAIN)]]
lab_hc <- hc_baze[[as.character(K_MAIN)]]
palyg <- rbind(
  eilute("k-means",       sprintf("k=%d, euklidinis, nstart=%d", K_MAIN, NSTART), lab_km, D_euk),
  eilute("hierarchinis",  sprintf("k=%d, Ward, euklidinis", K_MAIN), lab_hc, D_euk),
  eilute("DBSCAN",        sprintf("minPts=%d, eps=%.2f", MINPTS_MAIN, eps_use), lab_db, D_euk))
write.csv(palyg, isv("4_4_algoritmu_palyginimas.csv"), row.names = FALSE)
print(palyg, row.names = FALSE)

labs <- list("k-means" = lab_km, "hierarchinis" = lab_hc, "DBSCAN" = lab_db)
ari_m <- outer(seq_along(labs), seq_along(labs),
               Vectorize(function(i, j) round(ari(labs[[i]], labs[[j]]), 3)))
dimnames(ari_m) <- list(names(labs), names(labs))
cat("\nARI tarp algoritmų (DBSCAN triukšmas laikomas atskira grupe):\n"); print(ari_m)
write.csv(ari_m, isv("4_4_ARI_tarp_algoritmu.csv"))

# Sutampančių taškų skaičius tarp k-means ir hierarchinio (po optimalaus klasterių suderinimo)
lent_kh <- table(kmeans = lab_km, hierarchinis = lab_hc)
cat("\nk-means vs hierarchinis (sankirtos lentelė):\n"); print(lent_kh)
write.csv(as.data.frame.matrix(lent_kh), isv("4_4_kmeans_vs_hierarchinis.csv"))

png_pca("4_4_palyginimas_pca.png", unname(labs),
        c(sprintf("k-means, k = %d", K_MAIN), sprintf("Hierarchinis (Ward), k = %d", K_MAIN),
          sprintf("DBSCAN (eps = %.2f, minPts = %d)", eps_use, MINPTS_MAIN)))

# Atskaitos suskirstymai tolesniems skyriams (dimensijos mažinimas, stabilumas, interpretacija)
write.csv(data.frame(id = ids, class = klase, kmeans = lab_km, hierarchinis = lab_hc, dbscan = lab_db),
          isv("4_atskaitos_suskirstymai.csv"), row.names = FALSE)

cat("\nBaigta. Rezultatai aplanke:", isv_dir, "\n")