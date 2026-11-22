# =====================================================================
# 3 užduotis. 3 skyrius: optimalaus klasterių skaičiaus nustatymas
# Metodai: empyrinis (taisyklė k ~ sqrt(n/2)), Elbow, vidutinis siluetas,
#          dendrograma (hierarchinis Ward, sujungimo aukščių tarpai)
# Duomenys: A25_normuotas.csv (jau normuoti RobustScaler, 1050 x 25)
# =====================================================================
library(cluster)

# ---------------- Parametrai ----------------
K_MAX  <- 15      # tiriamas klasterių skaičius iki K_MAX
NSTART <- 25      # k-means pradžių skaičius
SEED   <- 0       # užfiksuota sėkla (atkartojamumas)

dir.create("klasteriu_sk_rezultatai", showWarnings = FALSE)
isv <- function(f) file.path("klasteriu_sk_rezultatai", f)

# ---------------- Duomenys ----------------
n_df     <- read.csv("A25_normuotas.csv")
pozymiai <- setdiff(names(n_df), c("class", "id"))
X_A      <- as.matrix(n_df[, pozymiai])          # A: visi 25 požymiai
stopifnot(ncol(X_A) == 25, nrow(X_A) == 1050, !anyNA(X_A))

# B: be stipriai koreliuotų požymių (|r| > 0,9) - tas pats sąrašas kaip 2 darbe
izmesti_B <- c("LogOfAreas", "SigmoidOfAreas", "Sum_of_Luminosity",
               "X_Perimeter", "Y_Perimeter", "Log_Y_Index", "Log_X_Index",
               "Y_Maximum", "X_Maximum")
X_B <- X_A[, !colnames(X_A) %in% izmesti_B]       # 16 požymių
stopifnot(ncol(X_B) == 16)

# Klasė NENAUDOJAMA klasterizavimui, tik palyginimui (3 žinomos klasės)
klase <- factor(n_df$class)

# ---------------- Pagalbinė funkcija ----------------
# Elbow taškas: didžiausias atstumas nuo tiesės, jungiančios pirmą ir paskutinį tašką
elbow_k <- function(k, w) {
  x <- (k - min(k)) / (max(k) - min(k))
  y <- (w - min(w)) / (max(w) - min(w))
  d <- abs(x + y - 1) / sqrt(2)
  k[which.max(d)]
}

# Visi trys metodai vienai požymių aibei
analizuoti <- function(X, vardas) {
  n <- nrow(X)
  D <- dist(X)                                    # euklidinis atstumas
  ks <- 1:K_MAX
  
  # ---- 1. Empyrinis metodas ----
  k_emp <- round(sqrt(n / 2))
  cat(sprintf("[%s] Empyrinis k ~ sqrt(n/2) = sqrt(%d/2) = %.2f -> %d\n",
              vardas, n, sqrt(n / 2), k_emp))
  
  # ---- 2. Elbow ir 3. siluetas (k-means kiekvienam k) ----
  wcss <- numeric(length(ks)); sil <- rep(NA_real_, length(ks))
  for (k in ks) {
    set.seed(SEED)
    km <- kmeans(X, centers = k, nstart = NSTART, iter.max = 100)
    wcss[k] <- km$tot.withinss
    if (k >= 2) sil[k] <- mean(silhouette(km$cluster, D)[, "sil_width"])
  }
  wcss_krit <- c(NA, 100 * (1 - wcss[-1] / wcss[-length(wcss)]))  # % sumažėjimas
  
  k_elbow <- elbow_k(ks, wcss)
  k_sil   <- which.max(sil)
  
  lentele <- data.frame(k = ks, WCSS = wcss,
                        WCSS_sumazejimas_proc = wcss_krit, Silhouette = sil)
  santrauka <- data.frame(
    pozymiu_aibe = vardas, n_pozymiu = ncol(X),
    Empyrinis = k_emp, Elbow = k_elbow, Silhouette = k_sil,
    max_silhouette = round(max(sil, na.rm = TRUE), 3))
  
  list(lentele = lentele, santrauka = santrauka, D = D, X = X,
       k_emp = k_emp, k_elbow = k_elbow, k_sil = k_sil)
}

# ---------------- Skaičiavimai ----------------
cat("== Klasterių skaičiaus nustatymas ==\n")
rez_A <- analizuoti(X_A, "A (25 pozymiai)")
rez_B <- analizuoti(X_B, "B (16 pozymiu)")

# ---------------- Lentelės ----------------
write.csv(rez_A$lentele, isv("A_metodu_lentele.csv"), row.names = FALSE)
write.csv(rez_B$lentele, isv("B_metodu_lentele.csv"), row.names = FALSE)
santr <- rbind(rez_A$santrauka, rez_B$santrauka)
write.csv(santr, isv("optimalus_k_santrauka.csv"), row.names = FALSE)
cat("\n== Rekomenduojamas k pagal metodus ==\n")
print(santr, row.names = FALSE)

cat("\n== A aibė: metrikos pagal k (suapvalinta) ==\n")
print(round(rez_A$lentele, 3), row.names = FALSE)
cat("\n== B aibė: metrikos pagal k (suapvalinta) ==\n")
print(round(rez_B$lentele, 3), row.names = FALSE)

# ---------------- Grafikai (1 x 3: Elbow, WCSS sumažėjimas, siluetas) ----------------
braizyti_metodus <- function(r, failas) {
  t <- r$lentele
  png(isv(failas), width = 2400, height = 800, res = 200)
  par(mfrow = c(1, 3), mar = c(4, 4.5, 3, 1))
  
  # Elbow
  plot(t$k, t$WCSS, type = "b", pch = 16, col = "steelblue",
       xlab = "Klasterių skaičius k", ylab = "WCSS", main = "Elbow metodas")
  abline(v = r$k_elbow, lty = 2, col = "firebrick")
  legend("topright", legend = paste("alkūnė: k =", r$k_elbow),
         lty = 2, col = "firebrick", bty = "n", cex = 0.9)
  
  # WCSS procentinis sumažėjimas
  barplot(t$WCSS_sumazejimas_proc[-1], names.arg = t$k[-1], col = "steelblue",
          xlab = "Klasterių skaičius k", ylab = "WCSS sumažėjimas (%)",
          main = "WCSS sumažėjimas didinant k")
  
  # Siluetas
  plot(t$k[-1], t$Silhouette[-1], type = "b", pch = 16, col = "darkgreen",
       xlab = "Klasterių skaičius k", ylab = "Vidutinis siluetas",
       main = "Vidutinio silueto metodas")
  abline(v = r$k_sil, lty = 2, col = "firebrick")
  legend("bottomright", legend = paste("max: k =", r$k_sil),
         lty = 2, col = "firebrick", bty = "n", cex = 0.9)
  
  dev.off()
}
braizyti_metodus(rez_A, "k1_A_metodai.png")
braizyti_metodus(rez_B, "k2_B_metodai.png")

# ---------------- Dendrograma (hierarchinis, Ward, euklidinis) ----------------
# Klasterių skaičius pagal tarpus tarp sujungimo aukščių: didelis tarpas = aiški riba
K_DENDR <- 5   # k, kuriam braižoma kirpimo linija ir stačiakampiai
pal <- c("steelblue", "firebrick", "goldenrod", "darkgreen", "purple",
         "orange3", "cyan4", "brown", "deeppink", "olivedrab")

dendrograma <- function(r, vardas) {
  h  <- hclust(r$D, method = "ward.D2")
  hh <- rev(tail(h$height, 15))                                   # k = 2..16
  kritimas <- c(NA, round(100 * (1 - hh[-1] / hh[-length(hh)]), 1)) # aukščio kritimas nuo k-1 iki k
  tarpas   <- c(kritimas[-1], NA)                                  # tarpas nuo k iki k+1 (kirpimas ties k)
  tab <- data.frame(k = 2:16, sujungimo_aukstis = round(hh, 2),
                    kritimas_proc_nuo_k_minus_1 = kritimas,
                    tarpas_proc_iki_k_plius_1 = tarpas)
  write.csv(tab, isv(paste0(vardas, "_dendrogramos_aukstis.csv")), row.names = FALSE)
  cat("\nDendrogramos sujungimo aukščiai, ", vardas, " aibė (k = klasterių skaičius):\n", sep = "")
  print(tab, row.names = FALSE)
  
  png(isv(paste0("k4_", vardas, "_dendrograma.png")), width = 3000, height = 1200, res = 200)
  layout(matrix(1:2, 1), widths = c(2, 1))
  par(mar = c(2, 4.5, 3, 1))
  plot(h, labels = FALSE, hang = -1, xlab = "", sub = "",
       ylab = "Atstumas (Ward, euklidinis)", main = "Dendrograma")
  rect.hclust(h, k = K_DENDR, border = pal[1:K_DENDR])
  abline(h = mean(sort(h$height, decreasing = TRUE)[(K_DENDR - 1):K_DENDR]), lty = 2)
  par(mar = c(4, 4.5, 3, 1))
  barplot(hh, names.arg = 2:16, col = "steelblue", xlab = "Klasterių skaičius k",
          ylab = "Sujungimo aukštis", main = "Paskutiniai 15 sujungimų")
  dev.off()
  tab
}
dend_A <- dendrograma(rez_A, "A")
dend_B <- dendrograma(rez_B, "B")

# ---------------- Silueto diagramos (A ir B aibės) ----------------
# Kiekvienai aibei: k = 3, silueto optimumas ir (B aibei) dendrogramos k = 4.
# B aibės silueto optimumas (k = 13) plokščias, todėl pridedamas ir k = 4 bei k = 5.
sil_diagramos <- function(r, vardas, k_papildomi = NULL) {
  ks_sil <- sort(unique(c(3, r$k_sil, k_papildomi)))
  sil_obj <- lapply(ks_sil, function(k) {
    set.seed(SEED)
    km <- kmeans(r$X, centers = k, nstart = NSTART, iter.max = 100)
    list(k = k, km = km, s = silhouette(km$cluster, r$D))
  })
  
  # Visos diagramos vienoje nuotraukoje (viena šalia kitos)
  png(isv(paste0("k3_", vardas, "_siluetos.png")),
      width = 1800 * length(sil_obj), height = 1400, res = 200)
  par(mfrow = c(1, length(sil_obj)))
  for (o in sil_obj) {
    plot(o$s, col = rainbow(o$k),
         main = paste0("Silueto diagrama, ", vardas, " aibė, k = ", o$k), border = NA)
  }
  dev.off()
  
  # Lentelės pagal klasterius
  for (o in sil_obj) {
    lent <- data.frame(
      klasteris = seq_len(o$k), dydis = as.integer(table(o$km$cluster)),
      vid_siluetas = round(tapply(o$s[, "sil_width"], o$km$cluster, mean), 3),
      neigiamu_siluetu_proc = round(100 * tapply(o$s[, "sil_width"] < 0,
                                                 o$km$cluster, mean), 1))
    cat("\n", vardas, " aibė, k = ", o$k, " - siluetas pagal klasterius:\n", sep = "")
    print(lent, row.names = FALSE)
    write.csv(lent, isv(paste0(vardas, "_siluetas_klasteriai_k", o$k, ".csv")),
              row.names = FALSE)
  }
  invisible(sil_obj)
}
sil_A <- sil_diagramos(rez_A, "A")                       # k = 3 ir 5
sil_B <- sil_diagramos(rez_B, "B", k_papildomi = c(4, 5)) # k = 3, 4, 5 ir silueto max (13)

# ---------------- k = 3 klasteriai vs tikros klasės (tik interpretacijai) ----------------
set.seed(SEED)
km3 <- kmeans(X_A, centers = 3, nstart = NSTART)
cat("\nA aibė, k = 3: klasterių ir tikrų klasių sankirta (tik interpretacijai):\n")
print(table(klasteris = km3$cluster, klase = klase))

set.seed(SEED)
km3B <- kmeans(X_B, centers = 3, nstart = NSTART)
cat("\nB aibė, k = 3: klasterių ir tikrų klasių sankirta (tik interpretacijai):\n")
print(table(klasteris = km3B$cluster, klase = klase))

cat("\nBaigta. Rezultatai aplanke: klasteriu_sk_rezultatai/\n")