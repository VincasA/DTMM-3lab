# DTMM 3-ias laboratorinis darbas

A25 duomenų dimensijos mažinimas, klasterizavimas ir klasterių kokybės bei stabilumo vertinimas.

## Reikalavimai

- Python 3 su `numpy` ir `pillow` bibliotekomis.
- R su terminale pasiekiama `Rscript` komanda. Su R paketais `cluster` ir `dbscan`:

```r
install.packages(c("cluster", "dbscan"))
```

- `A25_clean.csv` ir `A25_normuotas.csv` failai tame pačiame aplanke kaip skriptai. Juos reikia pridėti atskirai.

Atidarykite terminalą repozitorijos aplanke. Jei Python aplinka dar nesukurta, paleiskite šias komandas:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install numpy pillow
```

Windows PowerShell aplinkoje Python aplinka sukurta ir aktyvuota šiuo būdu:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install numpy pillow
```

## 3 skyrius. Optimalus klasterių skaičius

Skriptas `3 skyriaus kodas.R` nustato optimalų klasterių skaičių iki 15 dviejose požymių aibėse:

- A — visi 25 požymiai;
- B — 16 požymių po stipriai koreliuotų požymių pašalinimo.

Klasė `class` klasterizavimui nenaudojama, ji tik palyginama su gautais klasteriais interpretacijai. Skriptas taiko:

- empyrinę taisyklę $k \approx \sqrt{n/2}$;
- Elbow metodą;
- vidutinio silueto metodą;
- hierarchinę Ward dendrogramą ir sujungimo aukščius;
- silueto diagramas pasirinktoms klasterių reikšmėms.

Paleidimas:

```bash
Rscript "3 skyriaus kodas.R"
```

Rezultatai išsaugomi aplanke `klasteriu_sk_rezultatai`:

- `A_metodu_lentele.csv` ir `B_metodu_lentele.csv` — Elbow ir silueto metrikos pagal $k$;
- `optimalus_k_santrauka.csv` — metodais gautos $k$ reikšmės;
- `A_dendrogramos_aukstis.csv` ir `B_dendrogramos_aukstis.csv` — dendrogramos sujungimo aukščiai;
- `k1_A_metodai.png`, `k2_B_metodai.png` — Elbow ir silueto grafikai;
- `k3_A_siluetos.png`, `k3_B_siluetos.png` — silueto diagramos;
- `k4_A_dendrograma.png`, `k4_B_dendrograma.png` — dendrogramos;
- silueto lentelės ir k-means bei tikrų klasių sankirtos išvestys.

## 4 skyrius. Klasterizavimas originalioje erdvėje

Skriptas `4 skyriaus koas.R` klasterizuoja originalių požymių erdvėje ir palygina k-means, hierarchinį bei DBSCAN metodus.

Aibė pasirenkama parametro eilutėje:

```r
POZ_AIBE <- "A"  # "A" = 25 požymiai; "B" = 16 požymių
```

Pagrindinės nustatomos reikšmės:

- `K_MAIN <- 5` — pagrindinis klasterių skaičius pagal 3 skyrių;
- `K_ALT <- 3` — papildomas $k$, atitinkantis tris žinomas klases;
- k-means tikrinami su 1, 25 ir 100 pradžių bei su 10 skirtingų sėklų;
- hierarchiniam klasterizavimui tikrinami Euklido, Manheteno ir Maksimumo atstumai bei keli jungimo būdai;
- DBSCAN tinklelyje tikrinami `minPts = 5, 10, 25` ir `eps` nuo 1 iki 3.5 žingsniu 0.25;
- DBSCAN galutiniam variantui pagal nutylėjimą `eps` apskaičiuojamas pagal kelio metodą prie `minPts = 10`.

Paleidimas:

```bash
Rscript "4 skyriaus koas.R"
```

Rezultatai išsaugomi aplanke `klasterizavimo_rezultatai_A` arba `klasterizavimo_rezultatai_B`, priklausomai nuo `POZ_AIBE`.

Pagrindinės išvestys:

- `4_1_kmeans_parametrai.csv` ir `4_1_kmeans_seklos.csv`;
- `4_1_kmeans_pca.png`;
- `4_2_hierarchinis_parametrai.csv` ir `4_2_hierarchinis_pca.png`;
- `4_3_dbscan_parametrai.csv`, `4_3_dbscan_kdist.png`, `4_3_dbscan_parametrai.png` ir `4_3_dbscan_pca.png`;
- `4_4_algoritmu_palyginimas.csv` — bendras k-means, hierarchinio metodo ir DBSCAN palyginimas;
- `4_4_ARI_tarp_algoritmu.csv` — ARI tarp algoritmų;
- `4_4_kmeans_vs_hierarchinis.csv` — k-means ir hierarchinio metodų sankirta;
- `4_4_palyginimas_pca.png` — trijų metodų PCA vizualizacija;
- `4_atskaitos_suskirstymai.csv` — priskyrimai tolimesniems skyriaus veiksmams.

PCA šiame skripte naudojama tik vizualizacijai. Klasterizavimas visais metodais vykdomas originalioje 25 arba 16 požymių erdvėje.

## 5 skyrius. Dimensijos mažinimas

Apskaičiuoja PCA bei t-SNE koordinates ir sukuria projekcijų grafiką. Rezultatai išsaugomi aplanke `rezultatai_5`.

```bash
python A25_5_dimensijos_mazinimas.py \
  A25_normuotas.csv rezultatai_5
```

## 6 skyrius. Klasterizavimas

Paleiskite po 5 skyriaus. Skriptas taiko „k-means“, Ward ir DBSCAN metodus pradinei bei sumažintoms erdvėms. Klasterių priskyrimai, suvestinės ir grafikas išsaugomi aplanke `rezultatai_6`.

```bash
Rscript A25_6_klasterizavimas.R \
  A25_normuotas.csv \
  rezultatai_5/A25_5_koordinates.csv \
  rezultatai_6
```


## 7 skyrius. Vertinimas

Paleiskite po 5 ir 6 skyrių. Skriptas apskaičiuoja klasterių kokybės matus ir įvertina stabilumą dešimtyje kartotinių imčių. Rezultatų lentelės išsaugomos aplanke `rezultatai_7`.

```bash
python A25_7_vertinimas.py \
  A25_clean.csv \
  A25_normuotas.csv \
  rezultatai_5/A25_5_koordinates.csv \
  rezultatai_6/A25_6_priskyrimai.csv \
  A25_5_dimensijos_mazinimas.py \
  A25_6_klasterizavimas.R \
  rezultatai_7
```

### 7.2 poskyris. ARI grafikas

Paleiskite po 7 skyriaus. Skriptas braižo vidutinio ARI priklausomybę nuo imčių skaičiaus. Pradinės ir PCA erdvių rezultatus perskaičiuoja iki 20 imčių, o t-SNE rezultatus skaito iš `rezultatai_7/A25_7_imciu_poros.csv`.

Reikia tos pačios Python aplinkos kaip 7 skyriuje ir `Rscript` su „Cairo“ palaikymu. Failai `A25_7_vertinimas.py` ir `A25_6_klasterizavimas.R` turi būti tame pačiame aplanke kaip grafiko skriptas.

```bash
python A25_7_ARI_grafikas.py \
  A25_clean.csv \
  rezultatai_7/A25_7_imciu_poros.csv \
  rezultatai_7/A25_7_ARI_pakankamumas.png
```

Grafikas išsaugomas faile `rezultatai_7/A25_7_ARI_pakankamumas.png`.

## 8 skyrius. Klasterių profiliai

Reikia Python 3. Papildomų bibliotekų diegti nereikia.

Paleiskite komandą iš repozitorijos aplanko po 6 skyriaus. Skriptui reikia `A25_clean.csv` ir 6 skyriuje sugeneruoto `rezultatai_6/A25_6_priskyrimai.csv`.

```bash
python3 A25_8_profiliai.py \
  A25_clean.csv \
  rezultatai_6/A25_6_priskyrimai.csv \
  rezultatai_8/A25_8_profiliai.csv
```

Skriptas apskaičiuoja PCA ir „k-means“ klasterių dydžius, defektų klasių pasiskirstymą bei požymių medianas ir kvartilius. Rezultatus išsaugo faile `rezultatai_8/A25_8_profiliai.csv`.
