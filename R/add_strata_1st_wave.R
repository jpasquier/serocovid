# Lecture des données de l'échantillon brut (premier échantillon OFS)
# Le fichier se trouve dans
# /UNISANTE_DESS/S_MC_COVID19/UNDERCOVER_CONF/Listes OFS/Originaux OFS/
smpl <- read.csv("COVID19_VD_V1_Total.csv", sep = ";")
# Conversion de la date de naissance : Du format chaîne de charactères au
# format Date
smpl$dateOfBirth <- as.Date(smpl$dateOfBirth)
# On extrait toutes les combinaisons `Date de naisssance` / `Strate`
tmp <- unique(smpl[!is.na(smpl$strate), c("dateOfBirth", "strate")])
# On complète la liste avec les dates de naissance intermédiaire pour avoir une
# liste plus complète (au cas où certaines dates de naissance auraient été
# modifiées par la suite)
# La variable `strate` est renommée `stratum`
tmp <- do.call(rbind, lapply(sort(unique(tmp$strate)), function(i) {
  z <- tmp[tmp$strate == i, "dateOfBirth"]
  data.frame(dateOfBirth = seq(min(z), max(z), by = 1), stratum = i)
}))
# On ajoute les strates aux données
#   * data : nom du data.frame contenant les données -> A ADAPTER
#   * DOB : nom de la variable `Date de naissance` dans `data` -> A ADAPTER
data <- merge(data, tmp, by.x = "DOB", by.y = "dateOfBirth",
              all.x = TRUE, sort = FALSE)
# On vérifie qu'une strate ait été attribuée à chaque participant. Si ce n'est
# pas le cas, il faut chercher d'où vient le problème.
if (any(is.na(data$stratum))) {
  warning("missing strata")
  print(data[is.na(data$stratum), ])
}
# On fait le ménage
rm(smpl, tmp)
