#    ____                              __     __ _  ____    _  _   
#   / ___|   ___  _ __  ___    ___  ___\ \   / /(_)|  _ \  | || |  
#   \___ \  / _ \| '__|/ _ \  / __|/ _ \\ \ / / | || | | | | || |_ 
#    ___) ||  __/| |  | (_) || (__| (_) |\ V /  | || |_| | |__   _|
#   |____/  \___||_|   \___/  \___|\___/  \_/   |_||____/     |_|  

# Initialisation du générateur de nombres aléatoires
set.seed(666)

# Les données pour le tirage au sort doivent contenir les IDs est les startes
# Exemple : 5 strates avec 5 IDs chacune
fake_data <- data.frame(
  hid = paste0("id", sprintf("%02d", 1:25)),
  stratum = rep(1:5, 5)
)

# Taille de l'échantillon par strate
smpl_size <- read.table(header = TRUE, text = "
  stratum  n
        1  2
        2  2
        3  2
        4  2
        5  2
")

# Sélection (par strate) des IDs à l'aide de la fonction sample
smpl <- do.call(rbind, lapply(smpl_size$stratum, function(s) {
  smpl <- sample(
    fake_data[fake_data$stratum == s, "hid"],
    smpl_size[smpl_size$stratum == s, "n"],
    replace = FALSE
  )
  data.frame(stratum = s, hid = smpl)
}))


