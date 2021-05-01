library(seqinr)
library(ggplot2)
library(viridis)
library(Biostrings)
library(ape)
library(gridExtra)
library('DECIPHER')
library(ggtree)


calc_GCpec <- function(dna){
  
  size <- length(dna)
  
  g_count <- seqinr::count(dna,1, alphabet = "g") * 100 / size
  c_count <- seqinr::count(dna,1, alphabet = "c") * 100 / size
  
  return(g_count + c_count)
}


coronavirus <- c("MW133981","MT577009", "MT835383", "MT890462",
                 "MW056032", "MT470219", "MT594401", "MW030193",
                 "MT810758", "MW041156", "MT324062", "MT994849",
                 "MT670013", "MT940481", "MW134558","MT876433")


virus_table <- data.frame(name = 1:16, Longitud = 1:16, Adenina = 1:16,
                          Citosina= 1:16, Guanina = 1:16, Timina= 1:16, ContenidoGC=1:16)

virus <- read.GenBank(coronavirus, as.character = TRUE, species.names = FALSE)


for (i in 1:length(coronavirus)){
  virus_code <- coronavirus[i]
  
  flashed_m <- paste("Cargando Virus", i, "de", length(coronavirus),
                     "Codigo de virus:", virus_code, collapse = " ")
  print(flashed_m)
  
  virus_table[i, 1] <- virus_code
  virus_table[i, 2] <- length(virus[[i]])
  virus_table[i, 3:6] <- seqinr::count(virus[[i]], 1)
  virus_table[i, 7] <- calc_GCpec(virus[[i]])
 
  if (i == length(coronavirus)){
    print("Cargado Finalizado")
  }
}

virus_table


len_plot <- ggplot(virus_table, aes(name, Longitud, fill=name))+
  geom_col()+
  theme_classic() +
  labs(fill="Virus Color Code",x="Virus", y="Longitud", title = "Comparación Longitud de Genomas", subtitle="Grafica que compara las longitudes de distintos virus"
       ,caption = "Datos Obtenidos de NCBI database")

len_plot

new_data <- data.frame(nuc = 1:4, freq = 1:4)
new_data[1:4, 1] <- c("a", "c", "g", "t")

  
# 1
new_data[1:4,2] <- c(virus_table[1,3], virus_table[1,4], virus_table[1,5], virus_table[1,6])
MW133981_plot <- ggplot(new_data, aes(nuc, freq, fill=nuc))+
  geom_col()+
  theme_bw()+
  labs(title = paste("Contenido Nucleotido de:", virus_table[1,1], sep =" "), subtitle="Grafica que muestra el contenido nucleotido"
       ,caption = "Datos Obtenidos de NCBI database", x ="Nucleotidos", y="Cantidad")

# 2
new_data[1:4,2] <- c(virus_table[2,3], virus_table[2,4], virus_table[2,5], virus_table[2,6])
MT577009_plot <- ggplot(new_data, aes(nuc, freq, fill=nuc))+
  geom_col()+
  theme_bw()+
  labs(title = paste("Contenido Nucleotido de:", virus_table[2,1], sep =" "), subtitle="Grafica que muestra el contenido nucleotido"
       ,caption = "Datos Obtenidos de NCBI database", x ="Nucleotidos", y="Cantidad")

# 3
new_data[1:4,2] <- c(virus_table[3,3], virus_table[3,4], virus_table[3,5], virus_table[3,6])
MT835383_plot <- ggplot(new_data, aes(nuc, freq, fill=nuc))+
  geom_col()+
  theme_bw()+
  labs(title = paste("Contenido Nucleotido de:", virus_table[3,1], sep =" "), subtitle="Grafica que muestra el contenido nucleotido"
       ,caption = "Datos Obtenidos de NCBI database", x ="Nucleotidos", y="Cantidad")

# 4
new_data[1:4,2] <- c(virus_table[4,3], virus_table[4,4], virus_table[4,5], virus_table[4,6])
MT890462_plot <- ggplot(new_data, aes(nuc, freq, fill=nuc))+
  geom_col()+
  theme_bw()+
  labs(title = paste("Contenido Nucleotido de:", virus_table[4,1], sep =" "), subtitle="Grafica que muestra el contenido nucleotido"
       ,caption = "Datos Obtenidos de NCBI database", x ="Nucleotidos", y="Cantidad")

# 5
new_data[1:4,2] <- c(virus_table[5,3], virus_table[6,4], virus_table[6,5], virus_table[6,6])
MW056032_plot <- ggplot(new_data, aes(nuc, freq, fill=nuc))+
  geom_col()+
  theme_bw()+
  labs(title = paste("Contenido Nucleotido de:", virus_table[5,1], sep =" "), subtitle="Grafica que muestra el contenido nucleotido"
       ,caption = "Datos Obtenidos de NCBI database", x ="Nucleotidos", y="Cantidad")


# 6
new_data[1:4,2] <- c(virus_table[6,3], virus_table[6,4], virus_table[6,5], virus_table[6,6])
MT470219_plot <- ggplot(new_data, aes(nuc, freq, fill=nuc))+
  geom_col()+
  theme_bw()+
  labs(title = paste("Contenido Nucleotido de:", virus_table[6,1], sep =" "), subtitle="Grafica que muestra el contenido nucleotido"
       ,caption = "Datos Obtenidos de NCBI database", x ="Nucleotidos", y="Cantidad")

# 7
new_data[1:4,2] <- c(virus_table[7,3], virus_table[7,4], virus_table[7,5], virus_table[7,6])
MT594401_plot <- ggplot(new_data, aes(nuc, freq, fill=nuc))+
  geom_col()+
  theme_bw()+
  labs(title = paste("Contenido Nucleotido de:", virus_table[7,1], sep =" "), subtitle="Grafica que muestra el contenido nucleotido"
       ,caption = "Datos Obtenidos de NCBI database", x ="Nucleotidos", y="Cantidad")

# 8
new_data[1:4,2] <- c(virus_table[8,3], virus_table[8,4], virus_table[8,5], virus_table[8,6])
MW030193_plot <- ggplot(new_data, aes(nuc, freq, fill=nuc))+
  geom_col()+
  theme_bw()+
  labs(title = paste("Contenido Nucleotido de:", virus_table[8,1], sep =" "), subtitle="Grafica que muestra el contenido nucleotido"
       ,caption = "Datos Obtenidos de NCBI database", x ="Nucleotidos", y="Cantidad")

# 9
new_data[1:4,2] <- c(virus_table[9,3], virus_table[9,4], virus_table[9,5], virus_table[9,6])
MT810758_plot <- ggplot(new_data, aes(nuc, freq, fill=nuc))+
  geom_col()+
  theme_bw()+
  labs(title = paste("Contenido Nucleotido de:", virus_table[9,1], sep =" "), subtitle="Grafica que muestra el contenido nucleotido"
       ,caption = "Datos Obtenidos de NCBI database", x ="Nucleotidos", y="Cantidad")

# 10
new_data[1:4,2] <- c(virus_table[10,3], virus_table[10,4], virus_table[10,5], virus_table[10,6])
MW041156_plot <- ggplot(new_data, aes(nuc, freq, fill=nuc))+
  geom_col()+
  theme_bw()+
  labs(title = paste("Contenido Nucleotido de:", virus_table[10,1], sep =" "), subtitle="Grafica que muestra el contenido nucleotido"
       ,caption = "Datos Obtenidos de NCBI database", x ="Nucleotidos", y="Cantidad")

# 11
new_data[1:4,2] <- c(virus_table[11,3], virus_table[11,4], virus_table[11,5], virus_table[11,6])
MT324062_plot <- ggplot(new_data, aes(nuc, freq, fill=nuc))+
  geom_col()+
  theme_bw()+
  labs(title = paste("Contenido Nucleotido de:", virus_table[11,1], sep =" "), subtitle="Grafica que muestra el contenido nucleotido"
       ,caption = "Datos Obtenidos de NCBI database", x ="Nucleotidos", y="Cantidad")

# 12
new_data[1:4,2] <- c(virus_table[12,3], virus_table[12,4], virus_table[12,5], virus_table[12,6])
MT994849_plot <- ggplot(new_data, aes(nuc, freq, fill=nuc))+
  geom_col()+
  theme_bw()+
  labs(title = paste("Contenido Nucleotido de:", virus_table[12,1], sep =" "), subtitle="Grafica que muestra el contenido nucleotido"
       ,caption = "Datos Obtenidos de NCBI database", x ="Nucleotidos", y="Cantidad")

# 13
new_data[1:4,2] <- c(virus_table[13,3], virus_table[13,4], virus_table[13,5], virus_table[13,6])
MT670013_plot <- ggplot(new_data, aes(nuc, freq, fill=nuc))+
  geom_col()+
  theme_bw()+
  labs(title = paste("Contenido Nucleotido de:", virus_table[13,1], sep =" "), subtitle="Grafica que muestra el contenido nucleotido"
       ,caption = "Datos Obtenidos de NCBI database", x ="Nucleotidos", y="Cantidad")

# 14
new_data[1:4,2] <- c(virus_table[14,3], virus_table[14,4], virus_table[14,5], virus_table[14,6])
MT940481_plot <- ggplot(new_data, aes(nuc, freq, fill=nuc))+
  geom_col()+
  theme_bw()+
  labs(title = paste("Contenido Nucleotido de:", virus_table[14,1], sep =" "), subtitle="Grafica que muestra el contenido nucleotido"
       ,caption = "Datos Obtenidos de NCBI database", x ="Nucleotidos", y="Cantidad")

# 15
new_data[1:4,2] <- c(virus_table[15,3], virus_table[15,4], virus_table[15,5], virus_table[15,6])
MW134558_plot <- ggplot(new_data, aes(nuc, freq, fill=nuc))+
  geom_col()+
  theme_bw()+
  labs(title = paste("Contenido Nucleotido de:", virus_table[15,1], sep =" "), subtitle="Grafica que muestra el contenido nucleotido"
       ,caption = "Datos Obtenidos de NCBI database", x ="Nucleotidos", y="Cantidad")

# 16
new_data[1:4,2] <- c(virus_table[16,3], virus_table[16,4], virus_table[16,5], virus_table[16,6])
MT876433_plot <- ggplot(new_data, aes(nuc, freq, fill=nuc))+
  geom_col()+
  theme_bw()+
  labs(title = paste("Contenido Nucleotido de:", virus_table[16,1], sep =" "), subtitle="Grafica que muestra el contenido nucleotido"
       ,caption = "Datos Obtenidos de NCBI database", x ="Nucleotidos", y="Cantidad")

# Primeros 3
grid.arrange(MW133981_plot,MT577009_plot,MT835383_plot, ncol=3)
    
# 4-6
grid.arrange(MT890462_plot,MW056032_plot, MT470219_plot, ncol=3)

# 7-9
grid.arrange(MT594401_plot, MW030193_plot, MT810758_plot, ncol=3)

# 10-12
grid.arrange(MW041156_plot, MT324062_plot, MT994849_plot, ncol=3)

# 13-16
grid.arrange(MT670013_plot, MT940481_plot, MW134558_plot,MT876433_plot, ncol=2)
    
virus_table

gc_plot <- ggplot(virus_table, aes(name, ContenidoGC, fill=name)) +
  geom_point()+
  theme_light()+
  ylim(37.6, 38.2)+
  labs(fill="Virus Color Code",x="Virus", y="Porcentaje", title = "Contenido de GC por virus",
       subtitle="Grafica que muestra el contenido porcentual de GC"
       ,caption = "Datos Obtenidos de NCBI database")
gc_plot

#####################
# Parte 2 ###########
#####################

write.dna(virus, "coronavirus_seqs.fasta", format = "fasta")

viruses_not_align <- readDNAStringSet("coronavirus_seqs.fasta", format = "fasta")
viruses_not_align <- OrientNucleotides(viruses_not_align)

viruses_align <- AlignSeqs(viruses_not_align)
BrowseSeqs(viruses_align)

writeXStringSet(viruses_align, file = "viruses_align.fasta")  
virus_align <- read.alignment("viruses_align.fasta", format = "fasta")  

matriz_distancia <- dist.alignment(virus_align, matrix = "similarity")  

matriz_distancia

temp <- as.data.frame(as.matrix(matriz_distancia))  

ade4::table.paint(temp, cleg= 0,row.labels = row.names(temp), col.labels = names(temp),
                  clabel.row=1, clabel.col=1,) + scale_color_viridis()

virus_tree <- nj(matriz_distancia)
virus_tree <- ladderize(virus_tree)
plot(virus_tree)

ggtree(virus_tree, branch.length = 'none', layout = 'circular',) + geom_tiplab()

plot_virus_filogenia <- ggtree(virus_tree) + geom_tiplab() + ggtitle("Phylogenetic analysis of SARS COV genomes")
plot_virus_filogenia








