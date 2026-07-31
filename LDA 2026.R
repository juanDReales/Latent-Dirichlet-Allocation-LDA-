######################## Script LDA######################################### 

# Lista de paquetes necesarios

paquetes <- c(
  "tm", "quanteda", "textstem", "textmineR", "LDAvis", "ggplot2", 
  "SnowballC", "readxl", "openxlsx", "parallel", 
  "dplyr", "beepr", "tidytext", "tidyr", "broom", "gridExtra"
)

# Función para verificar e instalar paquetes
verificar_instalar <- function(paquete) {
  if (!requireNamespace(paquete, quietly = TRUE)) {
    cat(paste0("Instalando el paquete: ", paquete, "\n"))
    install.packages(paquete, dependencies = TRUE)
  } else {
    cat(paste0("El paquete ", paquete, " ya está instalado.\n"))
  }
}

# Aplicar la función a cada paquete
lapply(paquetes, verificar_instalar)

### Cargar librerias
library(tm)
library(quanteda)
library(textstem)
library(textmineR)
library(LDAvis)
library(ggplot2)
library(SnowballC)
library(readxl)
library(openxlsx)
library(parallel)
library(dplyr)
library(beepr)
library(tidytext)
library(tidyr)
library(broom)
library(gridExtra)

######################################################################
# Configurar el Working Directory (La carpeta que contiene los archivos de trabajo)
######################################################################

### Lectura de la base de datos
library(readxl)
bd <- read_excel("Bibliometrix depurada final.xlsx")

#####################################################################
## Preprocesamiento de LDA
StopWord <- read_excel("StopWord.xlsx")

# Configura el número de CPUs a utilizar como el número de núcleos del equipo menos dos
cpus <- parallel::detectCores() - 6
cpus

# Crear la matriz de términos por documento (dtm) 
dtm <- textmineR::CreateDtm(
  doc_vec = bd$AB,  # Vector de texto de los abstracts de los documentos.
  doc_names = bd$TI,   #    de los documentos para identificar cada uno.
  ngram_window = c(1,2),  # Incluye unigramas y bigramas en la matriz.
  stopword_vec = c(stopwords::stopwords("en"),
                   StopWord$StopWord), 
  lower = TRUE,                   # Convierte todos los términos a minúsculas.
  remove_punctuation = TRUE,      # Elimina la puntuación del texto.
  remove_numbers = TRUE,          # Remueve los números.
  stem_lemma_function = function(x) SnowballC::wordStem(x, "porter"), 
  # Aplica stemming para reducir palabras a su raíz.
  cpus = cpus                     # Utiliza el número de CPUs asignado
)

beepr::beep(5)
dim(dtm)

zdtm <- quanteda::as.dfm(dtm)
CONVERT <- quanteda::convert(zdtm, to = "topicmodels")
zdtmt <- removeSparseTerms(CONVERT, sparse = 0.995)
dim(zdtmt)

library(Matrix)
library(slam)
dtm_triplet <- as.simple_triplet_matrix(zdtmt)
zdtmF <- sparseMatrix(i = dtm_triplet$i,
                      j = dtm_triplet$j,
                      x = dtm_triplet$v,
                      dims = dim(dtm_triplet),
                      dimnames = dimnames(dtm_triplet))

class(zdtmF)
dim(zdtmF)

library(text2vec)
tf_mat <- textmineR::TermDocFreq(dtm = zdtmF)
write.xlsx(tf_mat, file = "tf_mat.xlsx")

## Inferencia 
# Determinación del número optimo de tópicos(k)
k_list <- seq(5, 40, by = 1) 

system.time(
  model_list <- TmParallelApply(
    X = k_list, 
    FUN = function(k) {
      m <- FitLdaModel(
        dtm = zdtmF, 
        k = k, 
        iterations = 500, 
        burnin = 50,
        alpha = 0.1,
        beta = colSums(zdtmF) / sum(zdtmF) * 100,
        optimize_alpha = TRUE,
        calc_likelihood = TRUE,
        calc_coherence = TRUE,
        calc_r2 = FALSE,
        cpus = cpus
      )
      m$k <- k
      m
    }, 
    export = ls(), 
    cpus = cpus    
  )
)
beepr::beep(3)

coherence_mat <- data.frame(
  k = sapply(model_list, function(x) nrow(x$phi)),             
  coherence = sapply(model_list, function(x) mean(x$coherence)), 
  stringsAsFactors = FALSE
)

write.xlsx(coherence_mat,file = "coherence_mat.xlsx")

k_max <- coherence_mat$k[which.max(coherence_mat$coherence)]
coherence_max <- max(coherence_mat$coherence)

coherence_plot <- ggplot(coherence_mat, aes(x = k, y = coherence)) +
  geom_line() +
  geom_point() +
  geom_point(aes(x = k_max, y = coherence_max), color = "red", size = 2)+  
  theme_bw()+
  theme(
    text = element_text(size = 14),               
    axis.title = element_text(size = 14),          
    axis.text = element_text(size = 12)           
  )   

print(coherence_plot)

ggsave(
  filename = "coherence_plot.png",
  plot = coherence_plot,
  dpi = 500,
  width = 8,   
  height = 6   
)

library(ggplot2)

coherence_plot2 <- ggplot(coherence_mat, aes(x = k, y = coherence)) +
  geom_line(color = "#4C9BE8", linewidth = 0.9) +
  geom_point(color = "#4C9BE8", size = 3, fill = "#4C9BE8", shape = 21, stroke = 0.5) +
  geom_point(
    aes(x = k_max, y = coherence_max, color = paste0("k=", k_max)),
    size = 5, shape = 16
  ) +
  scale_color_manual(values = setNames("red", paste0("k=", k_max))) +
  scale_x_continuous(breaks = seq(5, 40, by = 5)) +
  labs(
    title = "Topic coherence by k",
    x     = "Number of topics (k)",
    y     = "Coherence",
    color = NULL
  ) +
  theme(
    # Background
    plot.background   = element_rect(fill = "#f5f5f5", color = NA),
    panel.background  = element_rect(fill = "#f5f5f5", color = NA),
    # Grid
    panel.grid.major  = element_line(color = "white", linetype = "dashed", linewidth = 0.6),
    panel.grid.minor  = element_blank(),
    # Border
    panel.border      = element_rect(color = "#cccccc", fill = NA, linewidth = 0.8),
    # Axes
    axis.title        = element_text(size = 13, color = "#333333", family = "sans"),
    axis.text         = element_text(size = 11, color = "#333333"),
    axis.ticks        = element_line(color = "#aaaaaa"),
    # Title
    plot.title        = element_text(size = 15, face = "bold", color = "#222222",
                                     hjust = 0.5, margin = margin(b = 10)),
    # Legend
    legend.position        = c(0.92, 0.20),
    legend.background      = element_rect(fill = "white", color = "#cccccc", linewidth = 0.5),
    legend.key             = element_rect(fill = "white"),
    legend.text            = element_text(size = 11),
    legend.margin          = margin(4, 8, 4, 8),
    # Margins
    plot.margin = margin(15, 20, 10, 10)
  ) +
  guides(color = guide_legend(override.aes = list(size = 5, shape = 16)))

print(coherence_plot2)

ggsave(
  filename = "coherence_plot2.png",
  plot = coherence_plot2,
  dpi = 500,
  width = 8,   
  height = 6   
)
################
k <- k_max# escribir el número de tópicos optimo encontrado

system.time(
  model <- FitLdaModel(
    dtm = zdtmF,        
    k = k,              
    iterations = 500,    
    burnin = 50,        
    alpha = 50 / k,     
    beta = colSums(zdtmF) / sum(zdtmF) * 100, 
    optimize_alpha = TRUE, 
    calc_likelihood = TRUE,
    calc_coherence = TRUE, 
    calc_r2 = FALSE,        
    cpus = cpus            
  )
)
beepr::beep(3)
saveRDS(model, file = "model.rds")
model <- readRDS("model.rds")

########## Visualización con LDAvis ###########################################
tm2LDAvis <- function(model) {
  if ("lda_topic_model" %in% class(model)) {
    if (ncol(model$theta) < 3) stop("The model must contain > 2 topics")
    ldavis <- LDAvis::createJSON(
      phi = model$phi,      
      theta = model$theta,  
      vocab = colnames(model$phi), 
      doc.length = slam::row_sums(model$data, na.rm = TRUE), 
      term.frequency = slam::col_sums(model$data, na.rm = TRUE) 
    )
  } else if ("LDA_VEM" %in% class(model)) { 
    post <- topicmodels::posterior(model) 
    if (ncol(post[["topics"]]) < 3) stop("The model must contain > 2 topics")
    mat <- model@wordassignments
    ldavis <- LDAvis::createJSON(
      phi = post[["terms"]],  
      theta = post[["topics"]], 
      vocab = colnames(post[["terms"]]),  
      doc.length = slam::row_sums(mat, na.rm = TRUE), 
      term.frequency = slam::col_sums(mat, na.rm = TRUE) 
    )
  } else {
    stop("Class of topic model not recognised. Must be of class LDA_VEM or lda_topic_model")
  }
  return(ldavis)
}

LDAvis::serVis(tm2LDAvis(model))
##########################################################################

model$top_terms <- GetTopTerms(phi = model$phi, M = 20) 
model$prevalence <- colSums(model$theta) / sum(model$theta) * 100

model$labels <- LabelTopics(assignments = model$theta > 0.15, 
                            dtm = zdtmF,                     
                            M = 2)                           

model$summary <- data.frame(
  topic = rownames(model$phi),   
  coherence = round(model$coherence, 3),
  prevalence = round(model$prevalence, 3),
  top_terms = apply(model$top_terms, 2, function(x) {  
    paste(x, collapse = ", ")
  }),
  stringsAsFactors = FALSE
)

write.xlsx(model$summary, file = "summary.xlsx")
label_df <- data.frame(topic = rownames(model$labels), 
                       label = model$labels, stringsAsFactors = FALSE)

top_topics <- apply(model$theta, 1, 
                    function(x) names(x)[which.max(x)][1])

top_topics <- data.frame(document = names(top_topics), top_topic = top_topics, stringsAsFactors = FALSE)

write.xlsx(top_topics, file = "toptopic.xlsx")

## Tidy Matrix Creation 
library(tidyr)
tidy_beta <- data.frame(
  topic = as.integer(stringr::str_replace_all(rownames(model$phi), "t_", "")),  
  model$phi,                                                                    
  stringsAsFactors = FALSE
) |> 
  gather(term, beta, -topic) |>  
  tibble::as_tibble()

write.xlsx(tidy_beta, file = "tidy_beta.xlsx")

tidy_gamma <- data.frame(
  document = rownames(model$theta),    
  model$theta,                         
  stringsAsFactors = FALSE
) |> 
  gather(topic, gamma, -document) |>   
  tibble::as_tibble()

## Clasificación de documentos
classifications <- tidy_gamma |> 
  dplyr::group_by(topic) |> 
  dplyr::top_n(5, gamma) |>   
  dplyr::ungroup()

write.xlsx(classifications, file = "clasificacion.xlsx")

# -------------------------------------------------------------------------
# CORRECCIÓN DE LEFT JOIN (Reemplazo del cbind que generaba error)
# -------------------------------------------------------------------------

# Seleccionamos las columnas únicas desde la base original (evitando duplicados)
metadata <- bd |> 
  dplyr::select(TI, PY, SO, C1) |> 
  dplyr::distinct(TI, .keep_all = TRUE)

# Realizamos el cruce por el nombre del documento (TI)
theta <- tidy_gamma |> 
  dplyr::left_join(metadata, by = c("document" = "TI")) |> 
  dplyr::rename(Year = PY, source = SO, countries = C1)

write.xlsx(theta, file = "theta.xlsx")
# -------------------------------------------------------------------------

## Obtención de matrices para análisis posteriores (por ejemplo HJ-Biplot)
library(dplyr)
df_source_topic <- theta |> 
  group_by(source, topic) |> 
  summarise(value = mean(gamma), .groups = "drop")
write.xlsx(df_source_topic, file = "df_source_topic.xlsx")

df_year_topic <- theta |> 
  group_by(Year, topic) |> 
  summarise(value = mean(gamma), .groups = "drop")
write.xlsx(df_year_topic, file = "df_year_topic.xlsx")

df_countries_topic <- theta |> 
  group_by(countries, topic) |> 
  summarise(value = mean(gamma), .groups = "drop")
write.xlsx(df_countries_topic, file = "df_countries_topic.xlsx")
###########################################################################

# ----------------- LIBRERÍAS NECESARIAS ---------------------------
library(dplyr)
library(ggplot2)
library(broom)
library(patchwork)
library(openxlsx)

# ----------------- DATOS Y PREPARACIÓN ----------------------------
df.h <- data.frame( 
  s = as.factor(df_year_topic$topic),  
  x = df_year_topic$Year,
  y = df_year_topic$value
)

# Ajuste de regresiones lineales por grupo
dfreg <- df.h |>
  group_by(s) |>
  do(fit = tidy(lm(y ~ x, data = .))) |>
  unnest(fit)

write.xlsx(dfreg, file = "dfCogwr.xlsx")

# Clasificación según significancia y dirección
A <- filter(dfreg, term == "x")
B <- filter(A, p.value >= 0.05)
C <- filter(A, estimate <= 0, p.value < 0.05)
D <- filter(A, estimate > 0, p.value < 0.05)

hot   <- D$s
cold  <- C$s
equal <- B$s

.dfhot   <- filter(df.h, s %in% hot)
.dfcold  <- filter(df.h, s %in% cold)
.dfequal <- filter(df.h, s %in% equal)

# ----------------- AJUSTES PARA GRÁFICOS --------------------------
df.h$s <- factor(df.h$s,
                 levels = paste0("t_", sort(as.numeric(gsub("t_", "", unique(as.character(df.h$s))))))
)

df.h <- df.h %>% 
  mutate(
    category = case_when(
      s %in% hot   ~ "hot",
      s %in% cold  ~ "cold",
      s %in% equal ~ "equal",
      TRUE         ~ NA_character_
    )
  )

library(patchwork)

# ----------------- GENERACIÓN DE GRÁFICOS POR TÓPICO --------------
plots_all <- list()

for (topic_id in levels(df.h$s)) {
  data_topic <- filter(df.h, s == topic_id)
  if (nrow(data_topic) == 0) next
  
  color <- case_when(
    topic_id %in% hot   ~ "red",
    topic_id %in% cold  ~ "blue",
    topic_id %in% equal ~ "black",
    TRUE                ~ "grey70"
  )
  
  p <- ggplot(data_topic, aes(x = x, y = y)) +
    geom_point(size = 1.2, colour = color) +
    geom_smooth(method = "lm", se = FALSE, colour = color, size = 1.2) +
    scale_x_continuous(breaks = seq(min(data_topic$x), max(data_topic$x), 5)) +
    labs(title = topic_id, x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          panel.spacing = unit(0.4, "lines"))
  
  plots_all[[topic_id]] <- p
}

# ----------------- COMPROBACIÓN Y GUARDADO ------------------------
nplots <- length(plots_all)
message("Número de gráficos generados: ", nplots)

if (nplots > 0) {
  ggsave("grafico_prueb.png", plots_all[[1]], dpi = 300)
  ncol_deseada <- 5
  panel_final <- wrap_plots(plots_all, ncol = ncol_deseada)
  
  ggsave("tendencias_pane.png",
         plot   = panel_final,
         width  = 3.5 * ncol_deseada,
         height = 3.5 * ceiling(nplots / ncol_deseada),
         dpi    = 300)
} else {
  warning("No se generaron gráficos. Revisa los datos de entrada.")
}



