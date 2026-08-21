############################################################
# Brooks
# Time Lag between Deforestation and Bird Extinction
# in Tropical Forest Fragments
############################################################

# Pacotes ---------------------------------------------------

# Instale caso ainda não tenha:
# install.packages(c("tidyverse", "broom"))

library(tidyverse)
library(broom)


############################################################
# 1. DADOS DO ESTUDO
############################################################

# Valores gerais utilizados pelos autores

A_total <- 25000       # área original da floresta (ha)
S_total <- 73          # pool total de espécies
z_original <- 0.15     # z para estimar riqueza original
z_fragment <- 0.25     # z para estimar riqueza no equilíbrio


# Dados da Tabela 2 de Brooks et al. (1999)

dados <- tibble(
  
  fragmento = c(
    "Malava",
    "Kisere",
    "Ikuywa",
    "Yala",
    "Kakamega"
  ),
  
  area = c(
    100,
    400,
    1450,
    1500,
    8600
  ),
  
  S_now = c(
    19,
    32,
    44,
    44,
    59
  ),
  
  # anos desde o isolamento
  tempo = c(
    101,
    63,
    20,
    24,
    82
  )
)


############################################################
# 2. CALCULAR S_original
############################################################

dados <- dados %>%
  mutate(
    
    S_original_calculado =
      S_total * (area / A_total)^z_original
    
  )


############################################################
# 3. CALCULAR S_fragment
############################################################

dados <- dados %>%
  mutate(
    
    S_fragment_calculado =
      S_total * (area / A_total)^z_fragment
    
  )


############################################################
# 4. CALCULAR O ÍNDICE DE RELAXAMENTO
############################################################

dados <- dados %>%
  mutate(
    
    I =
      (S_now - S_fragment_calculado) /
      (S_original_calculado - S_fragment_calculado)
    
  )


############################################################
# 5. EXIBIR RESULTADOS
############################################################

dados %>%
  select(
    fragmento,
    area,
    S_original_calculado,
    S_now,
    S_fragment_calculado,
    I,
    tempo
  )


############################################################
# 6. AJUSTAR O MODELO EXPONENCIAL
############################################################

dados_modelo <- dados %>%
  filter(
    !is.na(tempo),
    I > 0,
    I <= 1
  )


# Modelo:
#
# I = exp(-k*t)
#
# Tomando log:
#
# log(I) = -k*t

modelo <- lm(
  log(I) ~ tempo,
  data = dados_modelo
)

summary(modelo)


############################################################
# 8. ESTIMAR A MEIA-VIDA
############################################################

# Se:
#
# I = exp(-k*t)
#
# então:
#
# T_half = log(2)/k

k <- -coef(modelo)[["tempo"]]

meia_vida <- log(2) / k

meia_vida


############################################################
# 9. GERAR VALORES PREDITOS
############################################################

dados_modelo <- dados_modelo %>%
  mutate(
    
    I_predito = predict(
      modelo,
      newdata = dados_modelo
    ) %>% exp(),
    
    meia_vida_modelo = meia_vida
    
  )

dados_modelo


############################################################
# 10. FIGURA SEMELHANTE À FIGURA 4 DO ARTIGO
############################################################

# Criar uma sequência de tempo

tempo_seq <- seq(
  0,
  150,
  length.out = 500
)


# Curvas exponenciais utilizando as meias-vidas
# do menor e maior fragmento

curva_malava <- tibble(
  tempo = tempo_seq,
  I = 2^(-tempo / 23),
  curva = "Malava: meia-vida = 23 anos"
)

curva_kakamega <- tibble(
  tempo = tempo_seq,
  I = 2^(-tempo / 80),
  curva = "Kakamega: meia-vida = 80 anos"
)


# Gráfico

ggplot() +
  
  geom_point(
    data = dados_modelo,
    aes(
      x = tempo,
      y = I,
      label = fragmento
    ),
    size = 3
  ) +
  
  geom_text(
    data = dados_modelo,
    aes(
      x = tempo,
      y = I,
      label = fragmento
    ),
    nudge_y = 0.04
  ) +
  
  geom_line(
    data = curva_malava,
    aes(
      x = tempo,
      y = I
    ),
    linewidth = 1
  ) +
  
  geom_line(
    data = curva_kakamega,
    aes(
      x = tempo,
      y = I
    ),
    linewidth = 1
  ) +
  
  labs(
    x = "Tempo desde o isolamento (anos)",
    y = "Proporção de espécies esperadas que permanece",
    title = "Relaxamento da comunidade de aves",
    subtitle = "Reconstrução baseada em Brooks et al. (1999)"
  ) +
  
  theme_classic()


############################################################
# 11. MODELO COM A MEIA-VIDA DIRETAMENTE
############################################################

# Outra maneira de ajustar o modelo é escrever:

modelo2 <- nls(
  
  I ~ 2^(-tempo / half_life),
  
  data = dados_modelo,
  
  start = list(
    half_life = 50
  )
  
)

summary(modelo2)


############################################################
# 12. PREDIÇÕES DO MODELO
############################################################

dados_modelo <- dados_modelo %>%
  mutate(
    
    I_nls = predict(modelo2)
    
  )

dados_modelo


############################################################
# 13. RELAÇÃO ENTRE TAMANHO DO FRAGMENTO
#     E MEIA-VIDA
############################################################

# Valores publicados no artigo

meias <- tabela_publicada %>%
  select(
    fragmento,
    meia_vida_publicada
  ) %>%
  left_join(
    dados %>% select(fragmento, area),
    by = "fragmento"
  )


# Modelo de potência:
#
# meia-vida = a * area^b
#
# log(meia-vida) = log(a) + b*log(area)

modelo_area <- lm(
  
  log(meia_vida_publicada) ~ log(area),
  
  data = meias
  
)

summary(modelo_area)


# Gráfico

ggplot(
  meias,
  aes(
    x = area,
    y = meia_vida_publicada
  )
) +
  
  geom_point(size = 3) +
  
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE
  ) +
  
  scale_x_log10() +
  scale_y_log10() +
  
  labs(
    x = "Área do fragmento (ha)",
    y = "Meia-vida (anos)",
    title = "Relação entre área do fragmento e meia-vida"
  ) +
  
  theme_classic()


############################################################
# FIM
############################################################