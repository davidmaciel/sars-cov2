library(tidyverse)
library(furrr)
library(cowplot)
library(lubridate)
library(propagate)
# dados -------------------------------------------------------------------
# pars <- read_rds("parametros.rds")
cov <- read_csv("https://covid.ourworldindata.org/data/owid-covid-data.csv", locale = 
                  locale(encoding = "ASCII"))


# funções -----------------------------------------------------------------

gompertz <- function(time, a, mu, lambda){
  a*exp(-exp(mu*exp(1)/a*(lambda-time)+1))
  
}

residuos <- function(time, obs, alfa, mu, lambda){
  curva <- gompertz(time, alfa, mu, lambda)
sum((obs-curva)^2)
}

sequencia <- function(var, tamanho) {
  seq(min(var), max(var), diff(range(var))/tamanho)
}

analises <- function(pars){
  min_res <- min(pars$res)
  melhor_curva <- pars %>% filter(res == min_res)
  
  #valores previstos x valores observas
  obs_prev <- tibble(previsto = gompertz(time, melhor_curva$alfa, 
                                         melhor_curva$mu, melhor_curva$lambda),
                     days = time,
                     observado = obs) %>% 
    gather(key = "tipo", value = "mortes", -days)
  
  #projeção
  dias_ate <- (today() - ymd("2021-03-18")) %>% as.integer() %>% abs()
  obs_df <- obs_prev %>% filter(tipo == "observado")
  proj <- tibble(days = 1:dias_ate,
                 tipo = "previsto",
                 mortes = gompertz(1:dias_ate, melhor_curva$alf, 
                                   melhor_curva$mu, melhor_curva$lambda))
  proj <- bind_rows(obs_df, proj)
  list(melhor_curva,
       obs_prev,
       proj)
}




plot_casos_diarios <- function(pais,  str_pais = ""){
    observa_pais(pais, vars = "new_deaths") %>% 
    ggplot(aes(x = dias, y = new_deaths)) +
    geom_col(col = "white", fill = "black") +
    xlab("Dias desde a primeira morte confirmada") +
    ylab("Novas mortes diárias") +
    ggtitle(str_pais) +
    theme_cowplot(font_size = 10)
}


plota_grade <- function(pais, str_pais = ""){
  a<-map2(pais,str_pais, plot_casos_diarios) 
  plot_grid(plotlist = a)
}

observa_pais <- function(pais, vars = "total_deaths"){
  vars <- enquo(vars)
 x<-cov %>% filter(total_deaths >= 1 & location %in% pais) %>% 
  dplyr::select(location, date, !!vars) %>% #seleção dos países
   group_by(location) %>% 
  mutate(dias = rev((lubridate::today() - date)) %>% as.integer())  
if(min(x$dias) == 2){
  x$dias <- x$dias-1
  }
 x
}  


plota_totais <- function(paises, str_paises = sort(paises)){
  map(paises, observa_pais) %>% bind_rows() %>% 
    ggplot(aes(x = dias, y = total_deaths, col = location)) +
    geom_line() + xlab("Dias deste a primeira morte confirmada") +
    scale_color_brewer("País", type = "qual", palette = "Dark2", labels = str_paises) +
    ylab("Mortes acumuladas") +
    theme_minimal_hgrid(font_size = 12)
}


pega_lambda <- function(pais){
  observa_pais(pais) %>%  
    mutate(diff = total_deaths- lag(total_deaths)) %>% 
    summarise(lambda = which.max(diff), 
              max_diff = max(diff, na.rm = T))
   
}
plota_max_assint <- function(pais, lambda, max_dia){
  observa_pais(pais) %>% 
    filter(between(dias,lambda, max_dia)) %>%
    ggplot(aes(x = dias, y = total_deaths)) +
    geom_line() + geom_smooth(method = "lm", se = F) +
    xlab("Dias desde a primeira\nmorte confirmada") +
    ylab("Mortes acumuladas") +
    theme_cowplot(font_size = 10)
    
}


acha_mu <- function(pais){
  lambda <- pega_lambda(pais)
  lambda <- lambda$lambda
  x <- observa_pais(pais) %>% 
    filter(between(dias,lambda, lambda + 20))
  coefs <- lm(total_deaths ~ dias, data = x) %>% coef()
  coefs["dias"]
}


mu_br <- function(pais = "Brazil"){
  cov %>% 
    filter(location == pais &
             total_deaths >= 1) %>% 
    select(new_deaths) %>%
    mutate(days = 1:nrow(.)) %>% tail() %>% 
    summarise(media = mean(new_deaths)) %>% 
    unlist()
}

expande_grade <- function(lambda,mu,alfa){
  expand.grid(lambda,mu,alfa) %>% 
    rename("lambda" = Var1, "mu" = Var2, "alfa" = Var3)
}

predicao <- function(m, min_day, max_day, alfa= 0.05){
  predictNLS(m, 
             newdata = data.frame(t = min_day:max_day),
             interval = "prediction", alpha = alfa)
}

constroi_obs_prev_df <- function(obs,prevs){
  obs <- tibble(
    mortes = obs,
    dias = 1: length(obs),
    tipo = "Observado"
  )
  prevs <- tibble(
    mortes = prevs,
    dias = 1:length(prevs),
    tipo = "Previsto"
  )
  bind_rows(obs, prevs)
}

plota_obs_e_proj <- function(obs_preds,
                             int = NULL,
                             tipo = c("obs", "proj"),
                             max_day = 65) {
  tipo <- match.arg(tipo)
  g <- ggplot(obs_preds, aes(x = dias, y = mortes)) +
    geom_line() +
    scale_y_continuous(labels = (function(x)
      scales::comma(
        x,
        big.mark = ".",
        decimal.mark = ","
      ))) +
    geom_ribbon(
      data = int,
      aes(x = dias,
          ymin = lower, ymax = upper),
      alpha = 0.2,
      inherit.aes = F
    ) +
    xlab("Dias desde a primeira morte confirmada") +
    ylab("Mortes acumuladas") +
    theme_minimal_hgrid(font_size = 11)
  if (tipo == "obs") {
    g <- ggplot(obs_preds, aes(x = dias, y = mortes, col = tipo)) +
      geom_line(size = 1, alpha = 0.5) +
      scale_color_brewer("", type = "qual", palette = "Set1") +
      scale_y_continuous(labels = (
        function(x)
          scales::comma(x, big.mark = ".", decimal.mark = ",")
      )) +
      geom_ribbon(
        data = int,
        aes(x = dias,
            ymin = lower, ymax = upper),
        alpha = 0.1,
        inherit.aes = F
      ) +
      xlab("Dias desde a primeira morte confirmada") +
      ylab("Mortes acumuladas") + coord_cartesian(xlim = c(0, max_day), ylim = c(0, 21000)) +
      theme_minimal_hgrid(font_size = 11)
  }
  g
}

constro_int <- function(prev){
  prev %>% 
    dplyr::select("upper" = `Prop.97.5%`,
                  "lower" = `Prop.2.5%`) %>% 
    mutate(dias = 1:nrow(.))
}

tabula_prev <- function(prev, 
                        max_day = "2021-03-18", 
                        today = "2020-05-21") {
  prev %>% 
  dplyr::select("previsao" = Prop.Mean.2,
                "lim.inf" = `Prop.2.5%`,
                "lim.sup" = `Prop.97.5%`) %>% 
  mutate(lim.inf = if_else(lim.inf < 0, 0, lim.inf),
         dias = seq(ymd("2020-03-19"), ymd(max_day), by = 1)-1) %>% 
  filter(dias > ymd(today)) %>% 
    mutate_at(c("previsao","lim.inf","lim.sup"), round) %>% 
    left_join(observa_pais("Brazil"), by = c("dias" = "date"))
}

# analises ----------------------------------------------------------------

#1) figura 2----
f2 <- plota_grade(c("Brazil", "Italy", "United States"),
            c("Brasil", "Itália", "Estados Unidos"))
save_plot("figura_2.jpeg", f2)

#2) figura 3----
f3 <- plota_totais(c("Brazil", "Italy", "United States"),
             c("Brasil", "Itália", "Estados Unidos"))
save_plot("figura_3.jpeg", f3)

#3)lambdas----
lambda_br <- pega_lambda("Brazil")

lambda_ita <- pega_lambda("Italy")
lambda_usa <- pega_lambda("United States")

#4)assintotas máximas----
ass_usa <- plota_max_assint("United States", lambda_usa$lambda, lambda_usa$lambda + 20)

ass_ita <- plota_max_assint("Italy", lambda_ita$lambda, lambda_ita$lambda + 20)
grid_ass <- plot_grid(ass_usa, ass_ita,
                      labels = c("Estados Unidos","Itália"),
                      label_size = 12,align = "h", 
                      hjust = c(-0.55, -1.9))
save_plot("figura_4.jpeg", grid_ass)


mu_ita <- acha_mu("Italy")
mu_usa <- acha_mu("United States")

#5) expansão de grade----
lambda <- seq(lambda_br$lambda-20, lambda_br$lambda + 20, by = 1)
mu <- seq(mu_ita, mu_usa, by = 10)
alfa <- seq(44000, 1150000, length.out = 10000) %>% round()

grade <- expande_grade(lambda, mu, alfa)



#6) cálculo paralelizado dos resídos----
plan("cluster")
system.time(grade <- grade %>% 
  mutate(
    res = future_pmap_dbl(list(
      mu = mu,
      lambda = lambda,
      alfa = alfa), residuos, time = 1:length(obs), obs = obs)))

#7) melhor curva----
melhor_curva <- grade %>% filter(
  res == min(res)
) %>% dplyr::select(-res) %>% unlist()

names(melhor_curva) <- c("lambda", "mu", "a")

obs <- observa_pais("Brazil")
obs <- obs$total_deaths

prev_obs <- constroi_obs_prev_df(obs,
                                 prev = gompertz(
                                   1:length(obs),
                                   a = melhor_curva["a"],
                                   lambda = melhor_curva["lambda"],
                                   mu = melhor_curva["mu"]
                                 ))
f5 <- ggplot(prev_obs, aes(x = dias, y = mortes, col = tipo)) +
  geom_line(size = 1, alpha = 0.5) +
  scale_color_brewer("", type = "qual", palette = "Set1") +
  xlab("Dias desde a primeira morte confirmada") +
  ylab("Mortes acumuladas") +
  theme_minimal_hgrid(font_size = 12)
  
save_plot("figura_5.jpeg", f5)

#8)modelo não linear e predição com intervalo----
rm(grade, f2,f3,f5)
m_nls <- nls(formula("y~a*exp(-exp(mu*exp(1)/a*(lambda-t)+1))"),
         data.frame(y = obs, t = 1:length(obs)),
         melhor_curva)


preds <- predicao(m_nls, 365, 366)



obs_x_preds <- constroi_obs_prev_df(
  obs = obs,
  prevs = preds$summary$Prop.Mean.2 %>% round()
)

int <- constro_int(preds$summary)


f5 <- plota_obs_e_proj(obs_x_preds, int, tipo = "obs") 

f6 <- plota_obs_e_proj(obs_x_preds, int, tipo = "proj") 

save_plot("figura_6.jpeg", f6)


tab <- tabula_prev(preds$summary)


#9)tabela com os valores previstos para os próximos dias

tab <- tabula_prev(preds$summary,today = "2020-03-18" ) %>% 
  dplyr::select("Data" = dias,
         "Lim.inf" = lim.inf,
         "Lim.sup" = lim.sup,
         "Total_de_óbitos_previstos" = previsao,
         "Total_de_óbitos_observados" = total_deaths) 

write.csv2(tab, file = "mortalidade.csv")


