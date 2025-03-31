library(tidyverse)
library(ggtext)
library(ggpubr)
library(lme4)
library(scales)
library(broom)

df<-read.csv('biolog_pm_21_25_2_temps.csv', row.names = 1)



#### rozkład danych


hist(df$values, breaks = 100, main = "Histogram absorbancji", col = "lightblue")

qqnorm(df$values)
qqline(df$values, col = 'red')

### rozkład nienormalny, prawoskośne


### dodanie czynnika 'wzrost'

df_growth <- df %>% 
  group_by(Sample,strain, replication, substrate, concentration, temperature) %>% 
  summarise(growth = values[hour==96] - values[hour == 0])

hist(df_growth$growth, breaks = 50, main = "Histogram absorbancji", col = "lightblue")


### przeskaluj dane do modelu, bo są wartości ujemne
df_growth$scaled_growth <- rescale(df_growth$growth, to = c(0,1))

df_growth <- df_growth %>%
  mutate(scaled_growth = scaled_growth + 1e-6)



# Fit the model
model <- glmer(scaled_growth ~ temperature + strain + concentration + replication + (1 | Sample),
               family = Gamma(link='log'), data = df_growth)

summary(model)

tidy(model) %>%
  arrange(p.value) %>%
  filter(p.value < 0.05) 


#check model
hist(residuals(model))
qqnorm(residuals(model))
qqline(residuals(model))

plot(fitted(model), residuals(model))
abline(h = 0, col = "red")


library(effects)
plot(allEffects(model))



########### check substances

full_model <- lmerTest::lmer(scaled_growth ~ substrate + (1|strain) + (1|temperature) + (1|concentration), 
                   data = df_growth)

summary(model)

# Ekstrakcja istotnych substancji
t_substances <- broom.mixed::tidy(full_model) %>% 
  mutate(substrate = gsub("substrate", "", term))

significant_substances <- broom.mixed::tidy(full_model) %>%
  filter(p.value < 0.05 ) %>%
  mutate(substrate = gsub("substrate", "", term)) %>%
  select(substrate, estimate, p.value) %>%
  arrange(desc(estimate))  


# Wykres pudełkowy
box<-significant_substances %>% 
  filter(estimate < 0) %>% 
  ggplot(aes(x = reorder(substrate, -estimate), y = estimate, )) +
  geom_bar(stat = "identity", fill = 'lightblue') +
  coord_flip()+
  labs(
    x = "Substance",
    y = "Estimate (Effect Size)") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.margin = margin(10,10,10,10))
#axis.text.x = element_text(angle = 90, hjust = 1, size = 10))



volcano <- t_substances %>%
  filter(term != '(Intercept)') %>% 
  ggplot(aes(x = estimate, y = -log10(p.value), color = ifelse(p.value < 0.05, "sig", "ns"))) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("sig" = "lightblue", "ns" = "lightgray")) +
  labs(x = "Effect Size (Estimate)", y = "-log10(p-value)") +
  theme(legend.position = 'none')+
  theme_minimal()+
  guides(colour=F)

ggarrange(volcano,box, widths = c(1,2), labels='AUTO')

ggsave('plots/manuscript2.png', width = 11, height = 6, bg='white')