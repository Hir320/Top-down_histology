rm(list = ls())
# ライブラリを読み込む
library(tidyverse)
library(plotrix)



# ディレクトリ設定 ----------------------------------------------------------------

# Define the directory where you want to search for the files
setwd('/home/rstudio/')

# ワーキングディレクトリを設定
Path = "./Data/Histology/Retinotopic_quantification/APgradation/"
setwd(Path)


# 関数 ----------------------------------------------------------------------

plot_intensity_LGN <- function(data) {
  p <- ggplot(data, aes(x = Altitude, y = mean_intensity)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, fill = 'black') +
    geom_errorbar(aes(ymax = mean_intensity + SEM, ymin = mean_intensity - SEM),
                  color = 'black',
                  position = position_dodge(0.8),
                  width = 0.3,
                  size = 0.1) +
    theme_classic() +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), expand = expansion(mult = c(0.0, 0.1))) +
    theme(
      axis.title = element_blank(),
      #legend.title = element_text(size = 6),
      #legend.text = element_text(size = 6),
      axis.text.y = element_text(size = 6, family = "Arial"),
      axis.text.x = element_text(size = 6, family = "Arial"),
      # プロット全体を透過。color = NAで枠線も消す
      plot.background = element_rect(fill = "transparent", color = NA),
      # パネル部分を透過
      panel.background = element_rect(fill = "transparent")
    )
  
  return(p)
}

scale_intensity <- function(csv_data) {
  result <- csv_data %>% 
    group_by(Region, ID) %>%
    mutate(
      Index = row_number(),
      Altitude = ifelse(Index <= 10, seq(-25, 25, by = 5)[Index], NA),
      Altitude = as.factor(Altitude),
      Mean_Scaled = (Mean - min(Mean, na.rm = TRUE)) / (max(Mean, na.rm = TRUE) - min(Mean, na.rm = TRUE))
    ) %>%
    ungroup() %>% 
    filter(Index < 11)
  
  return(result)
}



# csv読み込み -----------------------------------------------------------------

# List all CSV files in the directory
csv_files <- list.files(path = "./csvfiles/", pattern = "*.csv", full.names = TRUE)
csv_files_Ai34LGN <- list.files(path = "./LGN-GFP_Ai34data/", pattern = "*.csv", full.names = TRUE)

csv_data <- read_csv(csv_files, id = "name") %>% 
  mutate(
    name = str_extract(name, "\\d{2,3}_[A-Za-z]+(?:[A-Za-z]+)?")) %>% 
  separate(name, into = c("ID", "Region"), sep = "_", remove = FALSE) %>% 
  select(name, ID, Region, Mean)
#csv_data$Region <- fct_rev(csv_data$Region)

csv_data_Ai34LGN <- read_csv(csv_files_Ai34LGN, id = "name") %>% 
  mutate(name = str_extract(name, "Ai34_(\\d{1,3})"),
         Region = "LGN") %>% 
  separate(name, into = c("Genotype", "ID"), sep = "_", remove = FALSE) %>% 
  select(name, ID, Region, Mean)

MinMaxscaled_intensity <- scale_intensity(csv_data) 

MinMaxscaled_intensity_mean <- MinMaxscaled_intensity %>% 
  group_by(Region, Altitude) %>% 
  summarise(mean_intensity = mean(Mean_Scaled, na.rm = TRUE), 
            SEM = std.error(Mean_Scaled, na.rm = TRUE), 
            n = n())

MinMaxscaled_intensity_Ai34LGN <- scale_intensity(csv_data_Ai34LGN)
MinMaxscaled_intensity_mean_Ai34LGN <- MinMaxscaled_intensity_Ai34LGN %>% 
  group_by(Altitude) %>% 
  summarise(mean_intensity = mean(Mean_Scaled, na.rm = TRUE), 
            SEM = std.error(Mean_Scaled, na.rm = TRUE), 
            n = n())

MinMaxscaled_intensity_mean_RSCACC <- filter(MinMaxscaled_intensity_mean, Region %in% c("RSC", "ACC"))

ggplot(MinMaxscaled_intensity_mean_RSCACC, aes(x = Altitude, y = mean_intensity, fill = Region)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymax = mean_intensity + SEM, ymin = mean_intensity - SEM, color = Region), # geom_errorbarでエラーバー追加
                position = position_dodge(0.8), #geom_barと同じ値にする
                width = 0.3, # 横棒の長さ
                size = 0.1)+ # 線の太さ
  theme_classic() +
  #labs(title = "Azimuth 20~25˚", y = "Scaled Intensity", x = "Altitude (degree)") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), expand = expansion(mult = c(0.0, 0.1))) +
  scale_fill_manual(values = c("RSC" = "#ff4b00", "ACC" = "#03af7a")) +
  scale_color_manual(values = c("RSC" = "#ff4b00", "ACC" = "#03af7a")) +
  theme(#axis.title = element_text(size = 12, family = "Arial"),
        axis.title = element_blank(),
        #legend.title = element_text(size = 6),
        #legend.text = element_text(size = 6),
        legend.position = "None",
        axis.text.y = element_text(size = 6, family = "Arial"),
        axis.text.x = element_text(size = 6, family = "Arial"),
        # プロット全体を透過。color = NAで枠線も消す
        plot.background = element_rect(fill = "transparent", color = NA),
        # パネル部分を透過
        panel.background = element_rect(fill = "transparent"))
ggsave("./plots/ACCvsRSC.svg", width = 50, height = 32, units = "mm")
  
MinMaxscaled_intensity_mean_LGN <- filter(MinMaxscaled_intensity_mean, Region == "LGN")

mCerulean_LGN <- plot_intensity_LGN(MinMaxscaled_intensity_mean_LGN)
plot(mCerulean_LGN)
ggsave("./plots/LGNlayer1.svg", width = 50, height = 32, units = "mm")
Ai34LGN <- plot_intensity_LGN(MinMaxscaled_intensity_mean_Ai34LGN)
plot(Ai34LGN)
ggsave("./plots/LGNlayer1_Ai34LGN.png", width = 90, height = 57, units = "mm")

MinMaxscaled_intensity_mean_ACCplusRSC <- filter(MinMaxscaled_intensity_mean, Region == "ACCplusRSC")



# ggplot(MinMaxscaled_intensity_mean_ACCplusRSC, aes(x = Altitude, y = mean_intensity)) +
#   geom_col(position = position_dodge(width = 0.8), width = 0.7) +
#   geom_errorbar(aes(ymax = mean_intensity + SEM, ymin = mean_intensity - SEM), # geom_errorbarでエラーバー追加
#                 position = position_dodge(0.8), #geom_barと同じ値にする
#                 width = 0.2, # 横棒の長さ
#                 size = 0.5)+ # 線の太さ
#   theme_classic() +
#   labs(title = "Azimuth 20~25˚", y = "Scaled Intensity", x = "Altitude (degree)") +
#   scale_y_continuous(expand = expansion(mult = c(0.0, 0.2))) +
#   theme(axis.title = element_text(size = 12, family = "Arial"),
#         axis.text.y = element_text(size = 8, family = "Arial"),
#         axis.text.x = element_text(size = 8, family = "Arial"))
# ggsave("./plots/ACCplusRSC.png", width = 60, height = 40, units = "mm")
