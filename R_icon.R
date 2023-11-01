library(hexSticker)
library(ggplot2)
# Helper theme for ggplot icon
theme_icon <- function () {
  theme_void() +
    theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA)
    )
}
library(hexSticker)
p <- ggplot(iris, aes(Species, Sepal.Length)) +
  geom_boxplot(color = "white", fill = "transparent") +
  theme_icon()
ggsave(
  filename = "img/boxplot-icon_72px.png", p,
  dpi=72, width = 1, height = 1, bg = "transparent"
)
p2 <- ggplot()+ theme_icon()
p.sticker <- sticker(
  p2, package=" ", p_size=3,
  s_x=1, s_y=1.1, s_width=1.3, s_height=1.5,
  h_color = "#478bca", h_fill = "#478bca",
  filename="img/icon.png"
)
p.sticker
