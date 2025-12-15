if (Sys.info()[["user"]] %in% c("florale")) {
  redir <- "/Users/florale/Library/CloudStorage/OneDrive-Personal/github/projects/ukbiobank/"
  outputdir <- "/Users/florale/Library/CloudStorage/OneDrive-Personal/monash/papers/ukbiobank/ukb-cancer-24h/output/"
} else if (Sys.info()[["user"]] %in% c("flee0016")) {
  redir <- "/Users/flee0016/Library/CloudStorage/OneDrive-Personal/github/projects/ukbiobank/"
  outputdir <- "/Users/flee0016/Library/CloudStorage/OneDrive-Personal/monash/papers/ukbiobank/ukb-cancer-24h/output/"
}


# colour pallete
c("#978787", "#BEACA2", "#EFE3E0", "#A1B2C2", "#8399AE")
# pal_time <- c("#4682b4", "#9c6755", "#f5c98e", "#dfa398")


# pal_type <- wesanderson::wes_palette("Cavalcanti1", 15, type = "continuous")
# pal_type <- c("#666666", "#f6e0d2", "#9c6755", "#659794", "#dfa398", "#f5c98e", "#d65b5a", "#586085",
#               "#6ca9c3", "#d1a391", "#7b906f", "#c0a5aa", "#7083a4", "#ad616c", "#4d3944")
pal <- c("#708885", "#A9A9A9", "#ba6c6e")
pal_time <- c("#708885", "#A9A9A9", "#D2A7A7", "#ba6c6e", "#4E2F26")
pal_combined <- c("#708885", "#A9A9A9", "#4E2F26", "#ba6c6e", "#CA8F90", "#E4C7C7")
pal_combined <- c("#E4C7C7", "#CA8F90", "#ba6c6e", "#4E2F26", "#A9A9A9", "#708885")

pal_sen <- c(

  # "#EFE3E0", "#A1B2C2", "#8399AE", "#7083a4",
  "#763C3C", "#944C4C",  "#ba6c6e", "#CA8F90", "#d2a0a1", "#E4C7C7", 
  "#EAD3BF", "#EEAB96", "#CF9C81", "#c48462", "#978787", "#4E2F26", 
  "#708885"
)
# scales::show_col(wesanderson::wes_palette("Zissou1", 12, type = "continuous"))
pal_sen <- 
c(
  "#3A9AB2", "#6AAFBF", "#8AB8B8", "#9FBFA8",
  "#B4C58D", "#CDC965", "#DFC131", "#E4AB0C",
  "#E88F05", "#EC7304", "#EF5102", "#F11B00",
  "#708885"

)
pal_type <- c(
  "#666666", 
  "#586085", "#9D93B9", "#3b3960",
  "#4682b4", "#8CAACB", "#4F7375", "#f5c98e",
  "#EAD3BF", "#9c6755", "#dfa398", "#E2C3C3",
  "#ba6c6e", "#944C4C", "#4E2F26", "#A9A9A9",
  "#708885"
)

pal_type <- c(
  "#666666", 
  "#4E2F26", "#944C4C", 
  "#ba6c6e", "#E2C3C3",
  "#9c6755", "#dfa398", 
  "#f5c98e", "#EAD3BF", 
  "#C1D1E3", 
  "#8CAACB", "#6A84A7", 
  "#9D93B9", 
  "#586085", "#3b3960",
  "#A9A9A9",
  "#708885"
)

# pal <- c("#a8a8a8", "#A9A9A9", "#ba6c6e")
# pal_time <- c("#a8a8a8", "#A9A9A9", "#D2A7A7", "#B87474", "#4E2F26")
# pal_combined <- c("#a8a8a8", "#A9A9A9", "#4E2F26", "#ba6c6e", "#CA8F90", "#E4C7C7")
# pal_combined <- c("#E4C7C7", "#CA8F90", "#ba6c6e", "#4E2F26", "#A9A9A9", "#a8a8a8")
#
# pal_type <- c("#586085",
#               "#9D93B9", "#3b3960", "#4682b4", "#8CAACB", "#4F7375", "#7b906f", "#f5c98e", "#ea967c",
#               "#9c6755",
#               "#dfa398", "#efcbcb", "#d65b5a", "#944C4C",
#
#               "#4E2F26", "#A9A9A9",
#               "#a8a8a8"
# )

# pal_type <- c("#708885", "#A9A9A9",
#               "#4E2F26", "#944C4C", "#d65b5a", "#efcbcb", "#dfa398", "#9c6755",  "#ea967c", "#f5c98e",
#               "#4F7375", "#8CAACB", "#4682b4", "#3b3960", "#9D93B9",
#               "#586085"
# )

# pal_type <- c("#1e2142",
#               "#395b85", "#7083a4", "#8CAACB", "#7b906f", "#659794",
#               "#f5c98e", "#ea967c", "#9c6755", "#dfa398", "#efcbcb", "#c24841", "#ba6c6e",
#               "#3d251e"
# )

# pal_type_quantile <- c("#666666",
#                        "#d65b5a", "#ba6c6e", "#efcbcb", "#dfa398", "#9c6755", "#ea967c", "#f5c98e",
#                        "#659794", "#7b906f", "#8CAACB", "#7083a4", "#395b85", "#635761",
#                        "#ad616c",
#                        "#666666",
#                        "#d65b5a", "#ba6c6e", "#efcbcb", "#dfa398", "#9c6755", "#ea967c", "#f5c98e",
#                        "#659794", "#7b906f", "#8CAACB", "#7083a4", "#395b85", "#635761",
#                        "#ad616c",
#                        "#666666",
#                        "#d65b5a", "#ba6c6e", "#efcbcb", "#dfa398", "#9c6755", "#ea967c", "#f5c98e",
#                        "#659794", "#7b906f", "#8CAACB", "#7083a4", "#395b85", "#635761",
#                        "#ad616c"
# )

pal_type_quantile <- c(
  "#666666", "#635761", "#395b85", "#7083a4", "#8CAACB", "#7b906f", "#659794",
  "#f5c98e", "#ea967c", "#9c6755", "#dfa398", "#efcbcb", "#ba6c6e", 
  # "#d65b5a",
  "#944C4C", 
  "#4E2F26", 

  "#9D93B9", 
  "#586085", "#3b3960",
  "#708885",


  "#666666",
  "#ad616c",
  "#635761", "#395b85", "#7083a4", "#8CAACB", "#7b906f", "#659794",
  "#f5c98e", "#ea967c", "#9c6755", "#dfa398", "#efcbcb", "#ba6c6e", "#d65b5a",
  "#666666",
  "#ad616c",
  "#635761", "#395b85", "#7083a4", "#8CAACB", "#7b906f", "#659794",
  "#f5c98e", "#ea967c", "#9c6755", "#dfa398", "#efcbcb", "#ba6c6e", "#d65b5a",
  "#708885"
)

pal_time_quantile <- c(
  "#666666", "#4682b4", "#f5c98e", "#9c6755",
  "#666666", "#4682b4", "#f5c98e", "#9c6755",
  "#666666", "#4682b4", "#f5c98e", "#9c6755",
  "#666666", "#4682b4", "#f5c98e", "#9c6755"
)

# pal_type_quantile <- c(
#   alpha("#666666", .5),
#   alpha("#d65b5a", .5),
#   alpha("#ba6c6e", .5),
#   alpha("#efcbcb", .5),
#   alpha("#dfa398", .5),
#   alpha("#9c6755", .5),
#   alpha("#ea967c", .5),
#   alpha("#f5c98e", .5),
#   alpha("#659794", .5),
#   alpha("#7b906f", .5),
#   alpha("#8CAACB", .5),
#   alpha("#7083a4", .5),
#   alpha("#395b85", .5),
#   alpha("#635761", .5),
#   alpha("#ad616c", .5),
#
#   alpha("#666666", .7),
#   alpha("#d65b5a", .7),
#   alpha("#ba6c6e", .7),
#   alpha("#efcbcb", .7),
#   alpha("#dfa398", .7),
#   alpha("#9c6755", .7),
#   alpha("#ea967c", .7),
#   alpha("#f5c98e", .7),
#   alpha("#659794", .7),
#   alpha("#7b906f", .7),
#   alpha("#8CAACB", .7),
#   alpha("#7083a4", .7),
#   alpha("#395b85", .7),
#   alpha("#635761", .7),
#   alpha("#ad616c", .7),
#   "#666666" ,
#   "#d65b5a",
#   "#ba6c6e",
#   "#efcbcb",
#   "#dfa398",
#   "#9c6755",
#   "#ea967c",
#   "#f5c98e",
#   "#659794",
#   "#7b906f",
#   "#8CAACB",
#   "#7083a4",
#   "#395b85",
#   "#635761",
#   "#ad616c",
#   "#666666",
#   "#d65b5a",
#   "#ba6c6e",
#   "#efcbcb",
#   "#dfa398",
#   "#9c6755",
#   "#ea967c",
#   "#f5c98e",
#   "#659794",
#   "#7b906f",
#   "#8CAACB",
#   "#7083a4",
#   "#395b85",
#   "#635761",
#   "#ad616c"
# )

# ggsci::pal_jco()(10)
# scales::show_col(tvthemes:::brooklyn99_palette$Dark)
#
# scales::show_col(tvthemes:::hilda_palette$Day)
# scales::show_col(tvthemes:::hilda_palette$Dusk)
# scales::show_col(tvthemes:::hilda_palette$Night)
#
# scales::show_col(tvthemes:::theLastAirbender_palette$FireNation)
# scales::show_col(tvthemes:::theLastAirbender_palette$AirNomads)
# scales::show_col(tvthemes:::theLastAirbender_palette$EarthKingdom)
# scales::show_col(tvthemes:::theLastAirbender_palette$WaterTribe)
#
# scales::show_col(tvthemes:::westeros_palette$Stark)
# scales::show_col(tvthemes:::westeros_palette$Stannis)
# scales::show_col(tvthemes:::westeros_palette$Lannister)
# scales::show_col(tvthemes:::westeros_palette$Tyrell)
# scales::show_col(tvthemes:::westeros_palette$Targaryen)
# scales::show_col(tvthemes:::westeros_palette$Martell)
# scales::show_col(tvthemes:::westeros_palette$Tully)
# scales::show_col(tvthemes:::westeros_palette$Greyjoy)
# scales::show_col(tvthemes:::westeros_palette$Manderly)
# scales::show_col(tvthemes:::westeros_palette$Arryn)
#
