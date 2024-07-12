
redir <- "/Users/florale/Library/CloudStorage/OneDrive-MonashUniversity/GitHub/projects/ukbiobank/"
outputdir <- "/Users/florale/Library/CloudStorage/OneDrive-MonashUniversity/PhD/projects/ukbiobank/ukb-cancer-24h/output/"

# colour pallete
c("#978787", "#BEACA2", "#EFE3E0", "#A1B2C2", "#8399AE")
# pal_time <- c("#4682b4", "#9c6755", "#f5c98e", "#dfa398")
pal_time <- c("#EAD3BF", "#B49797", "#708885", "#5A6367")

# pal_type <- wes_palette("Cavalcanti1", 15, type = "continuous")
# pal_type <- wes_palette("Zissou1", 15, type = "continuous")
# pal_type <- c("#666666", "#f6e0d2", "#9c6755", "#659794", "#dfa398", "#f5c98e", "#d65b5a", "#586085", 
#               "#6ca9c3", "#d1a391", "#7b906f", "#c0a5aa", "#7083a4", "#ad616c", "#4d3944")

# pal_type <- c("#666666", 
#               "#d65b5a", "#ba6c6e", "#efcbcb", "#dfa398", "#9c6755", "#ea967c", "#f5c98e", 
#               "#659794", "#7b906f", "#8CAACB", "#7083a4", "#395b85", "#635761",
#               "#ad616c"
# )

pal_type <- c("#ad616c",
              "#635761", "#395b85", "#7083a4", "#8CAACB", "#7b906f", "#659794", 
              "#f5c98e", "#ea967c", "#9c6755", "#dfa398", "#efcbcb", "#ba6c6e", "#d65b5a", 
              "#666666"
)

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

pal_type_quantile <- c("#ad616c",
                       "#635761", "#395b85", "#7083a4", "#8CAACB", "#7b906f", "#659794", 
                       "#f5c98e", "#ea967c", "#9c6755", "#dfa398", "#efcbcb", "#ba6c6e", "#d65b5a", 
                       "#666666",
                       "#ad616c",
                       "#635761", "#395b85", "#7083a4", "#8CAACB", "#7b906f", "#659794", 
                       "#f5c98e", "#ea967c", "#9c6755", "#dfa398", "#efcbcb", "#ba6c6e", "#d65b5a", 
                       "#666666",
                       "#ad616c",
                       "#635761", "#395b85", "#7083a4", "#8CAACB", "#7b906f", "#659794", 
                       "#f5c98e", "#ea967c", "#9c6755", "#dfa398", "#efcbcb", "#ba6c6e", "#d65b5a", 
                       "#666666"
)

pal_time_quantile <- c("#666666", "#4682b4", "#f5c98e", "#9c6755",
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

# pal_jco()(10)
scales::show_col(dutchmasters::dutchmasters_pal("milkmaid")(15))

scales::show_col(tvthemes:::brooklyn99_palette$Dark)

scales::show_col(tvthemes:::hilda_palette$Day)
scales::show_col(tvthemes:::hilda_palette$Dusk)
scales::show_col(tvthemes:::hilda_palette$Night)

scales::show_col(tvthemes:::theLastAirbender_palette$FireNation)
scales::show_col(tvthemes:::theLastAirbender_palette$AirNomads)
scales::show_col(tvthemes:::theLastAirbender_palette$EarthKingdom)
scales::show_col(tvthemes:::theLastAirbender_palette$WaterTribe)

scales::show_col(tvthemes:::westeros_palette$Stark)
scales::show_col(tvthemes:::westeros_palette$Stannis)
scales::show_col(tvthemes:::westeros_palette$Lannister)
scales::show_col(tvthemes:::westeros_palette$Tyrell)
scales::show_col(tvthemes:::westeros_palette$Targaryen)
scales::show_col(tvthemes:::westeros_palette$Martell)
scales::show_col(tvthemes:::westeros_palette$Tully)
scales::show_col(tvthemes:::westeros_palette$Greyjoy)
scales::show_col(tvthemes:::westeros_palette$Manderly)
scales::show_col(tvthemes:::westeros_palette$Arryn)

