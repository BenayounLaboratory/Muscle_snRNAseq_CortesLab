# Venn Diagrams

install.packages("ggvenn")

library(ggvenn)
library(readxl)

setwd("C:/Users/Tia/CortesLab/")


#Myo IIB 
# Read data from Excel
F_runners <- read_excel("IIB GO terms for analysis.xlsx", sheet = 1)
F_runners <- as.character(F_runners[[1]])

F_TFEB <- read_excel("IIB GO terms for analysis.xlsx", sheet = 2)
F_TFEB <- as.character(F_TFEB[[1]])

M_runners <- read_excel("IIB GO terms for analysis.xlsx", sheet = 3)
M_runners <- as.character(M_runners[[1]])

M_TFEB <- read_excel("IIB GO terms for analysis.xlsx", sheet = 4)
M_TFEB <- as.character(M_TFEB[[1]])

# Combine into a named list for ggvenn
go_lists <- list(
  "Male Runner"   = M_runners,
  "Female TFEB"   = F_TFEB,
  "Male TFEB"     = M_TFEB,
  "Female Runner" = F_runners
)

# Define colors for petals (in order)
petal_colors <- c("#87CEFA","#993399","#1E90FF","#FF69B4")
# (Female runner = bright pink, Female TFEB = light pink,
#  Male runner = bright blue, Male TFEB = light blue)

# Plot the Venn diagram
ggvenn(
  go_lists,
  fill_color = petal_colors,
  fill_alpha = 0.6,
  stroke_color = "white",
  stroke_size = 1.2,
  text_size = 5,
  set_name_size = 5,
  show_percentage = FALSE
) +
  ggtitle("Myonuclei IIB") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )


#Myo IIX
# Read data from Excel
F_runners <- read_excel("IIX GO terms for analysis.xlsx", sheet = 1)
F_runners <- as.character(F_runners[[1]])

F_TFEB <- read_excel("IIX GO terms for analysis.xlsx", sheet = 2)
F_TFEB <- as.character(F_TFEB[[1]])

M_runners <- read_excel("IIX GO terms for analysis.xlsx", sheet = 3)
M_runners <- as.character(M_runners[[1]])

M_TFEB <- read_excel("IIX GO terms for analysis.xlsx", sheet = 4)
M_TFEB <- as.character(M_TFEB[[1]])

# Combine into a named list for ggvenn
go_lists <- list(
  "Male Runner"   = M_runners,
  "Female TFEB"   = F_TFEB,
  "Male TFEB"     = M_TFEB,
  "Female Runner" = F_runners
)

# Define colors for petals (in order)
petal_colors <- c("#87CEFA","#993399","#1E90FF","#FF69B4")
# (Female runner = bright pink, Female TFEB = light pink,
#  Male runner = bright blue, Male TFEB = light blue)

# Plot the Venn diagram
ggvenn(
  go_lists,
  fill_color = petal_colors,
  fill_alpha = 0.6,
  stroke_color = "white",
  stroke_size = 1.2,
  text_size = 5,
  set_name_size = 5,
  show_percentage = FALSE
) +
  ggtitle("Myonuclei IIX") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )


#Endothelial 
# Read data from Excel
F_runners <- read_excel("Endothelial GO terms for analysis.xlsx", sheet = 1)
F_runners <- as.character(F_runners[[1]])

F_TFEB <- read_excel("Endothelial GO terms for analysis.xlsx", sheet = 2)
F_TFEB <- as.character(F_TFEB[[1]])

M_runners <- read_excel("Endothelial GO terms for analysis.xlsx", sheet = 3)
M_runners <- as.character(M_runners[[1]])

M_TFEB <- read_excel("Endothelial GO terms for analysis.xlsx", sheet = 4)
M_TFEB <- as.character(M_TFEB[[1]])

# Combine into a named list for ggvenn
go_lists <- list(
  "Male Runner"   = M_runners,
  "Female TFEB"   = F_TFEB,
  "Male TFEB"     = M_TFEB,
  "Female Runner" = F_runners
)

# Define colors for petals (in order)
petal_colors <- c("#87CEFA","#993399","#1E90FF","#FF69B4")
# (Female runner = bright pink, Female TFEB = light pink,
#  Male runner = bright blue, Male TFEB = light blue)

# Plot the Venn diagram
ggvenn(
  go_lists,
  fill_color = petal_colors,
  fill_alpha = 0.6,
  stroke_color = "white",
  stroke_size = 1.2,
  text_size = 5,
  set_name_size = 5,
  show_percentage = FALSE
) +
  ggtitle("Endothelial") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )


#SMC 
# Read data from Excel
F_runners <- read_excel("SMC GO terms for analysis.xlsx", sheet = 1)
F_runners <- as.character(F_runners[[1]])

F_TFEB <- read_excel("SMC GO terms for analysis.xlsx", sheet = 2)
F_TFEB <- as.character(F_TFEB[[1]])

M_runners <- read_excel("SMC GO terms for analysis.xlsx", sheet = 3)
M_runners <- as.character(M_runners[[1]])

M_TFEB <- read_excel("SMC GO terms for analysis.xlsx", sheet = 4)
M_TFEB <- as.character(M_TFEB[[1]])

# Combine into a named list for ggvenn
go_lists <- list(
  "Male Runner"   = M_runners,
  "Female TFEB"   = F_TFEB,
  "Male TFEB"     = M_TFEB,
  "Female Runner" = F_runners
)

# Define colors for petals (in order)
petal_colors <- c("#87CEFA","#993399","#1E90FF","#FF69B4")
# (Female runner = bright pink, Female TFEB = light pink,
#  Male runner = bright blue, Male TFEB = light blue)

# Plot the Venn diagram
ggvenn(
  go_lists,
  fill_color = petal_colors,
  fill_alpha = 0.6,
  stroke_color = "white",
  stroke_size = 1.2,
  text_size = 5,
  set_name_size = 5,
  show_percentage = FALSE
) +
  ggtitle("SMC") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

#FAPS
# Read data from Excel
F_runners <- read_excel("FAPS GO_BP terms for analysis.xlsx", sheet = 1)
F_runners <- as.character(F_runners[[1]])

F_TFEB <- read_excel("FAPS GO_BP terms for analysis.xlsx", sheet = 2)
F_TFEB <- as.character(F_TFEB[[1]])

M_runners <- read_excel("FAPS GO_BP terms for analysis.xlsx", sheet = 3)
M_runners <- as.character(M_runners[[1]])

M_TFEB <- read_excel("FAPS GO_BP terms for analysis.xlsx", sheet = 4)
M_TFEB <- as.character(M_TFEB[[1]])

# Combine into a named list for ggvenn
go_lists <- list(
  "Male Runner"   = M_runners,
  "Female TFEB"   = F_TFEB,
  "Male TFEB"     = M_TFEB,
  "Female Runner" = F_runners
)

# Define colors for petals (in order)
petal_colors <- c("#87CEFA","#993399","#1E90FF","#FF69B4")
# (Female runner = bright pink, Female TFEB = light pink,
#  Male runner = bright blue, Male TFEB = light blue)

# Plot the Venn diagram
ggvenn(
  go_lists,
  fill_color = petal_colors,
  fill_alpha = 0.6,
  stroke_color = "white",
  stroke_size = 1.2,
  text_size = 3.8,
  set_name_size = 5,
) +
  ggtitle("FAPS") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )


#Myo IIA 
# Read data from Excel
F_runners <- read_excel("IIA GO terms for analysis.xlsx", sheet = 1)
F_runners <- as.character(F_runners[[1]])

F_TFEB <- read_excel("IIA GO terms for analysis.xlsx", sheet = 2)
F_TFEB <- as.character(F_TFEB[[1]])

M_runners <- read_excel("IIA GO terms for analysis.xlsx", sheet = 3)
M_runners <- as.character(M_runners[[1]])

M_TFEB <- read_excel("IIA GO terms for analysis.xlsx", sheet = 4)
M_TFEB <- as.character(M_TFEB[[1]])

# Combine into a named list for ggvenn
go_lists <- list(
  "Male Runner"   = M_runners,
  "Female TFEB"   = F_TFEB,
  "Male TFEB"     = M_TFEB,
  "Female Runner" = F_runners
)

# Define colors for petals (in order)
petal_colors <- c("#87CEFA","#993399","#1E90FF","#FF69B4")
# (Female runner = bright pink, Female TFEB = light pink,
#  Male runner = bright blue, Male TFEB = light blue)

# Plot the Venn diagram
ggvenn(
  go_lists,
  fill_color = petal_colors,
  fill_alpha = 0.6,
  stroke_color = "white",
  stroke_size = 1.2,
  text_size = 3.8,
  set_name_size = 5,
) +
  ggtitle("FAPS") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )






# Only GO BP: 
#Myo IIB 

# Read data from Excel
F_runners <- read_excel("IIB GO terms for analysis.xlsx", sheet = 1)
F_runners <- as.character(F_runners[[1]])

F_TFEB <- read_excel("IIB GO terms for analysis.xlsx", sheet = 2)
F_TFEB <- as.character(F_TFEB[[1]])

M_runners <- read_excel("IIB GO terms for analysis.xlsx", sheet = 3)
M_runners <- as.character(M_runners[[1]])

M_TFEB <- read_excel("IIB GO terms for analysis.xlsx", sheet = 4)
M_TFEB <- as.character(M_TFEB[[1]])

F_runners <- F_runners[grepl("^GOBP_", F_runners)]
F_TFEB   <- F_TFEB[grepl("^GOBP_", F_TFEB)]
M_runners <- M_runners[grepl("^GOBP_", M_runners)]
M_TFEB   <- M_TFEB[grepl("^GOBP_", M_TFEB)]

# Combine into a named list for ggvenn
go_lists <- list(
  "Male Runner"   = M_runners,
  "Female TFEB"   = F_TFEB,
  "Male TFEB"     = M_TFEB,
  "Female Runner" = F_runners
)

# Define colors for petals
petal_colors <- c("#87CEFA","#993399","#1E90FF","#FF69B4")

# Plot the Venn diagram
ggvenn(
  go_lists,
  fill_color = petal_colors,
  fill_alpha = 0.6,
  stroke_color = "white",
  stroke_size = 1.2,
  text_size = 3.7,
  set_name_size = 5,
) +
  ggtitle("Myonuclei IIB") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

print(F_runners)


# --- Find the 4-way overlap ---
ggvenn_myoIIB_4way <- Reduce(intersect, list(F_runners, F_TFEB, M_runners, M_TFEB))

# --- Convert to data frame for Excel export ---
ggvenn_myoIIB_4way_df <- data.frame(GOBP_Term = ggvenn_myoIIB_4way)

library(writexl)

# --- Save to Excel file ---
write_xlsx(ggvenn_myoIIB_4way_df, "ggvenn_myoIIB_4way.xlsx")




#IIX  

# Read data from Excel
F_runners <- read_excel("IIX GO terms for analysis.xlsx", sheet = 1)
F_runners <- as.character(F_runners[[1]])

F_TFEB <- read_excel("IIX GO terms for analysis.xlsx", sheet = 2)
F_TFEB <- as.character(F_TFEB[[1]])

M_runners <- read_excel("IIX GO terms for analysis.xlsx", sheet = 3)
M_runners <- as.character(M_runners[[1]])

M_TFEB <- read_excel("IIX GO terms for analysis.xlsx", sheet = 4)
M_TFEB <- as.character(M_TFEB[[1]])

F_runners <- F_runners[grepl("^GOBP_", F_runners)]
F_TFEB   <- F_TFEB[grepl("^GOBP_", F_TFEB)]
M_runners <- M_runners[grepl("^GOBP_", M_runners)]
M_TFEB   <- M_TFEB[grepl("^GOBP_", M_TFEB)]

# Combine into a named list for ggvenn
go_lists <- list(
  "Male Runner"   = M_runners,
  "Female TFEB"   = F_TFEB,
  "Male TFEB"     = M_TFEB,
  "Female Runner" = F_runners
)

# Define colors for petals
petal_colors <- c("#87CEFA","#993399","#1E90FF","#FF69B4")

# Plot the Venn diagram
ggvenn(
  go_lists,
  fill_color = petal_colors,
  fill_alpha = 0.6,
  stroke_color = "white",
  stroke_size = 1.2,
  text_size = 3.7,
  set_name_size = 5,
) +
  ggtitle("Myonuclei IIX") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

print(F_runners)


# --- Find the 4-way overlap ---
ggvenn_myoIIX_4way <- Reduce(intersect, list(F_runners, F_TFEB, M_runners, M_TFEB))

# --- Convert to data frame for Excel export ---
ggvenn_myoIIX_4way_df <- data.frame(GOBP_Term = ggvenn_myoIIX_4way)

# --- Save to Excel file ---
write_xlsx(ggvenn_myoIIX_4way_df, "ggvenn_myoIIX_4way.xlsx")



#Endothelial 
# Read data from Excel
F_runners <- read_excel("Endothelial GO terms for analysis.xlsx", sheet = 1)
F_runners <- as.character(F_runners[[1]])

F_TFEB <- read_excel("Endothelial GO terms for analysis.xlsx", sheet = 2)
F_TFEB <- as.character(F_TFEB[[1]])

M_runners <- read_excel("Endothelial GO terms for analysis.xlsx", sheet = 3)
M_runners <- as.character(M_runners[[1]])

M_TFEB <- read_excel("Endothelial GO terms for analysis.xlsx", sheet = 4)
M_TFEB <- as.character(M_TFEB[[1]])

F_runners <- F_runners[grepl("^GOBP_", F_runners)]
F_TFEB   <- F_TFEB[grepl("^GOBP_", F_TFEB)]
M_runners <- M_runners[grepl("^GOBP_", M_runners)]
M_TFEB   <- M_TFEB[grepl("^GOBP_", M_TFEB)]                   

# Combine into a named list for ggvenn
go_lists <- list(
  "Male Runner"   = M_runners,
  "Female TFEB"   = F_TFEB,
  "Male TFEB"     = M_TFEB,
  "Female Runner" = F_runners
)

# Define colors for petals (in order)
petal_colors <- c("#87CEFA","#993399","#1E90FF","#FF69B4")
# (Female runner = bright pink, Female TFEB = light pink,
#  Male runner = bright blue, Male TFEB = light blue)

# Plot the Venn diagram
ggvenn(
  go_lists,
  fill_color = petal_colors,
  fill_alpha = 0.6,
  stroke_color = "white",
  stroke_size = 1.2,
  text_size = 3.7,
  set_name_size = 5,
) +
  ggtitle("Endothelial") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

# --- Find the 4-way overlap ---
ggvenn_endo_4way <- Reduce(intersect, list(F_runners, F_TFEB, M_runners, M_TFEB))

# --- Convert to data frame for Excel export ---
ggvenn_endo_4way_df <- data.frame(GOBP_Term = ggvenn_endo_4way)

# --- Save to Excel file ---
write_xlsx(ggvenn_endo_4way_df, "ggvenn_endo_4way.xlsx")

#SMC 
# Read data from Excel
F_runners <- read_excel("SMC GO terms for analysis.xlsx", sheet = 1)
F_runners <- as.character(F_runners[[1]])

F_TFEB <- read_excel("SMC GO terms for analysis.xlsx", sheet = 2)
F_TFEB <- as.character(F_TFEB[[1]])

M_runners <- read_excel("SMC GO terms for analysis.xlsx", sheet = 3)
M_runners <- as.character(M_runners[[1]])

M_TFEB <- read_excel("SMC GO terms for analysis.xlsx", sheet = 4)
M_TFEB <- as.character(M_TFEB[[1]])

F_runners <- F_runners[grepl("^GOBP_", F_runners)]
F_TFEB   <- F_TFEB[grepl("^GOBP_", F_TFEB)]
M_runners <- M_runners[grepl("^GOBP_", M_runners)]
M_TFEB   <- M_TFEB[grepl("^GOBP_", M_TFEB)]                   

# Combine into a named list for ggvenn
go_lists <- list(
  "Male Runner"   = M_runners,
  "Female TFEB"   = F_TFEB,
  "Male TFEB"     = M_TFEB,
  "Female Runner" = F_runners
)

# Define colors for petals (in order)
petal_colors <- c("#87CEFA","#993399","#1E90FF","#FF69B4")
# (Female runner = bright pink, Female TFEB = light pink,
#  Male runner = bright blue, Male TFEB = light blue)

# Plot the Venn diagram
ggvenn(
  go_lists,
  fill_color = petal_colors,
  fill_alpha = 0.6,
  stroke_color = "white",
  stroke_size = 1.2,
  text_size = 3.7,
  set_name_size = 5,
) +
  ggtitle("SMC") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

# --- Find the 4-way overlap ---
ggvenn_SMC_4way <- Reduce(intersect, list(F_runners, F_TFEB, M_runners, M_TFEB))

# --- Convert to data frame for Excel export ---
ggvenn_SMC_4way_df <- data.frame(GOBP_Term = ggvenn_SMC_4way)

# --- Save to Excel file ---
write_xlsx(ggvenn_SMC_4way_df, "ggvenn_SMC_4way.xlsx")





#FAPS 
# Read data from Excel
F_runners <- read_excel("FAPS GO_BP terms for analysis.xlsx", sheet = 1)
F_runners <- as.character(F_runners[[1]])

F_TFEB <- read_excel("FAPS GO_BP terms for analysis.xlsx", sheet = 2)
F_TFEB <- as.character(F_TFEB[[1]])

M_runners <- read_excel("FAPS GO_BP terms for analysis.xlsx", sheet = 3)
M_runners <- as.character(M_runners[[1]])

M_TFEB <- read_excel("FAPS GO_BP terms for analysis.xlsx", sheet = 4)
M_TFEB <- as.character(M_TFEB[[1]])

F_runners <- F_runners[grepl("^GOBP_", F_runners)]
F_TFEB   <- F_TFEB[grepl("^GOBP_", F_TFEB)]
M_runners <- M_runners[grepl("^GOBP_", M_runners)]
M_TFEB   <- M_TFEB[grepl("^GOBP_", M_TFEB)]                   

# Combine into a named list for ggvenn
go_lists <- list(
  "Male Runner"   = M_runners,
  "Female TFEB"   = F_TFEB,
  "Male TFEB"     = M_TFEB,
  "Female Runner" = F_runners
)

# Define colors for petals (in order)
petal_colors <- c("#87CEFA","#993399","#1E90FF","#FF69B4")
# (Female runner = bright pink, Female TFEB = light pink,
#  Male runner = bright blue, Male TFEB = light blue)

# Plot the Venn diagram
ggvenn(
  go_lists,
  fill_color = petal_colors,
  fill_alpha = 0.6,
  stroke_color = "white",
  stroke_size = 1.2,
  text_size = 3.7,
  set_name_size = 5,
) +
  ggtitle("FAPs") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

# --- Find the 4-way overlap ---
ggvenn_FAPS_4way <- Reduce(intersect, list(F_runners, F_TFEB, M_runners, M_TFEB))

# --- Convert to data frame for Excel export ---
ggvenn_FAPS_4way_df <- data.frame(GOBP_Term = ggvenn_FAPS_4way)

# --- Save to Excel file ---
write_xlsx(ggvenn_FAPS_4way_df, "ggvenn_FAPS_4way.xlsx")





#Myo IIA 
# Read data from Excel
F_runners <- read_excel("IIA GO terms for analysis.xlsx", sheet = 5)
F_runners <- as.character(F_runners[[1]])

F_TFEB <- read_excel("IIA GO terms for analysis.xlsx", sheet = 5)
F_TFEB <- as.character(F_TFEB[[2]])

M_runners <- read_excel("IIA GO terms for analysis.xlsx", sheet = 9)
M_runners <- as.character(M_runners[[1]])

M_TFEB <- read_excel("IIA GO terms for analysis.xlsx", sheet = 9)
M_TFEB <- as.character(M_TFEB[[2]])

F_runners <- F_runners[grepl("^GOBP_", F_runners)]
F_TFEB   <- F_TFEB[grepl("^GOBP_", F_TFEB)]
M_runners <- M_runners[grepl("^GOBP_", M_runners)]
M_TFEB   <- M_TFEB[grepl("^GOBP_", M_TFEB)]                   

# Combine into a named list for ggvenn
go_lists <- list(
  "Male Runner"   = M_runners,
  "Female TFEB"   = F_TFEB,
  "Male TFEB"     = M_TFEB,
  "Female Runner" = F_runners
)

# Define colors for petals (in order)
petal_colors <- c("#87CEFA","#993399","#1E90FF","#FF69B4")
# (Female runner = bright pink, Female TFEB = light pink,
#  Male runner = bright blue, Male TFEB = light blue)

# Plot the Venn diagram
ggvenn(
  go_lists,
  fill_color = petal_colors,
  fill_alpha = 0.6,
  stroke_color = "white",
  stroke_size = 1.2,
  text_size = 3.7,
  set_name_size = 5,
) +
  ggtitle("Myonuclei IIA") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

# --- Find the 4-way overlap ---
ggvenn_MyoIIA_4way <- Reduce(intersect, list(F_runners, F_TFEB, M_runners, M_TFEB))

# --- Convert to data frame for Excel export ---
ggvenn_MyoIIA_4way_df <- data.frame(GOBP_Term = ggvenn_MyoIIA_4way)

# --- Save to Excel file ---
write_xlsx(ggvenn_MyoIIA_4way_df, "ggvenn_MyoIIA_4way.xlsx")

