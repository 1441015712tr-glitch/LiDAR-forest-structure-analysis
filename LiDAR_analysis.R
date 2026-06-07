# ================================================
# LiDAR Forest Structure Analysis
# Tianran Yin | AUT | 2026
# ================================================

# ---- 1. Load packages ----
library(lidR)
library(sf)
library(dplyr)
library(ggplot2)
library(patchwork)
library(terra)

# ---- 2. Field measurement data ----
field_data <- data.frame(
  tree_id = c("3058","10818","10850","3112","10801","3111","3104",
              "3383","3375","739",
              "1","2","3","4","5","6","7","8"),
  plot = c(rep("TeMuri", 7),
           rep("Pourewa_young", 3),
           rep("Pourewa_mature", 8)),
  circumference_cm = c(63, 70, 80, 60, 40, 65, 80,
                       50, 35, 50,
                       70, 40, 40, 100, 40, 60, 60, 50),
  measure_position_cm = c(10, 18, 12, 15, 7, 14, 20,
                          10, 7, 6,
                          32, 13, 11, 25, 22, 10, 15, 29)
)
field_data$diameter_cm <- field_data$circumference_cm / pi

# ---- 3. UAV point cloud processing ----

# 3.1 Te Muri
ctg_temuri <- readLAScatalog("C:/AUT/毕业论文/数据/rawdata/L2-TeMuri")

temuri_corners <- data.frame(
  x = c(1752750.932, 1752742.377, 1752740.688, 1752751.369),
  y = c(5957028.764, 5957040.396, 5957030.032, 5957039.863)
)
temuri_poly <- st_as_sf(temuri_corners, coords = c("x", "y"), crs = 2193) |>
  st_combine() |> st_convex_hull()

las_temuri <- clip_roi(ctg_temuri, temuri_poly)
las_temuri <- classify_ground(las_temuri, algorithm = csf())
las_temuri_norm <- normalize_height(las_temuri, knnidw())
chm_temuri <- rasterize_canopy(las_temuri_norm, res = 0.5, algorithm = p2r())
ttops_temuri <- locate_trees(las_temuri_norm, lmf(ws = 2))
las_temuri_seg <- segment_trees(las_temuri_norm, dalponte2016(chm_temuri, ttops_temuri))
trees_temuri <- las_temuri_seg@data %>%
  filter(!is.na(treeID)) %>%
  group_by(treeID) %>%
  summarise(height_m = round(max(Z), 2))

# 3.2 Pourewa Plot 6
ctg_pourewa <- readLAScatalog("C:/AUT/毕业论文/数据/rawdata/L2-Pourewa")

pourewa6_corners <- data.frame(
  x = c(1762522.263, 1762547.825, 1762549.061, 1762522.347),
  y = c(5918671.488, 5918674.725, 5918649.737, 5918648.537)
)
pourewa6_poly <- st_as_sf(pourewa6_corners, coords = c("x", "y"), crs = 2193) |>
  st_combine() |> st_convex_hull()

las_pourewa6 <- clip_roi(ctg_pourewa, pourewa6_poly)
las_pourewa6 <- classify_ground(las_pourewa6, algorithm = csf())
las_pourewa6_norm <- normalize_height(las_pourewa6, knnidw())
chm_pourewa6 <- rasterize_canopy(las_pourewa6_norm, res = 0.5, algorithm = p2r())
ttops_pourewa6 <- locate_trees(las_pourewa6_norm, lmf(ws = 2))
las_pourewa6_seg <- segment_trees(las_pourewa6_norm, dalponte2016(chm_pourewa6, ttops_pourewa6))
trees_pourewa6 <- las_pourewa6_seg@data %>%
  filter(!is.na(treeID)) %>%
  group_by(treeID) %>%
  summarise(height_m = round(max(Z), 2))

# 3.3 Pourewa Kepa
kepa_corners <- data.frame(
  x = c(1762870.118, 1762891.079, 1762888.726, 1762866.930),
  y = c(5918738.584, 5918738.795, 5918715.434, 5918713.343)
)
kepa_poly <- st_as_sf(kepa_corners, coords = c("x", "y"), crs = 2193) |>
  st_combine() |> st_convex_hull()

las_kepa <- clip_roi(ctg_pourewa, kepa_poly)
las_kepa <- classify_ground(las_kepa, algorithm = csf())
las_kepa_norm <- normalize_height(las_kepa, knnidw())
chm_kepa <- rasterize_canopy(las_kepa_norm, res = 0.5, algorithm = p2r())
ttops_kepa <- locate_trees(las_kepa_norm, lmf(ws = 2))
las_kepa_seg <- segment_trees(las_kepa_norm, dalponte2016(chm_kepa, ttops_kepa))
trees_kepa <- las_kepa_seg@data %>%
  filter(!is.na(treeID)) %>%
  group_by(treeID) %>%
  summarise(height_m = round(max(Z), 2))

# ---- 4. Handheld LiDAR processing — Pourewa Plot 6 ----

# 4.1 Read and subsample
ctg_h_pourewa6 <- readLAScatalog("C:/AUT/毕业论文/数据/rawdata/LiGrip-Pourewa-Plots6")
opt_filter(ctg_h_pourewa6) <- "-keep_random_fraction 0.02"
las_h_pourewa6 <- readLAS(ctg_h_pourewa6)

# 4.2 Coordinate registration to NZTM2000
local_pts <- matrix(c(
  9.503,  14.866,
  32.237, 12.501,
  2.370,  -7.192
), ncol=2, byrow=TRUE)

nztm_pts <- matrix(c(
  1762522, 5918671,
  1762548, 5918675,
  1762522, 5918649
), ncol=2, byrow=TRUE)

A <- cbind(local_pts, 1)
Tx <- lm.fit(A, nztm_pts[,1])$coefficients
Ty <- lm.fit(A, nztm_pts[,2])$coefficients

X_new <- Tx[1]*las_h_pourewa6$X + Tx[2]*las_h_pourewa6$Y + Tx[3]
Y_new <- Ty[1]*las_h_pourewa6$X + Ty[2]*las_h_pourewa6$Y + Ty[3]
las_h_pourewa6@header$`X offset` <- floor(min(X_new))
las_h_pourewa6@header$`Y offset` <- floor(min(Y_new))
las_h_pourewa6$X <- X_new
las_h_pourewa6$Y <- Y_new
lidR::projection(las_h_pourewa6) <- 2193

# 4.3 Clip to plot boundary
las_h_p6_nztm <- clip_roi(las_h_pourewa6, pourewa6_poly)

# 4.4 Ground classification and normalisation
las_h_p6_nztm <- classify_ground(las_h_p6_nztm, algorithm = csf())
las_h_p6_nztm_norm <- normalize_height(las_h_p6_nztm, tin())

# 4.5 CHM and tree segmentation
chm_h_p6_nztm <- rasterize_canopy(las_h_p6_nztm_norm, res = 0.5, algorithm = p2r())
ttops_h_p6_nztm <- locate_trees(las_h_p6_nztm_norm, lmf(ws = 2))
las_h_p6_nztm_seg <- segment_trees(las_h_p6_nztm_norm,
                                   dalponte2016(chm_h_p6_nztm, ttops_h_p6_nztm))
trees_h_p6_nztm <- las_h_p6_nztm_seg@data %>%
  filter(!is.na(treeID)) %>%
  group_by(treeID) %>%
  summarise(height_m = round(max(Z), 2))

# ---- 5. Figure generation ----

# Figure 1: UAV CHM comparison
df_temuri <- as.data.frame(chm_temuri, xy=TRUE); names(df_temuri)[3] <- "height"
df_p6 <- as.data.frame(chm_pourewa6, xy=TRUE); names(df_p6)[3] <- "height"
df_kepa <- as.data.frame(chm_kepa, xy=TRUE); names(df_kepa)[3] <- "height"

p1 <- ggplot(df_temuri, aes(x=x, y=y, fill=height)) + geom_raster() +
  scale_fill_gradientn(colours=terrain.colors(50), name="Height (m)",
                       limits=c(0,10), na.value="white") +
  labs(title="Te Muri (UAV)") + theme_bw() +
  theme(plot.title=element_text(hjust=0.5, face="bold"), axis.title=element_blank())

p2 <- ggplot(df_p6, aes(x=x, y=y, fill=height)) + geom_raster() +
  scale_fill_gradientn(colours=terrain.colors(50), name="Height (m)",
                       limits=c(0,10), na.value="white") +
  labs(title="Pourewa Plot 6 (UAV)") + theme_bw() +
  theme(plot.title=element_text(hjust=0.5, face="bold"), axis.title=element_blank())

p3 <- ggplot(df_kepa, aes(x=x, y=y, fill=height)) + geom_raster() +
  scale_fill_gradientn(colours=terrain.colors(50), name="Height (m)",
                       limits=c(0,10), na.value="white") +
  labs(title="Pourewa Kepa (UAV)") + theme_bw() +
  theme(plot.title=element_text(hjust=0.5, face="bold"), axis.title=element_blank())

ggsave("Figure1_CHM_UAV.png", p1+p2+p3, width=12, height=5, dpi=300)

# Figure 2: UAV vs handheld CHM
df_h_p6 <- as.data.frame(chm_h_p6_nztm, xy=TRUE); names(df_h_p6)[3] <- "height"

p4 <- ggplot(df_p6, aes(x=x, y=y, fill=height)) + geom_raster() +
  scale_fill_gradientn(colours=terrain.colors(50), name="Height (m)",
                       limits=c(0,6), na.value="white") +
  labs(title="Pourewa Plot 6 - UAV LiDAR") + theme_bw() +
  theme(plot.title=element_text(hjust=0.5, face="bold"), axis.title=element_blank())

p5 <- ggplot(df_h_p6, aes(x=x, y=y, fill=height)) + geom_raster() +
  scale_fill_gradientn(colours=terrain.colors(50), name="Height (m)",
                       limits=c(0,6), na.value="white") +
  labs(title="Pourewa Plot 6 - Handheld LiDAR") + theme_bw() +
  theme(plot.title=element_text(hjust=0.5, face="bold"), axis.title=element_blank())

ggsave("Figure2_CHM_comparison.png", p4+p5, width=10, height=5, dpi=300)

# Figure 3: Tree height boxplot
trees_temuri$plot <- "Te Muri"; trees_temuri$method <- "UAV"
trees_pourewa6$plot <- "Pourewa Plot 6"; trees_pourewa6$method <- "UAV"
trees_kepa$plot <- "Pourewa Kepa"; trees_kepa$method <- "UAV"
trees_h_p6_nztm$plot <- "Pourewa Plot 6"; trees_h_p6_nztm$method <- "Handheld"

all_trees <- bind_rows(trees_temuri, trees_pourewa6, trees_kepa, trees_h_p6_nztm)

fig3 <- ggplot(all_trees, aes(x=plot, y=height_m, fill=method)) +
  geom_boxplot(position=position_dodge(0.8), width=0.6) +
  scale_fill_manual(values=c("UAV"="#2196F3", "Handheld"="#FF9800")) +
  labs(title="Tree Height Distribution by Plot and Method",
       x="Plot", y="Tree Height (m)", fill="Method") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, face="bold"),
        axis.text.x=element_text(angle=15, hjust=1))

ggsave("Figure3_height_comparison.png", fig3, width=10, height=6, dpi=300)

# Figure 4: Tree detection count
tree_count <- data.frame(
  plot = c("Te Muri", "Pourewa Plot 6", "Pourewa Plot 6", "Pourewa Kepa"),
  method = c("UAV", "UAV", "Handheld", "UAV"),
  count = c(13, 116, 113, 69)
)

fig4 <- ggplot(tree_count, aes(x=plot, y=count, fill=method)) +
  geom_bar(stat="identity", position=position_dodge(0.8), width=0.6) +
  scale_fill_manual(values=c("UAV"="#2196F3", "Handheld"="#FF9800")) +
  labs(title="Tree Detection Count by Plot and Method",
       x="Plot", y="Number of Trees Detected", fill="Method") +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5, face="bold"),
        axis.text.x=element_text(angle=15, hjust=1))

ggsave("Figure4_tree_count.png", fig4, width=10, height=6, dpi=300)