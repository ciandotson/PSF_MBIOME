#### Individualized Histograms ####
## Make each categorical variable into a factor ##
# Soil #
decompose_ps(soil.ps, 'soil')
soil$met$Plants <- factor(soil$met$Plant, levels = c("S. helvola", "C. fasciculata", "D. canadense", "A. bracteata", "T. repens", "M. truncatula"))
soil$met$Comps <- factor(soil$met$Compartment, levels = c("Source Community", "Rhizosphere", "Nodule"))
soil$met$Soils <- factor(soil$met$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

# Root #
decompose_ps(root.ps, 'root')
root$met$Plants <- factor(root$met$Plant.Species, levels = c("S. helvola", "C. fasciculata", "D. canadense", "A. bracteata", "T. repens", "M. truncatula"))
root$met$Comps <- factor(root$met$Compartment, levels = c("Root Endosphere", "Nodule"))
root$met$roots <- factor(root$met$Soil.Origin, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

# Save the new metadata table to the original phyloseq object #
sample_data(soil.ps) <- soil$met
sample_data(root.ps) <- root$met

# Create phyloseq objects that have samples corresponding to their compartment #
bulk.ps <- subset_samples(soil.ps, Compartment == "Source Community")
rhiz.ps <- subset_samples(soil.ps, Compartment == "Rhizosphere")
root.ps <- subset_samples(root.ps, Compartment == "Root Endosphere")

### Create phyloseq objects breaking down each each compartment by plant species ### 
## Source Community ##
# Fuzzy Bean #
fb_bulk.ps <- subset_samples(bulk.ps, Plant == "S. helvola")
fb_bulk.ps <- merge_phyloseq(fb_bulk.ps, subset_samples(bulk.ps, Plant == "C. fasciculata" & Soil_Treatment == "Non-PSF Soil"))
fb_bulk.ps <- subset_taxa(fb_bulk.ps, taxa_sums(fb_bulk.ps) > 0)

# Chamaecrista #
cc_bulk.ps <- subset_samples(bulk.ps, Plant == "C. fasciculata")
cc_bulk.ps <- merge_phyloseq(cc_bulk.ps, subset_samples(bulk.ps, Soil_Treatment == "Common Soil"))
cc_bulk.ps <- subset_taxa(cc_bulk.ps, taxa_sums(cc_bulk.ps) > 0)

# Desmodium #
ds_bulk.ps <- subset_samples(bulk.ps, Plant == "D. canadense")
ds_bulk.ps <- merge_phyloseq(ds_bulk.ps, subset_samples(bulk.ps, Soil_Treatment == "Common Soil"))
ds_bulk.ps <- subset_taxa(ds_bulk.ps, taxa_sums(ds_bulk.ps) > 0)

# Hog Peanut #
hp_bulk.ps <- subset_samples(bulk.ps, Plant == "A. bracteata")
hp_bulk.ps <- merge_phyloseq(hp_bulk.ps, subset_samples(bulk.ps, Soil_Treatment == "Common Soil"))
hp_bulk.ps <- merge_phyloseq(hp_bulk.ps, subset_samples(bulk.ps, Plant == "D. canadense" & Soil_Treatment == "Non-PSF Soil"))
hp_bulk.ps <- subset_taxa(hp_bulk.ps, taxa_sums(hp_bulk.ps) > 0)

# Clover #
cl_bulk.ps <- subset_samples(bulk.ps, Plant == "T. repens")
cl_bulk.ps <- merge_phyloseq(cl_bulk.ps, subset_samples(bulk.ps, Soil_Treatment == "Common Soil"))
cl_bulk.ps <- subset_taxa(cl_bulk.ps, taxa_sums(cl_bulk.ps) > 0)

# Medicago #
md_bulk.ps <- subset_samples(bulk.ps, Plant == "M. truncatula")
md_bulk.ps <- merge_phyloseq(md_bulk.ps, subset_samples(bulk.ps, Soil_Treatment == "Common Soil"))
md_bulk.ps <- merge_phyloseq(md_bulk.ps, subset_samples(bulk.ps, Plant == "T. repens" & Soil_Treatment == "Non-PSF Soil"))
md_bulk.ps <- subset_taxa(md_bulk.ps, taxa_sums(md_bulk.ps) > 0)

## Rhizosphere ##
# Fuzzy Bean #
fb_rhiz.ps <- subset_samples(rhiz.ps, Plant == "S. helvola")
fb_rhiz.ps <- subset_taxa(fb_rhiz.ps, taxa_sums(fb_rhiz.ps) > 0)

# Chamaecrista #
cc_rhiz.ps <- subset_samples(rhiz.ps, Plant == "C. fasciculata")
cc_rhiz.ps <- subset_taxa(cc_rhiz.ps, taxa_sums(cc_rhiz.ps) > 0)

# Desmodium #
ds_rhiz.ps <- subset_samples(rhiz.ps, Plant == "D. canadense")
ds_rhiz.ps <- subset_taxa(ds_rhiz.ps, taxa_sums(ds_rhiz.ps) > 0)

# Hog Peanut #
hp_rhiz.ps <- subset_samples(rhiz.ps, Plant == "A. bracteata")
hp_rhiz.ps <- subset_taxa(hp_rhiz.ps, taxa_sums(hp_rhiz.ps) > 0)

# Clover #
cl_rhiz.ps <- subset_samples(rhiz.ps, Plant == "T. repens")
cl_rhiz.ps <- subset_taxa(cl_rhiz.ps, taxa_sums(cl_rhiz.ps) > 0)

# Medicago #
md_rhiz.ps <- subset_samples(rhiz.ps, Plant == "M. truncatula")
md_rhiz.ps <- subset_taxa(md_rhiz.ps, taxa_sums(md_rhiz.ps) > 0)

## Root Endosphere ##
# Fuzzy Bean #
fb_root.ps <- subset_samples(root.ps, Plant.Species == "S. helvola")
fb_root.ps <- subset_taxa(fb_root.ps, taxa_sums(fb_root.ps) > 0)

# Chamaecrista #
cc_root.ps <- subset_samples(root.ps, Plant.Species == "C. fasciculata")
cc_root.ps <- subset_taxa(cc_root.ps, taxa_sums(cc_root.ps) > 0)

# Desmodium #
ds_root.ps <- subset_samples(root.ps, Plant.Species == "D. canadense")
ds_root.ps <- subset_taxa(ds_root.ps, taxa_sums(ds_root.ps) > 0)

# Hog Peanut #
hp_root.ps <- subset_samples(root.ps, Plant.Species == "A. bracteata")
hp_root.ps <- subset_taxa(hp_root.ps, taxa_sums(hp_root.ps) > 0)

# Clover #
cl_root.ps <- subset_samples(root.ps, Plant.Species == "T. repens")
cl_root.ps <- subset_taxa(cl_root.ps, taxa_sums(cl_root.ps) > 0)

# Medicago #
md_root.ps <- subset_samples(root.ps, Plant.Species == "M. truncatula")
md_root.ps <- subset_taxa(md_root.ps, taxa_sums(md_root.ps) > 0)

### Separate each plant species such that it has a phyloseq object just containing PSF Samples ###
## Fuzzy Bean ##
# Source Community #
fb_pbulk.ps <- subset_samples(fb_bulk.ps, Soil_Treatment == "PSF Soil")
fb_pbulk.ps <- subset_taxa(fb_pbulk.ps, taxa_sums(fb_pbulk.ps) > 0)

# Rhizosphere #
fb_prhiz.ps <- subset_samples(fb_rhiz.ps, Soil_Treatment == "PSF Soil")
fb_prhiz.ps<- subset_taxa(fb_prhiz.ps, taxa_sums(fb_prhiz.ps) > 0)

# Root Endosphere #
fb_proot.ps <- subset_samples(fb_root.ps, Soil.Origin == "PSF Soil")
fb_proot.ps <- subset_taxa(fb_proot.ps, taxa_sums(fb_proot.ps) > 0)

## Chamaecrista ##
# Source Community #
cc_pbulk.ps <- subset_samples(cc_bulk.ps, Soil_Treatment == "PSF Soil")
cc_pbulk.ps <- subset_taxa(cc_pbulk.ps, taxa_sums(cc_pbulk.ps) > 0)

# Rhizosphere #
cc_prhiz.ps <- subset_samples(cc_rhiz.ps, Soil_Treatment == "PSF Soil")
cc_prhiz.ps<- subset_taxa(cc_prhiz.ps, taxa_sums(cc_prhiz.ps) > 0)

# Root Endosphere #
cc_proot.ps <- subset_samples(cc_root.ps, Soil.Origin == "PSF Soil")
cc_proot.ps <- subset_taxa(cc_proot.ps, taxa_sums(cc_proot.ps) > 0)

## Desmodium ##
# Source Community #
ds_pbulk.ps <- subset_samples(ds_bulk.ps, Soil_Treatment == "PSF Soil")
ds_pbulk.ps <- subset_taxa(ds_pbulk.ps, taxa_sums(ds_pbulk.ps) > 0)

# Rhizosphere #
ds_prhiz.ps <- subset_samples(ds_rhiz.ps, Soil_Treatment == "PSF Soil")
ds_prhiz.ps<- subset_taxa(ds_prhiz.ps, taxa_sums(ds_prhiz.ps) > 0)

# Root Endosphere #
ds_proot.ps <- subset_samples(ds_root.ps, Soil.Origin == "PSF Soil")
ds_proot.ps <- subset_taxa(ds_proot.ps, taxa_sums(ds_proot.ps) > 0)

## Hog Peanut ##
# Source Community #
hp_pbulk.ps <- subset_samples(hp_bulk.ps, Soil_Treatment == "PSF Soil")
hp_pbulk.ps <- subset_taxa(hp_pbulk.ps, taxa_sums(hp_pbulk.ps) > 0)

# Rhizosphere #
hp_prhiz.ps <- subset_samples(hp_rhiz.ps, Soil_Treatment == "PSF Soil")
hp_prhiz.ps<- subset_taxa(hp_prhiz.ps, taxa_sums(hp_prhiz.ps) > 0)

# Root Endosphere #
hp_proot.ps <- subset_samples(hp_root.ps, Soil.Origin == "PSF Soil")
hp_proot.ps <- subset_taxa(hp_proot.ps, taxa_sums(hp_proot.ps) > 0)

## Clover ##
# Source Community #
cl_pbulk.ps <- subset_samples(cl_bulk.ps, Soil_Treatment == "PSF Soil")
cl_pbulk.ps <- subset_taxa(cl_pbulk.ps, taxa_sums(cl_pbulk.ps) > 0)

# Rhizosphere #
cl_prhiz.ps <- subset_samples(cl_rhiz.ps, Soil_Treatment == "PSF Soil")
cl_prhiz.ps<- subset_taxa(cl_prhiz.ps, taxa_sums(cl_prhiz.ps) > 0)

# Root Endosphere #
cl_proot.ps <- subset_samples(cl_root.ps, Soil.Origin == "PSF Soil")
cl_proot.ps <- subset_taxa(cl_proot.ps, taxa_sums(cl_proot.ps) > 0)

## Medicago ##
# Source Community #
md_pbulk.ps <- subset_samples(md_bulk.ps, Soil_Treatment == "PSF Soil")
md_pbulk.ps <- subset_taxa(md_pbulk.ps, taxa_sums(md_pbulk.ps) > 0)

# Rhizosphere #
md_prhiz.ps <- subset_samples(md_rhiz.ps, Soil_Treatment == "PSF Soil")
md_prhiz.ps<- subset_taxa(md_prhiz.ps, taxa_sums(md_prhiz.ps) > 0)

# Root Endosphere #
md_proot.ps <- subset_samples(md_root.ps, Soil.Origin == "PSF Soil")
md_proot.ps <- subset_taxa(md_proot.ps, taxa_sums(md_proot.ps) > 0)

### Make the histograms for each Comp x Plant, PSF Interaction ###
## Create a palette for each genus represented in the phyloseq object ##
set.seed(248)
soil.colr <- createPalette(ntaxa(soil.ps), c("#ff0000", "#00ff00", "#0000ff"))
soil.colr <- as.data.frame(soil.colr)
rownames(soil.colr) <- sort(taxa_names(soil.ps))
soil.colr[ntaxa(soil.ps) + 1,] <- "#D4D4D4"
rownames(soil.colr)[ntaxa(soil.ps) + 1] <- "Other"

set.seed(496)
root.colr <- createPalette(ntaxa(root.ps), c("#ff0000", "#00ff00", "#0000ff"))
root.colr <- as.data.frame(root.colr)
rownames(root.colr) <- sort(taxa_names(root.ps))
root.colr[ntaxa(root.ps) + 1,] <- "#D4D4D4"
rownames(root.colr)[ntaxa(root.ps) + 1] <- "Other"

## Fuzzy Bean Source Community ##
# Create a phyloseq with the top 19 asvs and the rest converted to "Other
fb_pbulk_top.ps <- aggregate_top_taxa2(fb_pbulk.ps, 19, "ASV")
fb_pbulk_top.name <- names(sort(taxa_sums(fb_pbulk_top.ps), decreasing = TRUE))
fb_pbulk_top.colr <- soil.colr[fb_pbulk_top.name,]

# Create a data.frame with all the relevant information for plotting #
fb_pbulk_top.df <- psmelt(fb_pbulk_top.ps)
fb_pbulk_top.df$ASVs <- factor(fb_pbulk_top.df$ASV, levels = fb_pbulk_top.name)
fb_pbulk_top.df$Group <- factor(fb_pbulk_top.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

# Plot the figure #
fb_pbulk_top.plot <- ggplot(fb_pbulk_top.df, aes(x = Comps, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = fb_pbulk_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('S. helvola')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90))

## Fuzzy Bean Rhizosphere ##
# Create a phyloseq with the top 19 asvs and the rest converted to "Other
fb_prhiz_top.ps <- aggregate_top_taxa2(fb_prhiz.ps, 19, "ASV")
fb_prhiz_top.name <- names(sort(taxa_sums(fb_prhiz_top.ps), decreasing = TRUE))
fb_prhiz_top.colr <- soil.colr[fb_prhiz_top.name,]

# Create a data.frame with all the relevant information for plotting #
fb_prhiz_top.df <- psmelt(fb_prhiz_top.ps)
fb_prhiz_top.df$ASVs <- factor(fb_prhiz_top.df$ASV, levels = fb_prhiz_top.name)
fb_prhiz_top.df$Group <- factor(fb_prhiz_top.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

# Plot the figure #
fb_prhiz_top.plot <- ggplot(fb_prhiz_top.df, aes(x = Comps, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = fb_prhiz_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('S. helvola')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90))

## Fuzzy Bean Root Endosphere ##
# Create a phyloseq with the top 19 asvs and the rest converted to "Other
fb_proot_top.ps <- aggregate_top_taxa2(fb_proot.ps, 19, "ASV")
fb_proot_top.name <- names(sort(taxa_sums(fb_proot_top.ps), decreasing = TRUE))
fb_proot_top.colr <- root.colr[fb_proot_top.name,]

# Create a data.frame with all the relevant information for plotting #
fb_proot_top.df <- psmelt(fb_proot_top.ps)
fb_proot_top.df$ASVs <- factor(fb_proot_top.df$ASV, levels = fb_proot_top.name)
fb_proot_top.df$Group <- factor(fb_proot_top.df$Soil.Origin, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

# Plot the figure #
fb_proot_top.plot <- ggplot(fb_proot_top.df, aes(x = Comps, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = fb_proot_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('S. helvola')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90))

## Chamaecrista Source Community ##
# Create a phyloseq with the top 19 asvs and the rest converted to "Other
cc_pbulk_top.ps <- aggregate_top_taxa2(cc_pbulk.ps, 19, "ASV")
cc_pbulk_top.name <- names(sort(taxa_sums(cc_pbulk_top.ps), decreasing = TRUE))
cc_pbulk_top.colr <- soil.colr[cc_pbulk_top.name,]

# Create a data.frame with all the relevant information for plotting #
cc_pbulk_top.df <- psmelt(cc_pbulk_top.ps)
cc_pbulk_top.df$ASVs <- factor(cc_pbulk_top.df$ASV, levels = cc_pbulk_top.name)
cc_pbulk_top.df$Group <- factor(cc_pbulk_top.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

# Plot the figure #
cc_pbulk_top.plot <- ggplot(cc_pbulk_top.df, aes(x = Comps, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = cc_pbulk_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('C. fasciculata')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90))

## Chamaecrista Rhizosphere ##
# Create a phyloseq with the top 19 asvs and the rest converted to "Other
cc_prhiz_top.ps <- aggregate_top_taxa2(cc_prhiz.ps, 19, "ASV")
cc_prhiz_top.name <- names(sort(taxa_sums(cc_prhiz_top.ps), decreasing = TRUE))
cc_prhiz_top.name <- c("Other", cc_prhiz_top.name)
cc_prhiz_top.colr <- soil.colr[unique(cc_prhiz_top.name),]

# Create a data.frame with all the relevant information for plotting #
cc_prhiz_top.df <- psmelt(cc_prhiz_top.ps)
cc_prhiz_top.df$ASVs <- factor(cc_prhiz_top.df$ASV, levels = unique(cc_prhiz_top.name))
cc_prhiz_top.df$Group <- factor(cc_prhiz_top.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

# Plot the figure #
cc_prhiz_top.plot <- ggplot(cc_prhiz_top.df, aes(x = Comps, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = cc_prhiz_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('C. fasciulata')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90))

## Fuzzy Bean Root Endosphere ##
# Create a phyloseq with the top 19 asvs and the rest converted to "Other
cc_proot_top.ps <- aggregate_top_taxa2(cc_proot.ps, 19, "ASV")
cc_proot_top.name <- names(sort(taxa_sums(cc_proot_top.ps), decreasing = TRUE))
cc_proot_top.name <- c("Other", cc_proot_top.name)
cc_proot_top.colr <- root.colr[unique(cc_proot_top.name),]

# Create a data.frame with all the relevant information for plotting #
cc_proot_top.df <- psmelt(cc_proot_top.ps)
cc_proot_top.df$ASVs <- factor(cc_proot_top.df$ASV, levels = unique(cc_proot_top.name))
cc_proot_top.df$Group <- factor(cc_proot_top.df$Soil.Origin, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

# Plot the figure #
cc_proot_top.plot <- ggplot(cc_proot_top.df, aes(x = Comps, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = cc_proot_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('C. fasciculata')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90))

## Desmodium Source Community ##
# Create a phyloseq with the top 19 asvs and the rest converted to "Other
ds_pbulk_top.ps <- aggregate_top_taxa2(ds_pbulk.ps, 19, "ASV")
ds_pbulk_top.name <- names(sort(taxa_sums(ds_pbulk_top.ps), decreasing = TRUE))
ds_pbulk_top.colr <- soil.colr[ds_pbulk_top.name,]

# Create a data.frame with all the relevant information for plotting #
ds_pbulk_top.df <- psmelt(ds_pbulk_top.ps)
ds_pbulk_top.df$ASVs <- factor(ds_pbulk_top.df$ASV, levels = ds_pbulk_top.name)
ds_pbulk_top.df$Group <- factor(ds_pbulk_top.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

# Plot the figure #
ds_pbulk_top.plot <- ggplot(ds_pbulk_top.df, aes(x = Comps, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = ds_pbulk_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('D. canadense')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90))

## Desmodium Rhizosphere ##
# Create a phyloseq with the top 19 asvs and the rest converted to "Other
ds_prhiz_top.ps <- aggregate_top_taxa2(ds_prhiz.ps, 19, "ASV")
ds_prhiz_top.name <- names(sort(taxa_sums(ds_prhiz_top.ps), decreasing = TRUE))
ds_prhiz_top.name <- c("Other", ds_prhiz_top.name)
ds_prhiz_top.colr <- soil.colr[unique(ds_prhiz_top.name),]

# Create a data.frame with all the relevant information for plotting #
ds_prhiz_top.df <- psmelt(ds_prhiz_top.ps)
ds_prhiz_top.df$ASVs <- factor(ds_prhiz_top.df$ASV, levels = unique(ds_prhiz_top.name))
ds_prhiz_top.df$Group <- factor(ds_prhiz_top.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

# Plot the figure #
ds_prhiz_top.plot <- ggplot(ds_prhiz_top.df, aes(x = Comps, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = ds_prhiz_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('D. canadense')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90))

## Desmodium Root Endosphere ##
# Create a phyloseq with the top 19 asvs and the rest converted to "Other
ds_proot_top.ps <- aggregate_top_taxa2(ds_proot.ps, 19, "ASV")
ds_proot_top.name <- names(sort(taxa_sums(ds_proot_top.ps), decreasing = TRUE))
ds_proot_top.name <- c("Other", ds_proot_top.name)
ds_proot_top.colr <- root.colr[unique(ds_proot_top.name),]

# Create a data.frame with all the relevant information for plotting #
ds_proot_top.df <- psmelt(ds_proot_top.ps)
ds_proot_top.df$ASVs <- factor(ds_proot_top.df$ASV, levels = unique(ds_proot_top.name))
ds_proot_top.df$Group <- factor(ds_proot_top.df$Soil.Origin, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

# Plot the figure #
ds_proot_top.plot <- ggplot(ds_proot_top.df, aes(x = Comps, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = ds_proot_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('D. canadense')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90))

## Hog Peanut Source Community ##
# Create a phyloseq with the top 19 asvs and the rest converted to "Other
hp_pbulk_top.ps <- aggregate_top_taxa2(hp_pbulk.ps, 19, "ASV")
hp_pbulk_top.name <- names(sort(taxa_sums(hp_pbulk_top.ps), decreasing = TRUE))
hp_pbulk_top.colr <- soil.colr[hp_pbulk_top.name,]

# Create a data.frame with all the relevant information for plotting #
hp_pbulk_top.df <- psmelt(hp_pbulk_top.ps)
hp_pbulk_top.df$ASVs <- factor(hp_pbulk_top.df$ASV, levels = hp_pbulk_top.name)
hp_pbulk_top.df$Group <- factor(hp_pbulk_top.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

# Plot the figure #
hp_pbulk_top.plot <- ggplot(hp_pbulk_top.df, aes(x = Comps, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = hp_pbulk_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('A. bracteta')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90))

## Hog Peanut Rhizosphere ##
# Create a phyloseq with the top 19 asvs and the rest converted to "Other
hp_prhiz_top.ps <- aggregate_top_taxa2(hp_prhiz.ps, 19, "ASV")
hp_prhiz_top.name <- names(sort(taxa_sums(hp_prhiz_top.ps), decreasing = TRUE))
hp_prhiz_top.name <- c("Other", hp_prhiz_top.name)
hp_prhiz_top.colr <- soil.colr[unique(hp_prhiz_top.name),]

# Create a data.frame with all the relevant information for plotting #
hp_prhiz_top.df <- psmelt(hp_prhiz_top.ps)
hp_prhiz_top.df$ASVs <- factor(hp_prhiz_top.df$ASV, levels = unique(hp_prhiz_top.name))
hp_prhiz_top.df$Group <- factor(hp_prhiz_top.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

# Plot the figure #
hp_prhiz_top.plot <- ggplot(hp_prhiz_top.df, aes(x = Comps, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = hp_prhiz_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('A. bracteta')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90))

## Hog Peanut Root Endosphere ##
# Create a phyloseq with the top 19 asvs and the rest converted to "Other
hp_proot_top.ps <- aggregate_top_taxa2(hp_proot.ps, 19, "ASV")
hp_proot_top.name <- names(sort(taxa_sums(hp_proot_top.ps), decreasing = TRUE))
hp_proot_top.name <- c("Other", hp_proot_top.name)
hp_proot_top.colr <- root.colr[unique(hp_proot_top.name),]

# Create a data.frame with all the relevant information for plotting #
hp_proot_top.df <- psmelt(hp_proot_top.ps)
hp_proot_top.df$ASVs <- factor(hp_proot_top.df$ASV, levels = unique(hp_proot_top.name))
hp_proot_top.df$Group <- factor(hp_proot_top.df$Soil.Origin, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

# Plot the figure #
hp_proot_top.plot <- ggplot(hp_proot_top.df, aes(x = Comps, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = hp_proot_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('A. bracteata')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90))

## Clover Source Community ##
# Create a phyloseq with the top 19 asvs and the rest converted to "Other
cl_pbulk_top.ps <- aggregate_top_taxa2(cl_pbulk.ps, 19, "ASV")
cl_pbulk_top.name <- names(sort(taxa_sums(cl_pbulk_top.ps), decreasing = TRUE))
cl_pbulk_top.colr <- soil.colr[cl_pbulk_top.name,]

# Create a data.frame with all the relevant information for plotting #
cl_pbulk_top.df <- psmelt(cl_pbulk_top.ps)
cl_pbulk_top.df$ASVs <- factor(cl_pbulk_top.df$ASV, levels = cl_pbulk_top.name)
cl_pbulk_top.df$Group <- factor(cl_pbulk_top.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

# Plot the figure #
cl_pbulk_top.plot <- ggplot(cl_pbulk_top.df, aes(x = Comps, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = cl_pbulk_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('T. repens')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90))

## Clover Rhizosphere ##
# Create a phyloseq with the top 19 asvs and the rest converted to "Other
cl_prhiz_top.ps <- aggregate_top_taxa2(cl_prhiz.ps, 19, "ASV")
cl_prhiz_top.name <- names(sort(taxa_sums(cl_prhiz_top.ps), decreasing = TRUE))
cl_prhiz_top.name <- c("Other", cl_prhiz_top.name)
cl_prhiz_top.colr <- soil.colr[unique(cl_prhiz_top.name),]

# Create a data.frame with all the relevant information for plotting #
cl_prhiz_top.df <- psmelt(cl_prhiz_top.ps)
cl_prhiz_top.df$ASVs <- factor(cl_prhiz_top.df$ASV, levels = unique(cl_prhiz_top.name))
cl_prhiz_top.df$Group <- factor(cl_prhiz_top.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

# Plot the figure #
cl_prhiz_top.plot <- ggplot(cl_prhiz_top.df, aes(x = Comps, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = cl_prhiz_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('T. repens')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90))

## Clover Root Endosphere ##
# Create a phyloseq with the top 19 asvs and the rest converted to "Other
cl_proot_top.ps <- aggregate_top_taxa2(cl_proot.ps, 19, "ASV")
cl_proot_top.name <- names(sort(taxa_sums(cl_proot_top.ps), decreasing = TRUE))
cl_proot_top.name <- c("Other", cl_proot_top.name)
cl_proot_top.colr <- root.colr[unique(cl_proot_top.name),]

# Create a data.frame with all the relevant information for plotting #
cl_proot_top.df <- psmelt(cl_proot_top.ps)
cl_proot_top.df$ASVs <- factor(cl_proot_top.df$ASV, levels = unique(cl_proot_top.name))
cl_proot_top.df$Group <- factor(cl_proot_top.df$Soil.Origin, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

# Plot the figure #
cl_proot_top.plot <- ggplot(cl_proot_top.df, aes(x = Comps, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = cl_proot_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('T. repens')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90))

## Medicago Source Community ##
# Create a phyloseq with the top 19 asvs and the rest converted to "Other
md_pbulk_top.ps <- aggregate_top_taxa2(md_pbulk.ps, 19, "ASV")
md_pbulk_top.name <- names(sort(taxa_sums(md_pbulk_top.ps), decreasing = TRUE))
md_pbulk_top.colr <- soil.colr[md_pbulk_top.name,]

# Create a data.frame with all the relevant information for plotting #
md_pbulk_top.df <- psmelt(md_pbulk_top.ps)
md_pbulk_top.df$ASVs <- factor(md_pbulk_top.df$ASV, levels = md_pbulk_top.name)
md_pbulk_top.df$Group <- factor(md_pbulk_top.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

# Plot the figure #
md_pbulk_top.plot <- ggplot(md_pbulk_top.df, aes(x = Comps, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = md_pbulk_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('M. truncatula')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90))

## Medicago Rhizosphere ##
# Create a phyloseq with the top 19 asvs and the rest converted to "Other
md_prhiz_top.ps <- aggregate_top_taxa2(md_prhiz.ps, 19, "ASV")
md_prhiz_top.name <- names(sort(taxa_sums(md_prhiz_top.ps), decreasing = TRUE))
md_prhiz_top.name <- c("Other", md_prhiz_top.name)
md_prhiz_top.colr <- soil.colr[unique(md_prhiz_top.name),]

# Create a data.frame with all the relevant information for plotting #
md_prhiz_top.df <- psmelt(md_prhiz_top.ps)
md_prhiz_top.df$ASVs <- factor(md_prhiz_top.df$ASV, levels = unique(md_prhiz_top.name))
md_prhiz_top.df$Group <- factor(md_prhiz_top.df$Soil_Treatment, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

# Plot the figure #
md_prhiz_top.plot <- ggplot(md_prhiz_top.df, aes(x = Comps, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = md_prhiz_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('M. truncatula')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90))

## Medicago Root Endosphere ##
# Create a phyloseq with the top 19 asvs and the rest converted to "Other
md_proot_top.ps <- aggregate_top_taxa2(md_proot.ps, 19, "ASV")
md_proot_top.name <- names(sort(taxa_sums(md_proot_top.ps), decreasing = TRUE))
md_proot_top.name <- c("Other", md_proot_top.name)
md_proot_top.colr <- root.colr[unique(md_proot_top.name),]

# Create a data.frame with all the relevant information for plotting #
md_proot_top.df <- psmelt(md_proot_top.ps)
md_proot_top.df$ASVs <- factor(md_proot_top.df$ASV, levels = unique(md_proot_top.name))
md_proot_top.df$Group <- factor(md_proot_top.df$Soil.Origin, levels = c("PSF Soil", "Non-PSF Soil", "Common Soil"))

# Plot the figure #
md_proot_top.plot <- ggplot(md_proot_top.df, aes(x = Comps, y = Abundance, fill = ASVs)) +
  geom_bar(stat='identity', position = 'fill') +
  xlab('') +
  ylab('Relative Abundance') +
  scale_fill_manual(values = md_proot_top.colr) +
  scale_y_continuous(sec.axis = dup_axis(name = expression(italic('M. truncatula')))) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Liberation Sans", size = 18),
        axis.title = element_text(size = 22, family = "Liberation Sans"),
        strip.text = element_text(size = 18, family = "Liberation Sans"),
        legend.text = element_text(size = 18, family = "Liberation Sans"),
        legend.position = "right",
        legend.box.just = "bottom",
        legend.justification = "left",
        legend.margin = margin(t = 40, b = 20),
        plot.margin = margin(10, 40, 10, 10),
        legend.title = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.y.right = element_text(size = 18, face = 'bold', angle = -90))
