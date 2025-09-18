#### Individualized Histograms ####
# Make each categorical variable into a factor #
decompose_ps(soil.ps, 'soil')

soil$met$Plants <- factor(soil$met$Plant, levels = c("S. helvola", "C. fasciculata", "D. illinoense", "A. bracteata", "T. repens", ))