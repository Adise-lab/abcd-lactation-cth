#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# Script for plotting the reuslts of the paper using ABCD data release 5.1 entitled   #
# 'Breastfeeding duration is positively related to cortical thickness and cognition'  #
# by Jonatan Ottino-González et al. (2024). Last update: Nov 26th 2024                #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

#1# Load libraries 

pacman::p_load(tidyverse, data.table, ggeffects, cowplot, ggseg, ggsegDesterieux, grid)

#2# Set paths

drivedir <- '/Users/jongonzalez/Library/CloudStorage/GoogleDrive-jonatanottino@gmail.com/.shortcut-targets-by-id/1-H3wnyMZIfSE18QmegCkCkmyQpjcPqtZ/Adise_lab/'

## Set lobules

frontal_lobe <- c(
  "lh_G_and_S_frontomargin", "lh_G_front_inf-Opercular", "lh_G_front_inf-Orbital", 
  "lh_G_front_inf-Triangul", "lh_G_front_middle", "lh_G_front_sup", 
  "lh_G_precentral", "lh_S_front_inf", "lh_S_front_middle", "lh_S_front_sup", 
  "rh_G_and_S_frontomargin", "rh_G_front_inf-Opercular", "rh_G_front_inf-Orbital", 
  "rh_G_front_inf-Triangul", "rh_G_front_middle", "rh_G_front_sup", 
  "rh_G_precentral", "rh_S_front_inf", "rh_S_front_middle", "rh_S_front_sup",
  "lh_G_and_S_subcentral", "lh_G_and_S_transv_frontopol", "lh_G_orbital", 
  "lh_S_central", "lh_S_orbital_lateral", "lh_S_orbital_med-olfact", 
  "lh_S_orbital-H_Shaped", "lh_S_precentral-inf-part", "lh_S_precentral-sup-part", 
  "rh_G_and_S_subcentral", "rh_G_and_S_transv_frontopol", "rh_G_orbital", 
  "rh_S_central", "rh_S_orbital_lateral", "rh_S_orbital_med-olfact", 
  "rh_S_orbital-H_Shaped", "rh_S_precentral-inf-part", "rh_S_precentral-sup-part",
  "lh_S_suborbital", "rh_S_suborbital", 
  "lh_Lat_Fis-ant-Horizont", "lh_Lat_Fis-ant-Vertical", 
  "rh_Lat_Fis-ant-Horizont", "rh_Lat_Fis-ant-Vertical"
)

parietal_lobe <- c(
  "lh_G_and_S_paracentral", "lh_G_pariet_inf-Angular", "lh_G_pariet_inf-Supramar", 
  "lh_G_parietal_sup", "lh_G_postcentral", "lh_G_precuneus", 
  "lh_S_interm_prim-Jensen", "lh_S_intrapariet_and_P_trans", 
  "rh_G_and_S_paracentral", "rh_G_pariet_inf-Angular", "rh_G_pariet_inf-Supramar", 
  "rh_G_parietal_sup", "rh_G_postcentral", "rh_G_precuneus", 
  "rh_S_interm_prim-Jensen", "rh_S_intrapariet_and_P_trans",
  "lh_S_postcentral", "rh_S_postcentral", "lh_S_subparietal", "rh_S_subparietal",
  "rh_S_parieto_occipital", "lh_S_parieto_occipital"
)

temporal_lobe <- c(
  "lh_G_oc-temp_lat-fusifor", "lh_G_oc-temp_med-Lingual", "lh_G_oc-temp_med-Parahip", 
  "lh_G_temp_sup-G_T_transv", "lh_G_temp_sup-Lateral", "lh_G_temp_sup-Plan_polar", 
  "lh_G_temp_sup-Plan_tempo", "lh_G_temporal_inf", "lh_G_temporal_middle", 
  "lh_S_oc-temp_lat", "lh_S_oc-temp_med_and_Lingual", "lh_S_temporal_inf", 
  "lh_S_temporal_sup", "lh_S_temporal_transverse", 
  "rh_G_oc-temp_lat-fusifor", "rh_G_oc-temp_med-Lingual", "rh_G_oc-temp_med-Parahip", 
  "rh_G_temp_sup-G_T_transv", "rh_G_temp_sup-Lateral", "rh_G_temp_sup-Plan_polar", 
  "rh_G_temp_sup-Plan_tempo", "rh_G_temporal_inf", "rh_G_temporal_middle", 
  "rh_S_oc-temp_lat", "rh_S_oc-temp_med_and_Lingual", "rh_S_temporal_inf", 
  "rh_S_temporal_sup", "rh_S_temporal_transverse", 
  "lh_Pole_temporal", "rh_Pole_temporal",  "lh_S_collat_transv_ant", "lh_S_collat_transv_post", "rh_S_collat_transv_ant", 
  "rh_S_collat_transv_post", "lh_Lat_Fis-post",
  "rh_Lat_Fis-post")

occipital_lobe <- c(
  "lh_G_and_S_occipital_inf", "lh_G_occipital_middle", "lh_G_occipital_sup", 
  "lh_G_oc-temp_lat-fusifor", "lh_G_cuneus", "lh_S_oc_middle_and_Lunatus", 
  "lh_S_oc_sup_and_transversal", "lh_S_occipital_ant", 
  "rh_G_and_S_occipital_inf", "rh_G_occipital_middle", "rh_G_occipital_sup", 
  "rh_G_oc-temp_lat-fusifor", "rh_G_cuneus", "rh_S_oc_middle_and_Lunatus", 
  "rh_S_oc_sup_and_transversal", "rh_S_occipital_ant", 
  "lh_Pole_occipital", "rh_Pole_occipital", "lh_S_calcarine", "rh_S_calcarine"
  )

limbic_lobe <- c(
  "lh_G_and_S_cingul-Ant", "lh_G_and_S_cingul-Mid-Ant", "lh_G_and_S_cingul-Mid-Post", 
  "lh_G_cingul-Post-dorsal", "lh_G_cingul-Post-ventral", "lh_G_rectus", 
  "lh_G_subcallosal", "lh_S_cingul-Marginalis",
  "rh_G_and_S_cingul-Ant", "rh_G_and_S_cingul-Mid-Ant", "rh_G_and_S_cingul-Mid-Post", 
  "rh_G_cingul-Post-dorsal", "rh_G_cingul-Post-ventral", "rh_G_rectus", 
  "rh_G_subcallosal", "rh_S_cingul-Marginalis", "lh_S_pericallosal", "rh_S_pericallosal"
)

insula <- c(
  "lh_G_Ins_lg_and_S_cent_ins", "lh_G_insular_short", "lh_S_circular_insula_ant", 
  "lh_S_circular_insula_inf", "lh_S_circular_insula_sup", 
  "rh_G_Ins_lg_and_S_cent_ins", "rh_G_insular_short", "rh_S_circular_insula_ant", 
  "rh_S_circular_insula_inf", "rh_S_circular_insula_sup"
)

replacements <- c(
  "superior" = 'Sup',
  "inferior" = "Inf",
  "anterior" = "Ant",
  "posterior" = "Post",
  "middle" = "Mid",
  "medial" = "Med",
  "lateral" = "Lat",
  "vertical" = "Ver",
  "horizontal" = 'Hor.'
)

# Function to shorten text and format accordingly

replace_case_insensitive <- function(x, replacements) {
  sapply(x, function(element) {
    for (pattern in names(replacements)) {
      replacement <- replacements[pattern]
      
      # Replace at the start of the string (capitalize first letter)
      if (str_starts(element, regex(pattern, ignore_case = TRUE))) {
        element <- str_replace_all(
          element, 
          regex(paste0("^", pattern, "(?=\\b|\\s)"), ignore_case = TRUE), 
          function(match) {
            if (grepl("^.+-$", substring(element, 1, nchar(match) + 1))) {
              # Do not append a period if followed by a hyphen
              str_to_sentence(replacement)
            } else {
              # Append period normally
              paste0(str_to_sentence(replacement), ".")
            }
          }
        )
      }
      
      # Replace elsewhere in the string (lowercase)
      element <- str_replace_all(
        element, 
        regex(paste0("\\b", pattern, "\\b"), ignore_case = TRUE), 
        function(match) {
          if (grepl("-$", substring(element, regexpr(match, element), regexpr(match, element) + nchar(match)))) {
            # Do not append period if followed by a hyphen
            str_to_lower(replacement)
          } else {
            # Append period normally
            paste0(str_to_lower(replacement), ".")
          }
        }
      )
    }
    element
  }, USE.NAMES = FALSE)
}

# ••••••••••••••••••••••••••••••••••••••••••••••••••• #
# ~~~~~~~~~~~~~~~~~ F I G U R E  1 ~~~~~~~~~~~~~~~~~~ # CORTICAL THICKNESS
# ••••••••••••••••••••••••••••••••••••••••••••••••••• #

main.cth <- fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/CT_results.csv'))

new_rows <- data.frame(label = c("lh_Medial_wall", "rh_Medial_wall"))

main.cth <- bind_rows(main.cth, new_rows) %>% 
  mutate(n.order = seq(nrow(.)),
         lobule = case_when(label %in% frontal_lobe ~ 'Frontal',
                                   label %in% parietal_lobe ~ 'Parietal',
                                   label %in% occipital_lobe ~ 'Occipital',
                                   label %in% temporal_lobe ~ 'Temporal',
                                   label %in% limbic_lobe ~ 'Limbic',
                                   label %in% insula ~ 'Insula',
                                   TRUE ~ NA),
         col = case_when(lobule == 'Frontal' ~ 'gold',
                         lobule == 'Parietal' ~ 'skyblue1',
                         lobule == 'Temporal' ~ 'plum1',
                         lobule == 'Occipital' ~ 'olivedrab2',
                         lobule == 'Limbic' ~ 'darkorange2',
                         lobule == 'Insula' ~ 'slateblue',
                         TRUE ~ 'lightgray'))

(any.negative <- subset(main.cth$ROI, main.cth$beta < 0 & main.cth$fdr < 0.05))

(names <- fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/CT_results_fdr.csv'))$ROI)

# Forest plot

(lh.forest <- main.cth %>% 
  filter(!label %in% 'lh_Medial_wall' &
           str_starts(label, 'lh_') &
           ROI %in% names) %>% 
    full_join(., readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
                 rename(ROI = roi) %>% 
                mutate(ROI = paste0(ROI, '_cth'),
                       name = str_remove_all(name, '\\(|\\)|hemisphere|left') %>% 
                         str_squish() %>% 
                         str_to_sentence()) %>%
                 select(ROI, name) %>% 
                mutate(name = replace_case_insensitive(name, replacements)), by = 'ROI') %>% 
    drop_na() %>% 
    arrange(., lobule))

(a <-  lh.forest %>% 
    mutate(lobule_num = as.numeric(factor(lobule))) %>% 
    ggplot(., aes(x = beta,
                  y = reorder(name, lobule_num),
                  xmin = ci.min, 
                  xmax = ci.max)) +
    # Add colored rectangles using annotate
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = .35, ymax = 1.5, fill = 'gold', alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 5.5, fill = 'olivedrab2', alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 5.5, ymax = 11.5, fill = 'skyblue', alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 11.5, ymax = 13.5, fill = 'plum1', alpha = 0.2) +
    # Add horizontal dashed line
    geom_vline(xintercept = 0, 
               linetype = 3) +
    geom_pointrange(shape = 20) +
    labs(x = 'Standardized coefficient') +
    theme_classic() +
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 11,
                                     colour = 'black'),
          axis.title.x = element_text(size = 13,
                                      colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          axis.title.y = element_blank(),
          plot.margin = margin(0, 0, 0, .5, 'cm')) +
    scale_x_continuous(limits = c(0, 0.12),
                       breaks = c(0, 0.025, 0.05, 0.075, 0.10),
                       labels = c('0', '0.025', '0.05', '0.075', '0.10'))) 

(rh.forest <- main.cth %>% 
  filter(!label %in% 'rh_Medial_wall' &
           str_starts(label, 'rh_') &
           ROI %in% names) %>% 
    full_join(., readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
                rename(ROI = roi) %>% 
                mutate(ROI = paste0(ROI, '_cth'),
                       name = str_remove_all(name, '\\(|\\)|hemisphere|right') %>% 
                         str_squish() %>% 
                         str_to_sentence()) %>%
                select(ROI, name) %>% 
                mutate(name = replace_case_insensitive(name, replacements)), by = 'ROI') %>% 
    drop_na() %>% 
    arrange(., lobule))

(b <-  rh.forest %>% 
    mutate(lobule_num = as.numeric(factor(lobule))) %>% 
    ggplot(., aes(x = beta,
                  y = reorder(name, lobule_num),
                  xmin = ci.min, 
                  xmax = ci.max)) +
    # Add colored rectangles using annotate
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = .35, ymax = 1.5, fill = 'gold', alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2.5, fill = 'slateblue', alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2.5, ymax = 3.5, fill = 'darkorange', alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3.5, ymax = 8.5, fill = 'olivedrab2', alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 8.5, ymax = 12.5, fill = 'skyblue1', alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 12.5, ymax = 15.5, fill = 'plum1', alpha = 0.2) +
    geom_vline(xintercept = 0, linetype = 3) +
    geom_pointrange(shape = 20) +
    labs(x = 'Standardized coefficient') +
    theme_classic() +
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 11,
                                     colour = 'black'),
          axis.title.x = element_text(size = 13,
                                      colour = 'black'),
          axis.text.y = element_text(size = 10,
                                     colour = 'black'),
          axis.title.y = element_blank()) +
    scale_x_continuous(limits = c(0, 0.11),
                       breaks = c(0, 0.025, 0.05, 0.075, 0.10),
                       labels = c('0', '0.025', '0.05', '0.075', '0.10')))

# Brain

(c <- main.cth %>%
    filter(str_starts(label, 'lh_')) %>% 
    mutate(col2 = ifelse(fdr < 0.05 | is.na(fdr), col, 'lightgray')) %>% 
    ggplot() +
    geom_brain(atlas = ggsegDesterieux::desterieux,   
               hemi = 'left',
               colour = 'black',
               size = .25,
               aes(fill = col2)) +
    labs(title = 'Left hemisphere') +
    scale_fill_identity() +
    theme_brain2() +
    theme(plot.background = element_rect(fill = 'white', colour = NA),
          axis.title = element_blank(),
          plot.title = element_text(size = 15, colour = 'black', hjust = 0.5, family = 'Helvetica', face = 'bold')) +
    # ignore warning about expression 
    annotate("text", x = 40, y = 0, label = expression(italic('Lateral')), size = 3.5, colour = "black") +
    annotate("text", x = 650, y = 0, label = expression(italic('Medial')), size = 3.5, colour = "black"))

(d <- main.cth %>%
    filter(str_starts(label, 'rh_')) %>% 
    mutate(col2 = ifelse(fdr < 0.05 | is.na(fdr), col, 'lightgray')) %>% 
    ggplot() +
    geom_brain(atlas = ggsegDesterieux::desterieux,   
               hemi = 'right',
               position = position_brain(.~hemi + side),
               colour = 'black',
               size = .25,
               aes(fill = col2)) +
    labs(title = 'Right hemisphere') +
    scale_fill_identity() +
    theme_brain2() +
    theme(plot.background = element_rect(fill = 'white', colour = NA),
          axis.title = element_blank(),
          plot.title = element_text(size = 15, colour = 'black', hjust = 0.5, family = 'Helvetica', face = 'bold')) +
    # ignore warning about expression
    annotate("text", x = 50, y = 0, label = expression(italic('Lateral')), size = 3.5, colour = "black") +
    annotate("text", x = 625, y = 0, label = expression(italic('Medial')), size = 3.5, colour = "black"))

library(patchwork)

(combined <- (c / a))
(combined2 <- (d / b))

(p1 <- plot_grid(combined,
                 combined2,
                 ncol = 2,
                 rel_widths = c(.95, 1),
                 label_x = c(.095, 0),
                 label_y = .95,
                 label_size = 16) +
    theme(plot.background = element_rect(fill = 'white', color = NA)))

color_mapping <- c(
  'Frontal' = 'gold',
  'Occipital' = 'olivedrab2',
  'Parietal' = 'skyblue1',
  'Limbic' = 'darkorange2',
  'Insula' = 'slateblue2',
  'Temporal' = 'plum1')

custom_legend_order <- c('Temporal', 'Parietal', 'Occipital', 'Limbic', 'Insula', 'Frontal')

legend_plot <- ggplot(main.cth %>% filter(!is.na(lobule)), aes(x = lobule, fill = lobule)) + 
  geom_bar(show.legend = TRUE) +  # Ensure legend is included
  scale_fill_manual(values = color_mapping, breaks = custom_legend_order) +  # Use unique colors from your dataset
  theme(legend.position = "top", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 15), 
        legend.key.size = unit(1, "lines")) +  # Adjust legend key size
  guides(fill = guide_legend(nrow = 1))

legend <- ggpubr::get_legend(legend_plot)

# Function to center the legend to the final plot

(legend_centered <- plot_grid(
  NULL,          # Empty space to push the legend to the center
  legend,        # Your legend
  NULL,          # Empty space on the right
  ncol = 3,      # Split into 3 columns
  rel_widths = c(1.7, 1, 1)  # Equal widths for centering
))

(fig1 <- plot_grid(legend_centered,  
                   p1,
                   nrow = 2,
                   labels = 'Significant associations between breastfeeding duration and cortical thickness',
                   label_x = -.14,
                   label_y = 1.5,
                   label_size = 20,
                   rel_heights = c(0.1, 1)) +
  theme(plot.background = element_rect(fill = 'white', color = NA),
        plot.margin = margin(1.5, 1, 0, 0, 'cm')))

ggsave(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/figures/fig1.png'), dpi = 300, width = 5, height = 3, scale = 3, plot = fig1)

# ••••••••••••••••••••••••••••••••••••••••••••••••••• #
# ~~~~~~~~~~~~~~~~~ F I G U R E  2 ~~~~~~~~~~~~~~~~~~ #   SURFACE AREA
# ••••••••••••••••••••••••••••••••••••••••••••••••••• #

main.area <- fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/Area_results.csv')) %>% 
  mutate(n.order = seq(nrow(.)))

new_rows <- data.frame(label = c("lh_Medial_wall", "rh_Medial_wall"))

main.area <- bind_rows(main.area, new_rows) %>% 
    mutate(lobule = case_when(
      label %in% frontal_lobe ~ 'Frontal',
      label %in% parietal_lobe ~ 'Parietal',
      label %in% occipital_lobe ~ 'Occipital',
      label %in% temporal_lobe ~ 'Temporal',
      label %in% limbic_lobe ~ 'Limbic',
      label %in% insula ~ 'Insula',
      TRUE ~ NA_character_),
      col = case_when(lobule == 'Frontal' ~ 'gold',
                      lobule == 'Parietal' ~ 'skyblue1',
                      lobule == 'Temporal' ~ 'plum1',
                      lobule == 'Occipital' ~ 'olivedrab2',
                      lobule == 'Limbic' ~ 'darkorange2',
                      lobule == 'Insula' ~ 'slateblue',
                      TRUE ~ 'lightgray'))

(names <- fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/Area_results_fdr.csv')) %>% 
    filter(fdr < 0.05) %>% 
    pull(ROI))

# Forest plot

(lh.forest <- main.area %>% 
    filter(!label %in% 'lh_Medial_wall' &
             str_starts(label, 'lh_') &
             ROI %in% names) %>% 
    full_join(., readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
                rename(ROI = roi) %>% 
                mutate(ROI = paste0(ROI, '_area'),
                       name = str_remove_all(name, '\\(|\\)|hemisphere|left') %>% 
                         str_squish() %>% 
                         str_to_sentence()) %>%
                select(ROI, name) %>% 
                mutate(name = replace_case_insensitive(name, replacements)), by = 'ROI') %>% 
    drop_na() %>% 
    arrange(., lobule))

(a1 <-  lh.forest %>% 
    mutate(lobule_num = as.numeric(factor(lobule))) %>% 
    ggplot(., aes(x = beta,
                  y = reorder(name, lobule_num),
                  xmin = ci.min, 
                  xmax = ci.max)) +
    # Add colored rectangles using annotate
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.35, ymax = 6.5, fill = 'gold', alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 6.5, ymax = 9.5, fill = 'slateblue', alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 9.5, ymax = 14.5, fill = 'darkorange2', alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 14.5, ymax = 17.5, fill = 'olivedrab2', alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 17.5, ymax = 19.5, fill = 'skyblue', alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 19.5, ymax = 26.5, fill = 'plum1', alpha = 0.2) +
    # Add horizontal dashed line
    geom_vline(xintercept = 0, 
               linetype = 3) +
    geom_pointrange(shape = 20) +
    labs(x = 'Standardized coefficient') +
    theme_classic() +
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 11,
                                     colour = 'black'),
          axis.title.x = element_text(size = 13,
                                      colour = 'black'),
          axis.text.y = element_text(margin = margin(r = 2, unit = 'mm'),
                                     size = 10,
                                     colour = 'black'),
          axis.title.y = element_blank()) +
    scale_x_continuous(limits = c(0, 0.10),
                       breaks = c(0, 0.025, 0.05, 0.075),
                       labels = c('0', '0.025', '0.05', '0.075'))) 

(rh.forest <- main.area %>% 
    filter(!label %in% 'rh_Medial_wall' &
             str_starts(label, 'rh_') &
             ROI %in% names) %>% 
    full_join(., readxl::read_xlsx(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_ggseg.xlsx')) %>% 
                rename(ROI = roi) %>% 
                mutate(ROI = paste0(ROI, '_area'),
                       name = str_remove_all(name, '\\(|\\)|hemisphere|right') %>% 
                         str_squish() %>% 
                         str_to_sentence()) %>%
                select(ROI, name) %>% 
                mutate(name = replace_case_insensitive(name, replacements)), by = 'ROI') %>% 
    drop_na() %>% 
    arrange(., lobule))

(b1 <-  rh.forest %>% 
    mutate(lobule_num = as.numeric(factor(lobule))) %>% 
    ggplot(., aes(x = beta,
                  y = reorder(name, lobule_num),
                  xmin = ci.min, 
                  xmax = ci.max)) +
    # Add colored rectangles using annotate
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = .35, ymax = 7.5, fill = 'gold', alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 7.5, ymax = 9.5, fill = 'darkorange2', alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 9.5, ymax = 13.5, fill = 'olivedrab2', alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 13.5, ymax = 18.5, fill = 'skyblue', alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 18.5, ymax = 25.5, fill = 'plum1', alpha = 0.2) +
    geom_vline(xintercept = 0, linetype = 3) +
    geom_pointrange(shape = 20) +
    labs(x = 'Standardized coefficient') +
    theme_classic() +
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 11,
                                     colour = 'black'),
          axis.title.x = element_text(size = 13,
                                      colour = 'black'),
          axis.text.y = element_text(margin = margin(r = 2, unit = 'mm'),
                                     size = 10,
                                     colour = 'black'),
          axis.title.y = element_blank(),
          plot.margin = margin(0, 1, 0, 1, 'cm')) +
    scale_x_continuous(limits = c(0, 0.10),
                       breaks = c(0, 0.025, 0.05, 0.075),
                       labels = c('0', '0.025', '0.05', '0.075'))) 

# Brain

(c1 <- main.area %>%
    filter(str_starts(label, 'lh_')) %>% 
    mutate(col2 = ifelse(fdr < 0.05 | is.na(fdr), col, 'lightgray')) %>% 
    ggplot() +
    geom_brain(atlas = ggsegDesterieux::desterieux,   
               hemi = 'left',
               colour = 'black',
               size = .25,
               aes(fill = col2)) +
    labs(title = 'Left hemisphere') +
    scale_fill_identity() +
    theme_brain2() +
    theme(plot.background = element_rect(fill = 'white', colour = NA),
          axis.title = element_blank(),
          plot.title = element_text(size = 15, colour = 'black', hjust = 0.5, family = 'Helvetica', face = 'bold')) +
    annotate("text", x = 40, y = 0, label = expression(italic('Lateral')), size = 3.5, colour = "black") +
    annotate("text", x = 650, y = 0, label = expression(italic('Medial')), size = 3.5, colour = "black"))

(d1 <- main.area %>%
    filter(str_starts(label, 'rh_')) %>% 
    mutate(col2 = ifelse(fdr < 0.05 | is.na(fdr), col, 'lightgray')) %>% 
    ggplot() +
    geom_brain(atlas = ggsegDesterieux::desterieux,    # or dk
               hemi = 'right',
               position = position_brain(.~hemi + side),
               colour = 'black',
               size = .25,
               aes(fill = col2)) +
    labs(title = 'Right hemisphere') +
    scale_fill_identity() +
    theme_brain2() +
    theme(plot.background = element_rect(fill = 'white', colour = NA),
          axis.title = element_blank(),
          plot.title = element_text(size = 15, colour = 'black', hjust = 0.5, family = 'Helvetica', face = 'bold')) +
    annotate("text", x = 50, y = 0, label = expression(italic('Lateral')), size = 3.5, colour = "black") +
    annotate("text", x = 625, y = 0, label = expression(italic('Medial')), size = 3.5, colour = "black"))

library(patchwork) 

# merge ---------

(combined <- (c1 / a1))
(combined2 <- (d1 / b1))

(p2 <- plot_grid(combined,
                 combined2,
                 ncol = 2,
                 rel_widths = c(.96, 1),
                 label_x = c(.095, 0),
                 label_y = .95,
                 label_size = 16) +
    theme(plot.background = element_rect(fill = 'white', color = NA)))

(fig2 <- plot_grid(legend_centered,
                 p2,
                 nrow = 2,
                 labels = 'Significant associations between breastfeeding duration and surface area',
                 label_x = -.07,
                 label_y = 1.5,
                 label_size = 20,
                 rel_heights = c(0.1, 1)) +
    theme(plot.background = element_rect(fill = 'white', color = NA),
          plot.margin = margin(1.5, 1, 0, 0, 'cm')))

ggsave(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/figures/fig2.png'), dpi = 300, width = 5, height = 3, scale = 3, plot = fig2)

# ••••••••••••••••••••••••••••••••••••••••••••••••••• #
# ~~~~~~~~~~~~~~~~~ F I G U R E  3 ~~~~~~~~~~~~~~~~~~ #   CORTICAL MYELIN
# ••••••••••••••••••••••••••••••••••••••••••••••••••• #

(myel.cth <- 
  fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/Myel_CT_results.csv')) %>% 
    select(ROI, label, n, ends_with('main')))

(myel.cth2 <- 
    fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/Myel_CT_results.csv')) %>% 
    select(ROI, label, n, ends_with('int')))

(myel.area <- 
  fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/Myel_area_results.csv')) %>% 
    select(ROI, label, n, ends_with('main')))

(myel.area2 <- 
    fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/Myel_area_results.csv')) %>% 
    select(ROI, label, n, ends_with('int')))

(a2 <- main.cth %>%
    mutate(col2 = ifelse(fdr < 0.05 & !is.na(fdr), 'darkorchid', 'lightgray')) %>% 
    ggplot() +
    geom_brain(atlas = ggsegDesterieux::desterieux,   
               # position = position_brain(hemi + side ~.),
               # position = position_brain(hemi ~ side),
               position = position_brain(position = 'horizontal'),
               colour = 'black',
               size = .5,
               aes(fill = col2)) +
    scale_fill_identity() +
    theme_brain2() +
    theme(plot.background = element_rect(fill = 'white', colour = NA)))

(names.cth <- myel.cth %>% filter(`p-value.main` < 0.05) %>% pull(label))

(names.cth2 <- myel.cth2 %>% filter(`p-value.int` < 0.05) %>% pull(label))

(a3 <- main.cth %>%
    mutate(col2 = ifelse(label %in% names.cth, 'steelblue', 'lightgray')) %>% 
    ggplot() +
    geom_brain(atlas = ggsegDesterieux::desterieux,    # or dk
               colour = 'black',
               # position = position_brain(hemi + side ~.),
               # position = position_brain(hemi ~ side),
               position = position_brain(position = 'horizontal'),
               size = .5,
               aes(fill = col2)) +
    scale_fill_identity() +
    theme_brain2() +
    theme(plot.background = element_rect(fill = 'white', colour = NA),
          axis.title = element_blank()) +
    # Add arrow and text annotation
    annotate( "segment",x = 310, y = 10, xend = 280, yend = 35,  arrow = arrow(type = "closed", length = unit(.25, "cm")), colour = "indianred", size = 1))

(a4 <- main.cth %>%
    mutate(col2 = ifelse(label %in% names.cth2, 'cyan3', 'lightgray')) %>% 
    ggplot() +
    geom_brain(atlas = ggsegDesterieux::desterieux,    # or dk
               colour = 'black',
               # position = position_brain(hemi + side ~.),
               # position = position_brain(hemi ~ side),
               position = position_brain(position = 'horizontal'),
               size = .5,
               aes(fill = col2)) +
    scale_fill_identity() +
    theme_brain2() +
    theme(plot.background = element_rect(fill = 'white', colour = NA)))

color_mapping <- c('Cortical thickness (breastfeeding duration)' = 'darkorchid',
                   'Cortical myelin (breastfeeding duration)' = 'steelblue',
                   'Cortical myelin (breastfeeding-by-age)' = 'cyan3')

legend_plot <- ggplot(main.cth %>% 
                        mutate(mod = factor(ifelse(label %in% names.cth, 'Cortical myelin (breastfeeding duration)', 
                                                   ifelse(label %in% names.cth2, 'Cortical myelin (breastfeeding-by-age)', 'Cortical thickness (breastfeeding duration)')), 
                                            levels = c('Cortical thickness (breastfeeding duration)', 'Cortical myelin (breastfeeding duration)', 'Cortical myelin (breastfeeding-by-age)'))), 
                      aes(x = mod, fill = mod)) + 
  geom_bar(show.legend = TRUE) +  # Ensure legend is included
  scale_fill_manual(values = color_mapping) +  # Use unique colors from your dataset
  theme(legend.position = 'right', 
        legend.direction = 'vertical',
        legend.key.spacing.y = unit(.25, 'cm'),
        legend.title = element_blank(), 
        legend.text = element_text(size = 15), 
        legend.key.size = unit(2, "lines")) +  # Adjust legend key size
  guides(fill = guide_legend(nrow = 3))

legend_plot

legend1 <- ggpubr::get_legend(legend_plot)

(combined1 <- plot_grid(plot_grid(a2, 
                                  a3,
                                  a4,
                                  nrow = 3),
                        legend1, 
                        NULL,
                        ncol = 3, 
                        rel_widths = c(.8, 0.3, .04),
                        labels = c('LH', 'RH'),
                        label_size = 18,
                        label_x = c(.04, -1.35),
                        label_y = c(.95, .95)))

(c2 <- main.area %>%
    mutate(col2 = ifelse(fdr < 0.05 & !is.na(fdr), 'forestgreen', 'lightgray')) %>% 
    ggplot() +
    geom_brain(atlas = ggsegDesterieux::desterieux,   
               # position = position_brain(hemi + side ~.),
               # position = position_brain(hemi ~ side),
               position = position_brain(position = 'horizontal'),
               colour = 'black',
               size = .5,
               aes(fill = col2)) +
    scale_fill_identity() +
    theme_brain2() +
    theme(plot.background = element_rect(fill = 'white', colour = NA)))

(names.area <- myel.area %>% filter(`p-value.main` < 0.05 & !str_detect(label, 'polar')) %>% pull(label))      # exclude polar bc main and interaction effect showed only interaction

(names.area2 <- myel.area2 %>% filter(`p-value.int` < 0.05) %>% pull(label))

(c3 <- main.area %>%
    mutate(col2 = ifelse(label %in% names.area, 'steelblue', 'lightgray')) %>% 
    ggplot() +
    geom_brain(atlas = ggsegDesterieux::desterieux,    # or dk
               colour = 'black',
               # position = position_brain(hemi + side ~.),
               # position = position_brain(hemi ~ side),
               position = position_brain(position = 'horizontal'),
               size = .5,
               aes(fill = col2)) +
    scale_fill_identity() +
    theme_brain2() +
    theme(plot.background = element_rect(fill = 'white', colour = NA),
          axis.title = element_blank()) +
    # Add arrow and text annotation
    annotate( "segment",x = 1140, y = 230, xend = 1170, yend = 205,  arrow = arrow(type = "closed", length = unit(.25, "cm")), colour = "indianred", size = 1))

(c4 <- main.area %>%
    mutate(col2 = ifelse(label %in% names.area2, 'cyan3', 'lightgray')) %>% 
    ggplot() +
    geom_brain(atlas = ggsegDesterieux::desterieux,    # or dk
               colour = 'black',
               # position = position_brain(hemi + side ~.),
               # position = position_brain(hemi ~ side),
               position = position_brain(position = 'horizontal'),
               size = .5,
               aes(fill = col2)) +
    scale_fill_identity() +
    theme_brain2() +
    theme(plot.background = element_rect(fill = 'white', colour = NA)))

color_mapping <- c('Surface area (breastfeeding duration)' = 'forestgreen',
                   'Cortical myelin (breastfeeding duration)' = 'steelblue',
                   'Cortical myelin (breastfeeding-by-age)' = 'cyan3')

legend_plot2 <- ggplot(main.area %>% 
                        mutate(mod = factor(ifelse(label %in% names.area, 'Cortical myelin (breastfeeding duration)',
                                                   ifelse(label %in% names.area2, 'Cortical myelin (breastfeeding-by-age)', 'Surface area (breastfeeding duration)')), 
                                                   levels = c('Surface area (breastfeeding duration)', 'Cortical myelin (breastfeeding duration)', 'Cortical myelin (breastfeeding-by-age)'))), 
                      aes(x = mod, fill = mod)) + 
  geom_bar(show.legend = TRUE) +  # Ensure legend is included
  scale_fill_manual(values = color_mapping) +  # Use unique colors from your dataset
  theme(legend.position = 'right', 
        legend.direction = 'vertical',
        legend.key.spacing.y = unit(.25, 'cm'),
        legend.title = element_blank(), 
        legend.text = element_text(size = 15), 
        legend.key.size = unit(2, "lines")) +  # Adjust legend key size
  guides(fill = guide_legend(nrow = 3))

legend_plot2

legend2 <- ggpubr::get_legend(legend_plot2)

(combined2 <- plot_grid(plot_grid(c2, 
                                  c3,
                                  c4,
                                  nrow = 3),
                        legend2, 
                        NULL,
                        ncol = 3,
                        rel_widths = c(.8, .3, .05),
                        labels = c('LH', 'RH'),
                        label_size = 18,
                        label_x = c(.04, -1.35),
                        label_y = c(.95, .95)))

(fig3 <- plot_grid(combined1 + theme(plot.margin = margin(.5, 1, 0, -1, 'cm')),
                   NULL,
                   combined2 + theme(plot.margin = margin(.5, 1, 0, -1, 'cm')),
                   nrow = 3,
                   labels = c('A. Cortical thickness and cortical myelin overlapping regions', 
                              'B. Surface area and cortical myelin overlapping regions'),
                   label_x = c(-0.23, -0.2),
                   label_y = c(1.05, .55),
                   label_size = 20,
                   rel_heights = c(1, .1, 1)) +
    theme(plot.background = element_rect(fill = 'white', color = NA),
          plot.margin = margin(2, 0, 1, 1, 'cm'))) 

ggsave(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/figures/fig3.png'), dpi = 300, width = 6, height = 5, scale = 3, plot = fig3)

# For talk
# 
# (talk1 <- main.cth %>%
#     mutate(col2 = ifelse(label %in% names.cth, 'steelblue', 'lightgray')) %>% 
#     ggplot() +
#     geom_brain(atlas = ggsegDesterieux::desterieux,    # or dk
#                colour = 'black',
#                position = position_brain(hemi ~ side),
#                size = .5,
#                aes(fill = col2)) +
#     scale_fill_identity() +
#     theme_brain2() +
#     theme(plot.background = element_rect(fill = 'white', colour = NA)))
# 
# (talk2 <- main.area %>%
#     mutate(col2 = ifelse(label %in% names.area, 'steelblue', 'lightgray')) %>% 
#     ggplot() +
#     geom_brain(atlas = ggsegDesterieux::desterieux,    # or dk
#                colour = 'black',
#                position = position_brain(hemi ~ side),
#                size = .5,
#                aes(fill = col2)) +
#     scale_fill_identity() +
#     theme_brain2() +
#     theme(plot.background = element_rect(fill = 'white', colour = NA)))

# ••••••••••••••••••••••••••••••••••••••••••••••••••• #
# ~~~~~~~~~~~~~~~~~ F I G U R E  4 ~~~~~~~~~~~~~~~~~~ #   FLUID COGNITION -- THIS HAS TO BE ARRANGED IN BIO RENDER
# ••••••••••••••••••••••••••••••••••••••••••••••••••• #

(cog.cth <- fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/brain_cog_results.csv')) %>% 
    filter(str_ends(ROI, '_cth')) %>% 
    mutate(col = ifelse(`p-value` < 0.05, 'darkorchid', 'lightgray')))
  
(cog.area <- fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/brain_cog_results.csv')) %>% 
    filter(str_ends(ROI, '_area')) %>% 
    mutate(col = ifelse(`p-value` < 0.05, 'darkorange', 'lightgray')))

(cog.myel <- fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/tables/brain_cog_results.csv')) %>% 
    filter(str_ends(ROI, '_myelin')) %>% 
    mutate(col = ifelse(`p-value` < 0.05, 'steelblue', 'lightgray')))

(a3 <- cog.cth %>%
    mutate(col2 = ifelse(`p-value` < 0.05 &
                           beta > 0, 'darkorchid', 'lightgray')) %>% 
    ggplot() +
    geom_brain(atlas = ggsegDesterieux::desterieux, 
               colour = 'black',
               position = position_brain(. ~ hemi / side),
               size = .25,
               aes(fill = col2),
               show.legend = F) +
    theme_brain2() +
    theme(plot.background = element_rect(fill = 'white', colour = NA)) +
    labs(fill = expression(beta)) +
    scale_fill_manual(values = c('darkorchid', 'lightgray'), na.value = "lightgray"))

(b3 <- cog.area %>%
    mutate(col2 = ifelse(`p-value` < 0.05, 'forestgreen', 'lightgray')) %>% 
    ggplot() +
    geom_brain(atlas = ggsegDesterieux::desterieux, 
               colour = 'black',
               position = position_brain(. ~ hemi / side),
               size = .25,
               aes(fill = col2),
               show.legend = F) +
    theme_brain2() +
    theme(plot.background = element_rect(fill = 'white', colour = NA)) +
    labs(fill = expression(beta)) +
    scale_fill_manual(values = c('forestgreen', 'lightgray'), na.value = "lightgray"))

(c3 <- cog.myel %>%
    mutate(col2 = ifelse(`p-value` < 0.05, 'lightgray', 'steelblue')) %>% 
    ggplot() +
    geom_brain(atlas = ggsegDesterieux::desterieux,   
               position = position_brain(. ~ hemi / side),
               colour = 'black',
               size = .25,
               aes(fill = col),
               show.legend = F) +
    theme_brain2() +
    theme(plot.background = element_rect(fill = 'white', colour = NA)) +
    scale_fill_manual(values = c('lightgray', 'steelblue'), na.value = "lightgray"))

ggsave(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/figures/med1_cth.png'), dpi = 300, width = 3, height = 2, scale = 3, plot = a3)

ggsave(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/figures/med2_area.png'), dpi = 300, width = 3, height = 2, scale = 3, plot = b3)

ggsave(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/figures/med3_myelin.png'), dpi = 300, width = 3, height = 2, scale = 3, plot = c3)

