#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# Script for plotting the reuslts of the paper using ABCD data release 5.1 entitled   #
# 'Breastfeeding duration is positively related to cortical thickness and cognition'  #
# by Jonatan Ottino-González et al. (2024). Last update: Aug 14th 2024                #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

#1# Load libraries 

pacman::p_load(tidyverse, data.table, ggeffects, gridExtra, ggseg, ggsegDesterieux)

#2# Set paths

drivedir <- '/Users/jongonzalez/Library/CloudStorage/GoogleDrive-jonatanottino@gmail.com/.shortcut-targets-by-id/1-H3wnyMZIfSE18QmegCkCkmyQpjcPqtZ/Adise_lab/'

main <- fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/CT_results.csv'))
 
(a <- main %>%
    ggplot() +
    geom_brain(atlas = ggsegDesterieux::desterieux, 
               colour = 'black',
               position = position_brain(position = 'horizontal'),
               size = 1.5,
               aes(fill = beta)) +
    theme_brain2() +
    theme(plot.background = element_rect(fill = 'white', colour = NA)) +
    labs(fill = expression(beta)) +
    scale_fill_gradient2(low = 'blue',
                         mid = 'white',
                         high = 'red'))

(b <- main %>%
    mutate(beta = ifelse(fdr < 0.05, beta, NA)) %>%
    ggplot() +
    geom_brain(atlas = ggsegDesterieux::desterieux,    # or dk
               colour = 'black',
               position = position_brain(position = 'horizontal'),
               size = 1.5,
               aes(fill = beta)) +
    theme_brain2() +
    theme(plot.background = element_rect(fill = 'white', colour = NA)) +
    labs(fill = expression(beta)) +
    scale_fill_gradient2(low = 'blue',
                         mid = 'white',
                         high = 'red'))

myel <- fread(paste0(drivedir, 'Lab_member_files/Jonatan_files/abcd_bf_results/Myel_results.csv'))

(c <- myel %>%
    ggplot() +
    geom_brain(atlas = ggsegDesterieux::desterieux, 
               colour = 'black',
               position = position_brain(position = 'horizontal'),
               size = 1.5,
               aes(fill = beta)) +
    theme_brain2() +
    theme(plot.background = element_rect(fill = 'white', colour = NA)) +
    labs(fill = expression(beta)) +
    scale_fill_gradient2(low = 'blue',
                         mid = 'white',
                         high = 'red'))

(d <- myel %>%
    mutate(beta = ifelse(`p-value` < 0.05, beta, NA)) %>%
    ggplot() +
    geom_brain(atlas = ggsegDesterieux::desterieux,    # or dk
               colour = 'black',
               position = position_brain(position = 'horizontal'),
               size = 1.5,
               aes(fill = beta)) +
    theme_brain2() +
    theme(plot.background = element_rect(fill = 'white', colour = NA)) +
    labs(fill = expression(beta)) +
    scale_fill_gradient2(low = 'blue',
                         mid = 'white',
                         high = 'red'))
