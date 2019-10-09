###############################################################################
#                     __                      _      ___ ____                 #
#                    / _| ___  _ __ _ __ ___ (_) ___|_ _|  _ \                #
#                   | |_ / _ \| '__| '_ ` _ \| |/ __|| || | | |               #
#                   |  _| (_) | |  | | | | | | | (__ | || |_| |               #
#                   |_|  \___/|_|  |_| |_| |_|_|\___|___|____/                #
#                                                                             #
#                                Get species list                             #
#                                                                             #
###############################################################################
"""Description:
Script to create a csv file with the species that are imaged the most. As of now, the list is not entirely correct, as some specimens have extra images, like close-ups and these also count as images.
"""

# Packages
###############################################################################

# FormicID imports
from AntWeb.AW2_to_json import most_imaged_species_to_csv

most_imaged_species_to_csv("Get_list_test.csv", min_images=68)
