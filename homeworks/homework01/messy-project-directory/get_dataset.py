###############################################################################
#                     __                      _      ___ ____                 #
#                    / _| ___  _ __ _ __ ___ (_) ___|_ _|  _ \                #
#                   | |_ / _ \| '__| '_ ` _ \| |/ __|| || | | |               #
#                   |  _| (_) | |  | | | | | | | (__ | || |_| |               #
#                   |_|  \___/|_|  |_| |_| |_|_|\___|___|____/                #
#                                                                             #
#                             Downloading the dataset                         #
#                                                                             #
###############################################################################
"""Description:
Functions for downloading and handeling the datasets.
"""

# Packages
###############################################################################

# Standard library imports
import logging

# FormicID imports
from data_loader.data_input import remove_reproductives
from data_loader.data_input import split_in_directory
from data_scraper.scrape import get_dataset
from data_scraper.scrape import stitch_maker
from utils.load_config import process_config
from utils.utils import get_args

# Main code
###############################################################################
try:
    args = get_args()
    config = process_config(args.config)
except:
    logging.error("Missing or invalid arguments.")
    exit(0)

# Creating a dataset
###############################################################################
get_dataset(
    input="testgenusspecies.csv",
    n_jsonfiles=100,
    config=config,
    shottypes="dhp",
    quality="medium",
    update=False,
    offset_set=0,
    limit_set=99999,
    multi_only=True,
)

# stitch_maker(config=config)

# Split into training, validation and test
###############################################################################
# split_in_directory(config=config, bad="data/badspecimens_multi.csv")

# Remove reproductives from dataset
###############################################################################
# remove_reproductives(
#     csv="data/reproductives.csv",
#     dataset="top97species_Qmed_def_clean_wtest",
#     config=config,
# )
