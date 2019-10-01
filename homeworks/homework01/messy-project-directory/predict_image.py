###############################################################################
#                     __                      _      ___ ____                 #
#                    / _| ___  _ __ _ __ ___ (_) ___|_ _|  _ \                #
#                   | |_ / _ \| '__| '_ ` _ \| |/ __|| || | | |               #
#                   |  _| (_) | |  | | | | | | | (__ | || |_| |               #
#                   |_|  \___/|_|  |_| |_| |_|_|\___|___|____/                #
#                                                                             #
#                                 Predict image(s)                            #
#                                                                             #
###############################################################################
"""Description:
Load a model with trained weights and predict an image or a test set.
"""

# Packages
###############################################################################

# Standard library imports
import os

# FormicID imports
from models.models import compile_model
from models.models import load_model
from testers.tester import predict_image
from testers.tester import predictor
from utils.load_config import process_config
from utils.model_utils import weights_load
from utils.utils import get_args

# Parameters and settings
###############################################################################
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
# To disable the tf warning for compiling in SEE4.2
# 0 = all logs, 1 = info, 2 = warnings, 3 = error


# Arguments
###############################################################################
try:
    args = get_args()
    config = process_config(args.config)
except:
    logging.error("Missing or invalid arguments.")
    exit(0)

# Loading model
##############################################################################
model_formicID = load_model(config=config, num_species=97)
model_formicID = compile_model(model=model_formicID, config=config)
model_formicID = weights_load(
    model=model_formicID,
    weights="experiments/T97_CaAll_QuM_ShP_AugM_D05_LR0001_E200_I4_def_clean/checkpoint/weights_76-1.83.hdf5",
)

# predicting
##############################################################################
Y_true, Y_pred, labels, species_dict = predictor(
    model=model_formicID,
    config=config,
    species_json="data/species_dict.json",
    plot=True,
    n_img=10,
    n_cols=3,
)

# print(Y_true)
# print(Y_pred)

predict_image(
    model=model_formicID,
    species_dict=species_dict,
    image="data/statia2015_rmnh/images/profile/3-test/wasmannia_auropunctata/RMNH5084478_prof_4xMontage.jpg",
    url=None,
    show=False,
)
