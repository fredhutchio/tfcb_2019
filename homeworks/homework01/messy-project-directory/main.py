###############################################################################
#                     __                      _      ___ ____                 #
#                    / _| ___  _ __ _ __ ___ (_) ___|_ _|  _ \                #
#                   | |_ / _ \| '__| '_ ` _ \| |/ __|| || | | |               #
#                   |  _| (_) | |  | | | | | | | (__ | || |_| |               #
#                   |_|  \___/|_|  |_| |_| |_|_|\___|___|____/                #
#                                                                             #
#                                      main                                   #
#                                                                             #
###############################################################################
"""Description:
This is were it all happens. This file initializes the model and trains the
model, after which training is evaluated and the model is tested.
"""

# Packages
###############################################################################

# Standard library imports
import logging
import os

# Deeplearning tools imports
import tensorflow as tf
from keras import __version__ as keras_version
from keras import backend as K

# FormicID imports
from models.models import compile_model
from models.models import load_model
from testers.tester import evaluator
from testers.tester import plot_confusion_matrix
from testers.tester import predictor
from testers.tester import predictor_reports
from trainers.train import trainer_dir
from utils.load_config import process_config
from utils.logger import build_logger
from utils.logger import plot_history
from utils.model_utils import save_model
from utils.model_utils import weights_load
from utils.utils import create_dirs
from utils.utils import get_args

# Parameters and settings
###############################################################################
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
# To disable the tf warning for compiling in SEE4.2
# 0 = all logs, 1 = info, 2 = warnings, 3 = error


# Main
###############################################################################


def main():
    # Arguments
    ###########################################################################
    try:
        args = get_args()
        config = process_config(args.config)
    except:
        logging.error("Missing or invalid arguments.")
        exit(0)

    # Logging
    ###########################################################################
    logging.basicConfig(
        filename=os.path.join("logs", config.exp_name + ".log"),
        format="[%(asctime)s] - [%(levelname)s]: %(message)s",
        filemode="a",
        level=logging.DEBUG,
    )
    logging.info("Logging started.")
    logging.info("Keras version: {}".format(keras_version))

    # Session
    ###########################################################################
    sess = tf.Session()
    K.set_session(sess)

    # create experiment related directories
    ###########################################################################
    create_dirs([config.summary_dir, config.checkpoint_dir])

    # Initialize the model
    ###########################################################################
    model_formicID = load_model(config=config, num_species=97)
    model_formicID = compile_model(model=model_formicID, config=config)
    model_formicID = weights_load(
        model=model_formicID,
        weights="experiments/T97_CaAll_QuM_ShSti_AugM_D05_LR0001_E200_I4_def_clean/checkpoint/weights_55-1.76.hdf5",
    )

    # Training in batches with iterator
    ###########################################################################
    history = trainer_dir(
        model=model_formicID,
        config=config,
        callbacks=build_logger(config=config, model=model_formicID),
    )
    save_model(
        model=model_formicID, filename="final_weights.hdf5", config=config
    )

    # Evaluation
    ###########################################################################
    plot_history(history=history, config=config, theme="ggplot", save=None)
    evaluator(model=model_formicID, config=config, test_dir=None)

    # Testing
    ###########################################################################
    Y_true, Y_pred, labels, species_dict = predictor(
        model=model_formicID,
        config=config,
        # species_json="data/species_dict.json",
        plot=True,
        n_img=10,
        n_cols=3,
    )
    predictor_reports(
        Y_true=Y_true,
        Y_pred=Y_pred,
        config=config,
        species_dict=species_dict,
        target_names=labels,
        digits=5,
    )
    plot_confusion_matrix(
        Y_pred=Y_pred,
        Y_true=Y_true,
        config=config,
        target_names=labels,
        species_dict=species_dict,
        title=None,
        cmap="viridis",
        normalize=True,
        scores=True,
        score_size=8,
        save="confusion_matrix.png",
    )
    # Footer
    ###########################################################################
    K.clear_session()
    logging.info("Logging ended.")


if __name__ == "__main__":
    main()
