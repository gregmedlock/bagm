# resave_FBC.py - this script loads all draft models and resaves them.
# This is necessary because ModelSEED does not yet use the FBC extension
# for SBML, so many model components are in the wrong location (e.g., 
# metabolite charges are in the metabolite notes). Without this step,
# the SBML models will raise thousands of errors each time they are
# loaded.

import cobra
import os

# directory where GENREs are stored which need to be reloaded to make FBC-compliant
genre_dir = "../genres/"
# directory to save new SBML models
save_dir = "../genres/"
# new extension for resaved files
new_extension = ".fbc.sbml"

resaved = []
for model_file in os.listdir(genre_dir):
    # resave models if they don't already have the new extension.
    # This will overwrite models that are already in the new format.
    if model_file.endswith(".sbml") and not model_file.endswith(new_extension):
        model = cobra.io.read_sbml_model(genre_dir+model_file)
        cobra.io.write_sbml_model(model,save_dir+os.path.splitext(model_file)[0]+new_extension)
        resaved.append(model_file)
        
# print resaved models - need this at end, since many errors get raised for each model
for old_file in resaved:
    print("Resaved SBML model with FBC for: " + old_file)