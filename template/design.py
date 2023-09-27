import os
from swiftemulator.design import latin
from swiftemulator.io.swift import write_parameter_files
from swiftemulator import ModelSpecification

spec = ModelSpecification(
    number_of_parameters=4,
    parameter_names=[
        "SIMBAFeedback:FIRE_velocity_normalization",
        "SIMBAFeedback:FIRE_eta_normalization",
        "SIMBAAGN:torque_accretion_norm",
        "SIMBAAGN:jet_mass_min_Msun"
    ],
    parameter_printable_names=[
        "$v_{a}",
        "$\eta_{0}$",
        "$\varepsilon_{\\rm torque}$",
        "$\\log_{10}$ $M_{\\rm BH, lim}/M_{\\rm \odot}$"
    ],
    parameter_limits=[
        [0.5, 5.0],
        [5.0, 13.0],
        [0.00001, 0.008],
        [6.5, 8.0]
    ],
)

parameter_transforms = {"SIMBAAGN:jet_mass_min_Msun": lambda x: 10.0 ** x}

number_of_simulations = 32

model_parameters = latin.create_hypercube(
    model_specification=spec,
    number_of_samples=number_of_simulations,
)

base_parameter_file = "swimba_s25n256.yml"

output_path = "./original_ymls"

if not os.path.isdir(output_path):
    os.makedirs(output_path, mode = 0o755, exist_ok = True)

write_parameter_files(
    filenames={
        key: f"{output_path}/{key}.yml"
        for key in model_parameters.model_parameters.keys()
    },
    model_parameters=model_parameters,
    parameter_transforms=parameter_transforms,
    base_parameter_file=base_parameter_file,
)

