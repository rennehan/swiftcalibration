import os
from swiftemulator.design import latin
from swiftemulator.io.swift import write_parameter_files
from swiftemulator import ModelSpecification

spec = ModelSpecification(
    number_of_parameters=5,
    parameter_names=[
        "ObsidianAGN:torque_accretion_norm",
        "ObsidianAGN:jet_velocity_km_s",
        "ObsidianAGN:adaf_coupling",
        "ObsidianAGN:adaf_kick_factor",
        "ObsidianAGN:eddington_fraction_lower_boundary"
    ],
    parameter_printable_names=[
        "$\varepsilon_{\\rm torque}$",
        "$v_{\\rm jet}$",
        "$\varepsilon_{\\rm adaf}$",
        "$f_{\\rm adaf,kick}$",
        "$R_{\\rm lower}$"
    ],
    parameter_limits=[
        [0.001, 0.1],
        [2500, 10000],
        [0.00001, 0.01],
        [0., 1.],
        [0.01, 0.2]
    ],
)

parameter_transforms = {}

number_of_simulations = 32

model_parameters = latin.create_hypercube(
    model_specification=spec,
    number_of_samples=number_of_simulations,
)

base_parameter_file = "swimba_s37n256.yml"

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

