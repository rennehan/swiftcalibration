from swiftemulator.design import latin
from swiftemulator.io.swift import write_parameter_files
from swiftemulator import ModelSpecification

spec = ModelSpecification(
    number_of_parameters=4,
    parameter_names=[
        "SIMBAAGN:xray_f_gas_limit",
        "SIMBAAGN:f_accretion",
        "SIMBAAGN:torque_accretion_norm",
        "SIMBAAGN:wind_momentum_flux"
    ],
    parameter_printable_names=[
        "$f_{\\rm gas,xray}",
        "$f_{\\rm acc}$",
        "$\epsilon_{\\rm torque}$",
        "Momentum flux"
    ],
    parameter_limits=[
        [0.2, 1.0],
        [0.05, 0.4],
        [0.05, 0.5],
        [20, 100]
    ],
)

parameter_transforms = {}

number_of_simulations = 64

model_parameters = latin.create_hypercube(
    model_specification=spec,
    number_of_samples=number_of_simulations,
)

base_parameter_file = "template.yml"
output_path = "./original_ymls"

write_parameter_files(
    filenames={
        key: f"{output_path}/{key}.yml"
        for key in model_parameters.model_parameters.keys()
    },
    model_parameters=model_parameters,
    parameter_transforms=parameter_transforms,
    base_parameter_file=base_parameter_file,
)

