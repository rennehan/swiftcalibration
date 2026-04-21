import os
from swiftemulator.design import latin
from swiftemulator.io.swift import write_parameter_files
from swiftemulator import ModelSpecification
import numpy as np
import convert_paramfile as cpf
from pathlib import Path
import argparse
import re
import dill





print("\nDesigning calibration parameter files...\n")




def save_object(obj, filename):
    with open(filename, 'wb') as f:  # Overwrites any existing file.
        dill.dump(obj, f, dill.HIGHEST_PROTOCOL)



## Command line keyword arguments
# parser = argparse.ArgumentParser(description="Design calibration parameter files.")
# parser.add_argument('--base_param_file', type=str, required=True,
#                    help='Path to base template parameter file (.yml or .tex)')
# parser.add_argument('--output_path', type=str, required=True,
#                    help='Path to directory where output parameter files will be written')
# parser.add_argument('--convert_input_tex_to_yml', action=argparse.BooleanOptionalAction, default=False,
#                    help='Whether to convert input template .tex parameter file to .yml format')
# parser.add_argument('--convert_output_yml_to_tex', action=argparse.BooleanOptionalAction, default=False,
#                    help='Whether to convert output .yml parameter files to .tex format')
# parser.add_argument('--num_simulations', type=int, default=32,
#                    help='Number of simulation parameter files to generate')
# args = parser.parse_args()



## Specify calibration design parameters ##

num_simulations = 32
convert_input_tex_to_yml = True
convert_output_yml_to_tex = True

base_param_file = './sim_files/gizmo-simba-hyenasc/params-template.tex'
output_dir = '/path/to/output/parameter/files/'
model_info_filename = 'model_info.pkl'


## Define calibration model parameter space
## Current example is for gizmo-simba-hyenasc black hole feedback parameters; modify as needed for your calibration parameters and model
model_specification = ModelSpecification(
    number_of_parameters=3,
    parameter_names=[
        # 'Parameters:SeedBlackHoleMass',
        'Parameters:BlackHoleRadiativeEfficiency',
        'Parameters:BAL_f_accretion',
        # 'Parameters:EddingtonRatioForJet',
        # 'Parameters:BlackHoleDecoupleTime',
        'Parameters:BlackHoleJetVelocity',
        # 'Parameters:BlackHoleNgbFactor',
    ],
    parameter_printable_names=[
        # r'$M_{\mathrm{BH,seed}} \, [\mathrm{M_{\odot}/h}]$',
        r'$\eta$',
        r'$f_{\mathrm{acc}} \equiv \frac{\dot{M}_{\mathrm{BH}}}{\dot{M}_{\mathrm{acc}}}$',
        # r'$f_{\mathrm{Edd,jet}}$',
        # r'$t_{\mathrm{decouple}}/t_H(z)$',
        r'$v_{\mathrm{jet}} \, [\mathrm{km/s}]$',
        # r'$f_{\mathrm{BH,kernel}} \equiv N_{\mathrm{BH,kernel}}/64$',
    ],
    parameter_limits=[
        # [-8, -4],
        [-2, np.log10(0.8)],
        [-2, np.log10(0.8)],
        # [-3, 0],
        # [-6, -2],
        [5000, 12000],
        # [1, 64],
    ],
)

parameter_transforms = {
    # 'Parameters:SeedBlackHoleMass': lambda x: 10.0 ** x,
    'Parameters:BlackHoleRadiativeEfficiency': lambda x: 10.0 ** x,
    'Parameters:BAL_f_accretion': lambda x: 10.0 ** x,
    # 'Parameters:EddingtonRatioForJet': lambda x: 10.0 ** x,
    # 'Parameters:BlackHoleDecoupleTime': lambda x: 10.0 ** x,
}


## Generate calibration parameter hypercube
model_parameters = latin.create_hypercube(
    model_specification=model_specification,
    number_of_samples=num_simulations,
)


if not os.path.isdir(output_dir):
    print(f"\tOutput directory {output_dir} does not exist; creating it.")
    os.makedirs(output_dir, mode = 0o755, exist_ok = True)


## If input template parameter file is .tex, convert to .yml
## as this is what the swiftemulator functions expect; will convert back to .tex after writing output parameter files if specified
if convert_input_tex_to_yml:
    print(f"\tConverting input template parameter file {base_param_file} from .tex to .yml format...")
    base_yml_file = os.path.splitext(base_param_file)[0] + '.yml'

    cpf.convert_to_yaml(
        tex_file=Path(base_param_file),
        yaml_file=Path(base_yml_file),
    )
    
    base_param_file = base_yml_file


paramfile_type = os.path.splitext(base_param_file)[1]
print(f"Using base parameter file: {base_param_file}")
print(f"Base parameter file type: {paramfile_type}")
# if paramfile_type not in ['.yml', '.tex']:
#     raise ValueError("Base parameter file must be .yml or .tex format")


## Write calibration parameter files
print(f"\tWriting calibration parameter files to {output_dir}...")
write_parameter_files(
    filenames={
        key: f"{output_dir}/{key}.yml" for key in model_parameters.model_parameters.keys()
    },
    model_parameters=model_parameters,
    parameter_transforms=parameter_transforms,
    base_parameter_file=base_param_file,
)


## If output parameter files should be in .tex format, convert from .yml back to .tex
if convert_output_yml_to_tex:
    print(f"\tConverting output parameter files from .yml to .tex format...")
    for key in model_parameters.model_parameters.keys():
        yml_file = os.path.join(output_dir, f"{key}.yml")
        tex_file = os.path.join(output_dir, f"{key}.tex")
        # yml_file = f"{output_dir}/{key}.yml"
        # tex_file = f"{output_dir}/{key}.tex"
        cpf.convert_to_tex(
            yaml_file=Path(yml_file),
            tex_file=Path(tex_file),
        )



## Save information about calibration model
model_info = {
    'model_specification':model_specification,
    'model_parameters':model_parameters,
    'parameter_transforms':parameter_transforms,
}

model_info_path = os.path.join(output_dir, model_info_filename)
save_object(model_info, model_info_path)
print(f"\tSaved calibration model information to {model_info_path}")


print("\nFinished designing calibrations.\n")