import os
from swiftemulator.design import latin
from swiftemulator.io.swift import write_parameter_files
from swiftemulator import ModelSpecification
import convert_paramfile as cpf
from pathlib import Path
import argparse
import re


def write_parameter_files_generic(
    filenames: dict,
    model_parameters,
    parameter_transforms: dict,
    base_parameter_file: str,
):
    """
    Write parameter files by reading a base parameter file (any format) and 
    replacing parameter values with those from the model parameters.
    
    Args:
        filenames: Dictionary mapping simulation keys to output file paths
        model_parameters: ModelParameters object with parameter values
        parameter_transforms: Dictionary of parameter_name -> transform_function
        base_parameter_file: Path to base parameter file (any format)
    """
    # Read the base parameter file
    with open(base_parameter_file, 'r') as f:
        base_content = f.read()
    
    # For each simulation, create a new parameter file
    for key, output_path in filenames.items():
        print(f"Writing parameter file for simulation {key} to {output_path}")
        content = base_content
        
        # Get parameters for this simulation
        sim_params = model_parameters.model_parameters[key]
        print(sim_params)
        
        # Replace each parameter value in the content
        for param_name, param_value in sim_params.items():
            # Apply transformation if it exists
            if param_name in parameter_transforms:
                print('Transforming parameter:', param_name)
                param_value = parameter_transforms[param_name](param_value)
            print(f"  Setting {param_name} to {param_value}")
            
            # Try multiple replacement patterns for different file formats
            # Pattern 1: param_name = value or param_name: value
            pattern1 = f'^(\\s*{re.escape(param_name)}\\s*[=:])\\s*[^\n]*'
            if re.search(pattern1, content, re.MULTILINE):
                content = re.sub(pattern1, f'\\1 {param_value}', content, flags=re.MULTILINE)
            else:
                # Pattern 2: Look for the parameter name with word boundaries
                pattern2 = f'\\b{re.escape(param_name)}\\b\\s*[=:]\\s*[^\n]*'
                if re.search(pattern2, content):
                    content = re.sub(pattern2, f'{param_name} = {param_value}', content)
        
        # Write the output file
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w') as f:
            f.write(content)


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


num_simulations = 128
convert_input_tex_to_yml = True
convert_output_yml_to_tex = True

base_param_file = '/scratch/aspadawe/sims/HyenasC/L1/SimbaC_L1_Calibration/trillium/gizmo-hyenasc-l1-correct_jet-vary_params/params-template.tex'
output_path = f'/scratch/aspadawe/sims/HyenasC/L1/SimbaC_L1_Calibration/trillium/gizmo-hyenasc-l1-correct_jet-vary_params/cali_params/{num_simulations}_sims-6_params/'




spec = ModelSpecification(
    number_of_parameters=6,
    parameter_names=[
        # "Cosmology:h",
        # "Cosmology:Omega_cdm",
        'Parameters:SeedBlackHoleMass',
        'Parameters:BlackHoleRadiativeEfficiency',
        'Parameters:BAL_f_accretion',
        'Parameters:EddingtonRatioForJet',
        'Parameters:BlackHoleDecoupleTime',
        'Parameters:BlackHoleJetVelocity',
    ],
    parameter_printable_names=[
        # "h",
        # r'$\Omega_{\mathrm{cdm}}$',
        r'$M_{\mathrm{BH,seed}} \, [\mathrm{M_{\odot}/h}]$',
        r'\eta',
        r'$f_{\mathrm{acc}}$',
        r'$f_{\mathrm{Edd,jet}}$',
        r'$t_{\mathrm{decouple}}/t_H(z)$',
        r'$v_{\mathrm{jet}} \, [\mathrm{km/s}]$',
    ],
    parameter_limits=[
        # [0.5,1],
        # [0.2,0.4],
        [-8, -4],
        [-2, 0],
        [-3, 0],
        [-3, 0],
        [-6, -2],
        [1000, 20000],
    ],
)

parameter_transforms = {
    'Parameters:SeedBlackHoleMass': lambda x: 10.0 ** x,
    'Parameters:BlackHoleRadiativeEfficiency': lambda x: 10.0 ** x,
    'Parameters:BAL_f_accretion': lambda x: 10.0 ** x,
    'Parameters:EddingtonRatioForJet': lambda x: 10.0 ** x,
    'Parameters:BlackHoleDecoupleTime': lambda x: 10.0 ** x,
}

model_parameters = latin.create_hypercube(
    model_specification=spec,
    number_of_samples=num_simulations,
)

# base_parameter_file = '/project/b/babul/stlock/simulations/gizmo-mufasa/params_hyenasc_l0_BAL_f_01.tex'

# output_path = "./original_ymls"
# output_path = '/scratch/b/babul/stlock/swiftcalibration_files/BAL_f_accretion'

if not os.path.isdir(output_path):
    os.makedirs(output_path, mode = 0o755, exist_ok = True)

if convert_input_tex_to_yml:
    base_yml_file = os.path.splitext(base_param_file)[0] + '.yml'

    # if not os.path.isfile(base_yml_file):
    cpf.convert_to_yaml(
        tex_file=Path(base_param_file),
        yaml_file=Path(base_yml_file),
    )
    
    base_param_file = base_yml_file

print(base_param_file)

paramfile_type = os.path.splitext(base_param_file)[1]
print(paramfile_type)
# if paramfile_type not in ['.yml', '.tex']:
#     raise ValueError("Base parameter file must be .yml or .tex format")

write_parameter_files(
    filenames={
        key: f"{output_path}/{key}.yml" for key in model_parameters.model_parameters.keys()
    },
    model_parameters=model_parameters,
    parameter_transforms=parameter_transforms,
    base_parameter_file=base_param_file,
)

# write_parameter_files_generic(
#     filenames={
#         key: f"{output_path}/{key}{paramfile_type}" for key in model_parameters.model_parameters.keys()
#     },
#     model_parameters=model_parameters,
#     parameter_transforms=parameter_transforms,
#     base_parameter_file=base_param_file,
# )

if convert_output_yml_to_tex:
    for key in model_parameters.model_parameters.keys():
        yml_file = f"{output_path}/{key}.yml"
        tex_file = f"{output_path}/{key}.tex"
        cpf.convert_to_tex(
            yaml_file=Path(yml_file),
            tex_file=Path(tex_file),
        )