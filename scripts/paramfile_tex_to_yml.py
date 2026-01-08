import yaml
from pathlib import Path


runs = 10

for run in range(runs):
    input_tex_filepath = f"/project/b/babul/stlock/simulations/gizmo-mufasa/params_hyenasc_l0_BAL_f_{run+1:02}.tex"
    output_yml_filepath = f"/scratch/b/babul/stlock/snapshots/HyenasC/L0/halo_3224_finished_runs/BAL_variations/run{run+1:02}/params.yml"

    def parse_tex_paramfile(filepath: Path) -> dict:
        # Parses a GIZMO .tex-style parameter file and extracts parameters.
        params = {}
        with open(filepath, "r") as f:
            lines = f.readlines()

        for line in lines:
            line = line.strip()
            if not line or line.startswith("%"):
                continue  # Skip comments and empty lines

            # Split off comments
            line = line.split("%")[0].strip()
            if not line:
                continue

            parts = line.split()
            if len(parts) >= 2:
                key, value = parts[0], parts[1]
                try:
                    params[key] = float(value)
                except ValueError:
                    params[key] = value  # Keep as string if not a float
        return params

    def convert_to_yaml(tex_file: Path, yaml_file: Path, section: str = "Parameters"):
        """Converts the extracted parameters into a YAML file."""
        params = parse_tex_paramfile(tex_file)
        yaml_data = {section: params}

        with open(yaml_file, "w") as f:
            yaml.safe_dump(yaml_data, f, default_flow_style=False)

        print(f"Saved YAML to {yaml_file}")

    tex_path = Path(input_tex_filepath)
    yaml_path = Path(output_yml_filepath)

    convert_to_yaml(tex_path, yaml_path)