import numpy as np
import logging
import os
import subprocess
import re
from pathlib import Path
import math
from midas_data import __serpent2exe__
from scipy.io import loadmat

## Initialize logging for the present file
logger = logging.getLogger("MIDAS_logger")

def evaluate(solution, input):

    chromosome = solution.chromosome
    genome = input.genome
    template_dict = {}

    # Convert normalized values back to real values and store them in a dictionary for easy access
    for key, value in sorted(genome.items(), key = lambda item: item[1]['index']):
        if 'continuous_range' in value:
            chromosome[value['index']] = float(chromosome[value['index']]*(value['continuous_range'][1]-value['continuous_range'][0]) + value['continuous_range'][0])
        elif 'discrete_range' in value:
            indices = value['normalized_discrete_range'].index(chromosome[value['index']])
            chromosome[value['index']] = value['discrete_range'][indices]
        
        template_dict[key] = chromosome[value['index']]
    
    ## Create and move to unique directory for Serpent execution
    cwd = Path(os.getcwd())
    indv_dir = cwd.joinpath(input.results_dir_name / Path(solution.name))
    base_dir = indv_dir / "base"
    doppler_dir = indv_dir / "doppler"
    depletion_dir = indv_dir / "depletion"
    base_file = base_dir / "base_input"
    doppler_file = doppler_dir / "doppler_input"
    depletion_file = depletion_dir / "depletion_input"
    if not indv_dir.exists():
        logger.debug(f"Creating new results directory: {indv_dir}")
        os.mkdir(indv_dir)
    logger.debug(f"Changing to new working directory: {indv_dir}")
    os.chdir(indv_dir)
    #Create subdirectories for base, DTC, and depletion runs
    if not base_dir.exists():
        os.mkdir(base_dir)
    fill_template(input.template_file["loc"], base_file, template_dict)
    if not doppler_dir.exists() and "doppler_temperature_coefficient" in genome.keys():
        os.mkdir(doppler_dir)
    if "doppler_temperature_coefficient" in genome.keys():
        fill_template(input.template_file["loc"], doppler_file, template_dict)
        update_temp(doppler_file)
    if input.depletion_settings['apply'] and not depletion_dir.exists():
        os.mkdir(depletion_dir)
    if input.depletion_settings['apply']:
        fill_template(input.template_file["loc"], depletion_file, template_dict)
        with open(depletion_file, "a") as f:
            if input.depletion_settings['depletion_units'].lower() == 'days':
                f.write("\ndep daytot\n")
            else:
                f.write("\ndep butot\n")
            for step in input.depletion_settings['depletion_steps']:
                f.write(f"{step}\n")
    
    for file in [base_file, doppler_file, depletion_file]:
        if file.exists():
            with open(file, "a") as f:
                f.write(f"\nset pop {input.particles_per_history} {input.active_cycles} {input.inactive_cycles}\n")
                f.write(f"set acelib {input.xs_lib}\n")
                f.write(f"set dec {input.dec_lib}\n")
                f.write(f"set nfylib {input.fy_lib}\n")
    
    sss2cmd = __serpent2exe__
    if input.depletion_settings['apply']:
        os.chdir(depletion_dir)
        dep_cmd = ["mpirun","-np",f"{input.depletion_settings["mpi_ranks"]}",f"{sss2cmd}","-omp",f"{input.depletion_settings["omp_threads"]}","depletion_input"]
        dep_process = subprocess.Popen(dep_cmd)
    
    if doppler_dir.exists():
        os.chdir(doppler_dir)
        dop_cmd = ["mpirun","-np",f"{input.mpi_ranks}",f"{sss2cmd}","-omp",f"{input.omp_threads}","doppler_input"]
        dop_process = subprocess.Popen(dop_cmd)
        dop_process.wait()
        dop_exit_code = dop_process.returncode
        dop_results = get_serpent_results("doppler_input_res.m")
    
    os.chdir(base_dir)
    base_cmd = dop_cmd = ["mpirun","-np",f"{input.mpi_ranks}",f"{sss2cmd}","-omp",f"{input.omp_threads}","base_input"]
    base_process = subprocess.Popen(base_cmd)
    base_process.wait()
    base_exit_code = base_process.returncode
    base_results = get_serpent_results("base_input_res.m")
    base_det_results = get_serpent_results("base_input_det0.m")
    if depletion_dir.exists():
        dep_process.wait()
        dep_exit_code = dep_process.returncode
        dep_results = get_serpent_results("depletion_input_res.m")
        burn_days = dep_results["BURN_DAYS"]
        burn_keff = dep_results["ABS_KEFF"]
    
    


    raise NotImplementedError("Serpent evaluation function is not yet implemented.")

def fill_template(template_path, output_path, template_dict):
    """
    Fill a template file with values from template_dict and save to output_path.
    - Placeholders {var}, {var * 2}, {sin(var)}, etc. are supported.
    """

    # Safe math environment (only what you allow)
    safe_env = {k: getattr(math, k) for k in dir(math) if not k.startswith("__")}
    safe_env.update(template_dict)

    # Read template
    template_text = Path(template_path).read_text()

    # Regex: find everything inside {}
    pattern = re.compile(r"\{(.*?)\}")

    def replace_match(match):
        expr = match.group(1).strip()
        try:
            return str(eval(expr, {"__builtins__": {}}, safe_env))
        except Exception as e:
            raise ValueError(f"Error evaluating expression '{expr}': {e}")

    # Replace placeholders
    filled_text = pattern.sub(replace_match, template_text)

    # Save to output
    Path(output_path).write_text(filled_text)

def get_serpent_results(output_file):

    #Run .m file in matlab to convert into a python-readable .mat file
    subprocess.run(["matlab",
    "-nodisplay","-nosplash","-nodesktop",
    "-r",f"run('{output_file}'); save('{output_file}.mat'); exit"])

    data = loadmat(f"{output_file}.mat")

    return data

def update_temp(filename):
    tmp_pattern = re.compile(r"(tmp\s+)([-+]?\d*\.?\d+)", re.IGNORECASE)

    with open(filename, "r") as f:
        lines = f.readlines()

    updated_lines = []
    for line in lines:
        if "mat" in line and "fuel" in line:  # quick filter
            # Replace tmp number with number+50
            def repl(match):
                prefix = match.group(1)  # "tmp "
                number = float(match.group(2))
                return f"{prefix}{number + 50:g}"  # keep clean formatting

            new_line = tmp_pattern.sub(repl, line)
            updated_lines.append(new_line)
        else:
            updated_lines.append(line)

    with open(filename, "w") as f:
        f.writelines(updated_lines)
