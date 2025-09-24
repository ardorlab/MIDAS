import numpy as np
import logging
import os
import subprocess
import re
from pathlib import Path
import math
from midas_data import __serpent2exe__
from midas.utils import LWR_fuelcyclecost as fuelcyclecost
from scipy.io import loadmat
from scipy.interpolate import interp1d

## Initialize logging for the present file
logger = logging.getLogger("MIDAS_logger")

def evaluate(solution, input):

    #Predefine exit codes for optional runs so later logic doesn't fail
    dop_exit_code = 0
    dep_exit_code = 0
    results_dict = {}
    chromosome = solution.chromosome
    genome = input.genome
    template_dict = {}
    
    # Convert normalized values back to real values and store them in a dictionary for easy access
    for key, value in sorted(genome.items(), key = lambda item: item[1]['index']):
        if 'continuous_range' in value:
            chromosome[value['index']] = float(chromosome[value['index']]*(value['continuous_range'][1]-value['continuous_range'][0]) + value['continuous_range'][0])
        elif 'discrete_range' in value:
            #Account for floating point errors
            indices = np.where(np.isclose(value['normalized_discrete_range'], chromosome[value['index']],atol=1e-5))[0]
            chromosome[value['index']] = value['discrete_range'][indices[0]]
        
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
    fill_template(input.input_template["loc"], base_file, template_dict)
    if not doppler_dir.exists() and "doppler_temperature_coefficient" in genome.keys():
        os.mkdir(doppler_dir)
    if "doppler_temperature_coefficient" in genome.keys():
        fill_template(input.input_template["loc"], doppler_file, template_dict)
        update_temp(doppler_file)
    if input.depletion_settings['apply'] and not depletion_dir.exists():
        os.mkdir(depletion_dir)
    if input.depletion_settings['apply']:
        fill_template(input.input_template["loc"], depletion_file, template_dict)
        with open(depletion_file, "a") as f:
            f.write(f"\nset pop {input.depletion_settings["particles_per_history"]} {input.depletion_settings["active_cycles"]} {input.depletion_settings["inactive_cycles"]}\n")
            if input.depletion_settings['depletion_units'].lower() == 'days':
                f.write("\ndep daytot\n")
            else:
                f.write("\ndep butot\n")
            for step in input.depletion_settings['depletion_steps']:
                f.write(f"{step}\n")
        #Remove detector lines outside of base case to improve speed
        remove_detector_lines(depletion_file,input.power_peaking_detectors)
    
    #Add cross sections and population into serpent files
    for file in [base_file, doppler_file, depletion_file]:
        if file.exists():
            with open(file, "a") as f:
                if file is not depletion_file:
                    f.write(f"\nset pop {input.particles_per_history} {input.active_cycles} {input.inactive_cycles}\n")
                f.write(f"set acelib {input.xs_lib}\n")
                f.write(f"set dec {input.dec_lib}\n")
                f.write(f"set nfylib {input.fy_lib}\n")
    
    #Start depletion calc first since it takes the longest
    sss2cmd = __serpent2exe__
    if input.depletion_settings['apply']:
        os.chdir(depletion_dir)
        dep_cmd = ["mpirun","-np",f"{input.depletion_settings["mpi_ranks"]}",f"{sss2cmd}","-omp",f"{input.depletion_settings["omp_threads"]}","depletion_input"]
        dep_process = subprocess.Popen(dep_cmd)
    
    #Run base calc and get results
    os.chdir(base_dir)
    base_cmd = dop_cmd = ["mpirun","-np",f"{input.mpi_ranks}",f"{sss2cmd}","-omp",f"{input.omp_threads}","base_input"]
    base_process = subprocess.Popen(base_cmd)
    base_process.wait()
    base_exit_code = base_process.returncode
    #Unrealistically bad results to return if code fails to drive optimization away from this region
    fail_results = {'doppler_temperature_coefficient': 5, 'cycle_length': 0, 'cost_fuelcycle': 1000000000000,'total_mass': 1000000000,'fdeltah': 8,'max_doserate':50000}
    if base_exit_code == 0:
        base_results = get_serpent_results("base_input_res.m")
        base_det_results = get_serpent_results("base_input_det0.m")
        if input.power_peaking_detectors == 'ppw':
            peaking_results = base_results["PPW_POW"][:,1::2]

        else:
            peaking_results = [] #!TODO: Add parsing logic for going through detector results 

        peaking_results = peaking_results[peaking_results != 0]
        mean_pow = np.mean(peaking_results)
        peaking_factors = peaking_results / mean_pow

        results_dict["fdeltah"] = np.max(peaking_factors)

        mass_dict = get_masses(base_dir / "base_input.out")

        #Exclude materials not to be included in mass calculation specified by user
        results_dict["total_mass"] = 0
        for key in mass_dict:
            if key in input.mass_materials or input.mass_materials == 'all':
                #Store mass of each included material and convert from g to lb
                results_dict["total_mass"] += float(mass_dict[key])/453.6
    else:
        results_dict = fail_results

    
    

    if doppler_dir.exists():
        os.chdir(doppler_dir)
        #Remove lines in file for power peaking detectors to improve runtime
        remove_detector_lines(doppler_file,input.power_peaking_detectors)
        dop_cmd = ["mpirun","-np",f"{input.mpi_ranks}",f"{sss2cmd}","-omp",f"{input.omp_threads}","doppler_input"]
        dop_process = subprocess.Popen(dop_cmd)
        dop_process.wait()
        dop_exit_code = dop_process.returncode
        if dop_exit_code ==0 and base_exit_code == 0:
            dop_results = get_serpent_results("doppler_input_res.m")

            rho1 = (base_results["ABS_KEFF"][-1,0] - 1)/base_results["ABS_KEFF"][-1,0]
            rho2 = (dop_results["ABS_KEFF"][-1,0] - 1)/dop_results["ABS_KEFF"][-1,0]

            results_dict["doppler_temperature_coefficient"] = ((rho2 - rho1) / 50) * 10**5
    
    if depletion_dir.exists():
        dep_process.wait()
        dep_exit_code = dep_process.returncode
        if dep_exit_code == 0 and dop_exit_code == 0 and base_exit_code == 0:
            dep_results = get_serpent_results("depletion_input_res.m")
            burn_days = dep_results["BURN_DAYS"][:,0]
            burn_keff = dep_results["ABS_KEFF"][:,0]
    
            burn_days = np.unique(burn_days)
            burn_keff = np.unique(burn_keff)
            interpolate = interp1d(burn_days,burn_keff,kind='linear',fill_value='extrapolate')
            results_dict["cycle_length"] = interpolate(1.0)
    
    if 'cost_fuelcycle' in input.objectives():
        hm_frac = get_heavy_metal_percent(base_file)
        cycle_cost = fuelcyclecost.calc_fuelcost_triso(template_dict["enrichment"],(mass_dict['fuel']*453.6/1000)*hm_frac)
        results_dict['cost_fuelcycle'] = cycle_cost
    
    for key, value in results_dict.items():
        if key in input.objectives():
            solution.parameters[key]["value"] = value
        else:
            logger.info(f"Objective {key} is available in serpent but not currently used by the optimization")

    return solution

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

import re
from pathlib import Path

def get_masses(filepath):
    """
    Reads a Serpent output-like file and extracts the mass of all materials.

    Parameters
    ----------
    filepath : str or Path
        Path to the file to read.

    Returns
    -------
    dict
        Dictionary mapping material names to masses in grams, e.g.,
        { "moderator": 1.13e6, "fuel": 8.7e5 }
    """
    filepath = Path(filepath)
    masses = {}
    current_material = None
    mass_pattern = re.compile(r"- Mass\s+([0-9.E+-]+)\s+g")

    with filepath.open("r") as f:
        for line in f:
            line_strip = line.strip()
            # Detect start of material block
            if line_strip.startswith("Material "):
                # Extract material name in quotes
                match_name = re.match(r'Material\s+"(.+?)"', line_strip)
                if match_name:
                    current_material = match_name.group(1)
                else:
                    current_material = None
                continue

            # If inside a material block, look for mass
            if current_material:
                match_mass = mass_pattern.search(line_strip)
                if match_mass:
                    masses[current_material] = float(match_mass.group(1))
                    current_material = None  # done with this material

    return masses

def remove_detector_lines(filepath, detector):
    """
    Remove lines from a file based on detector type.

    Parameters
    ----------
    filepath : str or Path
        Path to the file to modify.
    detector : str or list of str
        If 'ppw', remove lines containing 'set adf' and 'set ppw'.
        If a list, remove lines containing 'det {detector[i]}' for each element.
    """
    filepath = Path(filepath)

    # Read all lines
    with filepath.open("r") as f:
        lines = f.readlines()

    # Determine lines to remove
    if detector == 'ppw':
        remove_keywords = ['set adf', 'set ppw']
    elif isinstance(detector, list):
        remove_keywords = [f"det {d}" for d in detector]

    # Filter lines
    new_lines = [line for line in lines if not any(keyword in line for keyword in remove_keywords)]

    # Write back
    with filepath.open("w") as f:
        f.writelines(new_lines)

def get_heavy_metal_percent(filepath):
    """
    Reads a Serpent material definition and calculates the weight percent
    of isotopes with Z = 92 (U) or 94 (Pu). Uses atomic masses from the 
    `periodictable` library when possible.
    """
    with open(filepath, "r") as f:
        lines = f.readlines()

    masses = {}
    mat_name = None
    fractions = {}

    # Regex for isotopes like 92235.06c
    iso_pattern = re.compile(r"^\s*(\d+)\.\d+\w*\s+(-?\d+\.?\d*([Ee][+-]?\d+)?)")

    for line in lines:
        if line.strip().startswith("mat") and "fuel" in line:
            mat_name = line.strip()
            continue

        match = iso_pattern.match(line)
        if match:
            zaid, frac, _ = match.groups()
            Z = int(zaid[:-3]) // 1000     # first digits = Z
            A = int(zaid[:-3]) % 1000      # last 3 = mass number

            frac = float(frac)

            # Get atomic mass from mass number
            try:
                from periodictable import elements
                mass = getattr(elements, Z).isotopes[A].mass
            except Exception:
                mass = float(A)

            # Store
            fractions[(Z, A)] = (frac, mass)

    if not fractions:
        raise ValueError("No valid isotopes found in file.")

    # Compute total mass contribution
    total_mass = sum(abs(frac) for frac, _ in fractions.values())

    # Contribution from U & Pu only
    upu_mass = sum(abs(frac) for (Z, A), (frac, _) in fractions.items() if Z in [92, 94])

    # Weight percent
    weight_percent = (upu_mass / total_mass) * 100.0

    return weight_percent