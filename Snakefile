import yaml
from utils.data_preprocessing import clean_region_name

# Load base configuration
with open('configs/config.yaml') as f:
    config = yaml.safe_load(f)

tech = config.get('tech', '')

# Select technology specific defaults
if tech in ['SolarPV', 'SolarPVTracking']:
    with open('configs/solar.yaml') as f:
        tech_cfg = yaml.safe_load(f)
else:
    with open('configs/wind.yaml') as f:
        tech_cfg = yaml.safe_load(f)

# Merge defaults where not provided in base config

def merge(dst, src):
    for k, v in src.items():
        if k not in dst:
            dst[k] = v
        elif isinstance(v, dict) and isinstance(dst[k], dict):
            merge(dst[k], v)
merge(config, tech_cfg)

region_folder = config['region_folder_name']
region_clean = clean_region_name(config['region_name'])

spatial_done = f"data/{region_folder}/spatial_data_prep.done"
exclusion_done = f"data/{region_folder}/land_exclusion.done"
weather_done = f"data/{region_folder}/weather_data_prep.done"
suitability_done = f"data/{region_folder}/suitability.done"
profile_csv = (
    f"data/{region_folder}/energy_profiles/"
    f"{config['tech']}_profile_{region_clean}.csv"
)

rule all:
    input: profile_csv

rule spatial_data_prep:
    output: spatial_done
    shell: 'python spatial_data_prep.py && touch {output}'

rule land_exclusion:
    input: spatial_done
    output: exclusion_done
    shell: 'python Exclusion.py && touch {output}'

rule weather_data_prep:
    input: exclusion_done
    output: weather_done
    shell: 'python weather_data_prep.py && touch {output}'

rule suitability:
    input: weather_done
    output: suitability_done
    shell: 'python suitability.py && touch {output}'

rule energy_profiles:
    input: suitability_done
    output: profile_csv
    shell: 'python energy_profiles.py'
