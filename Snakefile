regions = ["NeiMongol"]
techs = ["solar", "wind"]

rule all:
    input:
        expand("data/{region}/suitability.done", region=regions)

rule spatial_data_prep:
    output:
        touch("data/{region}/spatial_data_prep.done")
    params:
        region=lambda wc: wc.region
    script:
        "spatial_data_prep.py"

rule exclusion:
    input:
        "data/{region}/spatial_data_prep.done"
    output:
        touch("data/{region}/exclusion_{tech}.done")
    params:
        region=lambda wc: wc.region,
        tech=lambda wc: wc.tech
    script:
        "Exclusion.py"

rule suitability:
    input:
        expand("data/{{region}}/exclusion_{tech}.done", tech=techs)
    output:
        touch("data/{region}/suitability.done")
    params:
        region=lambda wc: wc.region
    script:
        "suitability.py"
