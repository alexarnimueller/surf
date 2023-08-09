from ord_schema.proto import reaction_pb2

doi_pattern = "10[.][0-9]{4,}(?:[.][0-9]+)*/(?:(?![\"&'<>])\S)+"
email_pattern = "[\w.+-]+@[\w-]+\.[\w.-]+"
orcid_pattern = "\d{4}-\d{4}-\d{4}-\d{3}[0-9X]"
rxn_smiles_left = "(startingmat|reagent|catalyst|ligand|additive)_[1-9]_smiles"
rxn_smiles_right = "product_[1-9]_smiles"
role_pattern = "(startingmat|reagent|catalyst|ligand|additive|product|standard)"
const_cols = [
    "rxn_id",
    "source_id",
    "source_type",
    "rxn_date",
    "rxn_type",
    "rxn_name",
    "rxn_tech",
    "temperature_deg_c",
    "time_h",
    "atmosphere",
    "stirring_shaking",
    "scale_mol",
    "concentration_mol_l",
]
last_cols = ["provenance", "comment", "procedure"]

reaction = reaction_pb2.Reaction()
mapping_atmo = {
    "AR": reaction.conditions.pressure.Atmosphere.ARGON,
    "ARGON": reaction.conditions.pressure.Atmosphere.ARGON,
    "N2": reaction.conditions.pressure.Atmosphere.NITROGEN,
    "NITROGEN": reaction.conditions.pressure.Atmosphere.NITROGEN,
    "AIR": reaction.conditions.pressure.Atmosphere.AIR,
    "NONE": reaction.conditions.pressure.Atmosphere.UNSPECIFIED,
    "NA": reaction.conditions.pressure.Atmosphere.UNSPECIFIED,
    "UNSPECIFIED": reaction.conditions.pressure.Atmosphere.UNSPECIFIED,
    "O3": reaction.conditions.pressure.Atmosphere.OZONE,
    "OZON": reaction.conditions.pressure.Atmosphere.OZONE,
    "OZONE": reaction.conditions.pressure.Atmosphere.OZONE,
    "O2": reaction.conditions.pressure.Atmosphere.OXYGEN,
    "OXYGEN": reaction.conditions.pressure.Atmosphere.OXYGEN,
    "CO2": reaction.conditions.pressure.Atmosphere.CARBON_DIOXIDE,
    "CARBON_DIOXIDE": reaction.conditions.pressure.Atmosphere.CARBON_DIOXIDE,
    "H2": reaction.conditions.pressure.Atmosphere.HYDROGEN,
    "HYDROGEN": reaction.conditions.pressure.Atmosphere.HYDROGEN,
    "CH4": reaction.conditions.pressure.Atmosphere.METHANE,
    "METHANE": reaction.conditions.pressure.Atmosphere.METHANE,
    "NH3": reaction.conditions.pressure.Atmosphere.AMMONIA,
    "AMMONIA": reaction.conditions.pressure.Atmosphere.AMMONIA,
    "AMONIA": reaction.conditions.pressure.Atmosphere.AMMONIA,
    "C2H4": reaction.conditions.pressure.Atmosphere.ETHYLENE,
    "ETHYLENE": reaction.conditions.pressure.Atmosphere.ETHYLENE,
    "C2H2": reaction.conditions.pressure.Atmosphere.ACETYLENE,
    "ACETYLENE": reaction.conditions.pressure.Atmosphere.ACETYLENE,
}

mapping_role = {
    "startingmat": "REACTANT",
    "reagent": "REAGENT",
    "additive": "REAGENT",
    "solvent": "SOLVENT",
    "catalyst": "CATALYST",
    "ligand": "REAGENT",
    "product": "PRODUCT",
}

mapping_role_ord = {
    1: "startingmat",
    2: "reagent",
    3: "solvent",
    4: "catalyst",
    6: "standard",
    7: "standard",
    8: "product",
    9: "product",
    10: "product",
}

mapping_analyses_ord = {
    1: "UNKNOWN",
    2: "LC",
    3: "GC",
    4: "IR",
    5: "NMR_1H",
    6: "NMR_13C",
    7: "NMR",
    8: "MP",
    9: "UV",
    10: "TLC",
    11: "MS",
    12: "HRMS",
    13: "MSMS",
    14: "WEIGHT",
    15: "LCMS",
    16: "GCMS",
    17: "ELSD",
    18: "CD",
    19: "SFC",
    20: "EPR",
    21: "XRD",
    22: "RAMAN",
    23: "ED",
}

mapping_stirring = {
    "UNSPECIFIED": 0,
    "CUSTOM": 1,
    "NONE": 2,
    "STIR_BAR": 3,
    "OVERHEAD_MIXER": 4,
    "AGITATION": 5,
    "BALL_MILLING": 6,
    "SONICATION": 7,
}

mapping_stirring_ord = {
    0: "UNSPECIFIED",
    1: "CUSTOM",
    2: "NONE",
    3: "STIR_BAR",
    4: "OVERHEAD_MIXER",
    5: "AGITATION",
    6: "BALL_MILLING",
    7: "SONICATION",
}

mapping_time = {0: 1, 1: 1, 2: 0.016666, 3: 0.000277, 4: 24}

mapping_vol_mol_mass = {1: 1000, 2: 1, 3: 0.001}
