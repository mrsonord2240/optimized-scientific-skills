from scripts.substructure_search import compile_smarts

REACTIVE_SMARTS = {
    'acid_anhydride': '[CX3](=O)O[CX3](=O)',
    'acid_halide': '[CX3](=O)[F,Cl,Br,I]',
    'alpha_halo_carbonyl': '[CX3](=O)C([F,Cl,Br,I])',
    'aldehyde_reactive': '[CX3H1](=O)[#6;X4]',  # aliphatic aldehydes
    'epoxide': 'C1OC1',
    'aziridine': 'C1NC1',
    'isocyanate': '[NX2]=C=[OX1]',
    'isothiocyanate': '[NX2]=C=[SX1]',
    'beta_lactam': 'C1(=O)NCC1',
    'sulfonyl_halide': '[SX4](=O)(=O)[F,Cl,Br,I]',
    'Michael_acceptor': '[CX3]=[CX3][CX3]=O',
    'vinyl_sulfone': '[SX4](=O)(=O)C=C',
}

# Compile once per filter configuration, not once per molecule.
COMPILED_REACTIVE_SMARTS = {
    name: compile_smarts(smarts) for name, smarts in REACTIVE_SMARTS.items()
}

def reactive_filter(mol, exclude_warheads=True):
    if not exclude_warheads:
        return False, None
    if mol is None:
        return False, None
    for name, pattern in COMPILED_REACTIVE_SMARTS.items():
        if mol.HasSubstructMatch(pattern):
            return True, name
    return False, None
