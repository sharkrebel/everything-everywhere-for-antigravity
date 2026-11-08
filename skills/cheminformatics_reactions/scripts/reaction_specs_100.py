# Explicit ChemDraw mechanism specifications for 100 elementary reactions
# Verified 100% against ElectronGraph chemical invariants & ReactionSMARTS

SPECS = {
    'rxn_001': {
        'add_h': {'r1': 7},
        'arrows': [
            {
                'desc': 'LDA deprotonation of alpha-H',
                'src': ('r2', 3),
                'tgt': ('r1', 8),
                'type': 'attack',
            },
            {
                'b1': (7, 8),
                'b2': (5, 7),
                'desc': 'C-H to C=C',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'bond': (5, 6),
                'desc': 'C=O to O',
                'mol': 'r1',
                'tgt': 6,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 270, 'r1': 270, 'r2': 0},
    },
    'rxn_002': {
        'add_h': {'r1': 1},
        'arrows': [
            {
                'desc': 'Methoxide deprotonation of alpha-H',
                'src': ('r2', 1),
                'tgt': ('r1', 8),
                'type': 'attack',
            },
            {
                'b1': (1, 8),
                'b2': (1, 2),
                'desc': 'C-H to C=C',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'bond': (2, 3),
                'desc': 'C=O to O',
                'mol': 'r1',
                'tgt': 3,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 180, 'r1': 180, 'r2': 180},
    },
    'rxn_003': {
        'add_h': {'r1': 3},
        'arrows': [
            {
                'desc': 'Hydroxide abstracts acetic acid proton',
                'src': ('r2', 0),
                'tgt': ('r1', 4),
                'type': 'attack',
            },
            {
                'bond': (3, 4),
                'desc': 'O-H bond cleaves onto oxygen',
                'mol': 'r1',
                'tgt': 3,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_004': {
        'add_h': {'r2': 0},
        'arrows': [
            {
                'desc': 'Triethylamine N attacks hydronium proton',
                'src': ('r1', 2),
                'tgt': ('r2', 1),
                'type': 'attack',
            },
            {
                'bond': (0, 1),
                'desc': 'H-O+ bond cleaves to water',
                'mol': 'r2',
                'tgt': 0,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_005': {
        'add_h': {'r1': 4},
        'arrows': [
            {
                'desc': 'Amide anion abstracts terminal alkyne proton',
                'src': ('r2', 0),
                'tgt': ('r1', 5),
                'type': 'attack',
            },
            {
                'bond': (4, 5),
                'desc': 'C-H bond cleaves to acetylide anion',
                'mol': 'r1',
                'tgt': 4,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_006': {
        'add_h': {'r2': 0},
        'arrows': [
            {
                'desc': 'Ethoxide O- abstracts water proton',
                'src': ('r1', 2),
                'tgt': ('r2', 1),
                'type': 'attack',
            },
            {
                'bond': (0, 1),
                'desc': 'H-O bond cleaves to hydroxide',
                'mol': 'r2',
                'tgt': 0,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_007': {
        'add_h': {'r1': 0},
        'arrows': [
            {
                'desc': 'Water abstracts alpha-proton',
                'src': ('r2', 0),
                'tgt': ('r1', 4),
                'type': 'attack',
            },
            {
                'b1': (0, 4),
                'b2': (0, 1),
                'desc': 'C-H bond forms C=C enol pi bond',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'bond': (1, 2),
                'desc': 'C=O+ pi bond cleaves to neutral OH',
                'mol': 'r1',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_008': {
        'add_h': {'r2': 0},
        'arrows': [
            {
                'atom': 2,
                'bond': (1, 2),
                'desc': 'Enolate O- reforms C=O',
                'mol': 'r1',
                'type': 'lone_pair_to_bond',
            },
            {
                'conn': 0,
                'desc': 'Enolate C=C attacks water proton',
                'src': ('r1', (0, 1)),
                'tgt': ('r2', 1),
                'type': 'attack',
            },
            {
                'bond': (0, 1),
                'desc': 'H-O bond cleaves to hydroxide',
                'mol': 'r2',
                'tgt': 0,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_009': {
        'add_h': {'r1': 5},
        'arrows': [
            {
                'desc': 'Ethoxide deprotonates malonate CH2',
                'src': ('r2', 2),
                'tgt': ('r1', 11),
                'type': 'attack',
            },
            {
                'bond': (5, 11),
                'desc': 'C-H bond cleaves to malonate carbanion',
                'mol': 'r1',
                'tgt': 5,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 270, 'r1': 270, 'r2': 180},
    },
    'rxn_010': {
        'add_h': {'r1': 0},
        'arrows': [
            {
                'desc': 'Hydroxide abstracts nitromethane proton',
                'src': ('r2', 0),
                'tgt': ('r1', 4),
                'type': 'attack',
            },
            {
                'b1': (0, 4),
                'b2': (0, 1),
                'desc': 'C-H bond forms C=N nitronate pi bond',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'bond': (1, 2),
                'desc': 'N=O pi bond cleaves to O-',
                'mol': 'r1',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 180, 'r1': 180, 'r2': 180},
    },
    'rxn_011': {
        'arrows': [
            {
                'conn': 0,
                'desc': 'MeMgBr attacks acetone C=O',
                'src': ('r2', (0, 1)),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 2),
                'desc': 'C=O pi bond onto O',
                'mol': 'r1',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_012': {
        'add_h': {'r2': 0},
        'arrows': [
            {
                'conn': 1,
                'desc': 'BH4- hydride attack on benzaldehyde',
                'src': ('r2', (0, 1)),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 0),
                'desc': 'C=O pi bond onto O',
                'mol': 'r1',
                'tgt': 0,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 60, 'r1': 60, 'r2': 180},
    },
    'rxn_013': {
        'arrows': [
            {
                'desc': 'Cyanide attack on cyclohexanone',
                'src': ('r2', 0),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 0),
                'desc': 'C=O pi bond onto O',
                'mol': 'r1',
                'tgt': 0,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 300, 'r1': 300, 'r2': 180},
    },
    'rxn_014': {
        'arrows': [
            {
                'desc': 'Methanol attack on protonated acetaldehyde',
                'src': ('r2', 1),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 2),
                'desc': 'C=O+ pi bond onto O',
                'mol': 'r1',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 180},
    },
    'rxn_015': {
        'arrows': [
            {
                'desc': 'Methylamine attack on cyclohexanone',
                'src': ('r2', 1),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 0),
                'desc': 'C=O pi bond onto O',
                'mol': 'r1',
                'tgt': 0,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 300, 'r1': 300, 'r2': 180},
    },
    'rxn_016': {
        'arrows': [
            {
                'desc': 'Pyrrolidine attack on cyclopentanone',
                'src': ('r2', 3),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 0),
                'desc': 'C=O pi bond onto O',
                'mol': 'r1',
                'tgt': 0,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 300, 'r1': 300, 'r2': 180},
    },
    'rxn_017': {
        'arrows': [
            {
                'desc': 'Bisulfite S attacks benzaldehyde',
                'src': ('r2', 1),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 0),
                'desc': 'C=O pi bond onto O',
                'mol': 'r1',
                'tgt': 0,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 60, 'r1': 60, 'r2': 180},
    },
    'rxn_018': {
        'arrows': [
            {
                'conn': 0,
                'desc': 'MeLi 1,2-addition to enone',
                'src': ('r2', (0, 1)),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 2),
                'desc': 'C=O pi bond onto O',
                'mol': 'r1',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_019': {
        'arrows': [
            {
                'conn': 5,
                'desc': 'Reformatsky reagent attacks acetophenone',
                'src': ('r2', (5, 6)),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 2),
                'desc': 'C=O pi bond onto O',
                'mol': 'r1',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 60, 'r1': 60, 'r2': 180},
    },
    'rxn_020': {
        'arrows': [
            {
                'conn': 0,
                'desc': 'Silyl enol ether attacks benzaldehyde',
                'src': ('r2', (0, 1)),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 0),
                'desc': 'C=O pi bond onto O',
                'mol': 'r1',
                'tgt': 0,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 60, 'r1': 60, 'r2': 180},
    },
    'rxn_021': {
        'arrows': [
            {
                'desc': 'Hydroxide attack on ester',
                'src': ('r2', 0),
                'tgt': ('r1', 3),
                'type': 'attack',
            },
            {
                'bond': (3, 4),
                'desc': 'C=O pi bond onto O',
                'mol': 'r1',
                'tgt': 4,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_022': {
        'arrows': [
            {
                'atom': 4,
                'bond': (3, 4),
                'desc': 'O- reforms C=O',
                'mol': 'r1',
                'type': 'lone_pair_to_bond',
            },
            {
                'bond': (2, 3),
                'desc': 'C-OEt cleaves to ethoxide',
                'mol': 'r1',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0},
    },
    'rxn_023': {
        'arrows': [
            {
                'desc': 'Dimethylamine attack on acetyl chloride',
                'src': ('r2', 1),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 2),
                'desc': 'C=O pi bond onto O',
                'mol': 'r1',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 180},
    },
    'rxn_024': {
        'arrows': [
            {
                'atom': 2,
                'bond': (1, 2),
                'desc': 'O- reforms C=O',
                'mol': 'r1',
                'type': 'lone_pair_to_bond',
            },
            {
                'bond': (1, 3),
                'desc': 'C-Cl cleaves to chloride',
                'mol': 'r1',
                'tgt': 3,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0},
    },
    'rxn_025': {
        'arrows': [
            {
                'desc': 'Ethanol attack on protonated acetic acid',
                'src': ('r2', 2),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 2),
                'desc': 'C=O+ pi bond onto O',
                'mol': 'r1',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 60, 'r1': 60, 'r2': 180},
    },
    'rxn_026': {
        'arrows': [
            {
                'desc': 'Methoxide attacks carbonyl carbon',
                'src': ('r2', 1),
                'tgt': ('r1', 2),
                'type': 'attack',
            },
            {
                'bond': (2, 3),
                'desc': 'C=O pi bond onto O',
                'mol': 'r1',
                'tgt': 3,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_027': {
        'arrows': [
            {
                'desc': 'Salicylic acid OH attacks anhydride carbonyl',
                'src': ('r2', 0),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 2),
                'desc': 'C=O pi bond onto O',
                'mol': 'r1',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 180},
    },
    'rxn_028': {
        'arrows': [
            {
                'desc': 'Hydroxide attack on acetamide',
                'src': ('r2', 0),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 2),
                'desc': 'C=O pi bond onto O',
                'mol': 'r1',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_029': {
        'arrows': [
            {
                'conn': 5,
                'desc': 'Phenyl Grignard attacks Weinreb carbonyl C',
                'src': ('r2', (5, 6)),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 2),
                'desc': 'C=O pi bond cleaves to tetrahedral O-',
                'mol': 'r1',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 180},
    },
    'rxn_030': {
        'arrows': [
            {
                'desc': 'Azide nucleophile attacks benzoyl chloride C=O',
                'src': ('r2', 0),
                'tgt': ('r1', 6),
                'type': 'attack',
            },
            {
                'bond': (6, 7),
                'desc': 'C=O pi bond cleaves to O-',
                'mol': 'r1',
                'tgt': 7,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 180},
    },
    'rxn_031': {
        'arrows': [
            {
                'desc': 'Cyanide backside attack on bromoethane',
                'src': ('r2', 0),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 2),
                'desc': 'C-Br bond cleaves to bromide',
                'mol': 'r1',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 240, 'r2': 180},
    },
    'rxn_032': {
        'arrows': [
            {
                'desc': 'Ethoxide backside attack on iodomethane',
                'src': ('r2', 2),
                'tgt': ('r1', 0),
                'type': 'attack',
            },
            {
                'bond': (0, 1),
                'desc': 'C-I bond cleaves to iodide',
                'mol': 'r1',
                'tgt': 1,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 180, 'r2': 180},
    },
    'rxn_033': {
        'arrows': [
            {
                'desc': 'Iodide backside attack on 1-chlorobutane',
                'src': ('r2', 0),
                'tgt': ('r1', 3),
                'type': 'attack',
            },
            {
                'bond': (3, 4),
                'desc': 'C-Cl bond cleaves to chloride',
                'mol': 'r1',
                'tgt': 4,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 240, 'r2': 0},
    },
    'rxn_034': {
        'arrows': [
            {
                'desc': 'Phthalimide N- attacks benzyl bromide benzylic carbon',
                'src': ('r2', 10),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (0, 1),
                'desc': 'C-Br bond cleaves to bromide',
                'mol': 'r1',
                'tgt': 0,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 240, 'r2': 0},
    },
    'rxn_035': {
        'arrows': [
            {
                'desc': 'Azide attack on benzyl bromide',
                'src': ('r2', 0),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (0, 1),
                'desc': 'C-Br bond cleaves to bromide',
                'mol': 'r1',
                'tgt': 0,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 180, 'r2': 180},
    },
    'rxn_036': {
        'arrows': [
            {
                'desc': 'Methoxide attack on oxirane carbon',
                'src': ('r2', 1),
                'tgt': ('r1', 0),
                'type': 'attack',
            },
            {
                'bond': (0, 2),
                'desc': 'C-O ring bond opens',
                'mol': 'r1',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 90, 'r2': 180},
    },
    'rxn_037': {
        'arrows': [
            {
                'desc': 'Methanol attack on protonated propylene oxide',
                'src': ('r2', 1),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 3),
                'desc': 'C-O+ bond opens',
                'mol': 'r1',
                'tgt': 3,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 90, 'r2': 180},
    },
    'rxn_038': {
        'arrows': [{
            'bond': (1, 4),
            'desc': 'C-Cl heterolysis to carbocation',
            'mol': 'r1',
            'tgt': 4,
            'type': 'bond_cleave',
        }],
        'rot': {'p1': 0, 'r1': 0},
    },
    'rxn_039': {
        'arrows': [{
            'desc': 'Water capture of tert-butyl cation',
            'src': ('r2', 0),
            'tgt': ('r1', 1),
            'type': 'attack',
        }],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_040': {
        'arrows': [
            {
                'desc': 'Alkoxide O- attacks C1',
                'src': ('r1', 4),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 2),
                'desc': 'C-Br bond cleaves to Br-',
                'mol': 'r1',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0},
    },
    'rxn_041': {
        'add_h': {'r1': 0},
        'arrows': [
            {
                'desc': 'Ethoxide abstracts beta-H',
                'src': ('r2', 2),
                'tgt': ('r1', 4),
                'type': 'attack',
            },
            {
                'b1': (0, 4),
                'b2': (0, 1),
                'desc': 'C-H bond collapses into C=C double bond',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'bond': (1, 2),
                'desc': 'C-Br bond cleaves to bromide leaving group',
                'mol': 'r1',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 180},
    },
    'rxn_042': {
        'add_h': {'r1': 3},
        'arrows': [
            {
                'desc': 't-BuO- abstracts methyl beta-H',
                'src': ('r2', 4),
                'tgt': ('r1', 6),
                'type': 'attack',
            },
            {
                'b1': (3, 6),
                'b2': (2, 3),
                'desc': 'C-H collapses into terminal C=C',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'bond': (2, 5),
                'desc': 'C-Br leaves as bromide',
                'mol': 'r1',
                'tgt': 5,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 120, 'r1': 120, 'r2': 0},
    },
    'rxn_043': {
        'add_h': {'r1': 1},
        'arrows': [
            {
                'desc': 'Ethoxide abstracts CH2 beta-H',
                'src': ('r2', 2),
                'tgt': ('r1', 5),
                'type': 'attack',
            },
            {
                'b1': (1, 5),
                'b2': (1, 2),
                'desc': 'C-H collapses into internal C=C (Zaitsev)',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'bond': (2, 4),
                'desc': 'C-Br leaves as bromide',
                'mol': 'r1',
                'tgt': 4,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 330, 'r1': 330, 'r2': 180},
    },
    'rxn_044': {
        'arrows': [
            {
                'atom': 3,
                'bond': (3, 4),
                'desc': 'Carbanion pushes into C=C',
                'mol': 'r1',
                'type': 'lone_pair_to_bond',
            },
            {
                'bond': (4, 5),
                'desc': 'Hydroxide leaves',
                'mol': 'r1',
                'tgt': 5,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0},
    },
    'rxn_045': {
        'arrows': [{
            'bond': (1, 4),
            'desc': 'C-Br heterolysis',
            'mol': 'r1',
            'tgt': 4,
            'type': 'bond_cleave',
        }],
        'rot': {'p1': 0, 'r1': 0},
    },
    'rxn_046': {
        'add_h': {'r1': 0},
        'arrows': [
            {
                'desc': 'Water abstracts beta-proton from carbocation',
                'src': ('r2', 0),
                'tgt': ('r1', 4),
                'type': 'attack',
            },
            {
                'b1': (0, 4),
                'b2': (0, 1),
                'desc': 'C-H sigma pair quenches carbocation to C=C',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 180},
    },
    'rxn_047': {
        'add_h': {'r1': 2},
        'arrows': [
            {
                'desc': 'Water abstracts beta-H',
                'src': ('r2', 0),
                'tgt': ('r1', 7),
                'type': 'attack',
            },
            {
                'b1': (2, 7),
                'b2': (2, 3),
                'desc': 'C-H collapses to cyclohexene C=C',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'bond': (3, 4),
                'desc': 'Water leaves as neutral leaving group',
                'mol': 'r1',
                'tgt': 4,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 210, 'r1': 210, 'r2': 180},
    },
    'rxn_048': {
        'add_h': {'r1': 1},
        'arrows': [
            {
                'desc': 'Thiocarbonyl S abstracts beta-H',
                'src': ('r1', 6),
                'tgt': ('r1', 9),
                'type': 'attack',
            },
            {
                'b1': (1, 9),
                'b2': (1, 2),
                'desc': 'C-H sigma pair forms C=C double bond',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'bond': (2, 4),
                'desc': 'C-O bond cleaves to xanthic acid',
                'mol': 'r1',
                'tgt': 4,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0},
    },
    'rxn_049': {
        'add_h': {'r1': 1},
        'arrows': [
            {
                'desc': 'Amine oxide O- abstracts beta-H',
                'src': ('r1', 7),
                'tgt': ('r1', 8),
                'type': 'attack',
            },
            {
                'b1': (1, 8),
                'b2': (1, 2),
                'desc': 'C-H sigma pair forms C=C double bond',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'bond': (2, 4),
                'desc': 'C-N bond cleaves to hydroxylamine',
                'mol': 'r1',
                'tgt': 4,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0},
    },
    'rxn_050': {
        'arrows': [
            {
                'desc': 'Carbanion forms C=C',
                'src': ('r1', 4),
                'tgt': ('r1', 5),
                'type': 'attack',
            },
            {
                'bond': (5, 6),
                'desc': 'C-O cleaves to acetate',
                'mol': 'r1',
                'tgt': 6,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0},
    },
    'rxn_051': {
        'arrows': [
            {
                'conn': 0,
                'desc': 'Ethylene C=C attacks Br',
                'src': ('r1', (0, 1)),
                'tgt': ('r2', 0),
                'type': 'attack',
            },
            {
                'bond': (0, 1),
                'desc': 'Br-Br bond cleaves to Br-',
                'mol': 'r2',
                'tgt': 1,
                'type': 'bond_cleave',
            },
            {
                'desc': 'Br attacks C1 to close bromonium',
                'src': ('r2', 0),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
        ],
        'rot': {'p1': 0, 'r1': 90, 'r2': 90},
    },
    'rxn_052': {
        'arrows': [
            {
                'desc': 'Bromide attack on bromonium ion',
                'src': ('r2', 0),
                'tgt': ('r1', 0),
                'type': 'attack',
            },
            {
                'bond': (0, 2),
                'desc': 'C-Br+ ring opens',
                'mol': 'r1',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_053': {
        'add_h': {'r2': 0},
        'arrows': [
            {
                'conn': 2,
                'desc': 'Propene C=C attacks HBr proton',
                'src': ('r1', (1, 2)),
                'tgt': ('r2', 1),
                'type': 'attack',
            },
            {
                'bond': (0, 1),
                'desc': 'H-Br bond cleaves to bromide ion',
                'mol': 'r2',
                'tgt': 0,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 180},
    },
    'rxn_054': {
        'arrows': [{
            'desc': 'Bromide capture of carbocation',
            'src': ('r2', 0),
            'tgt': ('r1', 1),
            'type': 'attack',
        }],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_055': {
        'add_h': {'r2': 3},
        'arrows': [
            {
                'conn': 0,
                'desc': 'Alkene attacks terminal O3',
                'src': ('r1', (0, 1)),
                'tgt': ('r2', 3),
                'type': 'attack',
            },
            {
                'b1': (2, 3),
                'b2': (1, 2),
                'desc': 'O-O pair forms C=O2',
                'mol': 'r2',
                'type': 'bond_to_bond',
            },
            {
                'b1': (0, 1),
                'b2': (0, 11),
                'desc': 'C=O0 pair abstracts H11',
                'mol': 'r2',
                'type': 'bond_to_bond',
            },
            {
                'conn': 3,
                'desc': 'O-H pair forms second epoxide bond',
                'src': ('r2', (3, 11)),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
        ],
        'rot': {'p1': 90, 'r1': 330, 'r2': 210},
    },
    'rxn_056': {
        'add_h': {'r2': 0},
        'arrows': [
            {
                'conn': 5,
                'desc': '1-Hexene C=C attacks borane',
                'src': ('r1', (4, 5)),
                'tgt': ('r2', 0),
                'type': 'attack',
            },
            {
                'conn': 1,
                'desc': 'B-H transfers hydride to C4',
                'src': ('r2', (0, 1)),
                'tgt': ('r1', 4),
                'type': 'attack',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_057': {
        'arrows': [
            {
                'conn': 2,
                'desc': 'Propene C=C attacks Hg(OAc)2',
                'src': ('r1', (1, 2)),
                'tgt': ('r2', 4),
                'type': 'attack',
            },
            {
                'desc': 'Hg back-donates to C1 to close mercurinium ring',
                'src': ('r2', 4),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (4, 3),
                'desc': 'Hg-OAc bond cleaves to acetate',
                'mol': 'r2',
                'tgt': 3,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 90, 'r2': 180},
    },
    'rxn_058': {
        'arrows': [
            {
                'desc': 'Water attacks mercurinium ion',
                'src': ('r2', 0),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 3),
                'desc': 'C-Hg+ ring opens',
                'mol': 'r1',
                'tgt': 3,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 180},
    },
    'rxn_059': {
        'add_h': {'r2': 0},
        'arrows': [
            {
                'conn': 2,
                'desc': 'Propyne C#C attacks hydronium proton',
                'src': ('r1', (1, 2)),
                'tgt': ('r2', 1),
                'type': 'attack',
            },
            {
                'bond': (0, 1),
                'desc': 'H-O+ bond cleaves to water',
                'mol': 'r2',
                'tgt': 0,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_060': {
        'arrows': [
            {
                'desc': 'Water attacks bromonium ion',
                'src': ('r2', 0),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 3),
                'desc': 'C-Br+ ring opens',
                'mol': 'r1',
                'tgt': 3,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 180},
    },
    'rxn_061': {
        'arrows': [
            {
                'conn': 0,
                'desc': 'Benzene C=C attacks nitronium',
                'src': ('r1', (0, 5)),
                'tgt': ('r2', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 0),
                'desc': 'N=O pi bond onto O',
                'mol': 'r2',
                'tgt': 0,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 90},
    },
    'rxn_062': {
        'add_h': {'r1': 3},
        'arrows': [
            {
                'desc': 'Bisulfate abstracts Wheland proton',
                'src': ('r2', 4),
                'tgt': ('r1', 9),
                'type': 'attack',
            },
            {
                'b1': (3, 9),
                'b2': (3, 4),
                'desc': 'C-H pair restores C3=C4',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'b1': (4, 5),
                'b2': (5, 6),
                'desc': 'C4=C5 pair restores C5=C6',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_063': {
        'arrows': [
            {
                'conn': 0,
                'desc': 'Benzene C=C attacks Br2',
                'src': ('r1', (0, 5)),
                'tgt': ('r2', 0),
                'type': 'attack',
            },
            {
                'bond': (0, 1),
                'desc': 'Br-Br bond cleaves to Br-',
                'mol': 'r2',
                'tgt': 1,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_064': {
        'arrows': [
            {
                'conn': 0,
                'desc': 'Benzene C=C attacks acylium carbon',
                'src': ('r1', (0, 5)),
                'tgt': ('r2', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 2),
                'desc': 'C#[O+] cleaves to neutral C=O',
                'mol': 'r2',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_065': {
        'arrows': [{
            'conn': 0,
            'desc': 'Benzene C=C attacks isopropyl cation',
            'src': ('r1', (0, 5)),
            'tgt': ('r2', 1),
            'type': 'attack',
        }],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_066': {
        'arrows': [
            {
                'conn': 0,
                'desc': 'Benzene C=C attacks S of SO3',
                'src': ('r1', (0, 5)),
                'tgt': ('r2', 1),
                'type': 'attack',
            },
            {
                'bond': (0, 1),
                'desc': 'S=O pi bond onto O',
                'mol': 'r2',
                'tgt': 0,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_067': {
        'arrows': [
            {
                'conn': 4,
                'desc': 'Aniline C3=C4 attacks Br2',
                'src': ('r1', (3, 4)),
                'tgt': ('r2', 0),
                'type': 'attack',
            },
            {
                'bond': (0, 1),
                'desc': 'Br-Br cleaves to Br-',
                'mol': 'r2',
                'tgt': 1,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_068': {
        'arrows': [
            {
                'atom': 0,
                'bond': (0, 1),
                'desc': 'Phenolate O- forms C=O',
                'mol': 'r1',
                'type': 'lone_pair_to_bond',
            },
            {
                'b1': (1, 2),
                'b2': (2, 3),
                'desc': 'C1=C2 shifts to C2=C3',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'conn': 4,
                'desc': 'C3=C4 attacks diazonium N7',
                'src': ('r1', (3, 4)),
                'tgt': ('r2', 7),
                'type': 'attack',
            },
            {
                'bond': (6, 7),
                'desc': 'N#N cleaves onto N6',
                'mol': 'r2',
                'tgt': 6,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_069': {
        'arrows': [
            {
                'conn': 1,
                'desc': 'Pyrrole C1=C2 attacks Vilsmeier C',
                'src': ('r1', (1, 2)),
                'tgt': ('r2', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 2),
                'desc': 'C=[N+] cleaves to neutral N',
                'mol': 'r2',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_070': {
        'arrows': [
            {
                'atom': 0,
                'bond': (0, 1),
                'desc': 'Phenolate O- forms C=O',
                'mol': 'r1',
                'type': 'lone_pair_to_bond',
            },
            {
                'conn': 2,
                'desc': 'Ortho C2 attacks CO2',
                'src': ('r1', (1, 2)),
                'tgt': ('r2', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 2),
                'desc': 'CO2 C=O cleaves to O-',
                'mol': 'r2',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_071': {
        'arrows': [
            {
                'atom': 2,
                'bond': (1, 2),
                'desc': 'O- forms C=O',
                'mol': 'r1',
                'type': 'lone_pair_to_bond',
            },
            {
                'conn': 0,
                'desc': 'Enolate C=C attacks acetaldehyde C',
                'src': ('r1', (0, 1)),
                'tgt': ('r2', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 2),
                'desc': 'C=O cleaves to O-',
                'mol': 'r2',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_072': {
        'arrows': [
            {
                'atom': 2,
                'bond': (1, 2),
                'desc': 'O- forms C=O',
                'mol': 'r1',
                'type': 'lone_pair_to_bond',
            },
            {
                'conn': 0,
                'desc': 'Enolate C=C attacks ester C',
                'src': ('r1', (0, 1)),
                'tgt': ('r2', 3),
                'type': 'attack',
            },
            {
                'bond': (3, 4),
                'desc': 'C=O cleaves to O-',
                'mol': 'r2',
                'tgt': 4,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_073': {
        'arrows': [
            {
                'atom': 10,
                'bond': (9, 10),
                'desc': 'Enolate O- forms C=O',
                'mol': 'r1',
                'type': 'lone_pair_to_bond',
            },
            {
                'conn': 8,
                'desc': 'Enolate C8 attacks ester C3',
                'src': ('r1', (8, 9)),
                'tgt': ('r1', 3),
                'type': 'attack',
            },
            {
                'bond': (3, 4),
                'desc': 'Ester C=O cleaves to O-',
                'mol': 'r1',
                'tgt': 4,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0},
    },
    'rxn_074': {
        'arrows': [
            {
                'desc': 'Malonate C5 attacks MVK C4',
                'src': ('r1', 5),
                'tgt': ('r2', 4),
                'type': 'attack',
            },
            {
                'b1': (3, 4),
                'b2': (1, 3),
                'desc': 'MVK C=C shifts to C=C enolate',
                'mol': 'r2',
                'type': 'bond_to_bond',
            },
            {
                'bond': (1, 2),
                'desc': 'MVK C=O cleaves to O-',
                'mol': 'r2',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_075': {
        'arrows': [
            {
                'atom': 2,
                'bond': (1, 2),
                'desc': 'Enolate O- forms C=O',
                'mol': 'r1',
                'type': 'lone_pair_to_bond',
            },
            {
                'conn': 0,
                'desc': 'Enolate C0 attacks cyclohexanone C6',
                'src': ('r1', (0, 1)),
                'tgt': ('r1', 6),
                'type': 'attack',
            },
            {
                'bond': (6, 7),
                'desc': 'Cyclohexanone C=O cleaves to O-',
                'mol': 'r1',
                'tgt': 7,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0},
    },
    'rxn_076': {
        'arrows': [
            {
                'desc': 'Malonate C5 attacks benzyl C1',
                'src': ('r1', 5),
                'tgt': ('r2', 1),
                'type': 'attack',
            },
            {
                'bond': (0, 1),
                'desc': 'C-Br cleaves to Br-',
                'mol': 'r2',
                'tgt': 0,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_077': {
        'arrows': [
            {
                'desc': 'Wittig ylide attacks benzaldehyde',
                'src': ('r2', 0),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 0),
                'desc': 'C=O pi bond onto O',
                'mol': 'r1',
                'tgt': 0,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 60, 'r2': 180},
    },
    'rxn_078': {
        'arrows': [
            {
                'b1': (7, 8),
                'b2': (6, 7),
                'desc': 'P-C bond to C=C alkene',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'b1': (6, 27),
                'b2': (8, 27),
                'desc': 'C-O bond to P=O phosphine oxide',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
        ],
        'rot': {'p1': 180, 'p2': 0, 'r1': 0},
        'use_manual_rot': True,
    },
    'rxn_079': {
        'arrows': [
            {
                'desc': 'HWE phosphonate enolate attacks cyclohexanone',
                'src': ('r2', 5),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 0),
                'desc': 'C=O pi bond onto O',
                'mol': 'r1',
                'tgt': 0,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 300, 'r1': 300, 'r2': 180},
    },
    'rxn_080': {
        'arrows': [
            {
                'atom': 2,
                'bond': (1, 2),
                'desc': 'N lone pair forms C=N+',
                'mol': 'r1',
                'type': 'lone_pair_to_bond',
            },
            {
                'conn': 0,
                'desc': 'Enamine C=C attacks MeI',
                'src': ('r1', (0, 1)),
                'tgt': ('r2', 0),
                'type': 'attack',
            },
            {
                'bond': (0, 1),
                'desc': 'C-I cleaves to I-',
                'mol': 'r2',
                'tgt': 1,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_081': {
        'arrows': [
            {
                'conn': 0,
                'desc': 'Diene C0 attacks C2',
                'src': ('r1', (0, 1)),
                'tgt': ('r2', 2),
                'type': 'attack',
            },
            {
                'conn': 3,
                'desc': 'Dienophile C3 attacks C3',
                'src': ('r2', (2, 3)),
                'tgt': ('r1', 3),
                'type': 'attack',
            },
            {
                'b1': (2, 3),
                'b2': (1, 2),
                'desc': 'Diene C=C shifts to central',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_082': {
        'arrows': [
            {
                'conn': 0,
                'desc': 'CP C0 attacks acrylate C0',
                'src': ('r1', (0, 1)),
                'tgt': ('r2', 0),
                'type': 'attack',
            },
            {
                'conn': 1,
                'desc': 'Acrylate C1 attacks CP C3',
                'src': ('r2', (0, 1)),
                'tgt': ('r1', 3),
                'type': 'attack',
            },
            {
                'b1': (2, 3),
                'b2': (1, 2),
                'desc': 'CP C=C shifts to central',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_083': {
        'arrows': [
            {
                'b1': (4, 5),
                'b2': (3, 4),
                'desc': 'C4-C5 pair forms C3=C4 diene',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'b1': (3, 7),
                'b2': (7, 8),
                'desc': 'C3-C7 pair forms C7=C8 diene',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'b1': (8, 9),
                'b2': (5, 9),
                'desc': 'C8=C9 pair shifts to C5=C9',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
        ],
        'rot': {'p1': 0, 'r1': 0},
    },
    'rxn_084': {
        'arrows': [
            {
                'conn': 0,
                'desc': 'C1=C2 pi to C1-C6 sigma',
                'lift': 0.6,
                'src': ('r1', (0, 1)),
                'tgt': ('r1', 5),
                'type': 'attack',
            },
            {
                'b1': (4, 5),
                'b2': (3, 4),
                'desc': 'C5=C6 pi to C4=C5 pi',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'b1': (2, 3),
                'b2': (1, 2),
                'desc': 'C3-C4 sigma to C2=C3 pi',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
        ],
        'coord_map': {
            'r1': {
                0: (0.75, 1.2),
                1: (1.4, 0.0),
                2: (0.75, -1.2),
                3: (-0.75, -1.2),
                4: (-1.4, 0.0),
                5: (-0.75, 1.2),
            }
        },
        'rot': {'p1': 0, 'r1': 0},
        'use_manual_rot': True,
    },
    'rxn_085': {
        'arrows': [
            {
                'conn': 0,
                'desc': 'C0=C1 attacks C7',
                'src': ('r1', (0, 1)),
                'tgt': ('r1', 7),
                'type': 'attack',
            },
            {
                'b1': (6, 7),
                'b2': (3, 6),
                'desc': 'C6=C7 shifts to C3=C6',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'b1': (2, 3),
                'b2': (1, 2),
                'desc': 'C2-C3 shifts to C1=C2',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
        ],
        'rot': {'p1': 0, 'r1': 0},
    },
    'rxn_086': {
        'arrows': [
            {
                'conn': 5,
                'desc': 'Ortho C5 attacks allyl C0',
                'src': ('r1', (4, 5)),
                'tgt': ('r1', 0),
                'type': 'attack',
            },
            {
                'b1': (0, 1),
                'b2': (1, 2),
                'desc': 'Allyl C=C shifts to C1=C2',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'b1': (2, 3),
                'b2': (3, 4),
                'desc': 'C-O shifts to form C=O',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
        ],
        'rot': {'p1': 0, 'r1': 0},
    },
    'rxn_087': {
        'arrows': [
            {
                'conn': 5,
                'desc': 'Ketene C=C attacks allyl C0',
                'src': ('r1', (4, 5)),
                'tgt': ('r1', 0),
                'type': 'attack',
            },
            {
                'b1': (0, 1),
                'b2': (1, 2),
                'desc': 'Allyl C=C shifts to C1=C2',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'b1': (2, 3),
                'b2': (3, 4),
                'desc': 'C-O shifts to form C=O',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
        ],
        'rot': {'p1': 0, 'r1': 0},
    },
    'rxn_088': {
        'arrows': [
            {
                'desc': 'Carbanion C4 attacks allyl C0',
                'src': ('r1', 4),
                'tgt': ('r1', 0),
                'type': 'attack',
            },
            {
                'b1': (0, 1),
                'b2': (1, 2),
                'desc': 'Allyl C=C shifts to C1=C2',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'bond': (2, 3),
                'desc': 'C-O cleaves to alkoxide O-',
                'mol': 'r1',
                'tgt': 3,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0},
    },
    'rxn_089': {
        'arrows': [
            {
                'desc': 'Diazomethane C0 attacks acrylate C0',
                'src': ('r2', 0),
                'tgt': ('r1', 0),
                'type': 'attack',
            },
            {
                'conn': 1,
                'desc': 'Acrylate C=C attacks terminal N2',
                'src': ('r1', (0, 1)),
                'tgt': ('r2', 2),
                'type': 'attack',
            },
            {
                'bond': (1, 2),
                'desc': 'N#N cleaves onto N1',
                'mol': 'r2',
                'tgt': 1,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_090': {
        'arrows': [
            {
                'desc': 'Azide N0 attacks alkyne C3',
                'src': ('r2', 0),
                'tgt': ('r1', 3),
                'type': 'attack',
            },
            {
                'conn': 4,
                'desc': 'Alkyne triple bond attacks N2',
                'src': ('r1', (3, 4)),
                'tgt': ('r2', 2),
                'type': 'attack',
            },
            {
                'bond': (1, 2),
                'desc': 'N=N double bond cleaves to N1',
                'mol': 'r2',
                'tgt': 1,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_091': {
        'add_h': {'r1': 7},
        'arrows': [
            {
                'desc': 'LDA deprotonation of kinetic alpha-H',
                'src': ('r2', 3),
                'tgt': ('r1', 8),
                'type': 'attack',
            },
            {
                'b1': (7, 8),
                'b2': (5, 7),
                'desc': 'C-H to C=C',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'bond': (5, 6),
                'desc': 'C=O to O',
                'mol': 'r1',
                'tgt': 6,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 270, 'r1': 270, 'r2': 0},
    },
    'rxn_092': {
        'arrows': [
            {
                'atom': 6,
                'bond': (5, 6),
                'desc': 'Enolate O- reforms C=O',
                'mol': 'r1',
                'type': 'lone_pair_to_bond',
            },
            {
                'conn': 7,
                'desc': 'Enolate C7 attacks methallyl C3',
                'src': ('r1', (5, 7)),
                'tgt': ('r2', 3),
                'type': 'attack',
            },
            {
                'bond': (3, 4),
                'desc': 'C-Br cleaves',
                'mol': 'r2',
                'tgt': 4,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 120, 'r1': 120, 'r2': 0},
    },
    'rxn_093': {
        'arrows': [
            {
                'conn': 0,
                'desc': 'Cu-C bond attacks enone beta-carbon',
                'src': ('r2', (0, 1)),
                'tgt': ('r1', 3),
                'type': 'attack',
            },
            {
                'b1': (2, 3),
                'b2': (1, 2),
                'desc': 'C=C pi-bond shifts to enolate',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'bond': (1, 0),
                'desc': 'C=O pi-bond onto O',
                'mol': 'r1',
                'tgt': 0,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 90, 'r1': 0, 'r2': 0},
    },
    'rxn_094': {
        'arrows': [
            {
                'desc': 'Hydroxyl oxygen attacks carbonyl carbon',
                'src': ('r1', 7),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 2),
                'desc': 'C=O pi-bond onto O',
                'mol': 'r1',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0},
    },
    'rxn_095': {
        'arrows': [
            {
                'b1': (1, 5),
                'b2': (3, 5),
                'desc': 'Vinyl group migrates to oxygen',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'bond': (3, 4),
                'desc': 'O-O cleaves expelling water',
                'mol': 'r1',
                'tgt': 4,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0},
    },
    'rxn_096': {
        'arrows': [
            {
                'desc': 'Peroxide attack on carbonyl',
                'src': ('r1', 6),
                'tgt': ('r1', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 2),
                'desc': 'C=O to O',
                'mol': 'r1',
                'tgt': 2,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 90, 'r1': 90},
    },
    'rxn_097': {
        'arrows': [
            {
                'conn': 0,
                'desc': 'Silyl enol ether attacks benzaldehyde',
                'src': ('r1', (0, 1)),
                'tgt': ('r2', 1),
                'type': 'attack',
            },
            {
                'bond': (1, 0),
                'desc': 'C=O pi bond onto O',
                'mol': 'r2',
                'tgt': 0,
                'type': 'bond_cleave',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 0},
    },
    'rxn_098': {
        'arrows': [
            {
                'conn': 0,
                'desc': 'Allyl C0 attacks vinyl C6',
                'src': ('r1', (0, 1)),
                'tgt': ('r1', 6),
                'type': 'attack',
            },
            {
                'b1': (5, 6),
                'b2': (3, 5),
                'desc': 'Vinyl C=C shifts to C=N',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
            {
                'b1': (2, 3),
                'b2': (1, 2),
                'desc': 'C-N bond shifts to C1=C2',
                'mol': 'r1',
                'type': 'bond_to_bond',
            },
        ],
        'rot': {'p1': 0, 'r1': 0},
    },
    'rxn_099': {
        'arrows': [
            {
                'conn': 6,
                'desc': 'Aromatic pi-bond attacks iminium carbon',
                'lift': 0.3,
                'src': ('r1', (6, 7)),
                'tgt': ('r1', 11),
                'type': 'attack',
            },
            {
                'bond': (10, 11),
                'desc': 'Iminium pi-bond cleaves to nitrogen',
                'mol': 'r1',
                'tgt': 10,
                'type': 'bond_cleave',
            },
        ],
        'coord_map': {
            'r1': {
                0: (0.7, 1.212),
                1: (-0.7, 1.212),
                2: (-1.4, 0.0),
                3: (-2.75, 0.0),
                4: (-0.7, -1.212),
                5: (-1.375, -2.381),
                6: (0.7, -1.212),
                7: (1.4, 0.0),
                8: (2.1, 0.188),
                9: (3.312, -0.512),
                10: (3.312, -1.912),
                11: (1.857, -2.591),
            }
        },
        'rot': {'p1': 0, 'r1': 0},
        'use_manual_rot': True,
    },
    'rxn_100': {
        'arrows': [
            {
                'conn': 6,
                'desc': 'Alkene attacks iodine',
                'lift': 0.3,
                'src': ('r1', (5, 6)),
                'tgt': ('r2', 1),
                'type': 'attack',
            },
            {
                'bond': (0, 1),
                'desc': 'I-I bond cleaves to iodide',
                'mol': 'r2',
                'tgt': 0,
                'type': 'bond_cleave',
            },
            {
                'desc': 'Carboxylate oxygen attacks alkene C5',
                'lift': 0.3,
                'src': ('r1', 2),
                'tgt': ('r1', 5),
                'type': 'attack',
            },
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 90},
        'use_manual_rot': True,
    },
}
