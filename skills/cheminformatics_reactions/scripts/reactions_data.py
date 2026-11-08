# Auto-generated reactions database with verified cubic Bezier ChemDraw arrows

REACTIONS = [
    {
        'id': 'rxn_001',
        'title': 'Kinetic Deprotonation of Ketone by LDA',
        'category': 'Acid-Base & Enolization',
        'desc': (
            'Sterically hindered LDA selectively abstracts the less hindered'
            ' alpha-proton at -78 C to afford the kinetic enolate.'
        ),
        'r1': {
            'smi': 'CC1CCCC(=O)C1',
            'name': '(3R)-3-methylcyclohexan-1-one',
            'rot': 120,
            'w': 4.5,
            'add_h': 7,
        },
        'r2': {'smi': 'CC(C)[N-]C(C)C', 'name': 'LDA', 'rot': 0, 'w': 3.8},
        'cond': {
            'top': '$\\text{THF}, -78^\\circ\\text{C}$',
            'bottom': '$-\\text{HN}(i\\text{-Pr})_2$',
        },
        'p1': {
            'smi': 'CC1CCCC([O-])=C1',
            'name': 'Lithium Enolate',
            'rot': 120,
            'w': 4.5,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (5.35, 0.99) to[out=140, in=30] (1.65, 0.63);',
            '\\draw[chemarrow] (0.82, 0.33) to[out=135, in=55] (0.03, 0.35);',
            (
                '\\draw[chemarrow] (-0.68, 1.03) to[out=120, in=205] (-0.78,'
                ' 1.65);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (4.522, 0.196) .. controls (3.560, 1.943)'
                ' and (2.976, 2.019) .. (1.602, 0.572);'
            ),
            (
                '\\draw[chemarrow] (0.855, 0.266) .. controls (0.788, 0.596)'
                ' and (0.074, 0.569) .. (-0.031, 0.249);'
            ),
            (
                '\\draw[chemarrow] (-0.702, 0.970) .. controls (-0.871, 1.045)'
                ' and (-0.891, 1.421) .. (-0.722, 1.496);'
            ),
        ],
        'rot': {'r1': 240, 'p1': 240, 'r2': 270},
    },
    {
        'id': 'rxn_002',
        'title': 'Thermodynamic Enolate Formation (NaOMe / MeOH)',
        'category': 'Acid-Base & Enolization',
        'desc': (
            'Under equilibrating conditions (protic solvent, room temp),'
            ' deprotonation yields the more substituted, thermodynamically'
            ' stable enolate.'
        ),
        'r1': {
            'smi': 'CC1C(=O)CCCC1',
            'name': '2-methylcyclohexan-1-one',
            'rot': 180,
            'w': 4.5,
            'add_h': 1,
        },
        'r2': {'smi': 'C[O-]', 'name': 'Methoxide', 'rot': 180, 'w': 2.5},
        'cond': {
            'top': '$\\text{MeOH}, 25^\\circ\\text{C}$',
            'bottom': '$-\\text{MeOH}$',
        },
        'p1': {
            'smi': 'CC1=C([O-])CCCC1',
            'name': 'Thermodynamic Enolate',
            'rot': 60,
            'w': 4.5,
        },
        'p2': None,
        'arrows': [
            (
                '\\draw[chemarrow] (3.75, 0.32) to[out=140, in=-132] (-1.64,'
                ' -1.81);'
            ),
            '\\draw[chemarrow] (0.90, 1.42) to[out=120, in=205] (1.43, 1.81);',
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.504, 0.187) .. controls (2.963, 1.092)'
                ' and (2.660, 1.114) .. (1.991, 0.299);'
            ),
            (
                '\\draw[chemarrow] (1.108, 0.286) .. controls (1.124, 0.004)'
                ' and (0.761, -0.211) .. (0.491, -0.128);'
            ),
            (
                '\\draw[chemarrow] (0.754, -1.271) .. controls (0.710, -1.493)'
                ' and (1.086, -1.762) .. (1.281, -1.646);'
            ),
        ],
        'rot': {'r1': 180, 'p1': 150, 'r2': 180},
    },
    {
        'id': 'rxn_003',
        'title': 'Deprotonation of Carboxylic Acid by Hydroxide',
        'category': 'Acid-Base & Enolization',
        'desc': (
            'Fast diffusion-controlled proton transfer from acetic acid to'
            ' hydroxide ion generating acetate and water.'
        ),
        'r1': {'smi': 'CC(=O)O', 'name': 'Acetic acid', 'rot': 120, 'w': 3.6},
        'r2': {'smi': '[OH-]', 'name': 'Hydroxide', 'rot': 0, 'w': 2.0},
        'cond': {
            'top': '$\\text{H}_2\\text{O}, 25^\\circ\\text{C}$',
            'bottom': '$\\text{fast}$',
        },
        'p1': {
            'smi': 'CC(=O)[O-]',
            'name': 'Acetate anion',
            'rot': 120,
            'w': 3.6,
        },
        'p2': {'smi': 'O', 'name': 'Water', 'rot': 0, 'w': 2.0},
        'arrows': [
            '\\draw[chemarrow] (4.04, 0.32) to[out=140, in=90] (0.07, 1.49);',
            (
                '\\draw[chemarrow] (-0.07, 0.39) to[out=120, in=205] (-0.17,'
                ' 1.24);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.797, 0.192) .. controls (3.097, 1.536)'
                ' and (2.682, 1.579) .. (1.723, 0.405);'
            ),
            (
                '\\draw[chemarrow] (0.991, 0.501) .. controls (0.825, 0.381)'
                ' and (0.418, 0.497) .. (0.440, 0.700);'
            ),
        ],
        'rot': {'r1': 0, 'p1': 210, 'r2': 0, 'p2': 0},
    },
    {
        'id': 'rxn_004',
        'title': 'Amine Protonation by Hydronium Ion',
        'category': 'Acid-Base & Enolization',
        'desc': (
            'Nucleophilic lone pair on triethylamine abstracts a proton from'
            ' hydronium yielding triethylammonium.'
        ),
        'r1': {
            'smi': 'CCN(CC)CC',
            'name': 'Triethylamine',
            'rot': 90,
            'w': 4.2,
        },
        'r2': {'smi': '[OH3+]', 'name': 'Hydronium', 'rot': 0, 'w': 2.2},
        'cond': {
            'top': '$\\text{H}_2\\text{O}$',
            'bottom': '$-\\text{H}_2\\text{O}$',
        },
        'p1': {
            'smi': 'CC[NH+](CC)CC',
            'name': 'Triethylammonium',
            'rot': 90,
            'w': 4.2,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (0.29, 0.32) to[out=40, in=158] (4.13, 0.19);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (0.678, 0.133) .. controls (1.778, 0.838)'
                ' and (2.333, 0.846) .. (3.453, 0.172);'
            ),
            (
                '\\draw[chemarrow] (4.196, -0.035) .. controls (4.276, -0.216)'
                ' and (4.656, -0.373) .. (4.736, -0.193);'
            ),
        ],
        'rot': {'r1': 300, 'p1': 300, 'r2': 180},
    },
    {
        'id': 'rxn_005',
        'title': 'Terminal Alkyne Deprotonation by Sodium Amide',
        'category': 'Acid-Base & Enolization',
        'desc': (
            'Strongly basic amide anion abstracts the sp-hybridized terminal'
            ' alkyne proton (pKa ~ 25) generating an acetylide.'
        ),
        'r1': {'smi': 'CCCC#C', 'name': '1-Pentyne', 'rot': 0, 'w': 4.5},
        'r2': {'smi': '[NH2-]', 'name': 'Amide ion', 'rot': 0, 'w': 2.2},
        'cond': {
            'top': '$\\text{liq. NH}_3, -33^\\circ\\text{C}$',
            'bottom': '$-\\text{NH}_3$',
        },
        'p1': {
            'smi': 'CCCC#[C-]',
            'name': 'Pentynide anion',
            'rot': 0,
            'w': 4.5,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (4.41, 0.36) to[out=140, in=17] (2.21, 0.30);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (4.155, 0.223) .. controls (3.418, 1.528)'
                ' and (3.017, 1.548) .. (2.150, 0.326);'
            ),
            (
                '\\draw[chemarrow] (1.566, 0.051) .. controls (1.536, -0.095)'
                ' and (1.244, -0.188) .. (1.156, -0.068);'
            ),
        ],
        'rot': {'r1': 0, 'p1': 0, 'r2': 0},
    },
    {
        'id': 'rxn_006',
        'title': 'Alkoxide Protonation by Water',
        'category': 'Acid-Base & Enolization',
        'desc': (
            'Quenching of sodium ethoxide with water yields ethanol and'
            ' regenerates hydroxide.'
        ),
        'r1': {'smi': 'CC[O-]', 'name': 'Ethoxide', 'rot': 0, 'w': 3.0},
        'r2': {'smi': 'O', 'name': 'Water', 'rot': 180, 'w': 2.2},
        'cond': {
            'top': '$\\text{H}_2\\text{O}$',
            'bottom': '$\\text{quenching}$',
        },
        'p1': {'smi': 'CCO', 'name': 'Ethanol', 'rot': 0, 'w': 3.0},
        'p2': {'smi': '[OH-]', 'name': 'Hydroxide', 'rot': 0, 'w': 2.0},
        'arrows': [
            '\\draw[chemarrow] (1.06, 0.05) to[out=40, in=158] (3.64, 0.16);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (1.294, -0.079) .. controls (1.806, 0.850)'
                ' and (2.112, 0.886) .. (2.825, 0.101);'
            ),
            (
                '\\draw[chemarrow] (3.645, -0.079) .. controls (3.737, -0.286)'
                ' and (4.194, -0.389) .. (4.286, -0.181);'
            ),
        ],
        'rot': {'r1': 0, 'p1': 0, 'r2': 180, 'p2': 180},
    },
    {
        'id': 'rxn_007',
        'title': 'Acid-Catalyzed Enolization (Proton Loss)',
        'category': 'Acid-Base & Enolization',
        'desc': (
            'Protonated ketone loses an alpha-proton to water, establishing'
            ' keto-enol equilibrium.'
        ),
        'r1': {
            'smi': 'CC(=[OH+])C',
            'name': 'Protonated acetone',
            'rot': 180,
            'w': 3.8,
        },
        'r2': {'smi': 'O', 'name': 'Water', 'rot': 180, 'w': 2.2},
        'cond': {
            'top': '$\\text{H}_3\\text{O}^+$',
            'bottom': '$-\\text{H}_3\\text{O}^+$',
        },
        'p1': {
            'smi': 'CC(O)=C',
            'name': 'Prop-1-en-2-ol (enol)',
            'rot': 180,
            'w': 3.8,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (4.34, 0.36) to[out=140, in=30] (1.75, 1.53);',
            (
                '\\draw[chemarrow] (0.14, -0.36) to[out=120, in=205] (0.30,'
                ' -1.34);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (4.109, 0.140) .. controls (2.921, 1.581)'
                ' and (2.453, 1.507) .. (1.769, -0.232);'
            ),
            (
                '\\draw[chemarrow] (0.993, -0.640) .. controls (0.915, -0.259)'
                ' and (0.091, -0.287) .. (-0.030, -0.657);'
            ),
            (
                '\\draw[chemarrow] (-0.781, 0.175) .. controls (-0.971, 0.259)'
                ' and (-1.034, 0.679) .. (-0.844, 0.764);'
            ),
        ],
        'rot': {'r1': 0, 'p1': 0, 'r2': 0},
    },
    {
        'id': 'rxn_008',
        'title': 'Base-Catalyzed Enolate C-Protonation',
        'category': 'Acid-Base & Enolization',
        'desc': (
            'Enolate resonance form protonates at carbon from solvent to'
            ' regenerate the neutral ketone.'
        ),
        'r1': {
            'smi': 'C=C([O-])C',
            'name': 'Acetone enolate',
            'rot': 0,
            'w': 3.8,
        },
        'r2': {'smi': 'O', 'name': 'Water', 'rot': 180, 'w': 2.2},
        'cond': {'top': '$\\text{H}_2\\text{O}$', 'bottom': '$-\\text{OH}^-$'},
        'p1': {'smi': 'CC(=O)C', 'name': 'Acetone', 'rot': 0, 'w': 3.8},
        'p2': None,
        'arrows': [
            (
                '\\draw[chemarrow] (-0.29, 1.27) to[out=120, in=190] (-0.11,'
                ' 0.41);'
            ),
            '\\draw[chemarrow] (-1.45, -0.99) to[out=40, in=158] (4.04, 0.16);',
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (1.496, -0.926) .. controls (1.652, -0.656)'
                ' and (0.791, -0.437) .. (0.636, -0.706);'
            ),
            (
                '\\draw[chemarrow] (-0.063, 0.543) .. controls (1.491, 2.332)'
                ' and (2.156, 2.255) .. (3.261, 0.160);'
            ),
            (
                '\\draw[chemarrow] (4.045, -0.079) .. controls (4.137, -0.286)'
                ' and (4.594, -0.389) .. (4.686, -0.181);'
            ),
        ],
        'rot': {'r1': 120, 'p1': 120, 'r2': 180},
    },
    {
        'id': 'rxn_009',
        'title': 'Diethyl Malonate Deprotonation by Ethoxide',
        'category': 'Acid-Base & Enolization',
        'desc': (
            'Ethoxide deprotonates the doubly-activated methylene position (pKa'
            ' ~ 13) to form a resonance-stabilized carbanion.'
        ),
        'r1': {
            'smi': 'CCOC(=O)CC(=O)OCC',
            'name': 'Diethyl malonate',
            'rot': 180,
            'w': 5.2,
        },
        'r2': {
            'smi': 'CC[O-]',
            'name': 'Sodium ethoxide',
            'rot': 180,
            'w': 3.0,
        },
        'cond': {'top': '$\\text{EtOH}$', 'bottom': '$-\\text{EtOH}$'},
        'p1': {
            'smi': 'CCOC(=O)[CH-]C(=O)OCC',
            'name': 'Malonate carbanion',
            'rot': 180,
            'w': 5.2,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (4.10, 0.62) to[out=140, in=90] (0.00, 0.75);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.865, 0.502) .. controls (3.218, 1.875)'
                ' and (2.803, 1.936) .. (1.790, 0.805);'
            ),
            (
                '\\draw[chemarrow] (1.233, 0.580) .. controls (1.199, 0.455)'
                ' and (0.943, 0.386) .. (0.873, 0.495);'
            ),
        ],
        'rot': {'r1': 240, 'p1': 240, 'r2': 180},
    },
    {
        'id': 'rxn_010',
        'title': 'Nitroalkane Deprotonation (Nitronate Formation)',
        'category': 'Acid-Base & Enolization',
        'desc': (
            'Deprotonation of nitromethane by hydroxide generates the'
            ' resonance-delocalized nitronate ion.'
        ),
        'r1': {
            'smi': 'C[N+](=O)[O-]',
            'name': 'Nitromethane',
            'rot': 180,
            'w': 3.6,
        },
        'r2': {'smi': '[OH-]', 'name': 'Hydroxide', 'rot': 0, 'w': 2.0},
        'cond': {
            'top': '$\\text{H}_2\\text{O}$',
            'bottom': '$-\\text{H}_2\\text{O}$',
        },
        'p1': {
            'smi': 'C=[N+]([O-])[O-]',
            'name': 'Nitronate anion',
            'rot': 180,
            'w': 3.6,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (4.04, 0.32) to[out=140, in=30] (1.62, 1.29);',
            (
                '\\draw[chemarrow] (0.21, -0.39) to[out=120, in=205] (0.31,'
                ' -1.23);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.796, 0.191) .. controls (3.087, 1.531)'
                ' and (2.672, 1.571) .. (1.718, 0.393);'
            ),
            (
                '\\draw[chemarrow] (0.958, 0.504) .. controls (0.888, 0.162)'
                ' and (0.148, 0.189) .. (0.040, 0.521);'
            ),
            (
                '\\draw[chemarrow] (-1.012, 0.400) .. controls (-1.167, 0.288)'
                ' and (-1.514, 0.469) .. (-1.494, 0.659);'
            ),
        ],
        'rot': {'r1': 180, 'p1': 180, 'r2': 0},
    },
    {
        'id': 'rxn_011',
        'title': 'Grignard Addition of Methylmagnesium Bromide to Acetone',
        'category': 'Carbonyl Addition',
        'desc': (
            'Nucleophilic carbanion of MeMgBr attacks the carbonyl carbon,'
            ' forming a tetrahedral magnesium alkoxide.'
        ),
        'r1': {'smi': 'CC(=O)C', 'name': 'Acetone', 'rot': 0, 'w': 3.6},
        'r2': {
            'smi': 'C[Mg]Br',
            'name': 'Methylmagnesium bromide',
            'rot': 0,
            'w': 3.6,
        },
        'cond': {
            'top': '$\\text{Et}_2\\text{O}, 0^\\circ\\text{C}$',
            'bottom': '$\\text{addition}$',
        },
        'p1': {
            'smi': 'CC(C)(C)[O-]',
            'name': 'tert-Butoxide adduct',
            'rot': 0,
            'w': 3.8,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (3.23, 0.37) to[out=140, in=30] (0.21, -0.38);',
            (
                '\\draw[chemarrow] (-0.14, 0.35) to[out=120, in=205] (-0.24,'
                ' 1.27);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.875, 0.097) .. controls (2.565, 2.050)'
                ' and (1.912, 2.048) .. (0.611, 0.091);'
            ),
            (
                '\\draw[chemarrow] (-0.254, -0.200) .. controls (-0.368,'
                ' -0.456) and (-0.941, -0.439) .. (-1.054, -0.183);'
            ),
        ],
        'rot': {'r1': 270, 'p1': 300, 'r2': 0},
    },
    {
        'id': 'rxn_012',
        'title': 'Hydride Transfer from NaBH4 to Benzaldehyde',
        'category': 'Carbonyl Addition',
        'desc': (
            'Borohydride delivers a nucleophilic hydride to the carbonyl'
            ' carbon, generating a benzyl alkoxide intermediate.'
        ),
        'r1': {
            'smi': 'O=Cc1ccccc1',
            'name': 'Benzaldehyde',
            'rot': 60,
            'w': 4.2,
        },
        'r2': {'smi': '[BH4-]', 'name': 'Borohydride', 'rot': 0, 'w': 2.2},
        'cond': {
            'top': '$\\text{MeOH}, 0^\\circ\\text{C}$',
            'bottom': '$\\text{reduction}$',
        },
        'p1': {
            'smi': '[O-]Cc1ccccc1',
            'name': 'Alkoxide intermediate',
            'rot': 60,
            'w': 4.2,
        },
        'p2': None,
        'arrows': [
            (
                '\\draw[chemarrow] (4.26, 0.36) to[out=140, in=-47] (-0.07,'
                ' -1.49);'
            ),
            (
                '\\draw[chemarrow] (-0.70, -1.57) to[out=120, in=205] (-1.11,'
                ' -1.76);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (4.448, 0.077) .. controls (2.804, 1.767)'
                ' and (2.155, 1.650) .. (1.203, -0.507);'
            ),
            (
                '\\draw[chemarrow] (0.745, -0.954) .. controls (0.596, -0.971)'
                ' and (0.540, -1.273) .. (0.653, -1.372);'
            ),
        ],
        'rot': {'r1': 240, 'p1': 240, 'r2': 0},
    },
    {
        'id': 'rxn_013',
        'title': 'Cyanohydrin Formation (Cyanide Addition to Cyclohexanone)',
        'category': 'Carbonyl Addition',
        'desc': (
            'Cyanide ion attacks the ketone carbonyl carbon forming the'
            ' cyanohydrin tetrahedral alkoxide.'
        ),
        'r1': {
            'smi': 'O=C1CCCCC1',
            'name': 'Cyclohexanone',
            'rot': 90,
            'w': 4.2,
        },
        'r2': {'smi': '[C-]#N', 'name': 'Cyanide ion', 'rot': 180, 'w': 2.2},
        'cond': {
            'top': '$\\text{HCN} / \\text{KCN}$',
            'bottom': '$\\text{reversible}$',
        },
        'p1': {
            'smi': '[O-]C1(C#N)CCCCC1',
            'name': 'Cyanohydrin adduct',
            'rot': 90,
            'w': 4.5,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (3.59, 0.32) to[out=140, in=30] (0.21, 0.56);',
            '\\draw[chemarrow] (0.14, 0.96) to[out=120, in=205] (0.24, 1.40);',
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.345, 0.192) .. controls (2.263, 2.265)'
                ' and (1.625, 2.329) .. (0.154, 0.512);'
            ),
            (
                '\\draw[chemarrow] (-0.157, 0.905) .. controls (-0.312, 0.975)'
                ' and (-0.338, 1.322) .. (-0.183, 1.391);'
            ),
        ],
        'rot': {'r1': 270, 'p1': 240, 'r2': 180},
    },
    {
        'id': 'rxn_014',
        'title': 'Hemiacetal Formation: Methanol Addition to Acetaldehyde',
        'category': 'Carbonyl Addition',
        'desc': (
            'Nucleophilic oxygen of methanol attacks protonated acetaldehyde to'
            ' form the oxonium hemiacetal intermediate.'
        ),
        'r1': {
            'smi': 'CC=[OH+]',
            'name': 'Protonated acetaldehyde',
            'rot': 0,
            'w': 2.8,
        },
        'r2': {'smi': 'CO', 'name': 'Methanol', 'rot': 180, 'w': 2.2},
        'cond': {'top': '$\\text{H}^+$ cat.', 'bottom': '$\\text{addition}$'},
        'p1': {
            'smi': 'CC(O)[OH+]C',
            'name': 'Oxonium hemiacetal',
            'rot': 0,
            'w': 4.2,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (3.53, 0.32) to[out=140, in=90] (-0.27, 0.72);',
            '\\draw[chemarrow] (0.42, 0.24) to[out=120, in=205] (1.19, -0.02);',
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (2.898, 0.101) .. controls (1.902, 1.222)'
                ' and (1.505, 1.158) .. (0.911, -0.219);'
            ),
            (
                '\\draw[chemarrow] (0.288, -0.703) .. controls (0.312, -0.934)'
                ' and (-0.004, -1.299) .. (-0.192, -1.162);'
            ),
        ],
        'rot': {'r1': 120, 'p1': 120, 'r2': 180},
    },
    {
        'id': 'rxn_015',
        'title': 'Imine Formation: Methylamine Addition to Cyclohexanone',
        'category': 'Carbonyl Addition',
        'desc': (
            'Primary amine nucleophile attacks cyclohexanone carbonyl yielding'
            ' the tetrahedral carbinolamine intermediate.'
        ),
        'r1': {
            'smi': 'O=C1CCCCC1',
            'name': 'Cyclohexanone',
            'rot': 90,
            'w': 4.2,
        },
        'r2': {'smi': 'CN', 'name': 'Methylamine', 'rot': 180, 'w': 2.5},
        'cond': {
            'top': '$\\text{pH } 4.5$',
            'bottom': '$\\text{carbinolamine}$',
        },
        'p1': {
            'smi': 'C[NH2+]C1([O-])CCCCC1',
            'name': 'Carbinolamine adduct',
            'rot': 90,
            'w': 4.5,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (4.02, 0.37) to[out=140, in=30] (0.21, 0.56);',
            '\\draw[chemarrow] (0.14, 0.96) to[out=120, in=205] (0.24, 1.40);',
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.819, 0.201) .. controls (2.520, 2.274)'
                ' and (1.786, 2.337) .. (0.152, 0.516);'
            ),
            (
                '\\draw[chemarrow] (-0.157, 0.905) .. controls (-0.312, 0.975)'
                ' and (-0.338, 1.322) .. (-0.183, 1.391);'
            ),
        ],
        'rot': {'r1': 270, 'p1': 270, 'r2': 180},
    },
    {
        'id': 'rxn_016',
        'title': 'Enamine Formation: Pyrrolidine Addition to Cyclopentanone',
        'category': 'Carbonyl Addition',
        'desc': (
            'Secondary amine adds to cyclopentanone giving the zwitterionic'
            ' carbinolamine intermediate.'
        ),
        'r1': {
            'smi': 'O=C1CCCC1',
            'name': 'Cyclopentanone',
            'rot': 90,
            'w': 4.0,
        },
        'r2': {'smi': 'C1CCNC1', 'name': 'Pyrrolidine', 'rot': 180, 'w': 3.2},
        'cond': {
            'top': '$p\\text{-TsOH}, \\text{benzene}$',
            'bottom': '$-\\text{H}_2\\text{O}$',
        },
        'p1': {
            'smi': '[O-]C1([NH+]2CCCC2)CCCC1',
            'name': 'Adduct intermediate',
            'rot': 90,
            'w': 4.6,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (5.78, 1.07) to[out=140, in=27] (0.21, 0.34);',
            '\\draw[chemarrow] (0.14, 0.80) to[out=120, in=205] (0.24, 1.31);',
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.595, 0.540) .. controls (2.091, 2.405)'
                ' and (1.400, 2.361) .. (0.140, 0.324);'
            ),
            (
                '\\draw[chemarrow] (-0.167, 0.746) .. controls (-0.343, 0.824)'
                ' and (-0.359, 1.219) .. (-0.183, 1.297);'
            ),
        ],
        'rot': {'r1': 270, 'p1': 270, 'r2': 60},
    },
    {
        'id': 'rxn_017',
        'title': 'Bisulfite Addition to Benzaldehyde',
        'category': 'Carbonyl Addition',
        'desc': (
            'Bisulfite sulfur attacks the aldehyde carbon to form a crystalline'
            ' bisulfite adduct.'
        ),
        'r1': {
            'smi': 'O=Cc1ccccc1',
            'name': 'Benzaldehyde',
            'rot': 60,
            'w': 4.2,
        },
        'r2': {
            'smi': 'O=S([O-])O',
            'name': 'Bisulfite ion',
            'rot': 180,
            'w': 3.0,
        },
        'cond': {
            'top': '$\\text{H}_2\\text{O}, 0^\\circ\\text{C}$',
            'bottom': '$\\text{precipitation}$',
        },
        'p1': {
            'smi': 'O=[S+]([O-])(O)C([O-])c1ccccc1',
            'name': 'Bisulfite adduct',
            'rot': 60,
            'w': 4.8,
        },
        'p2': None,
        'arrows': [
            (
                '\\draw[chemarrow] (4.80, 0.64) to[out=140, in=-47] (-0.07,'
                ' -1.49);'
            ),
            (
                '\\draw[chemarrow] (-0.70, -1.57) to[out=120, in=205] (-1.11,'
                ' -1.76);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (4.167, 0.090) .. controls (2.663, 1.429)'
                ' and (2.071, 1.310) .. (1.204, -0.508);'
            ),
            (
                '\\draw[chemarrow] (0.745, -0.954) .. controls (0.596, -0.971)'
                ' and (0.540, -1.273) .. (0.653, -1.372);'
            ),
        ],
        'rot': {'r1': 240, 'p1': 30, 'r2': 90},
    },
    {
        'id': 'rxn_018',
        'title': 'Organolithium 1,2-Addition to Enone',
        'category': 'Carbonyl Addition',
        'desc': (
            'Hard organolithium nucleophile (MeLi) attacks the carbonyl carbon'
            ' selectively over conjugate 1,4-addition.'
        ),
        'r1': {
            'smi': 'CC(=O)C=CC',
            'name': 'Pent-3-en-2-one',
            'rot': 0,
            'w': 4.5,
        },
        'r2': {'smi': 'C[Li]', 'name': 'Methyllithium', 'rot': 0, 'w': 2.2},
        'cond': {
            'top': '$\\text{THF}, -78^\\circ\\text{C}$',
            'bottom': '$1,2\\text{-addition}$',
        },
        'p1': {
            'smi': 'CC([O-])(C)C=CC',
            'name': '1,2-adduct alkoxide',
            'rot': 0,
            'w': 4.6,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (3.64, 0.32) to[out=140, in=-22] (1.32, 0.22);',
            (
                '\\draw[chemarrow] (1.32, -0.26) to[out=120, in=205] (1.51,'
                ' -0.90);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (4.421, 0.029) .. controls (2.471, 1.264)'
                ' and (1.662, 1.094) .. (0.376, -0.822);'
            ),
            (
                '\\draw[chemarrow] (-0.175, -1.159) .. controls (-0.220,'
                ' -1.320) and (-0.553, -1.394) .. (-0.642, -1.253);'
            ),
        ],
        'rot': {'r1': 90, 'p1': 90, 'r2': 0},
    },
    {
        'id': 'rxn_019',
        'title': 'Reformatsky Reagent Addition to Acetophenone',
        'category': 'Carbonyl Addition',
        'desc': (
            'Organozinc enolate (Reformatsky reagent) nucleophilically attacks'
            ' acetophenone to form a beta-hydroxy ester precursor.'
        ),
        'r1': {
            'smi': 'CC(=O)c1ccccc1',
            'name': 'Acetophenone',
            'rot': 60,
            'w': 4.2,
        },
        'r2': {
            'smi': 'CCOC(=O)C[Zn]Br',
            'name': 'Reformatsky reagent',
            'rot': 180,
            'w': 4.5,
        },
        'cond': {
            'top': '$\\text{THF}, \\text{reflux}$',
            'bottom': '$\\text{addition}$',
        },
        'p1': {
            'smi': 'CCOC(=O)CC(C)([O-])c1ccccc1',
            'name': 'Zinc alkoxide',
            'rot': 60,
            'w': 5.2,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (5.08, 0.96) to[out=140, in=60] (0.57, 0.92);',
            '\\draw[chemarrow] (0.35, 1.17) to[out=120, in=205] (0.17, 1.67);',
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (4.492, -0.104) .. controls (3.128, 2.876)'
                ' and (2.180, 3.071) .. (-0.248, 0.869);'
            ),
            (
                '\\draw[chemarrow] (-0.828, 0.660) .. controls (-0.890, 0.521)'
                ' and (-1.199, 0.488) .. (-1.261, 0.627);'
            ),
        ],
        'rot': {'r1': 240, 'p1': 90, 'r2': 150},
    },
    {
        'id': 'rxn_020',
        'title': 'Mukaiyama Aldol: Silyl Enol Ether Attack on Benzaldehyde',
        'category': 'Carbonyl Addition',
        'desc': (
            'Lewis-acid activated benzaldehyde is attacked by trimethylsilyl'
            ' enol ether generating the aldol precursor.'
        ),
        'r1': {
            'smi': 'O=Cc1ccccc1',
            'name': 'Benzaldehyde (BF3 adduct)',
            'rot': 60,
            'w': 4.2,
        },
        'r2': {
            'smi': 'C=C(O[Si](C)(C)C)c1ccccc1',
            'name': 'Silyl enol ether',
            'rot': 180,
            'w': 4.8,
        },
        'cond': {
            'top': '$\\text{BF}_3\\cdot\\text{OEt}_2, -78^\\circ\\text{C}$',
            'bottom': '$\\text{addition}$',
        },
        'p1': {
            'smi': 'C[Si](C)(C)O[C+](CC([O-])c1ccccc1)c1ccccc1',
            'name': 'Mukaiyama aldol adduct',
            'rot': 60,
            'w': 5.5,
        },
        'p2': None,
        'arrows': [
            (
                '\\draw[chemarrow] (5.35, 1.26) to[out=140, in=-47] (-0.07,'
                ' -1.49);'
            ),
            (
                '\\draw[chemarrow] (-0.70, -1.57) to[out=120, in=205] (-1.11,'
                ' -1.76);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (4.781, 0.114) .. controls (3.015, 1.791)'
                ' and (2.300, 1.667) .. (1.204, -0.508);'
            ),
            (
                '\\draw[chemarrow] (0.745, -0.954) .. controls (0.596, -0.971)'
                ' and (0.540, -1.273) .. (0.653, -1.372);'
            ),
        ],
        'rot': {'r1': 240, 'p1': 0, 'r2': 120},
    },
    {
        'id': 'rxn_021',
        'title': 'Tetrahedral Intermediate Formation in Ester Saponification',
        'category': 'Acyl Substitution',
        'desc': (
            'Hydroxide nucleophile attacks the carbonyl of ethyl acetate to'
            ' generate the tetrahedral intermediate.'
        ),
        'r1': {
            'smi': 'CCOC(=O)C',
            'name': 'Ethyl acetate',
            'rot': 60,
            'w': 4.0,
        },
        'r2': {'smi': '[OH-]', 'name': 'Hydroxide', 'rot': 0, 'w': 2.0},
        'cond': {
            'top': '$\\text{H}_2\\text{O}, 25^\\circ\\text{C}$',
            'bottom': '$\\text{step 1: addition}$',
        },
        'p1': {
            'smi': 'CCOC([O-])(O)C',
            'name': 'Tetrahedral intermediate',
            'rot': 60,
            'w': 4.2,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (4.24, 0.32) to[out=140, in=-38] (0.33, 0.64);',
            (
                '\\draw[chemarrow] (-0.16, 1.20) to[out=120, in=205] (-0.55,'
                ' 1.63);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (4.014, 0.213) .. controls (2.923, 2.406)'
                ' and (2.197, 2.542) .. (0.385, 0.892);'
            ),
            (
                '\\draw[chemarrow] (-0.196, 0.747) .. controls (-0.273, 0.624)'
                ' and (-0.571, 0.631) .. (-0.610, 0.770);'
            ),
        ],
        'rot': {'r1': 270, 'p1': 270, 'r2': 0},
    },
    {
        'id': 'rxn_022',
        'title': 'Tetrahedral Intermediate Collapse (Ethoxide Expulsion)',
        'category': 'Acyl Substitution',
        'desc': (
            'Alkoxide lone pair reforms the carbonyl pi-bond, expelling'
            ' ethoxide as leaving group.'
        ),
        'r1': {
            'smi': 'CCOC([O-])(O)C',
            'name': 'Tetrahedral intermediate',
            'rot': 0,
            'w': 4.2,
        },
        'r2': None,
        'cond': {'top': '$\\text{elimination}$', 'bottom': '$-\\text{EtO}^-$'},
        'p1': {'smi': 'CC(=O)O', 'name': 'Acetic acid', 'rot': 0, 'w': 3.6},
        'p2': {'smi': 'CC[O-]', 'name': 'Ethoxide', 'rot': 0, 'w': 2.8},
        'arrows': [
            '\\draw[chemarrow] (0.57, 0.94) to[out=120, in=190] (0.58, 0.43);',
            (
                '\\draw[chemarrow] (0.24, -0.38) to[out=120, in=205] (-0.10,'
                ' -0.62);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (-0.872, 0.656) .. controls (-0.937, 0.818)'
                ' and (-0.443, 0.759) .. (-0.378, 0.597);'
            ),
            (
                '\\draw[chemarrow] (0.278, 0.351) .. controls (0.422, 0.346)'
                ' and (0.622, 0.124) .. (0.522, 0.020);'
            ),
        ],
        'rot': {'r1': 270, 'p1': 270, 'p2': 0},
    },
    {
        'id': 'rxn_023',
        'title': 'Acid Chloride Aminolysis: Dimethylamine Attack',
        'category': 'Acyl Substitution',
        'desc': (
            'Dimethylamine nucleophilically attacks acetyl chloride to form a'
            ' tetrahedral intermediate.'
        ),
        'r1': {
            'smi': 'CC(=O)Cl',
            'name': 'Acetyl chloride',
            'rot': 0,
            'w': 3.8,
        },
        'r2': {'smi': 'CNC', 'name': 'Dimethylamine', 'rot': 180, 'w': 2.8},
        'cond': {
            'top': '$\\text{DCM}, 0^\\circ\\text{C}$',
            'bottom': '$\\text{step 1: addition}$',
        },
        'p1': {
            'smi': 'CC([O-])(Cl)[NH+](C)C',
            'name': 'Tetrahedral intermediate',
            'rot': 0,
            'w': 4.4,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (4.50, 0.14) to[out=140, in=30] (0.14, -0.29);',
            (
                '\\draw[chemarrow] (0.70, -0.70) to[out=120, in=205] (1.52,'
                ' -1.05);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.820, 0.493) .. controls (2.382, 2.451)'
                ' and (1.662, 2.452) .. (0.216, 0.501);'
            ),
            (
                '\\draw[chemarrow] (-0.666, 0.605) .. controls (-0.897, 0.437)'
                ' and (-1.392, 0.746) .. (-1.362, 1.030);'
            ),
        ],
        'rot': {'r1': 180, 'p1': 180, 'r2': 300},
    },
    {
        'id': 'rxn_024',
        'title': 'Acid Chloride Tetrahedral Collapse (Chloride Expulsion)',
        'category': 'Acyl Substitution',
        'desc': (
            'O- lone pair pushes down to reform C=O, expelling chloride to'
            ' yield the neutral amide.'
        ),
        'r1': {
            'smi': 'CC([O-])(Cl)N(C)C',
            'name': 'Tetrahedral intermediate',
            'rot': 0,
            'w': 4.4,
        },
        'r2': None,
        'cond': {'top': '$\\text{elimination}$', 'bottom': '$-\\text{Cl}^-$'},
        'p1': {
            'smi': 'CC(=O)N(C)C',
            'name': 'N,N-dimethylacetamide',
            'rot': 0,
            'w': 4.0,
        },
        'p2': {'smi': '[Cl-]', 'name': 'Chloride', 'rot': 0, 'w': 1.8},
        'arrows': [
            (
                '\\draw[chemarrow] (-0.69, -1.48) to[out=120, in=190] (-0.51,'
                ' -0.74);'
            ),
            (
                '\\draw[chemarrow] (-0.51, 0.74) to[out=120, in=205] (-0.67,'
                ' 1.56);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (-1.386, 0.659) .. controls (-1.386, 0.911)'
                ' and (-0.707, 0.671) .. (-0.707, 0.419);'
            ),
            (
                '\\draw[chemarrow] (0.583, 0.419) .. controls (0.677, 0.631)'
                ' and (1.130, 0.788) .. (1.224, 0.577);'
            ),
        ],
        'rot': {'r1': 90, 'p1': 270, 'p2': 0},
    },
    {
        'id': 'rxn_025',
        'title': 'Fischer Esterification: Ethanol Attack on Activated Acid',
        'category': 'Acyl Substitution',
        'desc': (
            'Ethanol attacks protonated acetic acid carbonyl in acid-catalyzed'
            ' ester formation.'
        ),
        'r1': {
            'smi': 'CC(=[OH+])O',
            'name': 'Protonated acetic acid',
            'rot': 60,
            'w': 3.8,
        },
        'r2': {'smi': 'CCO', 'name': 'Ethanol', 'rot': 180, 'w': 2.8},
        'cond': {
            'top': '$\\text{H}_2\\text{SO}_4, \\text{reflux}$',
            'bottom': '$\\text{addition}$',
        },
        'p1': {
            'smi': 'CC(O)(O)[OH+]CC',
            'name': 'Protonated adduct',
            'rot': 60,
            'w': 4.5,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (3.64, 0.56) to[out=140, in=-30] (0.21, 0.14);',
            (
                '\\draw[chemarrow] (-0.65, 0.47) to[out=120, in=205] (-1.39,'
                ' 0.72);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.421, 0.346) .. controls (2.081, 1.817)'
                ' and (1.499, 1.753) .. (0.513, 0.024);'
            ),
            (
                '\\draw[chemarrow] (-0.290, -0.269) .. controls (-0.392,'
                ' -0.497) and (-0.900, -0.543) .. (-1.001, -0.315);'
            ),
        ],
        'rot': {'r1': 270, 'p1': 210, 'r2': 180},
    },
    {
        'id': 'rxn_026',
        'title': 'Transesterification: Methoxide Attack on Ethyl Propionate',
        'category': 'Acyl Substitution',
        'desc': (
            'Base-catalyzed transesterification: methoxide attacks ethyl'
            ' propionate carbonyl.'
        ),
        'r1': {
            'smi': 'CCC(=O)OCC',
            'name': 'Ethyl propionate',
            'rot': 0,
            'w': 4.4,
        },
        'r2': {'smi': 'C[O-]', 'name': 'Sodium methoxide', 'rot': 0, 'w': 2.2},
        'cond': {
            'top': '$\\text{MeOH}$',
            'bottom': '$\\text{transesterification}$',
        },
        'p1': {
            'smi': 'CCC([O-])(OC)OCC',
            'name': 'Mixed tetrahedral intermediate',
            'rot': 0,
            'w': 4.8,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (-0.52, 0.11) to[out=120, in=190] (-0.86, 0.41);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.419, 0.189) .. controls (2.296, 1.545)'
                ' and (1.685, 1.593) .. (0.365, 0.427);'
            ),
            (
                '\\draw[chemarrow] (-0.150, 0.250) .. controls (-0.210, 0.143)'
                ' and (-0.460, 0.122) .. (-0.499, 0.238);'
            ),
        ],
        'rot': {'r1': 90, 'p1': 300, 'r2': 180},
    },
    {
        'id': 'rxn_027',
        'title': 'Aspirin Synthesis: Salicylic Acid Attack on Acetic Anhydride',
        'category': 'Acyl Substitution',
        'desc': (
            'Phenolic hydroxyl of salicylic acid attacks acetic anhydride'
            ' carbonyl.'
        ),
        'r1': {
            'smi': 'CC(=O)OC(=O)C',
            'name': 'Acetic anhydride',
            'rot': 0,
            'w': 4.5,
        },
        'r2': {
            'smi': 'Oc1ccccc1C(=O)O',
            'name': 'Salicylic acid',
            'rot': 0,
            'w': 4.2,
        },
        'cond': {
            'top': '$\\text{H}_3\\text{PO}_4, 85^\\circ\\text{C}$',
            'bottom': '',
        },
        'p1': {
            'smi': 'CC(=O)OC(C)([O-])[OH+]c1ccccc1C(=O)O',
            'name': 'Tetrahedral intermediate',
            'rot': 0,
            'w': 4.8,
        },
        'p2': None,
        'arrows': [
            (
                '\\draw[chemarrow] (5.66, -0.96) to[out=140, in=-30] (-0.80,'
                ' 0.17);'
            ),
            (
                '\\draw[chemarrow] (-1.15, -0.29) to[out=120, in=205] (-1.25,'
                ' -0.95);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (4.196, 0.879) .. controls (2.720, 3.208)'
                ' and (1.961, 3.222) .. (0.398, 0.952);'
            ),
            (
                '\\draw[chemarrow] (-0.200, 0.713) .. controls (-0.267, 0.563)'
                ' and (-0.603, 0.535) .. (-0.670, 0.685);'
            ),
        ],
        'rot': {'r1': 90, 'p1': 30, 'r2': 120},
    },
    {
        'id': 'rxn_028',
        'title': 'Amide Hydrolysis: Hydroxide Attack on Acetamide',
        'category': 'Acyl Substitution',
        'desc': (
            'Strong nucleophile hydroxide attacks acetamide carbonyl under'
            ' basic reflux.'
        ),
        'r1': {'smi': 'CC(=O)N', 'name': 'Acetamide', 'rot': 0, 'w': 3.5},
        'r2': {'smi': '[OH-]', 'name': 'Hydroxide', 'rot': 0, 'w': 2.0},
        'cond': {
            'top': '$\\text{NaOH}, \\text{reflux}$',
            'bottom': '$\\text{addition}$',
        },
        'p1': {
            'smi': 'CC([O-])(O)N',
            'name': 'Tetrahedral intermediate',
            'rot': 0,
            'w': 3.8,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (3.99, 0.32) to[out=140, in=30] (-0.07, -0.15);',
            (
                '\\draw[chemarrow] (-0.42, 0.41) to[out=120, in=205] (-0.52,'
                ' 1.17);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.742, 0.186) .. controls (2.555, 2.231)'
                ' and (1.891, 2.278) .. (0.425, 0.423);'
            ),
            (
                '\\draw[chemarrow] (0.090, -0.273) .. controls (-0.127, -0.370)'
                ' and (-0.127, -0.856) .. (0.091, -0.953);'
            ),
        ],
        'rot': {'r1': 180, 'p1': 210, 'r2': 0},
    },
    {
        'id': 'rxn_029',
        'title': 'Weinreb Ketone Synthesis: Grignard Addition',
        'category': 'Acyl Substitution',
        'desc': (
            'Organometallic reagent adds to Weinreb amide forming a stable'
            ' chelated tetrahedral intermediate.'
        ),
        'r1': {
            'smi': 'CC(=O)N(C)OC',
            'name': 'Weinreb amide',
            'rot': 60,
            'w': 4.5,
        },
        'r2': {
            'smi': 'c1ccccc1[Mg]Br',
            'name': 'Phenylmagnesium bromide',
            'rot': 180,
            'w': 4.2,
        },
        'cond': {
            'top': '$\\text{THF}, 0^\\circ\\text{C}$',
            'bottom': '$\\text{chelation}$',
        },
        'p1': {
            'smi': 'CC([O-])(c1ccccc1)N(C)OC',
            'name': 'Chelated magnesium adduct',
            'rot': 60,
            'w': 5.2,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (7.00, 1.08) to[out=140, in=17] (0.23, -0.13);',
            (
                '\\draw[chemarrow] (0.26, -0.73) to[out=120, in=205] (0.72,'
                ' -1.18);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (5.144, 0.354) .. controls (2.785, 1.121)'
                ' and (1.744, 0.888) .. (-0.062, -0.814);'
            ),
            (
                '\\draw[chemarrow] (-0.642, -1.253) .. controls (-0.665,'
                ' -1.451) and (-1.050, -1.587) .. (-1.182, -1.436);'
            ),
        ],
        'rot': {'r1': 120, 'p1': 30, 'r2': 240},
    },
    {
        'id': 'rxn_030',
        'title': 'Acyl Azide Formation from Acid Chloride',
        'category': 'Acyl Substitution',
        'desc': (
            'Sodium azide displaces chloride via addition-elimination to form'
            ' acyl azide (Curtius precursor).'
        ),
        'r1': {
            'smi': 'c1ccccc1C(=O)Cl',
            'name': 'Benzoyl chloride',
            'rot': 60,
            'w': 4.2,
        },
        'r2': {
            'smi': '[N-]=[N+]=[N-]',
            'name': 'Azide ion',
            'rot': 180,
            'w': 2.5,
        },
        'cond': {
            'top': '$\\text{acetone}, \\text{H}_2\\text{O}$',
            'bottom': '$-\\text{Cl}^-$',
        },
        'p1': {
            'smi': '[N-]=[N+]=NC([O-])(Cl)c1ccccc1',
            'name': 'Benzoyl azide',
            'rot': 60,
            'w': 4.6,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (3.54, 0.32) to[out=140, in=0] (-1.16, -0.84);',
            (
                '\\draw[chemarrow] (-1.30, -0.39) to[out=120, in=205] (-1.07,'
                ' 0.06);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.314, -0.206) .. controls (2.202, -2.361)'
                ' and (1.492, -2.474) .. (-0.232, -0.767);'
            ),
            (
                '\\draw[chemarrow] (-0.792, -0.558) .. controls (-0.851,'
                ' -0.426) and (-1.145, -0.391) .. (-1.203, -0.523);'
            ),
        ],
        'rot': {'r1': 120, 'p1': 120, 'r2': 0},
    },
    {
        'id': 'rxn_031',
        'title': 'Classic SN2: Cyanide Displacement of Bromoethane',
        'category': 'Nucleophilic Substitution',
        'desc': (
            'Concerted bimolecular nucleophilic substitution with backside'
            ' attack and inversion of configuration.'
        ),
        'r1': {'smi': 'CCBr', 'name': 'Bromoethane', 'rot': 240, 'w': 3.2},
        'r2': {'smi': '[C-]#N', 'name': 'Cyanide ion', 'rot': 180, 'w': 2.0},
        'cond': {
            'top': '$\\text{DMSO}, 25^\\circ\\text{C}$',
            'bottom': '$S_N2$',
        },
        'p1': {'smi': 'CCC#N', 'name': 'Propanenitrile', 'rot': 0, 'w': 3.4},
        'p2': {'smi': '[Br-]', 'name': 'Bromide', 'rot': 0, 'w': 1.8},
        'arrows': [
            '\\draw[chemarrow] (3.06, 0.32) to[out=140, in=-30] (0.99, -0.43);',
            (
                '\\draw[chemarrow] (0.20, -0.81) to[out=120, in=205] (-0.45,'
                ' -1.36);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (2.778, 0.124) .. controls (1.813, 1.239)'
                ' and (1.440, 1.173) .. (0.913, -0.204);'
            ),
            (
                '\\draw[chemarrow] (0.236, -0.719) .. controls (0.263, -0.976)'
                ' and (-0.101, -1.368) .. (-0.310, -1.216);'
            ),
        ],
        'rot': {'r1': 120, 'p1': 120, 'r2': 180, 'p2': 180},
    },
    {
        'id': 'rxn_032',
        'title': 'Williamson Ether Synthesis: Ethoxide on Iodomethane',
        'category': 'Nucleophilic Substitution',
        'desc': (
            'Ethoxide nucleophile attacks methyl iodide via SN2 displacement to'
            ' yield ethyl methyl ether.'
        ),
        'r1': {'smi': 'CI', 'name': 'Iodomethane', 'rot': 180, 'w': 2.6},
        'r2': {
            'smi': 'CC[O-]',
            'name': 'Sodium ethoxide',
            'rot': 180,
            'w': 2.8,
        },
        'cond': {
            'top': '$\\text{EtOH}, 25^\\circ\\text{C}$',
            'bottom': '$S_N2$',
        },
        'p1': {'smi': 'CCOC', 'name': 'Ethyl methyl ether', 'rot': 0, 'w': 3.4},
        'p2': {'smi': '[I-]', 'name': 'Iodide', 'rot': 0, 'w': 1.8},
        'arrows': [
            '\\draw[chemarrow] (2.78, 0.60) to[out=140, in=0] (1.40, -0.00);',
            (
                '\\draw[chemarrow] (0.01, -0.14) to[out=120, in=205] (-1.21,'
                ' -0.24);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (2.493, 0.395) .. controls (1.812, 1.121)'
                ' and (1.570, 1.065) .. (1.283, 0.112);'
            ),
            (
                '\\draw[chemarrow] (0.122, -0.080) .. controls (-0.047, -0.459)'
                ' and (-0.891, -0.562) .. (-1.060, -0.183);'
            ),
        ],
        'rot': {'r1': 180, 'p1': 180, 'r2': 180, 'p2': 180},
    },
    {
        'id': 'rxn_033',
        'title': 'Finkelstein Reaction: Iodide on 1-Chlorobutane',
        'category': 'Nucleophilic Substitution',
        'desc': (
            'Equilibrium driven by precipitation of NaCl in acetone driving SN2'
            ' substitution to 1-iodobutane.'
        ),
        'r1': {'smi': 'CCCCCl', 'name': '1-Chlorobutane', 'rot': 240, 'w': 3.8},
        'r2': {'smi': '[I-]', 'name': 'Sodium iodide', 'rot': 0, 'w': 2.0},
        'cond': {
            'top': '$\\text{acetone}, \\text{reflux}$',
            'bottom': '$-\\text{NaCl}\\downarrow$',
        },
        'p1': {'smi': 'CCCCI', 'name': '1-Iodobutane', 'rot': 0, 'w': 3.8},
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (4.06, 0.32) to[out=140, in=-30] (0.33, -0.98);',
            (
                '\\draw[chemarrow] (-0.21, -1.21) to[out=120, in=205] (-0.61,'
                ' -1.62);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.779, 0.129) .. controls (2.383, 1.556)'
                ' and (1.798, 1.473) .. (0.857, -0.288);'
            ),
            (
                '\\draw[chemarrow] (0.643, -0.810) .. controls (0.504, -0.872)'
                ' and (0.347, -1.154) .. (0.486, -1.216);'
            ),
        ],
        'rot': {'r1': 60, 'p1': 60, 'r2': 0},
    },
    {
        'id': 'rxn_034',
        'title': 'Gabriel Synthesis Step: Phthalimide Alkylation',
        'category': 'Nucleophilic Substitution',
        'desc': (
            'Potassium phthalimide anion attacks benzyl bromide in clean SN2'
            ' displacement.'
        ),
        'r1': {
            'smi': 'BrCc1ccccc1',
            'name': 'Benzyl bromide',
            'rot': 240,
            'w': 4.2,
        },
        'r2': {
            'smi': 'O=C1c2ccccc2C(=O)[N-]1',
            'name': 'Potassium phthalimide',
            'rot': 0,
            'w': 4.2,
        },
        'cond': {
            'top': '$\\text{DMF}, 80^\\circ\\text{C}$',
            'bottom': '$S_N2$',
        },
        'p1': {
            'smi': 'O=C1c2ccccc2C(=O)N1Cc1ccccc1',
            'name': 'N-Benzylphthalimide',
            'rot': 0,
            'w': 5.2,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (6.18, -0.37) to[out=140, in=73] (1.05, 0.66);',
            '\\draw[chemarrow] (1.42, 0.17) to[out=120, in=205] (1.80, -0.02);',
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.765, 0.107) .. controls (2.370, 1.433)'
                ' and (1.856, 1.313) .. (1.197, -0.495);'
            ),
            (
                '\\draw[chemarrow] (0.874, -0.983) .. controls (0.724, -1.000)'
                ' and (0.491, -1.219) .. (0.604, -1.318);'
            ),
        ],
        'rot': {'r1': 240, 'p1': 270, 'r2': 0},
    },
    {
        'id': 'rxn_035',
        'title': 'Azidation: Sodium Azide on Benzyl Bromide',
        'category': 'Nucleophilic Substitution',
        'desc': (
            'Azide ion attacks the benzylic position displacing bromide to'
            ' yield benzyl azide.'
        ),
        'r1': {
            'smi': 'BrCc1ccccc1',
            'name': 'Benzyl bromide',
            'rot': 180,
            'w': 4.2,
        },
        'r2': {
            'smi': '[N-]=[N+]=[N-]',
            'name': 'Sodium azide',
            'rot': 180,
            'w': 2.5,
        },
        'cond': {
            'top': '$\\text{DMSO}, 25^\\circ\\text{C}$',
            'bottom': '$S_N2$',
        },
        'p1': {
            'smi': '[N-]=[N+]=NCc1ccccc1',
            'name': 'Benzyl azide',
            'rot': 0,
            'w': 4.6,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (5.45, 0.32) to[out=140, in=73] (1.05, 0.66);',
            '\\draw[chemarrow] (1.42, 0.17) to[out=120, in=205] (1.80, -0.02);',
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.250, 0.093) .. controls (1.988, 1.408)'
                ' and (1.576, 1.292) .. (1.192, -0.491);'
            ),
            (
                '\\draw[chemarrow] (0.874, -0.983) .. controls (0.724, -1.000)'
                ' and (0.491, -1.219) .. (0.604, -1.318);'
            ),
        ],
        'rot': {'r1': 240, 'p1': 240, 'r2': 0},
    },
    {
        'id': 'rxn_036',
        'title': 'Base-Catalyzed Oxirane Opening: Methoxide on Ethylene Oxide',
        'category': 'Nucleophilic Substitution',
        'desc': (
            'Methoxide nucleophile attacks the less hindered carbon of ethylene'
            ' oxide with ring strain relief.'
        ),
        'r1': {'smi': 'C1CO1', 'name': 'Ethylene oxide', 'rot': 90, 'w': 3.0},
        'r2': {'smi': 'C[O-]', 'name': 'Methoxide', 'rot': 180, 'w': 2.2},
        'cond': {
            'top': '$\\text{MeOH}, 25^\\circ\\text{C}$',
            'bottom': '$S_N2$ opening',
        },
        'p1': {
            'smi': 'COCC[O-]',
            'name': '2-Methoxyethoxide',
            'rot': 0,
            'w': 3.6,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (2.96, 0.32) to[out=140, in=90] (-0.06, 1.29);',
            '\\draw[chemarrow] (0.64, 0.13) to[out=120, in=205] (1.26, -0.74);',
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (2.706, 0.171) .. controls (2.042, 1.024)'
                ' and (1.716, 1.019) .. (1.076, 0.147);'
            ),
            (
                '\\draw[chemarrow] (0.183, -0.460) .. controls (0.221, -0.819)'
                ' and (-0.363, -1.275) .. (-0.656, -1.063);'
            ),
        ],
        'rot': {'r1': 0, 'p1': 270, 'r2': 180},
    },
    {
        'id': 'rxn_037',
        'title': (
            'Acid-Catalyzed Oxirane Opening: Methanol on Protonated Propylene'
            ' Oxide'
        ),
        'category': 'Nucleophilic Substitution',
        'desc': (
            'Methanol attacks the more substituted carbon due to significant'
            ' partial positive charge build-up.'
        ),
        'r1': {
            'smi': 'CC1C[OH+]1',
            'name': 'Protonated propylene oxide',
            'rot': 90,
            'w': 3.4,
        },
        'r2': {'smi': 'CO', 'name': 'Methanol', 'rot': 180, 'w': 2.2},
        'cond': {
            'top': '$\\text{MeOH}, \\text{H}^+$',
            'bottom': '$\\text{regioselective}$',
        },
        'p1': {
            'smi': 'CC([OH+]C)CO',
            'name': 'Protonated ether-alcohol',
            'rot': 0,
            'w': 4.0,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (3.43, 0.32) to[out=140, in=15] (-0.05, 0.05);',
            (
                '\\draw[chemarrow] (0.17, -0.51) to[out=120, in=205] (0.59,'
                ' -0.92);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.210, 0.124) .. controls (2.087, 1.682)'
                ' and (1.565, 1.657) .. (0.597, -0.001);'
            ),
            (
                '\\draw[chemarrow] (-0.132, 0.153) .. controls (-0.327, 0.012)'
                ' and (-0.815, 0.111) .. (-0.789, 0.350);'
            ),
        ],
        'rot': {'r1': 60, 'p1': 240, 'r2': 180},
    },
    {
        'id': 'rxn_038',
        'title': 'SN1 Ionization: Heterolysis of tert-Butyl Chloride',
        'category': 'Nucleophilic Substitution',
        'desc': (
            'Rate-determining unimolecular heterolysis of C-Cl bond forming'
            ' planar tert-butyl carbocation.'
        ),
        'r1': {
            'smi': 'CC(C)(C)Cl',
            'name': 'tert-Butyl chloride',
            'rot': 0,
            'w': 3.5,
        },
        'r2': None,
        'cond': {
            'top': '$\\text{H}_2\\text{O} / \\text{acetone}$',
            'bottom': '$S_N1\\text{ (RDS)}$',
        },
        'p1': {
            'smi': 'C[C+](C)C',
            'name': 'tert-Butyl cation',
            'rot': 0,
            'w': 3.2,
        },
        'p2': {'smi': '[Cl-]', 'name': 'Chloride ion', 'rot': 0, 'w': 1.8},
        'arrows': [
            '\\draw[chemarrow] (-0.48, 0.47) to[out=120, in=205] (-1.02, 1.09);'
        ],
        'arrows_tikz': [
            '\\draw[chemarrow] (-0.080, 0.451) .. controls (-0.279, 0.539) and'
            ' (-0.436, 0.963) .. (-0.238, 1.052);'
        ],
        'rot': {'r1': 0, 'p1': 330, 'p2': 0},
    },
    {
        'id': 'rxn_039',
        'title': 'SN1 Solvent Capture: Water Attack on tert-Butyl Cation',
        'category': 'Nucleophilic Substitution',
        'desc': (
            'Fast nucleophilic capture of planar sp2 carbocation by water'
            ' solvent.'
        ),
        'r1': {
            'smi': 'C[C+](C)C',
            'name': 'tert-Butyl cation',
            'rot': 0,
            'w': 3.2,
        },
        'r2': {'smi': 'O', 'name': 'Water', 'rot': 180, 'w': 2.0},
        'cond': {'top': '$\\text{fast capture}$', 'bottom': '$S_N1$ step 2'},
        'p1': {
            'smi': 'CC(C)(C)[OH2+]',
            'name': 'Protonated tert-butanol',
            'rot': 0,
            'w': 3.5,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (3.93, 0.36) to[out=140, in=90] (0.48, 1.12);'
        ],
        'arrows_tikz': [
            '\\draw[chemarrow] (3.712, 0.168) .. controls (2.472, 2.121) and'
            ' (1.853, 2.121) .. (0.618, 0.164);'
        ],
        'rot': {'r1': 270, 'p1': 300, 'r2': 0},
    },
    {
        'id': 'rxn_040',
        'title': 'Intramolecular SN2: Bromohydrin Cyclization to Epoxide',
        'category': 'Nucleophilic Substitution',
        'desc': (
            'Base-generated alkoxide performs intramolecular backside'
            ' displacement of bromide.'
        ),
        'r1': {
            'smi': 'CC(Br)C[O-]',
            'name': 'Bromohydrin alkoxide',
            'rot': 0,
            'w': 3.8,
        },
        'r2': None,
        'cond': {
            'top': '$\\text{NaOH}, \\text{H}_2\\text{O}$',
            'bottom': '$\\text{intramolecular } S_N2$',
        },
        'p1': {'smi': 'CC1CO1', 'name': 'Propylene oxide', 'rot': 0, 'w': 3.0},
        'p2': {'smi': '[Br-]', 'name': 'Bromide', 'rot': 0, 'w': 1.8},
        'arrows': [
            (
                '\\draw[chemarrow] (0.37, -0.65) to[out=140, in=150] (-0.87,'
                ' -0.26);'
            ),
            (
                '\\draw[chemarrow] (-0.80, 0.22) to[out=120, in=205] (-0.96,'
                ' 0.91);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (1.142, -0.214) .. controls (0.440, 0.618)'
                ' and (0.109, 0.604) .. (-0.517, -0.287);'
            ),
            (
                '\\draw[chemarrow] (-0.742, 0.162) .. controls (-0.925, 0.243)'
                ' and (-1.083, 0.631) .. (-0.899, 0.713);'
            ),
        ],
        'rot': {'r1': 0, 'p1': 180, 'p2': 0},
    },
    {
        'id': 'rxn_041',
        'title': 'Classic E2: Ethoxide-Promoted Dehydrohalogenation',
        'category': 'Elimination',
        'desc': (
            'Concerted anti-periplanar abstraction of beta-proton, pi-bond'
            ' formation, and bromide departure.'
        ),
        'r1': {
            'smi': 'CC(Br)C',
            'name': '2-Bromopropane',
            'rot': 120,
            'w': 3.4,
        },
        'r2': {
            'smi': 'CC[O-]',
            'name': 'Sodium ethoxide',
            'rot': 180,
            'w': 2.8,
        },
        'cond': {
            'top': '$\\text{EtOH}, 55^\\circ\\text{C}$',
            'bottom': '$E2$ concerted',
        },
        'p1': {'smi': 'CC=C', 'name': 'Propene', 'rot': 0, 'w': 3.0},
        'p2': {'smi': '[Br-]', 'name': 'Bromide', 'rot': 0, 'w': 1.8},
        'arrows': [
            '\\draw[chemarrow] (3.18, 0.60) to[out=140, in=-30] (1.66, -1.19);',
            (
                '\\draw[chemarrow] (-0.45, -0.81) to[out=120, in=205] (-1.10,'
                ' -1.37);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (2.878, 0.344) .. controls (2.046, 0.949)'
                ' and (1.784, 0.849) .. (1.570, -0.157);'
            ),
            (
                '\\draw[chemarrow] (0.884, -0.564) .. controls (0.815, -0.226)'
                ' and (0.083, -0.254) .. (-0.024, -0.582);'
            ),
            (
                '\\draw[chemarrow] (-0.623, 0.157) .. controls (-0.795, 0.234)'
                ' and (-0.953, 0.596) .. (-0.781, 0.673);'
            ),
        ],
        'rot': {'r1': 0, 'p1': 0, 'r2': 180, 'p2': 180},
    },
    {
        'id': 'rxn_042',
        'title': 'Hofmann Elimination: Steric Anti-Zaitsev E2',
        'category': 'Elimination',
        'desc': (
            'Bulky potassium tert-butoxide abstracts the more accessible'
            ' primary proton to give 2-methylbut-1-ene.'
        ),
        'r1': {
            'smi': 'CCC(C)(C)Br',
            'name': '2-Bromo-2-methylbutane',
            'rot': 120,
            'w': 4.0,
        },
        'r2': {
            'smi': 'CC(C)(C)[O-]',
            'name': 'Potassium tert-butoxide',
            'rot': 180,
            'w': 3.4,
        },
        'cond': {
            'top': '$t\\text{-BuOH}, 75^\\circ\\text{C}$',
            'bottom': '$\\text{Hofmann } E2$',
        },
        'p1': {
            'smi': 'CCC(=C)C',
            'name': '2-Methylbut-1-ene',
            'rot': 0,
            'w': 3.8,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (5.60, -0.82) to[out=140, in=90] (0.00, 1.76);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.672, 1.259) .. controls (2.820, 2.448)'
                ' and (2.454, 2.413) .. (1.844, 1.084);'
            ),
            (
                '\\draw[chemarrow] (1.173, 0.581) .. controls (1.012, 0.901)'
                ' and (0.286, 0.678) .. (0.267, 0.320);'
            ),
            (
                '\\draw[chemarrow] (-0.801, 0.627) .. controls (-0.832, 0.821)'
                ' and (-1.156, 1.071) .. (-1.280, 0.918);'
            ),
        ],
        'rot': {'r1': 120, 'p1': 300, 'r2': 330},
    },
    {
        'id': 'rxn_043',
        'title': 'Zaitsev Elimination: E2 Formation of trans-2-Butene',
        'category': 'Elimination',
        'desc': (
            'Small base ethoxide abstracts secondary proton yielding the more'
            ' substituted trans-alkene predominantly.'
        ),
        'r1': {
            'smi': 'CCC(C)Br',
            'name': '2-Bromobutane',
            'rot': 120,
            'w': 3.8,
        },
        'r2': {
            'smi': 'CC[O-]',
            'name': 'Sodium ethoxide',
            'rot': 180,
            'w': 2.8,
        },
        'cond': {
            'top': '$\\text{EtOH}, 60^\\circ\\text{C}$',
            'bottom': '$\\text{Zaitsev } E2$',
        },
        'p1': {'smi': 'C/C=C/C', 'name': 'trans-2-Butene', 'rot': 0, 'w': 3.5},
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (3.38, 0.60) to[out=140, in=30] (0.07, -0.64);',
            '\\draw[chemarrow] (0.39, 1.07) to[out=120, in=205] (0.98, 1.59);',
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.136, 0.470) .. controls (2.650, 1.379)'
                ' and (2.369, 1.404) .. (1.731, 0.594);'
            ),
            (
                '\\draw[chemarrow] (0.941, 0.409) .. controls (1.050, 0.081)'
                ' and (0.403, -0.260) .. (0.147, -0.030);'
            ),
            (
                '\\draw[chemarrow] (-0.695, -0.569) .. controls (-0.772,'
                ' -0.741) and (-1.133, -0.898) .. (-1.210, -0.726);'
            ),
        ],
        'rot': {'r1': 300, 'p1': 300, 'r2': 180},
    },
    {
        'id': 'rxn_044',
        'title': 'E1cB Elimination: Enolate Beta-Elimination of Hydroxide',
        'category': 'Elimination',
        'desc': (
            'Enolate intermediate expels poor leaving group hydroxide in the'
            ' condensation step of aldol reaction.'
        ),
        'r1': {
            'smi': 'CC(=O)[CH-]C(O)c1ccccc1',
            'name': 'Aldol enolate intermediate',
            'rot': 0,
            'w': 4.8,
        },
        'r2': None,
        'cond': {'top': '$\\text{NaOH}, \\Delta$', 'bottom': '$E1cB$'},
        'p1': {
            'smi': 'CC(=O)/C=C/c1ccccc1',
            'name': 'alpha,beta-Unsaturated enone',
            'rot': 0,
            'w': 4.8,
        },
        'p2': {'smi': '[OH-]', 'name': 'Hydroxide', 'rot': 0, 'w': 1.8},
        'arrows': [
            '\\draw[chemarrow] (-0.64, 0.11) to[out=120, in=205] (-0.35,'
            ' -0.25);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (0.425, 0.772) .. controls (0.393, 0.926)'
                ' and (0.071, 0.614) .. (0.103, 0.460);'
            ),
            (
                '\\draw[chemarrow] (-0.419, 0.613) .. controls (-0.368, 0.725)'
                ' and (-0.448, 0.969) .. (-0.565, 0.932);'
            ),
        ],
        'rot': {'r1': 120, 'p1': 120, 'p2': 0},
    },
    {
        'id': 'rxn_045',
        'title': 'E1 Step 1: Ionization of 2-Bromo-2-methylpropane',
        'category': 'Elimination',
        'desc': (
            'Heterolysis of C-Br bond to form stable tertiary carbocation in'
            ' polar protic solvent.'
        ),
        'r1': {
            'smi': 'CC(C)(C)Br',
            'name': 'tert-Butyl bromide',
            'rot': 0,
            'w': 3.5,
        },
        'r2': None,
        'cond': {
            'top': '$\\text{EtOH} / \\text{H}_2\\text{O}$',
            'bottom': '$E1\\text{ (RDS)}$',
        },
        'p1': {
            'smi': 'C[C+](C)C',
            'name': 'tert-Butyl cation',
            'rot': 0,
            'w': 3.2,
        },
        'p2': {'smi': '[Br-]', 'name': 'Bromide', 'rot': 0, 'w': 1.8},
        'arrows': [
            '\\draw[chemarrow] (-0.48, 0.48) to[out=120, in=205] (-1.02, 1.11);'
        ],
        'arrows_tikz': [
            '\\draw[chemarrow] (-0.080, 0.459) .. controls (-0.280, 0.548) and'
            ' (-0.438, 0.975) .. (-0.238, 1.064);'
        ],
        'rot': {'r1': 0, 'p1': 330, 'p2': 0},
    },
    {
        'id': 'rxn_046',
        'title': 'E1 Step 2: Deprotonation of Carbocation to Isobutene',
        'category': 'Elimination',
        'desc': (
            'Weak base solvent abstracts adjacent beta-proton from carbocation'
            ' establishing alkene double bond.'
        ),
        'r1': {
            'smi': 'C[C+](C)C',
            'name': 'tert-Butyl cation',
            'rot': 0,
            'w': 3.2,
        },
        'r2': {'smi': 'O', 'name': 'Water', 'rot': 180, 'w': 2.0},
        'cond': {
            'top': '$\\text{fast}$',
            'bottom': '$-\\text{H}_3\\text{O}^+$',
        },
        'p1': {'smi': 'CC(=C)C', 'name': 'Isobutene', 'rot': 0, 'w': 3.2},
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (3.93, 0.36) to[out=140, in=150] (-1.65, 0.86);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.726, 0.188) .. controls (2.966, 1.524)'
                ' and (2.532, 1.562) .. (1.553, 0.377);'
            ),
            (
                '\\draw[chemarrow] (0.830, 0.462) .. controls (0.765, 0.146)'
                ' and (0.080, 0.172) .. (-0.021, 0.480);'
            ),
        ],
        'rot': {'r1': 0, 'p1': 180, 'r2': 0},
    },
    {
        'id': 'rxn_047',
        'title': 'Acid-Catalyzed Alcohol Dehydration: Water Elimination',
        'category': 'Elimination',
        'desc': (
            'Protonated cyclohexanol loses water with concerted/stepwise proton'
            ' abstraction to form cyclohexene.'
        ),
        'r1': {
            'smi': 'C1CCC([OH2+])CC1',
            'name': 'Protonated cyclohexanol',
            'rot': 90,
            'w': 4.0,
        },
        'r2': {'smi': 'O', 'name': 'Water', 'rot': 180, 'w': 2.0},
        'cond': {
            'top': '$\\text{H}_2\\text{SO}_4, 160^\\circ\\text{C}$',
            'bottom': '$-\\text{H}_3\\text{O}^+$',
        },
        'p1': {'smi': 'C1=CCCCC1', 'name': 'Cyclohexene', 'rot': 90, 'w': 3.5},
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (4.33, 0.36) to[out=140, in=-30] (1.02, 0.01);',
            (
                '\\draw[chemarrow] (0.14, -0.81) to[out=120, in=205] (0.30,'
                ' -1.36);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (4.145, 0.211) .. controls (3.525, 1.983)'
                ' and (3.053, 2.081) .. (1.781, 0.700);'
            ),
            (
                '\\draw[chemarrow] (1.011, 0.536) .. controls (1.114, 0.227)'
                ' and (0.503, -0.095) .. (0.260, 0.122);'
            ),
            (
                '\\draw[chemarrow] (0.112, -0.768) .. controls (0.008, -0.912)'
                ' and (0.040, -1.284) .. (0.217, -1.265);'
            ),
        ],
        'rot': {'r1': 210, 'p1': 330, 'r2': 0},
    },
    {
        'id': 'rxn_048',
        'title': 'Chugaev Elimination: Pyrolytic Xanthate Elimination',
        'category': 'Elimination',
        'desc': (
            'Thermal syn-elimination of methyl xanthate via 6-membered cyclic'
            ' transition state.'
        ),
        'r1': {
            'smi': 'CCC(C)OC(=S)SC',
            'name': 'sec-Butyl xanthate',
            'rot': 0,
            'w': 4.8,
        },
        'r2': None,
        'cond': {
            'top': '$200^\\circ\\text{C}$',
            'bottom': '$\\text{syn-elimination}$',
        },
        'p1': {'smi': 'C/C=C/C', 'name': 'trans-2-Butene', 'rot': 0, 'w': 3.5},
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (-0.54, 0.47) to[out=40, in=139] (1.80, 0.31);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (-0.875, -1.099) .. controls (-3.064, 1.263)'
                ' and (-2.808, 1.804) .. (0.407, 1.604);'
            ),
            (
                '\\draw[chemarrow] (0.373, 1.422) .. controls (0.573, 1.305)'
                ' and (0.392, 0.848) .. (0.161, 0.851);'
            ),
            (
                '\\draw[chemarrow] (-0.140, 0.283) .. controls (-0.256, 0.336)'
                ' and (-0.509, 0.252) .. (-0.471, 0.130);'
            ),
        ],
        'rot': {'r1': 270, 'p1': 270},
    },
    {
        'id': 'rxn_049',
        'title': 'Cope Elimination: Thermal Syn-Elimination of Amine Oxide',
        'category': 'Elimination',
        'desc': (
            'Pyrolysis of tertiary amine oxide gives alkene and hydroxylamine'
            ' via 5-membered cyclic transition state.'
        ),
        'r1': {
            'smi': 'CCC(C)[N+](C)(C)[O-]',
            'name': 'Amine oxide',
            'rot': 0,
            'w': 4.5,
        },
        'r2': None,
        'cond': {
            'top': '$120^\\circ\\text{C}$',
            'bottom': '$\\text{Cope syn-elimination}$',
        },
        'p1': {'smi': 'C/C=C/C', 'name': 'trans-2-Butene', 'rot': 0, 'w': 3.5},
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (-1.50, 1.87) to[out=40, in=125] (1.13, -0.25);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (-0.738, -1.304) .. controls (-2.316, 0.316)'
                ' and (-2.163, 0.864) .. (0.026, 1.434);'
            ),
            (
                '\\draw[chemarrow] (0.054, 1.164) .. controls (0.326, 1.114)'
                ' and (0.312, 0.528) .. (0.051, 0.437);'
            ),
            (
                '\\draw[chemarrow] (0.108, -0.389) .. controls (0.199, -0.509)'
                ' and (0.203, -0.829) .. (0.053, -0.816);'
            ),
        ],
        'rot': {'r1': 270, 'p1': 270},
    },
    {
        'id': 'rxn_050',
        'title': 'Julia Olefination: Reductive Elimination Step',
        'category': 'Elimination',
        'desc': (
            'Radical-mediated beta-elimination of acetoxy sulfone yields'
            ' (E)-alkene stereoselectively.'
        ),
        'r1': {
            'smi': 'CCCC[CH-]C(OC(=O)C)CCC',
            'name': 'Acetoxy sulfone',
            'rot': 0,
            'w': 5.5,
        },
        'r2': None,
        'cond': {
            'top': '$\\text{SmI}_2 \\text{ or } \\text{Na(Hg)}$',
            'bottom': '$(E)\\text{-selective}$',
        },
        'p1': {'smi': 'CCC/C=C/CCCC', 'name': '(E)-Alkene', 'rot': 0, 'w': 4.5},
        'p2': {'smi': 'CC(=O)[O-]', 'name': 'Acetate', 'rot': 0, 'w': 2.0},
        'arrows': [
            '\\draw[chemarrow] (0.74, -0.44) to[out=120, in=205] (0.96, -0.81);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (-0.302, 0.789) .. controls (-0.873, 1.448)'
                ' and (-0.928, 1.414) .. (-0.579, 0.614);'
            ),
            (
                '\\draw[chemarrow] (-0.667, 0.147) .. controls (-0.571, 0.080)'
                ' and (-0.511, -0.158) .. (-0.625, -0.185);'
            ),
        ],
        'rot': {'r1': 180, 'p1': 180, 'p2': 0},
    },
    {
        'id': 'rxn_051',
        'title': 'Bromination Step 1: Bromonium Ion Formation',
        'category': 'Electrophilic Addition',
        'desc': (
            'Alkene pi-electrons attack electrophilic bromine, with'
            ' back-donation of bromine lone pair to form cyclic bromonium ion.'
        ),
        'r1': {'smi': 'C=C', 'name': 'Ethene', 'rot': 90, 'w': 2.8},
        'r2': {'smi': 'BrBr', 'name': 'Bromine', 'rot': 90, 'w': 2.8},
        'cond': {
            'top': '$\\text{DCM}, 0^\\circ\\text{C}$',
            'bottom': '$\\text{bromonium ion}$',
        },
        'p1': {
            'smi': 'C1C[Br+]1',
            'name': 'Bromonium intermediate',
            'rot': 0,
            'w': 3.2,
        },
        'p2': {'smi': '[Br-]', 'name': 'Bromide', 'rot': 0, 'w': 1.8},
        'arrows': [
            '\\draw[chemarrow] (0.00, -0.76) to[out=40, in=-90] (3.94, -1.22);',
            '\\draw[chemarrow] (3.80, -0.00) to[out=120, in=205] (3.64, 0.93);',
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (0.076, 0.047) .. controls (1.153, 1.701)'
                ' and (1.724, 1.724) .. (2.930, 0.161);'
            ),
            (
                '\\draw[chemarrow] (3.933, -0.080) .. controls (4.048, -0.340)'
                ' and (4.615, -0.498) .. (4.731, -0.238);'
            ),
            (
                '\\draw[chemarrow] (2.930, 0.161) .. controls (1.270, 2.432)'
                ' and (0.461, 2.417) .. (-1.115, 0.088);'
            ),
        ],
        'rot': {'r1': 180, 'p1': 180, 'r2': 0, 'p2': 0},
    },
    {
        'id': 'rxn_052',
        'title': 'Bromination Step 2: Bromide Opening of Bromonium Ion',
        'category': 'Electrophilic Addition',
        'desc': (
            'Bromide attacks the bridged bromonium ion from the back face,'
            ' yielding trans-1,2-dibromoethane (anti-addition).'
        ),
        'r1': {'smi': 'C1C[Br+]1', 'name': 'Bromonium ion', 'rot': 0, 'w': 3.2},
        'r2': {'smi': '[Br-]', 'name': 'Bromide', 'rot': 0, 'w': 1.8},
        'cond': {
            'top': '$\\text{anti-addition}$',
            'bottom': '$S_N2\\text{ opening}$',
        },
        'p1': {
            'smi': 'BrCCBr',
            'name': '1,2-Dibromoethane',
            'rot': 0,
            'w': 3.5,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (3.73, 0.32) to[out=140, in=0] (1.31, 0.06);',
            (
                '\\draw[chemarrow] (0.20, -0.60) to[out=120, in=205] (-0.58,'
                ' -1.23);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.468, 0.168) .. controls (2.556, 1.425)'
                ' and (2.106, 1.422) .. (1.214, 0.150);'
            ),
            (
                '\\draw[chemarrow] (0.265, -0.491) .. controls (0.305, -0.872)'
                ' and (-0.285, -1.395) .. (-0.595, -1.169);'
            ),
        ],
        'rot': {'r1': 0, 'p1': 90, 'r2': 0},
    },
    {
        'id': 'rxn_053',
        'title': 'Markovnikov Hydrobromination: Alkene Protonation',
        'category': 'Electrophilic Addition',
        'desc': (
            'Propene pi-electrons attack HBr forming the more stable secondary'
            ' carbocation (Markovnikov rule).'
        ),
        'r1': {'smi': 'CC=C', 'name': 'Propene', 'rot': 90, 'w': 3.0},
        'r2': {'smi': 'Br', 'name': 'Hydrogen bromide', 'rot': 0, 'w': 2.2},
        'cond': {
            'top': '$0^\\circ\\text{C}$',
            'bottom': '$\\text{Markovnikov (RDS)}$',
        },
        'p1': {
            'smi': 'C[CH+]C',
            'name': 'Isopropyl cation',
            'rot': 0,
            'w': 3.2,
        },
        'p2': {'smi': '[Br-]', 'name': 'Bromide', 'rot': 0, 'w': 1.8},
        'arrows': [
            '\\draw[chemarrow] (0.25, 1.48) to[out=40, in=158] (3.60, 0.14);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (0.736, -0.327) .. controls (1.307, 1.081)'
                ' and (1.720, 1.165) .. (2.798, 0.094);'
            ),
            (
                '\\draw[chemarrow] (3.684, -0.080) .. controls (3.785, -0.308)'
                ' and (4.278, -0.466) .. (4.380, -0.238);'
            ),
        ],
        'rot': {'r1': 60, 'p1': 60, 'r2': 180, 'p2': 180},
    },
    {
        'id': 'rxn_054',
        'title': 'Hydrobromination Step 2: Bromide Capture of Carbocation',
        'category': 'Electrophilic Addition',
        'desc': (
            'Bromide ion captures the planar secondary carbocation yielding'
            ' 2-bromopropane.'
        ),
        'r1': {
            'smi': 'C[CH+]C',
            'name': 'Isopropyl cation',
            'rot': 0,
            'w': 3.2,
        },
        'r2': {'smi': '[Br-]', 'name': 'Bromide', 'rot': 0, 'w': 1.8},
        'cond': {
            'top': '$\\text{fast capture}$',
            'bottom': '$\\text{product}$',
        },
        'p1': {'smi': 'CC(Br)C', 'name': '2-Bromopropane', 'rot': 0, 'w': 3.4},
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (3.73, 0.32) to[out=140, in=90] (0.00, 0.73);'
        ],
        'arrows_tikz': [
            '\\draw[chemarrow] (3.487, 0.194) .. controls (2.651, 1.928) and'
            ' (2.127, 1.996) .. (0.872, 0.537);'
        ],
        'rot': {'r1': 60, 'p1': 60, 'r2': 0},
    },
    {
        'id': 'rxn_055',
        'title': 'Epoxidation of Cyclohexene by mCPBA (Concerted Butterfly)',
        'category': 'Electrophilic Addition',
        'desc': (
            'Concerted peroxy acid oxygen transfer to alkene double bond with'
            ' retention of stereochemistry.'
        ),
        'r1': {'smi': 'C1=CCCCC1', 'name': 'Cyclohexene', 'rot': 330, 'w': 3.8},
        'r2': {
            'smi': 'O=C(OO)c1cccc(Cl)c1',
            'name': 'mCPBA',
            'rot': 30,
            'w': 4.5,
        },
        'cond': {
            'top': '$\\text{DCM}, 0^\\circ\\text{C}$',
            'bottom': '$\\text{concerted syn-addition}$',
        },
        'p1': {
            'smi': 'C1CC2OC2CC1',
            'name': 'Cyclohexene oxide',
            'rot': 90,
            'w': 3.8,
        },
        'p2': {
            'smi': 'O=C(O)c1cccc(Cl)c1',
            'name': 'm-Chlorobenzoic acid',
            'rot': 0,
            'w': 4.2,
        },
        'arrows': [
            '\\draw[chemarrow] (-0.00, 1.79) to[out=40, in=-177] (4.92, -0.33);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (1.355, -0.758) .. controls (2.580, -1.774)'
                ' and (3.050, -1.666) .. (3.707, -0.217);'
            ),
            (
                '\\draw[chemarrow] (4.268, -0.009) .. controls (4.227, 0.212)'
                ' and (4.677, 0.370) .. (4.822, 0.198);'
            ),
            (
                '\\draw[chemarrow] (4.833, 0.765) .. controls (4.584, 0.552)'
                ' and (4.371, -0.021) .. (4.160, 0.229);'
            ),
            (
                '\\draw[chemarrow] (3.591, -0.410) .. controls (3.184, -2.327)'
                ' and (2.668, -2.544) .. (1.013, -1.494);'
            ),
        ],
        'rot': {'r1': 0, 'p1': 270, 'r2': 180, 'p2': 30},
    },
    {
        'id': 'rxn_056',
        'title': 'Hydroboration: Concerted Syn-Addition of Borane',
        'category': 'Electrophilic Addition',
        'desc': (
            'Concerted four-center addition of B-H across alkene with boron'
            ' adding to less substituted carbon.'
        ),
        'r1': {'smi': 'CCCCC=C', 'name': '1-Hexene', 'rot': 90, 'w': 4.2},
        'r2': {'smi': 'B', 'name': 'Borane (BH3)', 'rot': 0, 'w': 2.2},
        'cond': {
            'top': '$\\text{THF}, 0^\\circ\\text{C}$',
            'bottom': '$\\text{syn, anti-Markovnikov}$',
        },
        'p1': {'smi': 'CCCCCCB', 'name': 'Hexylborane', 'rot': 0, 'w': 4.5},
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (-0.10, 1.94) to[out=40, in=158] (3.96, 0.16);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (1.694, 0.801) .. controls (2.858, 1.794)'
                ' and (3.273, 1.680) .. (3.768, 0.233);'
            ),
            (
                '\\draw[chemarrow] (4.461, 0.098) .. controls (3.438, 1.827)'
                ' and (2.841, 1.890) .. (1.476, 0.417);'
            ),
        ],
        'rot': {'r1': 330, 'p1': 150, 'r2': 0},
    },
    {
        'id': 'rxn_057',
        'title': 'Oxymercuration: Mercurinium Ion Formation',
        'category': 'Electrophilic Addition',
        'desc': (
            'Electrophilic attack of mercuric acetate on propene forming a'
            ' bridged mercurinium ion.'
        ),
        'r1': {'smi': 'CC=C', 'name': 'Propene', 'rot': 90, 'w': 3.0},
        'r2': {
            'smi': 'CC(=O)O[Hg]OC(=O)C',
            'name': 'Mercuric acetate',
            'rot': 180,
            'w': 4.5,
        },
        'cond': {
            'top': '$\\text{THF} / \\text{H}_2\\text{O}$',
            'bottom': '$\\text{mercurinium ion}$',
        },
        'p1': {
            'smi': 'CC1C[Hg+]1OC(=O)C',
            'name': 'Mercurinium intermediate',
            'rot': 0,
            'w': 3.5,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (-0.42, 0.32) to[out=40, in=90] (4.95, -0.28);'
        ],
        'arrows_tikz': [
            '\\draw[chemarrow] (0.734, -0.325) .. controls (1.618, 1.890) and'
            ' (2.283, 2.046) .. (4.059, 0.454);'
        ],
        'rot': {'r1': 60, 'p1': 60, 'r2': 300},
    },
    {
        'id': 'rxn_058',
        'title': 'Oxymercuration: Water Capture of Mercurinium Ion',
        'category': 'Electrophilic Addition',
        'desc': (
            'Water attacks the more substituted carbon of mercurinium ion with'
            ' anti-stereochemistry.'
        ),
        'r1': {
            'smi': 'CC1C[Hg+]1OC(=O)C',
            'name': 'Mercurinium ion',
            'rot': 0,
            'w': 3.5,
        },
        'r2': {'smi': 'O', 'name': 'Water', 'rot': 180, 'w': 2.0},
        'cond': {
            'top': '$\\text{Markovnikov capture}$',
            'bottom': '$\\text{anti-addition}$',
        },
        'p1': {
            'smi': 'CC([OH2+])C[Hg]OC(=O)C',
            'name': 'Organomercury alcohol',
            'rot': 0,
            'w': 4.0,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (4.08, 0.36) to[out=140, in=75] (0.15, 0.35);',
            (
                '\\draw[chemarrow] (-0.45, -0.35) to[out=120, in=205] (-0.91,'
                ' -0.80);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.859, 0.162) .. controls (2.453, 2.049)'
                ' and (1.799, 2.016) .. (0.590, -0.003);'
            ),
            (
                '\\draw[chemarrow] (-0.155, 0.159) .. controls (-0.355, 0.014)'
                ' and (-0.857, 0.122) .. (-0.831, 0.368);'
            ),
        ],
        'rot': {'r1': 60, 'p1': 240, 'r2': 0},
    },
    {
        'id': 'rxn_059',
        'title': 'Alkyne Hydration: Protonation to Vinyl Cation',
        'category': 'Electrophilic Addition',
        'desc': (
            'Protonation of propyne by acid yields the secondary vinyl cation'
            ' intermediate.'
        ),
        'r1': {'smi': 'CC#C', 'name': 'Propyne', 'rot': 0, 'w': 3.0},
        'r2': {'smi': '[OH3+]', 'name': 'Hydronium', 'rot': 0, 'w': 2.2},
        'cond': {
            'top': '$\\text{H}_2\\text{SO}_4, \\text{HgSO}_4$',
            'bottom': '$\\text{vinyl cation}$',
        },
        'p1': {'smi': 'C[C+]=C', 'name': 'Vinyl cation', 'rot': 0, 'w': 3.2},
        'p2': {'smi': 'O', 'name': 'Water', 'rot': 0, 'w': 2.0},
        'arrows': [
            '\\draw[chemarrow] (1.35, 0.32) to[out=40, in=158] (3.53, 0.19);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (0.066, -0.517) .. controls (0.806, 1.309)'
                ' and (1.359, 1.438) .. (2.833, 0.130);'
            ),
            (
                '\\draw[chemarrow] (3.596, -0.035) .. controls (3.676, -0.216)'
                ' and (4.056, -0.373) .. (4.136, -0.193);'
            ),
        ],
        'rot': {'r1': 90, 'p1': 90, 'r2': 180, 'p2': 180},
    },
    {
        'id': 'rxn_060',
        'title': 'Halohydrin Formation: Water Attack on Bromonium Ion',
        'category': 'Electrophilic Addition',
        'desc': (
            'Water captures the more substituted position of bromonium ion with'
            ' inversion, yielding a bromohydrin.'
        ),
        'r1': {
            'smi': 'CC1C[Br+]1',
            'name': 'Propylene bromonium ion',
            'rot': 0,
            'w': 3.5,
        },
        'r2': {'smi': 'O', 'name': 'Water', 'rot': 180, 'w': 2.0},
        'cond': {
            'top': '$\\text{H}_2\\text{O}$',
            'bottom': '$\\text{halohydrin}$',
        },
        'p1': {
            'smi': 'CC([OH2+])CBr',
            'name': 'Protonated bromohydrin',
            'rot': 0,
            'w': 4.0,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (4.08, 0.36) to[out=140, in=75] (0.11, 0.30);',
            (
                '\\draw[chemarrow] (-0.52, -0.42) to[out=120, in=205] (-1.00,'
                ' -0.87);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.859, 0.162) .. controls (2.437, 2.049)'
                ' and (1.775, 2.016) .. (0.549, -0.003);'
            ),
            (
                '\\draw[chemarrow] (-0.196, 0.159) .. controls (-0.396, 0.014)'
                ' and (-0.898, 0.122) .. (-0.872, 0.368);'
            ),
        ],
        'rot': {'r1': 60, 'p1': 240, 'r2': 0},
    },
    {
        'id': 'rxn_061',
        'title': 'Benzene Nitration Step 1: Attack on Nitronium Ion',
        'category': 'Electrophilic Aromatic Substitution',
        'desc': (
            'Aromatic pi-system attacks strongly electrophilic nitronium ion'
            ' (NO2+) to form resonance-stabilized Wheland intermediate.'
        ),
        'r1': {'smi': 'c1ccccc1', 'name': 'Benzene', 'rot': 30, 'w': 3.2},
        'r2': {'smi': 'O=[N+]=O', 'name': 'Nitronium ion', 'rot': 0, 'w': 2.5},
        'cond': {
            'top': '$\\text{HNO}_3, \\text{H}_2\\text{SO}_4$',
            'bottom': '$\\text{Wheland RDS}$',
        },
        'p1': {
            'smi': 'O=[N+]([O-])C1C=CC=C[CH+]1',
            'name': 'Arenium (Wheland) ion',
            'rot': 30,
            'w': 3.8,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (1.07, -0.30) to[out=40, in=0] (5.20, 0.85);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (0.620, 0.951) .. controls (2.323, 2.566)'
                ' and (2.955, 2.421) .. (3.782, 0.225);'
            ),
            (
                '\\draw[chemarrow] (4.138, -0.377) .. controls (4.261, -0.432)'
                ' and (4.301, -0.705) .. (4.178, -0.759);'
            ),
        ],
        'rot': {'r1': 330, 'p1': 210, 'r2': 90},
    },
    {
        'id': 'rxn_062',
        'title': 'Benzene Nitration Step 2: Wheland Rearomatization',
        'category': 'Electrophilic Aromatic Substitution',
        'desc': (
            'Bisulfate base abstracts the sp3 proton restoring aromatic'
            ' stabilization and yielding nitrobenzene.'
        ),
        'r1': {
            'smi': 'O=[N+]([O-])C1C=C[CH+]C=C1',
            'name': 'Wheland intermediate',
            'rot': 30,
            'w': 3.8,
        },
        'r2': {
            'smi': 'OS(=O)(=O)[O-]',
            'name': 'Bisulfate ion',
            'rot': 180,
            'w': 2.8,
        },
        'cond': {
            'top': '$\\text{fast}$',
            'bottom': '$\\text{rearomatization}$',
        },
        'p1': {
            'smi': 'O=[N+]([O-])c1ccccc1',
            'name': 'Nitrobenzene',
            'rot': 30,
            'w': 3.6,
        },
        'p2': {
            'smi': 'OS(=O)(=O)O',
            'name': 'Sulfuric acid',
            'rot': 0,
            'w': 2.8,
        },
        'arrows': [
            '\\draw[chemarrow] (5.32, 0.86) to[out=140, in=-30] (0.14, -1.07);',
            (
                '\\draw[chemarrow] (-0.39, -1.30) to[out=120, in=205] (-0.80,'
                ' -1.56);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.179, 0.643) .. controls (2.156, 1.717)'
                ' and (1.778, 1.635) .. (1.286, 0.236);'
            ),
            (
                '\\draw[chemarrow] (0.679, -0.053) .. controls (0.748, -0.209)'
                ' and (0.584, -0.405) .. (0.413, -0.416);'
            ),
            (
                '\\draw[chemarrow] (0.017, -1.010) .. controls (-0.091, -0.797)'
                ' and (-0.531, -1.060) .. (-0.542, -1.298);'
            ),
        ],
        'rot': {'r1': 300, 'p1': 300, 'r2': 240, 'p2': 30},
    },
    {
        'id': 'rxn_063',
        'title': 'Benzene Bromination: Attack on Br2-FeBr3 Complex',
        'category': 'Electrophilic Aromatic Substitution',
        'desc': (
            'Lewis acid activates bromine for electrophilic attack by benzene'
            ' ring.'
        ),
        'r1': {'smi': 'c1ccccc1', 'name': 'Benzene', 'rot': 30, 'w': 3.2},
        'r2': {
            'smi': 'BrBr',
            'name': 'Bromine (FeBr3 complex)',
            'rot': 90,
            'w': 3.0,
        },
        'cond': {
            'top': '$\\text{FeBr}_3, 25^\\circ\\text{C}$',
            'bottom': '$S_E\\text{Ar}$',
        },
        'p1': {
            'smi': 'BrC1C=CC=C[CH+]1',
            'name': 'Bromonium arenium ion',
            'rot': 30,
            'w': 3.8,
        },
        'p2': {'smi': '[Br-]', 'name': 'Bromide', 'rot': 0, 'w': 1.8},
        'arrows': [
            '\\draw[chemarrow] (1.07, -0.30) to[out=40, in=-90] (4.23, -1.28);',
            '\\draw[chemarrow] (4.09, -0.00) to[out=120, in=205] (3.94, 0.99);',
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (1.155, 0.643) .. controls (2.234, 1.714)'
                ' and (2.647, 1.630) .. (3.221, 0.222);'
            ),
            (
                '\\draw[chemarrow] (4.228, -0.080) .. controls (4.353, -0.361)'
                ' and (4.966, -0.518) .. (5.091, -0.238);'
            ),
        ],
        'rot': {'r1': 0, 'p1': 240, 'r2': 0, 'p2': 0},
    },
    {
        'id': 'rxn_064',
        'title': 'Friedel-Crafts Acylation: Attack on Acylium Ion',
        'category': 'Electrophilic Aromatic Substitution',
        'desc': (
            'Benzene attacks resonance-stabilized acylium ion generated from'
            ' acetyl chloride and AlCl3.'
        ),
        'r1': {'smi': 'c1ccccc1', 'name': 'Benzene', 'rot': 30, 'w': 3.2},
        'r2': {'smi': 'CC#[O+]', 'name': 'Acylium ion', 'rot': 180, 'w': 2.5},
        'cond': {
            'top': '$\\text{AlCl}_3, 80^\\circ\\text{C}$',
            'bottom': '$\\text{acylation}$',
        },
        'p1': {
            'smi': 'CC(=O)C1C=CC=C[CH+]1',
            'name': 'Acyl arenium ion',
            'rot': 30,
            'w': 4.0,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (1.07, -0.30) to[out=40, in=90] (4.10, 0.24);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (1.153, 0.649) .. controls (2.526, 2.041)'
                ' and (3.073, 1.945) .. (3.887, 0.167);'
            ),
            (
                '\\draw[chemarrow] (4.226, -0.351) .. controls (4.351, -0.407)'
                ' and (4.374, -0.684) .. (4.249, -0.740);'
            ),
        ],
        'rot': {'r1': 0, 'p1': 60, 'r2': 90},
    },
    {
        'id': 'rxn_065',
        'title': 'Friedel-Crafts Alkylation: Attack on Isopropyl Cation',
        'category': 'Electrophilic Aromatic Substitution',
        'desc': (
            'Benzene attacks secondary carbocation generated from isopropyl'
            ' chloride and Lewis acid.'
        ),
        'r1': {'smi': 'c1ccccc1', 'name': 'Benzene', 'rot': 30, 'w': 3.2},
        'r2': {
            'smi': 'C[CH+]C',
            'name': 'Isopropyl cation',
            'rot': 0,
            'w': 3.0,
        },
        'cond': {'top': '$\\text{AlCl}_3$', 'bottom': '$\\text{alkylation}$'},
        'p1': {
            'smi': 'CC(C)C1C=CC=C[CH+]1',
            'name': 'Cumene arenium ion',
            'rot': 30,
            'w': 4.0,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (1.07, -0.30) to[out=40, in=90] (4.30, 0.71);'
        ],
        'arrows_tikz': [
            '\\draw[chemarrow] (1.149, 0.657) .. controls (2.110, 1.886) and'
            ' (2.566, 1.869) .. (3.430, 0.570);'
        ],
        'rot': {'r1': 0, 'p1': 60, 'r2': 300},
    },
    {
        'id': 'rxn_066',
        'title': 'Aromatic Sulfonation: Attack on Sulfur Trioxide',
        'category': 'Electrophilic Aromatic Substitution',
        'desc': (
            'Electrophilic attack of benzene on neutral sulfur trioxide in'
            ' fuming sulfuric acid.'
        ),
        'r1': {'smi': 'c1ccccc1', 'name': 'Benzene', 'rot': 30, 'w': 3.2},
        'r2': {
            'smi': 'O=S(=O)=O',
            'name': 'Sulfur trioxide',
            'rot': 0,
            'w': 2.5,
        },
        'cond': {
            'top': '$\\text{fuming } \\text{H}_2\\text{SO}_4$',
            'bottom': '$\\text{sulfonation}$',
        },
        'p1': {
            'smi': 'O=S(=O)([O-])C1C=CC=C[CH+]1',
            'name': 'Sulfonate arenium ion',
            'rot': 30,
            'w': 4.2,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (1.07, -0.30) to[out=40, in=150] (3.77, -0.11);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (1.141, 0.671) .. controls (1.721, 1.624)'
                ' and (2.072, 1.674) .. (2.894, 0.919);'
            ),
            (
                '\\draw[chemarrow] (3.457, 0.428) .. controls (3.438, 0.250)'
                ' and (3.743, 0.048) .. (3.888, 0.153);'
            ),
        ],
        'rot': {'r1': 0, 'p1': 240, 'r2': 60},
    },
    {
        'id': 'rxn_067',
        'title': 'Aniline Bromination: Activated Ortho/Para Attack',
        'category': 'Electrophilic Aromatic Substitution',
        'desc': (
            'Strongly activating amino group delocalizes lone pair into ring,'
            ' driving rapid ortho-attack on bromine without catalyst.'
        ),
        'r1': {'smi': 'Nc1ccccc1', 'name': 'Aniline', 'rot': 30, 'w': 3.5},
        'r2': {'smi': 'BrBr', 'name': 'Bromine', 'rot': 90, 'w': 2.8},
        'cond': {
            'top': '$\\text{H}_2\\text{O}, 25^\\circ\\text{C}$',
            'bottom': '$\\text{activated rapid}$',
        },
        'p1': {
            'smi': 'NC1=C[CH+]C(Br)C=C1',
            'name': 'Iminium arenium ion',
            'rot': 30,
            'w': 4.0,
        },
        'p2': {'smi': '[Br-]', 'name': 'Bromide', 'rot': 0, 'w': 1.8},
        'arrows': [
            '\\draw[chemarrow] (-1.57, -0.25) to[out=40, in=-90] (4.29, -1.22);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (0.460, -1.153) .. controls (2.101, -2.323)'
                ' and (2.679, -2.141) .. (3.349, -0.241);'
            ),
            (
                '\\draw[chemarrow] (4.283, 0.080) .. controls (4.398, 0.340)'
                ' and (4.965, 0.498) .. (5.081, 0.238);'
            ),
        ],
        'rot': {'r1': 270, 'p1': 150, 'r2': 0, 'p2': 0},
    },
    {
        'id': 'rxn_068',
        'title': 'Azo Coupling: Phenolate Attack on Diazonium Ion',
        'category': 'Electrophilic Aromatic Substitution',
        'desc': (
            'Activated phenolate anion attacks weakly electrophilic'
            ' benzenediazonium ion to form an azo dye.'
        ),
        'r1': {
            'smi': '[O-]c1ccccc1',
            'name': 'Phenolate anion',
            'rot': 30,
            'w': 3.5,
        },
        'r2': {
            'smi': 'c1ccccc1[N+]#N',
            'name': 'Benzenediazonium',
            'rot': 180,
            'w': 4.0,
        },
        'cond': {
            'top': '$\\text{pH } 9-10, 0^\\circ\\text{C}$',
            'bottom': '$\\text{azo coupling}$',
        },
        'p1': {
            'smi': 'O=C1C=CC(N=Nc2ccccc2)C=C1',
            'name': 'Azo intermediate',
            'rot': 30,
            'w': 5.2,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (-1.56, -0.30) to[out=40, in=180] (2.92, 0.00);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (-1.275, 1.360) .. controls (-1.166, 1.550)'
                ' and (-0.786, 1.053) .. (-0.895, 0.864);'
            ),
            (
                '\\draw[chemarrow] (0.192, 0.575) .. controls (0.257, 0.255)'
                ' and (0.933, 0.432) .. (1.035, 0.742);'
            ),
            (
                '\\draw[chemarrow] (1.598, 0.012) .. controls (1.726, 1.008)'
                ' and (2.009, 1.238) .. (3.010, 1.166);'
            ),
            (
                '\\draw[chemarrow] (3.767, 1.162) .. controls (3.894, 1.255)'
                ' and (4.175, 1.099) .. (4.159, 0.943);'
            ),
        ],
        'rot': {'r1': 210, 'p1': 270, 'r2': 210},
    },
    {
        'id': 'rxn_069',
        'title': 'Vilsmeier-Haack Formylation: Pyrrole Attack',
        'category': 'Electrophilic Aromatic Substitution',
        'desc': (
            'Electron-rich pyrrole attacks chloroiminium intermediate'
            ' (generated from DMF and POCl3).'
        ),
        'r1': {'smi': '[nH]1cccc1', 'name': 'Pyrrole', 'rot': 90, 'w': 3.2},
        'r2': {
            'smi': 'ClC=[N+](C)C',
            'name': 'Vilsmeier reagent',
            'rot': 180,
            'w': 3.8,
        },
        'cond': {
            'top': '$\\text{DMF}, \\text{POCl}_3$',
            'bottom': '$\\text{formylation}$',
        },
        'p1': {
            'smi': 'CN(C)C(Cl)C1[CH+]C=CN1',
            'name': 'Iminium adduct',
            'rot': 90,
            'w': 4.5,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (-1.00, 0.34) to[out=40, in=90] (4.33, -0.66);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (0.703, 0.752) .. controls (2.131, 2.165)'
                ' and (2.725, 2.076) .. (3.674, 0.305);'
            ),
            (
                '\\draw[chemarrow] (4.193, -0.124) .. controls (4.173, -0.307)'
                ' and (4.400, -0.620) .. (4.549, -0.512);'
            ),
        ],
        'rot': {'r1': 60, 'p1': 240, 'r2': 240},
    },
    {
        'id': 'rxn_070',
        'title': 'Kolbe-Schmitt Reaction: Phenolate Carboxylation',
        'category': 'Electrophilic Aromatic Substitution',
        'desc': (
            'Sodium phenolate attacks carbon dioxide under high pressure and'
            ' temperature to yield salicylate.'
        ),
        'r1': {
            'smi': '[O-]c1ccccc1',
            'name': 'Sodium phenolate',
            'rot': 30,
            'w': 3.5,
        },
        'r2': {'smi': 'O=C=O', 'name': 'Carbon dioxide', 'rot': 0, 'w': 2.2},
        'cond': {
            'top': '$125^\\circ\\text{C}, 100\\text{ atm}$',
            'bottom': '$\\text{salicylate}$',
        },
        'p1': {
            'smi': 'O=C([O-])C1C=CC=CC1=O',
            'name': 'Dienone intermediate',
            'rot': 30,
            'w': 4.2,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (0.34, -0.30) to[out=40, in=90] (4.05, 0.24);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (0.575, -1.365) .. controls (0.409, -1.461)'
                ' and (0.371, -0.915) .. (0.537, -0.819);'
            ),
            (
                '\\draw[chemarrow] (-0.051, -0.393) .. controls (1.731, -2.913)'
                ' and (2.524, -2.854) .. (3.913, -0.099);'
            ),
            (
                '\\draw[chemarrow] (4.185, 0.335) .. controls (4.291, 0.382)'
                ' and (4.338, 0.614) .. (4.233, 0.661);'
            ),
        ],
        'rot': {'r1': 60, 'p1': 30, 'r2': 90},
    },
    {
        'id': 'rxn_071',
        'title': 'Aldol Addition: Enolate Attack on Neutral Aldehyde',
        'category': 'Enolate Chemistry',
        'desc': (
            'Acetaldehyde enolate attacks neutral acetaldehyde generating'
            ' beta-hydroxy alkoxide intermediate.'
        ),
        'r1': {
            'smi': 'C=C[O-]',
            'name': 'Acetaldehyde enolate',
            'rot': 180,
            'w': 3.0,
        },
        'r2': {'smi': 'CC=O', 'name': 'Acetaldehyde', 'rot': 240, 'w': 3.0},
        'cond': {
            'top': '$\\text{NaOH}, \\text{H}_2\\text{O}, 0^\\circ\\text{C}$',
            'bottom': '$\\text{addition}$',
        },
        'p1': {
            'smi': 'CC([O-])CC=O',
            'name': 'Beta-hydroxy alkoxide',
            'rot': 0,
            'w': 4.2,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (1.22, 0.03) to[out=120, in=190] (0.54, 0.21);',
            '\\draw[chemarrow] (-1.28, 0.07) to[out=40, in=150] (3.27, 0.47);',
            '\\draw[chemarrow] (4.05, 0.84) to[out=120, in=205] (4.69, 1.33);',
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (-1.360, -0.009) .. controls (-1.494,'
                ' -0.240) and (-0.746, -0.395) .. (-0.612, -0.163);'
            ),
            (
                '\\draw[chemarrow] (0.752, -0.036) .. controls (1.502, 1.731)'
                ' and (2.017, 1.822) .. (3.325, 0.417);'
            ),
            (
                '\\draw[chemarrow] (4.204, 0.423) .. controls (4.411, 0.271)'
                ' and (4.800, 0.631) .. (4.773, 0.887);'
            ),
        ],
        'rot': {'r1': 180, 'p1': 180, 'r2': 300},
    },
    {
        'id': 'rxn_072',
        'title': 'Claisen Condensation: Ester Enolate Attack on Ester',
        'category': 'Enolate Chemistry',
        'desc': (
            'Ethyl acetate enolate attacks neutral ethyl acetate forming'
            ' beta-keto ester tetrahedral adduct.'
        ),
        'r1': {
            'smi': 'C=C([O-])OCC',
            'name': 'Ethyl acetate enolate',
            'rot': 0,
            'w': 4.2,
        },
        'r2': {
            'smi': 'CCOC(=O)C',
            'name': 'Ethyl acetate',
            'rot': 60,
            'w': 4.0,
        },
        'cond': {
            'top': '$\\text{NaOEt}, \\text{EtOH}$',
            'bottom': '$\\text{addition}$',
        },
        'p1': {
            'smi': 'CCOC(=O)CC([O-])(OCC)C',
            'name': 'Tetrahedral adduct',
            'rot': 0,
            'w': 5.2,
        },
        'p2': None,
        'arrows': [
            (
                '\\draw[chemarrow] (-1.46, -0.78) to[out=120, in=190] (-1.20,'
                ' -0.27);'
            ),
            '\\draw[chemarrow] (-1.84, 1.24) to[out=40, in=-158] (5.22, 0.70);',
            '\\draw[chemarrow] (5.14, 1.20) to[out=120, in=205] (4.75, 1.63);',
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (0.579, 1.417) .. controls (0.784, 1.500)'
                ' and (0.779, 0.867) .. (0.574, 0.784);'
            ),
            (
                '\\draw[chemarrow] (1.372, 0.178) .. controls (2.393, 1.855)'
                ' and (2.954, 1.891) .. (4.177, 0.354);'
            ),
            (
                '\\draw[chemarrow] (4.416, -0.217) .. controls (4.564, -0.310)'
                ' and (4.536, -0.665) .. (4.369, -0.712);'
            ),
        ],
        'rot': {'r1': 150, 'p1': 270, 'r2': 180},
    },
    {
        'id': 'rxn_073',
        'title': 'Dieckmann Cyclization: Intramolecular Enolate Attack',
        'category': 'Enolate Chemistry',
        'desc': (
            'Intramolecular Claisen condensation of diethyl adipate enolate to'
            ' form ethyl 2-oxocyclopentanecarboxylate.'
        ),
        'r1': {
            'smi': 'CCOC(=O)CCCC=C([O-])OCC',
            'name': 'Diethyl adipate enolate',
            'rot': 0,
            'w': 5.5,
        },
        'r2': None,
        'cond': {
            'top': '$\\text{NaOEt}, \\Delta$',
            'bottom': '$\\text{cyclization}$',
        },
        'p1': {
            'smi': 'CCOC(=O)C1CCCC1([O-])OCC',
            'name': 'Cyclic beta-keto ester',
            'rot': 0,
            'w': 4.5,
        },
        'p2': None,
        'arrows': [
            (
                '\\draw[chemarrow] (-1.18, -0.49) to[out=120, in=190] (-0.93,'
                ' -0.32);'
            ),
            '\\draw[chemarrow] (-0.03, -0.09) to[out=40, in=102] (0.94, 0.05);',
            (
                '\\draw[chemarrow] (1.20, -0.45) to[out=120, in=205] (1.37,'
                ' -0.80);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (2.307, 0.610) .. controls (2.398, 0.510)'
                ' and (2.012, 0.484) .. (1.921, 0.585);'
            ),
            (
                '\\draw[chemarrow] (1.270, 0.569) .. controls (0.255, 1.714)'
                ' and (-0.172, 1.659) .. (-0.865, 0.295);'
            ),
            (
                '\\draw[chemarrow] (-1.175, 0.427) .. controls (-1.264, 0.446)'
                ' and (-1.357, 0.609) .. (-1.284, 0.663);'
            ),
        ],
        'rot': {'r1': 180, 'p1': 150},
    },
    {
        'id': 'rxn_074',
        'title': 'Michael Addition: Malonate Enolate on MVK',
        'category': 'Enolate Chemistry',
        'desc': (
            'Resonance-stabilized soft carbanion adds in conjugate 1,4-fashion'
            ' to methyl vinyl ketone.'
        ),
        'r1': {
            'smi': 'CCOC(=O)[CH-]C(=O)OCC',
            'name': 'Diethyl malonate anion',
            'rot': 0,
            'w': 4.8,
        },
        'r2': {
            'smi': 'CC(=O)C=C',
            'name': 'Methyl vinyl ketone',
            'rot': 180,
            'w': 3.6,
        },
        'cond': {
            'top': '$\\text{NaOEt}$ cat.',
            'bottom': '$1,4\\text{-addition}$',
        },
        'p1': {
            'smi': 'CCOC(=O)C(CC=C([O-])C)C(=O)OCC',
            'name': 'Michael enolate adduct',
            'rot': 0,
            'w': 5.5,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (-0.54, 0.17) to[out=40, in=-150] (3.57, 0.27);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (0.615, 0.190) .. controls (1.676, 1.901)'
                ' and (2.278, 1.955) .. (3.626, 0.459);'
            ),
            (
                '\\draw[chemarrow] (4.454, 0.468) .. controls (4.525, 0.120)'
                ' and (5.260, 0.310) .. (5.370, 0.647);'
            ),
            (
                '\\draw[chemarrow] (6.114, -0.174) .. controls (6.307, -0.259)'
                ' and (6.316, -0.692) .. (6.123, -0.778);'
            ),
        ],
        'rot': {'r1': 270, 'p1': 60, 'r2': 180},
    },
    {
        'id': 'rxn_075',
        'title': 'Robinson Annulation: Intramolecular Aldol Ring Closure',
        'category': 'Enolate Chemistry',
        'desc': (
            'Enolate generated by Michael addition cyclizes onto second'
            ' carbonyl in 6-membered aldol closure.'
        ),
        'r1': {
            'smi': 'C=C([O-])CCC1C(=O)CCCC1',
            'name': '1,5-Diketone enolate',
            'rot': 0,
            'w': 4.8,
        },
        'r2': None,
        'cond': {
            'top': '$\\text{KOH}, \\text{EtOH}$',
            'bottom': '$\\text{cyclization}$',
        },
        'p1': {
            'smi': 'O=C1CCC2CCCCC2([O-])C1',
            'name': 'Aldol alkoxide adduct',
            'rot': 0,
            'w': 4.8,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (4.74, 0.32) to[out=140, in=133] (-2.41, 0.60);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (-0.733, 1.844) .. controls (-0.690, 1.986)'
                ' and (-0.426, 1.656) .. (-0.469, 1.514);'
            ),
            (
                '\\draw[chemarrow] (0.149, 1.503) .. controls (1.502, 0.833)'
                ' and (1.551, 0.424) .. (0.392, -0.544);'
            ),
            (
                '\\draw[chemarrow] (0.638, -0.579) .. controls (0.736, -0.621)'
                ' and (0.927, -0.510) .. (0.891, -0.409);'
            ),
        ],
        'rot': {'r1': 90, 'p1': 270},
    },
    {
        'id': 'rxn_076',
        'title': 'Malonic Ester Synthesis: Malonate Alkylation',
        'category': 'Enolate Chemistry',
        'desc': (
            'Malonate carbanion displaces bromide from benzyl bromide via clean'
            ' SN2 pathway.'
        ),
        'r1': {
            'smi': 'CCOC(=O)[CH-]C(=O)OCC',
            'name': 'Diethyl malonate anion',
            'rot': 0,
            'w': 4.8,
        },
        'r2': {
            'smi': 'BrCc1ccccc1',
            'name': 'Benzyl bromide',
            'rot': 180,
            'w': 4.2,
        },
        'cond': {
            'top': '$\\text{EtOH}, \\Delta$',
            'bottom': '$S_N2\\text{ alkylation}$',
        },
        'p1': {
            'smi': 'CCOC(=O)C(Cc1ccccc1)C(=O)OCC',
            'name': 'Alkylated malonate',
            'rot': 0,
            'w': 5.5,
        },
        'p2': {'smi': '[Br-]', 'name': 'Bromide', 'rot': 0, 'w': 1.8},
        'arrows': [
            '\\draw[chemarrow] (-0.54, 0.17) to[out=40, in=73] (6.75, 0.66);',
            '\\draw[chemarrow] (7.12, 0.17) to[out=120, in=205] (7.50, -0.02);',
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (0.293, 0.178) .. controls (1.701, 2.199)'
                ' and (2.450, 2.233) .. (4.036, 0.348);'
            ),
            (
                '\\draw[chemarrow] (4.146, 0.808) .. controls (4.306, 0.928)'
                ' and (4.375, 1.339) .. (4.179, 1.378);'
            ),
        ],
        'rot': {'r1': 300, 'p1': 0, 'r2': 30, 'p2': 30},
    },
    {
        'id': 'rxn_077',
        'title': 'Wittig Reaction Step 1: Ylide Addition to Aldehyde',
        'category': 'Enolate Chemistry',
        'desc': (
            'Phosphonium ylide carbanion attacks benzaldehyde carbonyl to form'
            ' the betaine intermediate.'
        ),
        'r1': {
            'smi': 'O=Cc1ccccc1',
            'name': 'Benzaldehyde',
            'rot': 60,
            'w': 4.2,
        },
        'r2': {
            'smi': '[CH2-][P+](c1ccccc1)(c1ccccc1)c1ccccc1',
            'name': 'Methylenetriphenylphosphorane',
            'rot': 180,
            'w': 4.8,
        },
        'cond': {
            'top': '$\\text{THF}, 0^\\circ\\text{C}$',
            'bottom': '$\\text{betaine formation}$',
        },
        'p1': {
            'smi': '[O-]C(c1ccccc1)C[P+](c1ccccc1)(c1ccccc1)c1ccccc1',
            'name': 'Betaine intermediate',
            'rot': 0,
            'w': 5.5,
        },
        'p2': None,
        'arrows': [
            (
                '\\draw[chemarrow] (5.70, 1.71) to[out=140, in=-47] (-0.07,'
                ' -1.49);'
            ),
            (
                '\\draw[chemarrow] (-0.70, -1.57) to[out=120, in=205] (-1.11,'
                ' -1.76);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (4.211, 0.112) .. controls (2.614, 1.778)'
                ' and (2.012, 1.654) .. (1.204, -0.509);'
            ),
            (
                '\\draw[chemarrow] (0.745, -0.954) .. controls (0.596, -0.971)'
                ' and (0.540, -1.273) .. (0.653, -1.372);'
            ),
        ],
        'rot': {'r1': 240, 'p1': 270, 'r2': 90},
    },
    {
        'id': 'rxn_078',
        'title': 'Wittig Reaction Step 2: Oxaphosphetane Collapse',
        'category': 'Enolate Chemistry',
        'desc': (
            'Four-membered oxaphosphetane cycloreversion driven by formation of'
            ' the strong P=O bond (535 kJ/mol).'
        ),
        'r1': {
            'smi': 'c1ccccc1C1CP(c2ccccc2)(c2ccccc2)(c2ccccc2)O1',
            'name': 'Oxaphosphetane',
            'rot': 0,
            'w': 5.0,
        },
        'r2': None,
        'cond': {
            'top': '$\\text{cycloreversion}$',
            'bottom': '$\\text{alkene formation}$',
        },
        'p1': {'smi': 'C=Cc1ccccc1', 'name': 'Styrene', 'rot': 0, 'w': 3.8},
        'p2': {
            'smi': 'O=P(c1ccccc1)(c1ccccc1)c1ccccc1',
            'name': 'Triphenylphosphine oxide',
            'rot': 0,
            'w': 4.5,
        },
        'arrows': [
            '\\draw[chemarrow] (0.42, 0.33) to[out=120, in=205] (0.53, 0.67);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (0.158, -0.158) .. controls (0.102, -0.046)'
                ' and (-0.104, -0.067) .. (-0.172, -0.172);'
            ),
            (
                '\\draw[chemarrow] (-0.158, 0.158) .. controls (-0.102, 0.046)'
                ' and (0.104, 0.067) .. (0.172, 0.172);'
            ),
        ],
        'rot': {'p1': 180, 'p2': 0, 'r1': 0},
    },
    {
        'id': 'rxn_079',
        'title': 'Horner-Wadsworth-Emmons (HWE) Reaction',
        'category': 'Enolate Chemistry',
        'desc': (
            'Phosphonate carbanion attacks ketone to yield'
            ' (E)-alpha,beta-unsaturated ester with water-soluble phosphate'
            ' byproduct.'
        ),
        'r1': {
            'smi': 'O=C1CCCCC1',
            'name': 'Cyclohexanone',
            'rot': 90,
            'w': 4.2,
        },
        'r2': {
            'smi': 'CCOC(=O)[CH-]P(=O)(OCC)OCC',
            'name': 'Triethyl phosphonoacetate anion',
            'rot': 180,
            'w': 4.8,
        },
        'cond': {
            'top': '$\\text{NaH}, \\text{THF}, 25^\\circ\\text{C}$',
            'bottom': '$(E)\\text{-selective}$',
        },
        'p1': {
            'smi': 'CCOC(=O)C(P(=O)(OCC)OCC)C1([O-])CCCCC1',
            'name': 'Phosphonate alkoxide intermediate',
            'rot': 0,
            'w': 4.8,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (5.17, 0.52) to[out=140, in=30] (0.21, 0.56);',
            '\\draw[chemarrow] (0.14, 0.96) to[out=120, in=205] (0.24, 1.40);',
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (4.226, 0.408) .. controls (2.660, 2.754)'
                ' and (1.845, 2.777) .. (0.148, 0.524);'
            ),
            (
                '\\draw[chemarrow] (-0.157, 0.905) .. controls (-0.312, 0.975)'
                ' and (-0.338, 1.322) .. (-0.183, 1.391);'
            ),
        ],
        'rot': {'r1': 270, 'p1': 300, 'r2': 150},
    },
    {
        'id': 'rxn_080',
        'title': 'Stork Enamine Alkylation: Attack on Methyl Iodide',
        'category': 'Enolate Chemistry',
        'desc': (
            'Enamine nitrogen lone pair pushes electrons to beta-carbon, which'
            ' attacks methyl iodide via clean SN2 displacement.'
        ),
        'r1': {
            'smi': 'C1=C(N2CCCC2)CCCC1',
            'name': 'Cyclohexanone pyrrolidine enamine',
            'rot': 30,
            'w': 4.5,
        },
        'r2': {'smi': 'CI', 'name': 'Methyl iodide', 'rot': 180, 'w': 2.5},
        'cond': {
            'top': '$\\text{dioxane}, \\Delta$',
            'bottom': '$\\text{iminium salt}$',
        },
        'p1': {
            'smi': 'CC1CCCCC1=[N+]1CCCC1',
            'name': 'Iminium intermediate',
            'rot': 0,
            'w': 4.5,
        },
        'p2': {'smi': '[I-]', 'name': 'Iodide', 'rot': 0, 'w': 1.8},
        'arrows': [
            '\\draw[chemarrow] (1.19, -0.07) to[out=40, in=0] (6.05, -0.00);',
            (
                '\\draw[chemarrow] (4.71, -0.14) to[out=120, in=205] (3.53,'
                ' -0.24);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (-0.556, -0.292) .. controls (-0.693,'
                ' -0.212) and (-0.330, 0.057) .. (-0.193, -0.022);'
            ),
            (
                '\\draw[chemarrow] (0.571, 0.244) .. controls (1.540, -1.484)'
                ' and (2.113, -1.547) .. (3.437, -0.074);'
            ),
            (
                '\\draw[chemarrow] (4.583, 0.080) .. controls (4.745, 0.444)'
                ' and (5.554, 0.547) .. (5.716, 0.183);'
            ),
        ],
        'rot': {'r1': 300, 'p1': 270, 'r2': 0, 'p2': 0},
    },
    {
        'id': 'rxn_081',
        'title': (
            'Diels-Alder [4+2] Cycloaddition: Butadiene + Maleic Anhydride'
        ),
        'category': 'Pericyclic Reactions',
        'desc': (
            'Concerted [4s + 2s] cycloaddition with suprafacial overlap of'
            ' HOMO(diene) and LUMO(dienophile).'
        ),
        'r1': {
            'smi': 'C=CC=C',
            'name': 's-cis-1,3-Butadiene',
            'rot': 90,
            'w': 3.5,
        },
        'r2': {
            'smi': 'O=C1C=CC(=O)O1',
            'name': 'Maleic anhydride',
            'rot': 90,
            'w': 3.5,
        },
        'cond': {
            'top': '$\\text{benzene}, 80^\\circ\\text{C}$',
            'bottom': '$\\text{concerted } [4+2]$',
        },
        'p1': {
            'smi': 'O=C1C2CC=CCC2C(=O)O1',
            'name': 'Tetrahydrophthalic anhydride',
            'rot': 0,
            'w': 4.5,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (0.09, -1.03) to[out=40, in=144] (3.91, 0.50);',
            (
                '\\draw[chemarrow] (4.11, -0.03) to[out=140, in=113] (-0.19,'
                ' 1.57);'
            ),
            '\\draw[chemarrow] (-0.12, 0.09) to[out=135, in=55] (0.27, 0.95);',
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (1.103, -0.079) .. controls (2.098, 1.191)'
                ' and (2.599, 1.193) .. (3.606, -0.067);'
            ),
            (
                '\\draw[chemarrow] (3.879, -0.455) .. controls (1.996, 2.695)'
                ' and (0.935, 2.780) .. (-1.425, -0.031);'
            ),
            (
                '\\draw[chemarrow] (-0.924, -0.114) .. controls (-0.897,'
                ' -0.451) and (-0.179, -0.355) .. (-0.036, -0.048);'
            ),
        ],
        'rot': {'r1': 180, 'p1': 300, 'r2': 240},
    },
    {
        'id': 'rxn_082',
        'title': (
            'Diels-Alder [4+2] Cycloaddition: Cyclopentadiene + Methyl Acrylate'
        ),
        'category': 'Pericyclic Reactions',
        'desc': (
            'Endo-selective Diels-Alder cycloaddition driven by secondary'
            ' orbital interactions.'
        ),
        'r1': {
            'smi': 'C1=CC=CC1',
            'name': 'Cyclopentadiene',
            'rot': 105,
            'w': 3.2,
        },
        'r2': {
            'smi': 'C=CC(=O)OC',
            'name': 'Methyl acrylate',
            'rot': 270,
            'w': 3.5,
        },
        'cond': {
            'top': '$\\text{toluene}, 25^\\circ\\text{C}$',
            'bottom': '$\\text{endo-selective}$',
        },
        'p1': {
            'smi': 'COC(=O)C1CC2C=CC1C2',
            'name': 'Endo-norbornene ester',
            'rot': 0,
            'w': 4.5,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (0.00, 1.55) to[out=40, in=-60] (4.48, -1.56);',
            '\\draw[chemarrow] (3.97, -0.36) to[out=140, in=18] (1.53, 0.37);',
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (0.809, 0.650) .. controls (1.814, 1.777)'
                ' and (2.222, 1.715) .. (2.848, 0.341);'
            ),
            (
                '\\draw[chemarrow] (3.307, 0.392) .. controls (2.197, -2.305)'
                ' and (1.333, -2.546) .. (-1.011, -0.811);'
            ),
            (
                '\\draw[chemarrow] (-0.865, 0.000) .. controls (-0.491, -0.125)'
                ' and (-0.233, 0.600) .. (-0.430, 0.941);'
            ),
        ],
        'rot': {'r1': 0, 'p1': 90, 'r2': 0},
    },
    {
        'id': 'rxn_083',
        'title': 'Retro-Diels-Alder of Dicyclopentadiene',
        'category': 'Pericyclic Reactions',
        'desc': (
            'Thermal cycloreversion cracking dicyclopentadiene into two'
            ' molecules of monomeric cyclopentadiene.'
        ),
        'r1': {
            'smi': 'C1C=CC2C1C3CC2C=C3',
            'name': 'Dicyclopentadiene',
            'rot': 0,
            'w': 4.5,
        },
        'r2': None,
        'cond': {
            'top': '$160-180^\\circ\\text{C}$',
            'bottom': '$\\text{retro-}[4+2]$',
        },
        'p1': {
            'smi': 'C1=CCC=C1.C1=CCC=C1',
            'name': 'Cyclopentadiene (2 equiv)',
            'rot': 0,
            'w': 3.5,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (0.73, -0.90) to[out=120, in=205] (0.21, -0.76);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (-0.584, 0.321) .. controls (-0.360, 0.487)'
                ' and (0.057, 0.150) .. (-0.000, -0.123);'
            ),
            (
                '\\draw[chemarrow] (0.584, 0.321) .. controls (0.236, 0.319)'
                ' and (0.257, 1.200) .. (0.603, 1.237);'
            ),
            (
                '\\draw[chemarrow] (0.000, 1.504) .. controls (0.079, 1.266)'
                ' and (-0.386, 1.113) .. (-0.603, 1.237);'
            ),
        ],
        'rot': {'r1': 90, 'p1': 90},
    },
    {
        'id': 'rxn_084',
        'title': 'Cope Rearrangement: [3,3]-Sigmatropic Shift',
        'category': 'Pericyclic Reactions',
        'desc': (
            'Thermal concerted [3,3]-sigmatropic rearrangement of 1,5-hexadiene'
            ' via chair-like transition state.'
        ),
        'r1': {'smi': 'C=CCCC=C', 'name': '1,5-Hexadiene', 'rot': 0, 'w': 4.2},
        'r2': None,
        'cond': {
            'top': '$300^\\circ\\text{C}$',
            'bottom': '$[3,3]\\text{-sigmatropic}$',
        },
        'p1': {
            'smi': 'C=CCCC=C',
            'name': 'Degenerate hexadiene',
            'rot': 0,
            'w': 4.2,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (-1.89, 0.20) to[out=40, in=27] (2.10, 0.23);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (1.393, 0.878) .. controls (0.724, 1.842)'
                ' and (0.276, 1.998) .. (-0.848, 1.659);'
            ),
            (
                '\\draw[chemarrow] (-1.133, 0.638) .. controls (-0.593, 0.542)'
                ' and (-0.873, -0.625) .. (-1.398, -0.781);'
            ),
            (
                '\\draw[chemarrow] (-0.000, -1.540) .. controls (-0.191,'
                ' -0.966) and (0.938, -0.390) .. (1.398, -0.781);'
            ),
        ],
        'rot': {'p1': 0, 'r1': 0},
    },
    {
        'id': 'rxn_085',
        'title': 'Anionic Oxy-Cope Rearrangement',
        'category': 'Pericyclic Reactions',
        'desc': (
            'Potassium alkoxide accelerates [3,3]-sigmatropic shift by 10^12 to'
            ' 10^17 rate enhancement at room temp.'
        ),
        'r1': {
            'smi': 'C=CCC(C)([O-])C=C',
            'name': 'Potassium 1,5-dien-3-olate',
            'rot': 0,
            'w': 4.8,
        },
        'r2': None,
        'cond': {
            'top': '$\\text{KH}, 18\\text{-crown-}6, 25^\\circ\\text{C}$',
            'bottom': '$\\text{rate accel. } 10^{14}$',
        },
        'p1': {
            'smi': 'C=CCCC=C(C)[O-]',
            'name': 'Enolate enone precursor',
            'rot': 0,
            'w': 4.8,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (2.16, 0.66) to[out=140, in=11] (-0.66, 0.68);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (-0.967, 1.399) .. controls (1.002, 1.909)'
                ' and (1.493, 1.528) .. (1.487, -0.506);'
            ),
            (
                '\\draw[chemarrow] (0.926, -0.695) .. controls (0.911, -1.015)'
                ' and (0.229, -0.945) .. (0.084, -0.660);'
            ),
            (
                '\\draw[chemarrow] (-0.689, -0.479) .. controls (-0.364,'
                ' -0.357) and (-0.506, 0.364) .. (-0.848, 0.419);'
            ),
        ],
        'rot': {'r1': 60, 'p1': 30},
    },
    {
        'id': 'rxn_086',
        'title': 'Claisen Rearrangement: Allyl Phenyl Ether',
        'category': 'Pericyclic Reactions',
        'desc': (
            'Concerted [3,3]-sigmatropic rearrangement of allyl phenyl ether'
            ' via chair transition state to ortho-allylphenol.'
        ),
        'r1': {
            'smi': 'C=CCOc1ccccc1',
            'name': 'Allyl phenyl ether',
            'rot': 0,
            'w': 4.5,
        },
        'r2': None,
        'cond': {
            'top': '$200^\\circ\\text{C}$',
            'bottom': '$[3,3]\\text{-shift}$',
        },
        'p1': {
            'smi': 'O=C1C=CC=CC1CC=C',
            'name': 'Dienone intermediate',
            'rot': 0,
            'w': 4.5,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (-2.02, 0.59) to[out=40, in=139] (0.54, 0.66);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (0.068, -0.478) .. controls (-2.547, 0.440)'
                ' and (-2.534, 0.855) .. (0.138, 1.596);'
            ),
            (
                '\\draw[chemarrow] (0.152, 1.371) .. controls (0.313, 1.303)'
                ' and (0.137, 0.972) .. (-0.037, 0.951);'
            ),
            (
                '\\draw[chemarrow] (-0.145, 0.322) .. controls (0.030, 0.248)'
                ' and (-0.071, -0.141) .. (-0.259, -0.163);'
            ),
        ],
        'rot': {'r1': 90, 'p1': 0},
    },
    {
        'id': 'rxn_087',
        'title': 'Ireland-Claisen Rearrangement of Silyl Ketenacetal',
        'category': 'Pericyclic Reactions',
        'desc': (
            'Diastereoselective [3,3]-sigmatropic shift of ester silyl'
            ' ketenacetal giving gamma,delta-unsaturated acid.'
        ),
        'r1': {
            'smi': 'C=CCOC(=C)O[Si](C)(C)C',
            'name': 'Silyl ketene acetal',
            'rot': 0,
            'w': 5.0,
        },
        'r2': None,
        'cond': {
            'top': '$\\text{THF}, \\text{reflux}$',
            'bottom': '$\\text{stereoselective}$',
        },
        'p1': {
            'smi': 'C=CCCC(=O)O[Si](C)(C)C',
            'name': 'Silyl ester adduct',
            'rot': 0,
            'w': 5.0,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (2.25, 1.00) to[out=140, in=-126] (-0.92,'
            ' -1.32);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (0.687, -0.377) .. controls (1.672, 1.240)'
                ' and (1.446, 1.684) .. (-0.440, 1.842);'
            ),
            (
                '\\draw[chemarrow] (-0.412, 1.576) .. controls (-0.564, 1.467)'
                ' and (-0.303, 1.164) .. (-0.118, 1.184);'
            ),
            (
                '\\draw[chemarrow] (0.142, 0.565) .. controls (-0.022, 0.447)'
                ' and (0.175, 0.066) .. (0.375, 0.088);'
            ),
        ],
        'rot': {'r1': 270, 'p1': 270},
    },
    {
        'id': 'rxn_088',
        'title': '[2,3]-Wittig Rearrangement of Allylic Ether',
        'category': 'Pericyclic Reactions',
        'desc': (
            'Base-generated carbanion undergoes concerted [2,3]-sigmatropic'
            ' shift to give homoallylic alkoxide.'
        ),
        'r1': {
            'smi': 'C=CCO[CH-]c1ccccc1',
            'name': 'Allylic ether carbanion',
            'rot': 0,
            'w': 4.8,
        },
        'r2': None,
        'cond': {
            'top': '$n\\text{-BuLi}, \\text{THF}, -78^\\circ\\text{C}$',
            'bottom': '$[2,3]\\text{-sigmatropic}$',
        },
        'p1': {
            'smi': 'C=CCC([O-])c1ccccc1',
            'name': 'Homoallylic alkoxide',
            'rot': 0,
            'w': 4.8,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (2.25, 0.66) to[out=140, in=102] (-0.24, 0.06);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (-0.500, 0.120) .. controls (-1.187, 0.961)'
                ' and (-1.111, 1.282) .. (-0.119, 1.726);'
            ),
            (
                '\\draw[chemarrow] (-0.067, 1.515) .. controls (-0.223, 1.513)'
                ' and (-0.186, 1.181) .. (-0.049, 1.105);'
            ),
            (
                '\\draw[chemarrow] (-0.034, 0.615) .. controls (0.058, 0.607)'
                ' and (0.211, 0.481) .. (0.143, 0.417);'
            ),
        ],
        'rot': {'r1': 90, 'p1': 90},
    },
    {
        'id': 'rxn_089',
        'title': '1,3-Dipolar Cycloaddition: Diazomethane on Methyl Acrylate',
        'category': 'Pericyclic Reactions',
        'desc': 'Concerted [3+2] cycloaddition yielding pyrazoline ring.',
        'r1': {
            'smi': 'C=CC(=O)OC',
            'name': 'Methyl acrylate',
            'rot': 90,
            'w': 3.8,
        },
        'r2': {
            'smi': '[CH2-][N+]#N',
            'name': 'Diazomethane',
            'rot': 180,
            'w': 2.5,
        },
        'cond': {
            'top': '$\\text{Et}_2\\text{O}, 0^\\circ\\text{C}$',
            'bottom': '$[3+2]\\text{ cycloaddition}$',
        },
        'p1': {
            'smi': 'COC(=O)C1CCN=N1',
            'name': '1-Pyrazoline',
            'rot': 0,
            'w': 4.2,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (5.26, 0.32) to[out=140, in=-60] (-0.09, -1.67);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.135, 0.797) .. controls (2.792, 1.741)'
                ' and (2.540, 1.790) .. (1.872, 1.040);'
            ),
            (
                '\\draw[chemarrow] (1.511, 0.521) .. controls (3.453, 2.051)'
                ' and (4.187, 1.867) .. (5.180, -0.396);'
            ),
            (
                '\\draw[chemarrow] (4.774, -0.538) .. controls (4.624, -0.647)'
                ' and (4.304, -0.443) .. (4.323, -0.259);'
            ),
        ],
        'rot': {'r1': 150, 'p1': 300, 'r2': 30},
    },
    {
        'id': 'rxn_090',
        'title': 'Click Chemistry: Azide-Alkyne 1,3-Dipolar Cycloaddition',
        'category': 'Pericyclic Reactions',
        'desc': (
            'Copper(I)-catalyzed azide-alkyne cycloaddition (CuAAC) forming'
            ' 1,4-disubstituted 1,2,3-triazole.'
        ),
        'r1': {'smi': 'CCCC#C', 'name': '1-Pentyne', 'rot': 0, 'w': 3.8},
        'r2': {
            'smi': '[N-]=[N+]=NCc1ccccc1',
            'name': 'Benzyl azide',
            'rot': 180,
            'w': 4.2,
        },
        'cond': {
            'top': '$\\text{CuSO}_4, \\text{Na ascorbate}$',
            'bottom': '$\\text{CuAAC 1,4-regioselective}$',
        },
        'p1': {
            'smi': 'CCCc1cn(Cc2ccccc2)nn1',
            'name': '1,4-Disubstituted triazole',
            'rot': 0,
            'w': 5.0,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (6.92, 0.39) to[out=140, in=17] (1.90, 0.27);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (3.087, 1.189) .. controls (1.767, 2.027)'
                ' and (1.339, 1.858) .. (0.950, 0.343);'
            ),
            (
                '\\draw[chemarrow] (1.293, 0.654) .. controls (2.526, 2.665)'
                ' and (3.186, 2.695) .. (4.592, 0.800);'
            ),
            (
                '\\draw[chemarrow] (4.447, 0.646) .. controls (4.372, 0.563)'
                ' and (4.139, 0.530) .. (4.127, 0.642);'
            ),
        ],
        'rot': {'r1': 330, 'p1': 330, 'r2': 30},
    },
    {
        'id': 'rxn_091',
        'title': 'Kinetic Enolization of 3-Substituted Cyclohexanone by LDA',
        'category': 'Total Synthesis Steps',
        'desc': (
            'Kinetic alpha-deprotonation of (3R)-3-methylcyclohexanone by LDA'
            ' generating conjugated lithium enolate.'
        ),
        'r1': {
            'smi': 'CC1CCCC(=O)C1',
            'name': '(3R)-3-methylcyclohexan-1-one',
            'rot': 120,
            'w': 4.5,
            'add_h': 7,
        },
        'r2': {'smi': 'CC(C)[N-]C(C)C', 'name': 'LDA', 'rot': 0, 'w': 3.8},
        'cond': {
            'top': '$\\text{THF}, -78^\\circ\\text{C}$',
            'bottom': '$-\\text{HN}(i\\text{-Pr})_2$',
        },
        'p1': {
            'smi': 'CC1CCCC([O-])=C1',
            'name': 'Lithium Enolate',
            'rot': 120,
            'w': 4.5,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (5.35, 0.99) to[out=140, in=30] (1.65, 0.63);',
            '\\draw[chemarrow] (0.82, 0.33) to[out=135, in=55] (0.03, 0.35);',
            (
                '\\draw[chemarrow] (-0.68, 1.03) to[out=120, in=205] (-0.78,'
                ' 1.65);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (4.522, 0.196) .. controls (3.560, 1.943)'
                ' and (2.976, 2.019) .. (1.602, 0.572);'
            ),
            (
                '\\draw[chemarrow] (0.855, 0.266) .. controls (0.788, 0.596)'
                ' and (0.074, 0.569) .. (-0.031, 0.249);'
            ),
            (
                '\\draw[chemarrow] (-0.702, 0.970) .. controls (-0.871, 1.045)'
                ' and (-0.891, 1.421) .. (-0.722, 1.496);'
            ),
        ],
        'rot': {'r1': 240, 'p1': 240, 'r2': 270},
    },
    {
        'id': 'rxn_092',
        'title': 'Diastereoselective Enolate C-Alkylation by Methallyl Bromide',
        'category': 'Total Synthesis Steps',
        'desc': (
            'Stereoselective equatorial alkylation by methallyl bromide anti to'
            ' the (3R)-methyl substituent (>15:1 dr).'
        ),
        'r1': {
            'smi': 'CC1CCCC([O-])=C1',
            'name': 'Lithium Enolate',
            'rot': 120,
            'w': 4.5,
        },
        'r2': {
            'smi': 'CC(=C)CBr',
            'name': 'Methallyl bromide',
            'rot': 0,
            'w': 4.0,
        },
        'cond': {
            'top': '$-78^\\circ\\text{C} \\rightarrow -20^\\circ\\text{C}$',
            'bottom': '$-\\text{LiBr}$',
        },
        'p1': {
            'smi': 'CC(=C)CC1C(=O)CCCC1C',
            'name': '2-Methallyl-3-methylcyclohexan-1-one',
            'rot': 120,
            'w': 4.8,
        },
        'p2': {'smi': '[Br-]', 'name': 'Bromide', 'rot': 0, 'w': 1.8},
        'arrows': [
            (
                '\\draw[chemarrow] (-0.80, 1.58) to[out=120, in=190] (-0.62,'
                ' 1.03);'
            ),
            '\\draw[chemarrow] (0.48, 0.24) to[out=40, in=90] (5.86, 1.19);',
            '\\draw[chemarrow] (6.47, 0.75) to[out=120, in=205] (7.18, 0.53);',
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (1.971, -1.359) .. controls (2.172, -1.242)'
                ' and (1.655, -0.826) .. (1.453, -0.942);'
            ),
            (
                '\\draw[chemarrow] (0.537, -0.479) .. controls (2.510, -3.203)'
                ' and (3.271, -3.087) .. (4.343, 0.101);'
            ),
            (
                '\\draw[chemarrow] (4.551, 0.737) .. controls (4.730, 0.816)'
                ' and (4.888, 1.194) .. (4.709, 1.274);'
            ),
        ],
        'rot': {'r1': 30, 'p1': 30, 'r2': 240, 'p2': 240},
    },
    {
        'id': 'rxn_093',
        'title': 'Gilman Cuprate 1,4-Conjugate Addition to Cyclohexenone',
        'category': 'Total Synthesis Steps',
        'desc': (
            '1,4-Conjugate addition of lithium dimethylcuprate Me2CuLi to'
            ' cyclohexenone establishing 3-methyl stereocenter.'
        ),
        'r1': {
            'smi': 'O=C1C=CCCC1',
            'name': 'Cyclohex-2-en-1-one',
            'rot': 90,
            'w': 4.0,
        },
        'r2': {'smi': 'C[Cu-]C', 'name': 'Me2CuLi cuprate', 'rot': 0, 'w': 3.0},
        'cond': {
            'top': '$\\text{Et}_2\\text{O}, -78^\\circ\\text{C}$',
            'bottom': '$1,4\\text{-addition}$',
        },
        'p1': {
            'smi': 'CC1CCCC([O-])=C1',
            'name': 'Lithium enolate intermediate',
            'rot': 90,
            'w': 4.5,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (3.35, -0.13) to[out=140, in=90] (-0.00, -1.30);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (4.171, 0.760) .. controls (3.225, 1.990)'
                ' and (2.777, 1.973) .. (1.928, 0.674);'
            ),
            (
                '\\draw[chemarrow] (1.078, 0.667) .. controls (1.003, 0.298)'
                ' and (0.221, 0.498) .. (0.104, 0.856);'
            ),
            (
                '\\draw[chemarrow] (-0.922, 1.030) .. controls (-0.898, 1.255)'
                ' and (-1.296, 1.488) .. (-1.478, 1.355);'
            ),
        ],
        'rot': {'r1': 210, 'p1': 90, 'r2': 30},
    },
    {
        'id': 'rxn_094',
        'title': 'Intramolecular Hydroxy-Ketone Hemiacetal Cyclization',
        'category': 'Total Synthesis Steps',
        'desc': (
            'Intramolecular attack of primary hydroxyl on ketone carbonyl'
            ' forming a six-membered cyclic hemiacetal.'
        ),
        'r1': {
            'smi': 'CC(=O)CCCCO',
            'name': '6-Hydroxyhexan-2-one',
            'rot': 0,
            'w': 4.8,
        },
        'r2': None,
        'cond': {
            'top': '$p\\text{-TsOH}, \\text{benzene}$',
            'bottom': '$\\text{lactol closure}$',
        },
        'p1': {
            'smi': 'CC1([O-])CCCC[OH+]1',
            'name': 'Cyclic hemiacetal',
            'rot': 0,
            'w': 3.8,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (-1.32, 0.73) to[out=40, in=93] (1.58, 0.44);',
            (
                '\\draw[chemarrow] (1.75, -0.12) to[out=120, in=205] (1.87,'
                ' -0.53);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (1.644, -0.172) .. controls (0.455, 1.453)'
                ' and (-0.149, 1.460) .. (-1.379, -0.134);'
            ),
            (
                '\\draw[chemarrow] (-1.695, 0.116) .. controls (-1.810, 0.156)'
                ' and (-1.875, 0.397) .. (-1.769, 0.456);'
            ),
        ],
        'rot': {'r1': 180, 'p1': 150},
    },
    {
        'id': 'rxn_095',
        'title': (
            'Hock Rearrangement and Acid-Catalyzed Cleavage of Allylic'
            ' Hydroperoxide'
        ),
        'category': 'Total Synthesis Steps',
        'desc': (
            'Acid-catalyzed rearrangement of allylic hydroperoxide with C-C'
            ' bond cleavage to oxocarbenium intermediate.'
        ),
        'r1': {
            'smi': 'CC(C)([O][OH2+])C=C',
            'name': 'Protonated allylic hydroperoxide',
            'rot': 0,
            'w': 4.2,
        },
        'r2': None,
        'cond': {
            'top': '$\\text{TFA}, \\text{DCM}, 0^\\circ\\text{C}$',
            'bottom': '$\\text{Hock rearrangement}$',
        },
        'p1': {
            'smi': 'C=CO[C+](C)C',
            'name': 'Oxocarbenium intermediate',
            'rot': 0,
            'w': 4.2,
        },
        'p2': {'smi': 'O', 'name': 'Water', 'rot': 0, 'w': 2.0},
        'arrows': [
            '\\draw[chemarrow] (0.67, 1.09) to[out=120, in=205] (1.33, 1.18);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (-0.307, -0.589) .. controls (-0.500,'
                ' -0.525) and (0.032, -0.032) .. (-0.167, -0.074);'
            ),
            (
                '\\draw[chemarrow] (0.135, 1.258) .. controls (0.117, 1.428)'
                ' and (0.317, 1.726) .. (0.455, 1.625);'
            ),
        ],
        'rot': {'r1': 270, 'p1': 120, 'p2': 0},
    },
    {
        'id': 'rxn_096',
        'title': 'Acid-Catalyzed Endoperoxide 1,2,4-Trioxane Cyclization',
        'category': 'Total Synthesis Steps',
        'desc': (
            'Acid-catalyzed intramolecular cyclization of peroxy-ester'
            ' precursor yielding a 1,2,4-trioxane ring core.'
        ),
        'r1': {
            'smi': 'CC(=O)OC(OO)c1ccccc1',
            'name': 'Peroxy-ester precursor',
            'rot': 90,
            'w': 4.8,
        },
        'r2': None,
        'cond': {
            'top': '$\\text{TFA}, 25^\\circ\\text{C}$',
            'bottom': '$\\text{trioxane cyclization}$',
        },
        'p1': {
            'smi': 'CC1([O-])OC(c2ccccc2)O[OH+]1',
            'name': '1,2,4-Trioxane core adduct',
            'rot': 90,
            'w': 4.2,
        },
        'p2': None,
        'arrows': [
            (
                '\\draw[chemarrow] (-0.64, -0.56) to[out=140, in=67] (-1.56,'
                ' 0.67);'
            ),
            (
                '\\draw[chemarrow] (-1.94, 0.12) to[out=120, in=205] (-2.22,'
                ' -0.30);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (-0.657, 1.413) .. controls (-0.377, 2.301)'
                ' and (-0.201, 2.316) .. (0.222, 1.486);'
            ),
            (
                '\\draw[chemarrow] (0.028, 1.416) .. controls (-0.052, 1.340)'
                ' and (-0.276, 1.380) .. (-0.279, 1.491);'
            ),
        ],
        'rot': {'r1': 90, 'p1': 90},
    },
    {
        'id': 'rxn_097',
        'title': 'Mukaiyama Aldol Reaction of Silyl Enol Ether',
        'category': 'Named Organic Reactions',
        'desc': (
            'Lewis acid-mediated Mukaiyama aldol coupling of a trimethylsilyl'
            ' enol ether with an aromatic aldehyde.'
        ),
        'r1': {
            'smi': 'C=C(O[Si](C)(C)C)C1CCCC1',
            'name': 'Trimethylsilyl enol ether',
            'rot': 0,
            'w': 4.8,
        },
        'r2': {
            'smi': 'O=Cc1ccccc1',
            'name': 'Benzaldehyde',
            'rot': 60,
            'w': 3.8,
        },
        'cond': {
            'top': '$\\text{BF}_3\\cdot\\text{OEt}_2, -78^\\circ\\text{C}$',
            'bottom': '$\\text{Mukaiyama aldol coupling}$',
        },
        'p1': {
            'smi': 'C[Si](C)(C)O[C+](CC([O-])c1ccccc1)C1CCCC1',
            'name': 'Mukaiyama aldol adduct',
            'rot': 0,
            'w': 5.2,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (0.21, -0.67) to[out=40, in=133] (5.13, -1.02);',
            (
                '\\draw[chemarrow] (4.87, -1.43) to[out=120, in=205] (4.51,'
                ' -1.61);'
            ),
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (0.858, -0.129) .. controls (1.948, 1.069)'
                ' and (2.565, 1.156) .. (3.944, 0.305);'
            ),
            (
                '\\draw[chemarrow] (4.228, 0.765) .. controls (4.373, 0.873)'
                ' and (4.234, 1.218) .. (4.056, 1.254);'
            ),
        ],
        'rot': {'r1': 300, 'p1': 210, 'r2': 30},
    },
    {
        'id': 'rxn_098',
        'title': 'Cationic Aza-Cope [3,3]-Sigmatropic Rearrangement',
        'category': 'Named Organic Reactions',
        'desc': (
            'Acid-catalyzed aza-Cope [3,3]-sigmatropic rearrangement of an'
            ' allylic enamine/iminium system forging an unsaturated aldimine.'
        ),
        'r1': {
            'smi': 'C=CCN(C)C=C',
            'name': 'Allylic amino-diene',
            'rot': 0,
            'w': 4.2,
        },
        'r2': None,
        'cond': {
            'top': '$\\text{CSA}, 80^\\circ\\text{C}$',
            'bottom': '$\\text{aza-Cope [3,3]-rearrangement}$',
        },
        'p1': {
            'smi': 'C=CCCC=NC',
            'name': 'N-Methyl aldimine',
            'rot': 0,
            'w': 4.2,
        },
        'p2': None,
        'arrows': [
            '\\draw[chemarrow] (0.58, 0.19) to[out=140, in=68] (-0.86, 0.70);'
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (-1.351, 0.472) .. controls (-0.371, 2.253)'
                ' and (0.242, 2.349) .. (1.715, 0.948);'
            ),
            (
                '\\draw[chemarrow] (1.431, 0.594) .. controls (1.549, 0.357)'
                ' and (1.028, 0.138) .. (0.808, 0.285);'
            ),
            (
                '\\draw[chemarrow] (0.094, -0.120) .. controls (0.146, 0.161)'
                ' and (-0.426, 0.364) .. (-0.610, 0.145);'
            ),
        ],
        'rot': {'r1': 0, 'p1': 30},
    },
    {
        'id': 'rxn_099',
        'title': 'Pictet-Spengler Tetrahydroisoquinoline Cyclization',
        'category': 'Named Organic Reactions',
        'desc': (
            'Intramolecular electrophilic aromatic cyclization of an'
            ' electron-rich phenethyliminium ion forging a'
            ' tetrahydroisoquinoline core.'
        ),
        'r1': {
            'smi': 'c1cc(O)c(O)cc1CC[NH+]=C',
            'name': 'Dopamine iminium ion',
            'rot': 0,
            'w': 4.5,
        },
        'r2': None,
        'cond': {
            'top': '$\\text{H}^+, 50^\\circ\\text{C}$',
            'bottom': '$\\text{Pictet-Spengler cyclization}$',
        },
        'p1': {
            'smi': 'OC1=C(O)C2CNCC[C+]2C=C1',
            'name': 'Tetrahydroisoquinoline core',
            'rot': 0,
            'w': 4.5,
        },
        'p2': None,
        'arrows': [],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (0.425, 0.016) .. controls (0.895, -0.179)'
                ' and (0.973, -0.349) .. (0.815, -0.832);'
            ),
            (
                '\\draw[chemarrow] (1.198, -0.957) .. controls (1.299, -1.045)'
                ' and (1.579, -1.028) .. (1.576, -0.894);'
            ),
        ],
        'rot': {'p1': 0, 'r1': 0},
    },
    {
        'id': 'rxn_100',
        'title': (
            'Halolactonization (Iodolactonization) of Unsaturated Carboxylate'
        ),
        'category': 'Named Organic Reactions',
        'desc': (
            'Intramolecular iodolactonization of cyclopentene acetate yielding'
            ' a bicyclic iodolactone.'
        ),
        'r1': {
            'smi': 'O=C([O-])CC1C=CCC1',
            'name': 'Cyclopent-2-ene-1-acetate',
            'rot': 0,
            'w': 4.5,
        },
        'r2': {'smi': 'II', 'name': 'Iodine (KI, NaHCO3)', 'rot': 90, 'w': 2.8},
        'cond': {
            'top': '$\\text{I}_2, \\text{KI}, \\text{NaHCO}_3$',
            'bottom': '$\\text{iodolactonization}$',
        },
        'p1': {
            'smi': 'O=C1CC2CCC(I)C2O1',
            'name': 'Bicyclic iodolactone',
            'rot': 0,
            'w': 4.2,
        },
        'p2': {'smi': '[I-]', 'name': 'Iodide', 'rot': 0, 'w': 2.0},
        'arrows': [
            '\\draw[chemarrow] (0.65, 0.33) to[out=40, in=-90] (4.85, -1.16);',
            '\\draw[chemarrow] (-0.97, -0.27) to[out=40, in=122] (1.53, 0.14);',
        ],
        'arrows_tikz': [
            (
                '\\draw[chemarrow] (1.742, -0.585) .. controls (2.917, -0.285)'
                ' and (3.498, -0.308) .. (4.647, -0.697);'
            ),
            (
                '\\draw[chemarrow] (4.930, -0.084) .. controls (5.202, 0.037)'
                ' and (5.305, 0.639) .. (5.033, 0.760);'
            ),
            (
                '\\draw[chemarrow] (-1.667, -0.645) .. controls (-0.583,'
                ' -1.121) and (-0.021, -1.188) .. (1.145, -0.979);'
            ),
        ],
        'rot': {'p1': 0, 'r1': 0, 'r2': 90},
    },
]
