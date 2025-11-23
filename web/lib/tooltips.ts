// Tooltip content definitions for InfoTooltip components throughout the app

export const TOOLTIPS = {
    drugScore: {
        title: 'Drug Score',
        content: 'Composite metric (0-10) evaluating drug-likeness. Combines Lipinski\'s Rule of Five (4 points), Veber\'s Rule (2 points), Blood-Brain Barrier permeability (2 points), and synthetic accessibility (2 points). Score 7-10 indicates excellent drug candidate with good oral bioavailability and ease of synthesis.'
    },
    vqc: {
        title: 'Variational Quantum Circuit (VQC)',
        content: 'Quantum machine learning model using 4 qubits to predict molecular properties. Processes molecular embeddings through quantum gates (RX for encoding, Rot for learning, CNOT for entanglement) to capture complex molecular patterns that classical ML may miss. Currently runs on PennyLane\'s quantum simulator.'
    },
    quantumCircuit: {
        title: 'Quantum Circuit',
        content: 'Visual representation of quantum operations performed during prediction. Shows qubit wires, RX gates (blue) for encoding molecular features, Rot gates (purple) for learning patterns, CNOT gates (pink) for creating quantum entanglement, and measurement gates (green) for obtaining predictions.'
    },
    lipinski: {
        title: 'Lipinski\'s Rule of Five',
        content: 'Famous drug-likeness criteria named after Christopher Lipinski. A molecule passes if: (1) Molecular weight ≤ 500 Da, (2) LogP ≤ 5, (3) H-bond donors ≤ 5, (4) H-bond acceptors ≤ 10. Molecules passing all four rules have better oral bioavailability potential.'
    },
    logp: {
        title: 'LogP (Partition Coefficient)',
        content: 'Measures how easily a molecule dissolves in fat vs water. Higher LogP = more lipophilic (fat-soluble). Optimal range for drug candidates is 0-5. Values too low mean poor membrane permeability, too high mean poor water solubility.'
    },
    molecularWeight: {
        title: 'Molecular Weight (MW)',
        content: 'Total mass of the molecule in Daltons (Da). Important for drug-likeness - should be ≤ 500 Da for oral drugs (Lipinski\'s rule). Larger molecules have difficulty crossing biological membranes and may have poor absorption.'
    },
    hDonors: {
        title: 'H-Bond Donors (HBD)',
        content: 'Number of atoms that can donate hydrogen bonds (typically O-H or N-H groups). Affects solubility and permeability. Should be ≤ 5 for drug-likeness (Lipinski\'s rule). Too many H-donors can reduce membrane permeability.'
    },
    hAcceptors: {
        title: 'H-Bond Acceptors (HBA)',
        content: 'Number of atoms that can accept hydrogen bonds (typically O or N atoms). Affects solubility. Should be ≤ 10 for drug-likeness (Lipinski\'s rule). Higher numbers improve water solubility but may reduce membrane permeability.'
    },
    rotatableBonds: {
        title: 'Rotatable Bonds',
        content: 'Number of single bonds that can rotate freely (excluding bonds in rings and terminal methyl groups). Affects molecular flexibility and oral bioavailability. Should be ≤ 10 for good oral absorption (Veber\'s rule).'
    },
    tpsa: {
        title: 'TPSA (Topological Polar Surface Area)',
        content: 'Sum of surface areas of polar atoms (O, N, S) in the molecule. Predicts oral bioavailability and permeability. Should be ≤ 140 Å² for good oral absorption (Veber\'s rule). Lower TPSA typically means better membrane permeability.'
    },
    admet: {
        title: 'ADMET Properties',
        content: 'Absorption, Distribution, Metabolism, Excretion, and Toxicity - the five key properties determining drug effectiveness and safety. Absorption: How well absorbed into bloodstream. Distribution: How drug spreads through body. Metabolism: How body breaks it down. Excretion: How drug leaves body. Toxicity: Potential harmful effects.'
    },
    bbbPermeability: {
        title: 'BBB Permeability',
        content: 'Blood-Brain Barrier permeability - predicts if a molecule can cross the blood-brain barrier to reach the central nervous system. Scored based on LogP (0-3 ideal), molecular weight (< 450 Da), and TPSA (< 90 Å²). Important for CNS-targeting drugs.'
    },
    syntheticAccessibility: {
        title: 'Synthetic Accessibility Score',
        content: 'Measures how easy it is to synthesize the molecule (scale 1-10, lower is easier). Based on molecular complexity, number of rotatable bonds, and ring structures. Scores ≤ 5 indicate easier synthesis, while > 7 suggest difficult synthesis.'
    },
    veber: {
        title: 'Veber\'s Rule',
        content: 'Rule for predicting oral bioavailability. A molecule passes if: (1) Rotatable bonds ≤ 10, (2) TPSA ≤ 140 Å². Molecules passing Veber\'s rule have better chances of good oral absorption. Part of the drug score calculation.'
    },
    xgboost: {
        title: 'XGBoost (Classical ML)',
        content: 'Gradient boosting machine learning algorithm for solubility prediction. Trained on molecular descriptors (MW, LogP, H-bonds, etc.) plus 1536-dimensional embeddings. Achieves R² = 0.849. Represents the classical (non-quantum) approach to molecular property prediction.'
    },
    rdkit: {
        title: 'RDKit',
        content: 'Open-source cheminformatics toolkit used to calculate molecular descriptors, generate 3D structures, and analyze molecular properties. Powers feature extraction including molecular weight, LogP, hydrogen bonds, TPSA, and other descriptors used in predictions.'
    },
    pennylane: {
        title: 'PennyLane',
        content: 'Quantum machine learning framework that powers the Variational Quantum Circuit (VQC). Enables construction and training of quantum circuits for molecular property prediction. Uses quantum simulators to process molecular data through quantum gates.'
    },
    smiles: {
        title: 'SMILES Notation',
        content: 'Simplified Molecular Input Line Entry System - a text-based way to represent chemical structures. Uses atoms (C, O, N), bonds (= double, # triple), parentheses for branches, and numbers for ring connections. Example: CCO = ethanol, c1ccccc1 = benzene.'
    },
    aromaticRings: {
        title: 'Aromatic Rings',
        content: 'Number of benzene-like rings (aromatic rings) in the molecule. Aromatic rings contribute to molecular stability and can affect drug-receptor interactions. Typically ranges from 0 to several rings in larger drug molecules.'
    },
    totalAtoms: {
        title: 'Total Atoms',
        content: 'Simple count of all atoms in the molecule. Used as a molecular descriptor to characterize molecular size. Larger molecules typically have more atoms, which can affect various properties like solubility and permeability.'
    },
    cyp450: {
        title: 'CYP450 (Cytochrome P450)',
        content: 'Family of enzymes responsible for drug metabolism in the liver. CYP450 Substrate indicates if the molecule is metabolized by these enzymes. CYP450 Inhibitor indicates if it blocks metabolism of other drugs, which can cause drug-drug interactions.'
    },
    hepatotoxicity: {
        title: 'Hepatotoxicity',
        content: 'Risk of liver damage. Predicted based on molecular properties like LogP. Low risk is preferred for drug candidates. High hepatotoxicity can cause serious side effects and is a common reason for drug candidate failure.'
    },
    amesMutagenicity: {
        title: 'Ames Mutagenicity',
        content: 'Risk of DNA mutations (cancer risk) as measured by the Ames test. Based on the presence of aromatic rings and other structural features. Low risk is critical for drug safety. Mutagenic compounds are typically eliminated early in drug development.'
    },
    rxGate: {
        title: 'RX Gate (Rotation)',
        content: 'Quantum rotation gate that encodes molecular features onto qubits. Rotates the qubit state around the X-axis based on input data. Part of the encoding layer that converts classical molecular data into quantum states.'
    },
    rotGate: {
        title: 'Rot Gate (Variational)',
        content: 'Three-parameter rotation gate that is trained during model optimization. Can rotate around X, Y, and Z axes. These gates form the learnable (variational) part of the quantum circuit that discovers patterns in molecular data.'
    },
    cnotGate: {
        title: 'CNOT Gate (Entanglement)',
        content: 'Controlled-NOT gate that creates quantum entanglement between qubits. The dot (•) is the control qubit, the ⊕ is the target qubit. Entanglement allows the circuit to capture complex correlations that enable powerful quantum pattern recognition.'
    },
    measurement: {
        title: 'Quantum Measurement',
        content: 'Final step where quantum states are measured to obtain a prediction. In this circuit, measures the expectation value of the Pauli-Z operator on qubit 0. The measurement converts quantum information into a classical number used for solubility prediction.'
    },
    expectationValue: {
        title: 'Expectation Value',
        content: 'The measured quantum value representing the average outcome of multiple quantum measurements. In this VQC, it measures ⟨Z⟩ (Pauli-Z expectation) on qubit 0. This value is converted to a solubility prediction by adding a trained bias term.'
    },
    inputEncoding: {
        title: 'Input Encoding',
        content: 'The molecular embedding values (4 dimensions) encoded onto the 4 qubits using RX rotation gates. These values are normalized to the range [-π, π] and represent the molecular features that the quantum circuit processes.'
    }
};

