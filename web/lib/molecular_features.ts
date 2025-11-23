// Pure TypeScript/JavaScript molecular feature extraction
// Lightweight alternative to RDKit for Vercel serverless

interface Atom {
    symbol: string;
    bonds: number;
    isAromatic?: boolean;
}

interface Bond {
    begin: number;
    end: number;
    type: number; // 1=single, 2=double, 3=triple, 4=aromatic
}

// Periodic table atomic weights
const ATOMIC_WEIGHTS: Record<string, number> = {
    'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998,
    'P': 30.974, 'S': 32.065, 'Cl': 35.453, 'Br': 79.904, 'I': 126.904,
    'B': 10.811, 'Si': 28.085, 'Se': 78.96, 'Te': 127.60
};

// Hydrogen bond donor atoms (O-H, N-H, S-H, etc.)
const H_DONOR_PATTERNS = /[ONSH]/i;

// Hydrogen bond acceptor atoms (O, N, F, Cl, etc.)
const H_ACCEPTOR_PATTERNS = /[ONFSCl]/i;

// Aromatic atoms
const AROMATIC_PATTERNS = /[cnops]/i;

export function parseSmilesBasic(smiles: string): { atoms: Atom[], bonds: Bond[] } | null {
    // Simple SMILES parser for basic features
    // This is a simplified version - for production, use a proper parser
    try {
        const atoms: Atom[] = [];
        const bonds: Bond[] = [];
        
        let i = 0;
        while (i < smiles.length) {
            const char = smiles[i];
            
            // Two-letter element
            if (char >= 'A' && char <= 'Z' && i + 1 < smiles.length && smiles[i + 1] >= 'a' && smiles[i + 1] <= 'z') {
                atoms.push({ symbol: char + smiles[i + 1], bonds: 0 });
                i += 2;
            }
            // Single letter element
            else if (char >= 'A' && char <= 'Z') {
                atoms.push({ symbol: char, bonds: 0 });
                i++;
            }
            // Lowercase aromatic
            else if (char >= 'a' && char <= 'z') {
                atoms.push({ symbol: char.toUpperCase(), bonds: 0, isAromatic: true });
                i++;
            }
            // Ring closures, branches, etc. - skip for now
            else {
                i++;
            }
        }
        
        return { atoms, bonds };
    } catch (e) {
        return null;
    }
}

export function calculateMolecularWeight(atoms: Atom[]): number {
    let mw = 0;
    for (const atom of atoms) {
        const weight = ATOMIC_WEIGHTS[atom.symbol] || 12.0; // Default to carbon
        mw += weight;
    }
    return mw;
}

export function calculateLogP(atoms: Atom[], bonds: Bond[]): number {
    // Simplified LogP calculation using fragment contributions
    // This is an approximation - real LogP uses complex fragment-based methods
    let logp = 0;
    
    for (const atom of atoms) {
        switch (atom.symbol) {
            case 'C':
                logp += 0.5; // Carbon contributes to lipophilicity
                break;
            case 'N':
                logp -= 1.0; // Nitrogen is hydrophilic
                break;
            case 'O':
                logp -= 0.5; // Oxygen is hydrophilic
                break;
            case 'F':
            case 'Cl':
            case 'Br':
            case 'I':
                logp += 0.5; // Halogens increase lipophilicity
                break;
        }
    }
    
    // Adjust for aromatic rings (slightly more lipophilic)
    const aromaticCount = atoms.filter(a => a.isAromatic).length;
    logp += aromaticCount * 0.2;
    
    return logp;
}

export function countHBondDonors(atoms: Atom[]): number {
    // Count atoms that can donate H-bonds (O, N, S with H)
    // Simplified: count O, N, S atoms (assuming they might have H)
    return atoms.filter(atom => 
        atom.symbol === 'O' || atom.symbol === 'N' || atom.symbol === 'S'
    ).length;
}

export function countHBondAcceptors(atoms: Atom[]): number {
    // Count atoms that can accept H-bonds (O, N, F, Cl, etc.)
    return atoms.filter(atom => 
        ['O', 'N', 'F', 'Cl', 'S'].includes(atom.symbol)
    ).length;
}

export function countRotatableBonds(bonds: Bond[], atoms: Atom[]): number {
    // Simplified: count single bonds (excluding ring bonds, terminal, etc.)
    // This is an approximation
    return Math.max(0, bonds.filter(b => b.type === 1).length - atoms.length / 2);
}

export function countAromaticRings(atoms: Atom[]): number {
    // Count aromatic atoms and estimate rings
    const aromaticAtoms = atoms.filter(a => a.isAromatic).length;
    // Rough estimate: ~6 atoms per ring
    return Math.floor(aromaticAtoms / 6);
}

export function calculateTPSA(atoms: Atom[]): number {
    // Topological Polar Surface Area approximation
    // Simplified using atomic contributions
    let tpsa = 0;
    
    for (const atom of atoms) {
        switch (atom.symbol) {
            case 'O':
                tpsa += 20.23; // Single O
                break;
            case 'N':
                tpsa += 17.07; // Single N
                break;
            case 'S':
                tpsa += 38.8;
                break;
            case 'P':
                tpsa += 13.59;
                break;
        }
    }
    
    return tpsa;
}

export interface MolecularFeatures {
    molecular_weight: number;
    logp: number;
    num_h_donors: number;
    num_h_acceptors: number;
    num_rotatable_bonds: number;
    num_aromatic_rings: number;
    tpsa: number;
    num_atoms: number;
}

export function extractMolecularFeatures(smiles: string): MolecularFeatures | null {
    try {
        const parsed = parseSmilesBasic(smiles);
        if (!parsed) {
            return null;
        }
        
        const { atoms, bonds } = parsed;
        
        // If parsing is too simple, use a fallback calculation
        if (atoms.length === 0) {
            // Fallback: very basic estimation
            const basicAtoms = smiles.match(/[A-Z][a-z]?/g) || [];
            const atomSymbols: Atom[] = basicAtoms.map(symbol => ({
                symbol: symbol.length === 1 ? symbol : symbol,
                bonds: 0
            }));
            
            return {
                molecular_weight: calculateMolecularWeight(atomSymbols),
                logp: calculateLogP(atomSymbols, []),
                num_h_donors: countHBondDonors(atomSymbols),
                num_h_acceptors: countHBondAcceptors(atomSymbols),
                num_rotatable_bonds: Math.max(0, atomSymbols.length / 3),
                num_aromatic_rings: countAromaticRings(atomSymbols),
                tpsa: calculateTPSA(atomSymbols),
                num_atoms: atomSymbols.length
            };
        }
        
        return {
            molecular_weight: calculateMolecularWeight(atoms),
            logp: calculateLogP(atoms, bonds),
            num_h_donors: countHBondDonors(atoms),
            num_h_acceptors: countHBondAcceptors(atoms),
            num_rotatable_bonds: countRotatableBonds(bonds, atoms),
            num_aromatic_rings: countAromaticRings(atoms),
            tpsa: calculateTPSA(atoms),
            num_atoms: atoms.length
        };
    } catch (e) {
        console.error('Error extracting molecular features:', e);
        return null;
    }
}

