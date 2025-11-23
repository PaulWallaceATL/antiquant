// Enhanced molecular feature extraction using better SMILES parsing
// This is a more robust implementation for Vercel serverless

// Atomic weights for common elements
const ATOMIC_WEIGHTS: Record<string, number> = {
    'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012,
    'B': 10.811, 'C': 12.011, 'N': 14.007, 'O': 15.999,
    'F': 18.998, 'Ne': 20.180, 'Na': 22.990, 'Mg': 24.305,
    'Al': 26.982, 'Si': 28.085, 'P': 30.974, 'S': 32.065,
    'Cl': 35.453, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078,
    'Br': 79.904, 'I': 126.904, 'Se': 78.96, 'Te': 127.60
};

// Parse SMILES more robustly
export function parseSmiles(smiles: string): { atoms: string[], bonds: any[], aromatic: boolean[] } | null {
    const atoms: string[] = [];
    const bonds: any[] = [];
    const aromatic: boolean[] = [];
    
    let i = 0;
    let atomIndex = 0;
    const ringClosures = new Map<number, number>();
    
    while (i < smiles.length) {
        const char = smiles[i];
        const isUpper = char >= 'A' && char <= 'Z';
        const isLower = char >= 'a' && char <= 'z';
        
        if (isUpper || isLower) {
            // Check for two-letter element
            if (i + 1 < smiles.length && isLower) {
                const twoChar = char + smiles[i + 1];
                if (ATOMIC_WEIGHTS[twoChar] || twoChar === 'Br' || twoChar === 'Cl') {
                    atoms.push(twoChar);
                    aromatic.push(false);
                    atomIndex++;
                    i += 2;
                    continue;
                }
            }
            
            // Single letter element
            if (isUpper) {
                atoms.push(char);
                aromatic.push(false);
            } else {
                // Lowercase = aromatic
                atoms.push(char.toUpperCase());
                aromatic.push(true);
            }
            atomIndex++;
            i++;
        } else if (char === '=') {
            // Double bond
            if (atoms.length >= 2) {
                bonds.push({ begin: atoms.length - 2, end: atoms.length - 1, type: 2 });
            }
            i++;
        } else if (char === '#') {
            // Triple bond
            if (atoms.length >= 2) {
                bonds.push({ begin: atoms.length - 2, end: atoms.length - 1, type: 3 });
            }
            i++;
        } else if (char >= '0' && char <= '9') {
            // Ring closure
            const ringNum = parseInt(char);
            if (ringClosures.has(ringNum)) {
                const prevIdx = ringClosures.get(ringNum)!;
                bonds.push({ begin: prevIdx, end: atoms.length - 1, type: 1 });
                ringClosures.delete(ringNum);
            } else {
                ringClosures.set(ringNum, atoms.length - 1);
            }
            i++;
        } else if (char === '(' || char === ')' || char === '[' || char === ']') {
            // Branch or atom specification - skip for now
            i++;
        } else {
            i++;
        }
    }
    
    return { atoms, bonds, aromatic };
}

export function calculateMW(atoms: string[]): number {
    return atoms.reduce((sum, atom) => sum + (ATOMIC_WEIGHTS[atom] || 12.0), 0);
}

export function estimateLogP(atoms: string[], aromatic: boolean[]): number {
    // Fragment-based LogP estimation
    let logp = 0;
    const atomCounts: Record<string, number> = {};
    
    atoms.forEach((atom, i) => {
        atomCounts[atom] = (atomCounts[atom] || 0) + 1;
        const isAro = aromatic[i];
        
        switch (atom) {
            case 'C':
                logp += isAro ? 1.68 : 0.5; // Aromatic C more lipophilic
                break;
            case 'N':
                logp -= 1.4;
                break;
            case 'O':
                logp -= 0.7;
                break;
            case 'F':
                logp += 0.13;
                break;
            case 'Cl':
                logp += 0.71;
                break;
            case 'Br':
                logp += 0.86;
                break;
            case 'I':
                logp += 1.12;
                break;
            case 'S':
                logp += 0.05;
                break;
        }
    });
    
    // Aromatic rings add lipophilicity
    const aromaticRings = Math.floor(aromatic.filter(a => a).length / 6);
    logp += aromaticRings * 0.4;
    
    return logp;
}

export function countHDonors(atoms: string[]): number {
    // O, N, S with potential H
    return atoms.filter(a => ['O', 'N', 'S'].includes(a)).length;
}

export function countHAcceptors(atoms: string[]): number {
    // O, N, F, Cl, S
    return atoms.filter(a => ['O', 'N', 'F', 'Cl', 'S'].includes(a)).length;
}

export function estimateRotatableBonds(atoms: string[], bonds: any[]): number {
    // Simplified: estimate based on number of atoms and bonds
    // Single bonds that aren't in rings
    const singleBonds = bonds.filter(b => b.type === 1).length;
    // Rough estimate: subtract ring bonds and terminal bonds
    return Math.max(0, Math.floor(singleBonds - atoms.length * 0.3));
}

export function countAromaticRings(aromatic: boolean[]): number {
    const aromaticCount = aromatic.filter(a => a).length;
    // Estimate rings: ~6 atoms per aromatic ring
    return Math.floor(aromaticCount / 6);
}

export function estimateTPSA(atoms: string[]): number {
    // Topological Polar Surface Area using atomic contributions
    let tpsa = 0;
    const atomCounts: Record<string, number> = {};
    
    atoms.forEach(atom => {
        atomCounts[atom] = (atomCounts[atom] || 0) + 1;
    });
    
    // TPSA contributions (approximate)
    tpsa += (atomCounts['O'] || 0) * 20.23;
    tpsa += (atomCounts['N'] || 0) * 17.07;
    tpsa += (atomCounts['S'] || 0) * 38.8;
    tpsa += (atomCounts['P'] || 0) * 13.59;
    
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
        const parsed = parseSmiles(smiles);
        if (!parsed || parsed.atoms.length === 0) {
            // Fallback: very basic estimation
            const basicAtoms = smiles.match(/[A-Z][a-z]?/g) || smiles.match(/[a-z]/g)?.map(s => s.toUpperCase()) || [];
            
            return {
                molecular_weight: basicAtoms.reduce((sum, a) => sum + (ATOMIC_WEIGHTS[a] || 12.0), 0),
                logp: basicAtoms.length * 0.3 - (basicAtoms.filter(a => ['O', 'N'].includes(a)).length * 0.5),
                num_h_donors: basicAtoms.filter(a => ['O', 'N', 'S'].includes(a)).length,
                num_h_acceptors: basicAtoms.filter(a => ['O', 'N', 'F', 'Cl', 'S'].includes(a)).length,
                num_rotatable_bonds: Math.max(0, Math.floor(basicAtoms.length / 3)),
                num_aromatic_rings: Math.floor((smiles.match(/[a-z]/g) || []).length / 6),
                tpsa: (basicAtoms.filter(a => a === 'O').length * 20.23) + (basicAtoms.filter(a => a === 'N').length * 17.07),
                num_atoms: basicAtoms.length
            };
        }
        
        const { atoms, bonds, aromatic } = parsed;
        
        return {
            molecular_weight: calculateMW(atoms),
            logp: estimateLogP(atoms, aromatic),
            num_h_donors: countHDonors(atoms),
            num_h_acceptors: countHAcceptors(atoms),
            num_rotatable_bonds: estimateRotatableBonds(atoms, bonds),
            num_aromatic_rings: countAromaticRings(aromatic),
            tpsa: estimateTPSA(atoms),
            num_atoms: atoms.length
        };
    } catch (e) {
        console.error('Error extracting molecular features:', e);
        // Return a very basic fallback
        return {
            molecular_weight: 100,
            logp: 1.0,
            num_h_donors: 1,
            num_h_acceptors: 1,
            num_rotatable_bonds: 2,
            num_aromatic_rings: 0,
            tpsa: 40,
            num_atoms: 10
        };
    }
}

