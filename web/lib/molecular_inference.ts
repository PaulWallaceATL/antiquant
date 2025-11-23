import { spawn } from 'child_process';
import path from 'path';

const INFERENCE_SCRIPT = path.join(process.cwd(), '../scripts/molecular_inference.py');

export interface QuantumCircuit {
    n_qubits: number;
    n_layers: number;
    operations: Array<{
        gate: string;
        qubit?: number;
        parameter?: number;
        parameters?: number[];
        control?: number;
        target?: number;
        layer: string;
    }>;
    expectation_value: number;
    final_state_probabilities: number[];
    input_encoding: number[];
}

export interface Structure3D {
    sdf: string;
    atoms: Array<{
        element: string;
        x: number;
        y: number;
        z: number;
        idx: number;
    }>;
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

export interface PredictionResult {
    prediction: string;
    property: string;
    model_used: string;
    confidence: number;
    quantum_circuit?: QuantumCircuit;
}

export interface ADMET {
    absorption: {
        caco2_permeability: string;
        hia: string;
    };
    distribution: {
        bbb_permeability: string;
        plasma_protein_binding: string;
    };
    metabolism: {
        cyp450_substrate: string;
        cyp450_inhibitor: string;
    };
    excretion: {
        renal_clearance: string;
    };
    toxicity: {
        ames_mutagenicity: string;
        hepatotoxicity: string;
        skin_sensitization: string;
    };
}

export interface DrugLikeness {
    lipinski: {
        pass: boolean;
        violations: number;
        details: {
            molecular_weight: { value: number; limit: number; pass: boolean };
            logp: { value: number; limit: number; pass: boolean };
            h_donors: { value: number; limit: number; pass: boolean };
            h_acceptors: { value: number; limit: number; pass: boolean };
        };
    };
    veber: {
        pass: boolean;
        rotatable_bonds: number;
        tpsa: number;
    };
    bbb_permeability: string;
    synthetic_accessibility: number;
    drug_score: number;
    admet: ADMET;
}

export interface CompetitorDrug {
    name: string;
    smiles: string;
    similarity: number;
    indication: string;
    status: string;
    company: string;
}

export interface PatentInfo {
    patentNumber: string;
    title: string;
    filingDate: string;
    status: string;
    assignee: string;
    relevance: string;
}

export interface MarketSize {
    therapeuticArea: string;
    estimatedMarketSize: string;
    currency: string;
    year: number;
    growthRate: string;
}

export interface MarketAnalysis {
    competitors: CompetitorDrug[];
    patents: PatentInfo[];
    marketSize: MarketSize;
    regulatoryStatus: string;
}

export interface MolecularAnalysisResult {
    features: MolecularFeatures;
    smiles: string;
    structure_3d?: Structure3D;
    drug_likeness?: DrugLikeness;
    marketAnalysis?: MarketAnalysis;
    prediction?: string;
    property?: string;
    model_used?: string;
    confidence?: number;
    quantum_circuit?: QuantumCircuit;
    classical?: PredictionResult;
    quantum?: PredictionResult;
    error?: string;
}

async function runInference(smiles: string, mode: 'classical' | 'quantum' | 'both'): Promise<MolecularAnalysisResult> {
    return new Promise((resolve, reject) => {
        const pythonProcess = spawn('python3', [INFERENCE_SCRIPT]);

        let outputData = '';
        let errorData = '';

        pythonProcess.stdout.on('data', (data) => {
            outputData += data.toString();
        });

        pythonProcess.stderr.on('data', (data) => {
            errorData += data.toString();
        });

        pythonProcess.on('close', (code) => {
            if (code !== 0) {
                console.error(`Python script exited with code ${code}`);
                console.error(`Stderr: ${errorData}`);
                reject(new Error(`Inference failed: ${errorData}`));
                return;
            }

            try {
                const result = JSON.parse(outputData);
                resolve(result);
            } catch (e) {
                reject(new Error(`Failed to parse output: ${outputData}`));
            }
        });

        // Send input
        const input = JSON.stringify({ smiles, mode });
        pythonProcess.stdin.write(input);
        pythonProcess.stdin.end();
    });
}

export async function predictClassical(smiles: string): Promise<MolecularAnalysisResult> {
    return runInference(smiles, 'classical');
}

export async function predictQuantum(smiles: string): Promise<MolecularAnalysisResult> {
    return runInference(smiles, 'quantum');
}

export async function predictBoth(smiles: string): Promise<MolecularAnalysisResult> {
    return runInference(smiles, 'both');
}
