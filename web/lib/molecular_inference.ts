import { spawn } from 'child_process';
import path from 'path';

const INFERENCE_SCRIPT = path.join(process.cwd(), 'scripts/molecular_inference.py');

export interface MolecularAnalysisResult {
    prediction: string;
    property: string;
    model_used: string;
    features: {
        molecular_weight: number;
        logp: number;
        num_h_donors: number;
        num_h_acceptors: number;
        num_rotatable_bonds: number;
        num_aromatic_rings: number;
        tpsa: number;
        num_atoms: number;
    };
    smiles: string;
    error?: string;
}

async function runInference(smiles: string, mode: 'classical' | 'quantum'): Promise<MolecularAnalysisResult> {
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
