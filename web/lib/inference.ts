import { spawn } from 'child_process';
import path from 'path';

const INFERENCE_SCRIPT = path.join(process.cwd(), '../scripts/inference.py');

export interface AnalysisResult {
    prediction: string;
    confidence: number;
    model_used: string;
    features: {
        lines_of_code: number;
        num_functions: number;
        cyclomatic_complexity: number;
        max_nesting_depth: number;
    };
    error?: string;
}

async function runInference(code: string, mode: 'classical' | 'quantum'): Promise<AnalysisResult> {
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
        const input = JSON.stringify({ code, mode });
        pythonProcess.stdin.write(input);
        pythonProcess.stdin.end();
    });
}

export async function predictClassical(code: string): Promise<AnalysisResult> {
    return runInference(code, 'classical');
}

export async function predictQuantum(code: string): Promise<AnalysisResult> {
    return runInference(code, 'quantum');
}
