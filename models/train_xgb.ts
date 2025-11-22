import { execSync } from 'child_process';
import path from 'path';

const SCRIPT_PATH = path.join(__dirname, 'train_xgb.py');

console.log("Starting XGBoost training via Python script...");

try {
    const output = execSync(`python3 "${SCRIPT_PATH}"`, { encoding: 'utf-8' });
    console.log(output);
} catch (e) {
    console.error("Training failed:", e);
    process.exit(1);
}
