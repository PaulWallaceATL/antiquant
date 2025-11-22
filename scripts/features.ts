import fs from 'fs';
console.log("Starting script...");
import path from 'path';
import { execSync } from 'child_process';
import * as csv from 'fast-csv';

const INPUT_CSV = path.join(__dirname, '../data/codenet_sample.csv');
const OUTPUT_CSV = path.join(__dirname, '../data/features.csv');

interface CodeNetRow {
    id: string;
    language: string;
    problem_id: string;
    code: string;
}

interface Features {
    id: string;
    lines_of_code: number;
    num_functions: number;
    cyclomatic_complexity: number;
    max_nesting_depth: number;
}

function calculateNestingDepth(code: string): number {
    let maxDepth = 0;
    let currentDepth = 0;
    const lines = code.split('\n');

    // Simple indentation-based depth for Python
    // This is an approximation.
    for (const line of lines) {
        const trimmed = line.trim();
        if (!trimmed || trimmed.startsWith('#')) continue;

        // Count leading spaces
        const leadingSpaces = line.match(/^ */)?.[0].length || 0;
        // Assuming 4 spaces per indent, or just tracking changes
        // A better heuristic for Python:
        // Depth is roughly leading_spaces / 4. 
        // But mixed tabs/spaces can be tricky. 
        // Let's assume standard 4-space indent for CodeNet or just count indentation levels.

        // Actually, let's just use the raw indentation level relative to the base.
        // But we need to handle the fact that the code might be wrapped.

        const depth = Math.floor(leadingSpaces / 4);
        if (depth > maxDepth) maxDepth = depth;
    }
    return maxDepth;
}

async function processRows() {
    const rows: CodeNetRow[] = [];

    // Read CSV
    await new Promise<void>((resolve, reject) => {
        fs.createReadStream(INPUT_CSV)
            .pipe(csv.parse({ headers: true }))
            .on('error', (error: any) => reject(error))
            .on('data', (row: any) => rows.push(row))
            .on('end', () => resolve());
    });

    console.log(`Loaded ${rows.length} rows.`);

    const featuresList: Features[] = [];
    let processed = 0;

    for (const row of rows) {
        const { id, code } = row;

        // Write code to temp file
        const tempFile = path.join(__dirname, `temp_${id}.py`);
        fs.writeFileSync(tempFile, code);

        let loc = 0;
        let numFunctions = 0;
        let ccn = 0;

        try {
            // Run lizard
            // lizard --csv returns: NLOC, CCN, token_count, param_count, length, location, file, function_name, ...
            // But if we run on a file, it gives a summary or per-function list.
            // We want aggregate metrics for the snippet.

            // lizard output format (default):
            // NLOC    CCN   token  PARAM  length  location  
            //       3      1     14      0       3 temp_s1.py@1-3@read_root
            // 1 file analyzed.
            // ==============================================================
            // NLOC    Avg.NLOC  AvgCCN  Avg.token  function_cnt    file
            //       3       3.0     1.0       14.0         1     temp_s1.py

            // We can parse the summary line at the bottom.
            const output = execSync(`python3 -m lizard "${tempFile}"`, { encoding: 'utf-8' });
            const lines = output.trim().split('\n');
            const summaryLine = lines[lines.length - 1];

            // Summary line format: NLOC Avg.NLOC AvgCCN Avg.token function_cnt file
            // We need to be careful parsing this.
            // Let's use regex to find the line ending with the filename.

            // Example: "      3       3.0     1.0       14.0         1     temp_s1.py"
            const parts = summaryLine.trim().split(/\s+/);

            // parts: [NLOC, Avg.NLOC, AvgCCN, Avg.token, function_cnt, file]
            // If function_cnt is 0, the format might be different?
            // If no functions, lizard might just show file stats.

            if (parts.length >= 6) {
                loc = parseInt(parts[0], 10);
                // AvgCCN is average per function. If we want total complexity, maybe sum?
                // Or just use AvgCCN as the metric? The prompt says "Cyclomatic complexity".
                // Usually for a file, it's either sum of functions or max.
                // Let's use AvgCCN * function_cnt or just AvgCCN.
                // Actually, let's parse the per-function lines to get max or sum.
                // But for simplicity, let's use the summary AvgCCN.
                // Wait, if there are 0 functions, AvgCCN might be 0.

                const avgCcn = parseFloat(parts[2]);
                numFunctions = parseInt(parts[4], 10);

                // If numFunctions > 0, Total CCN approx = AvgCCN * function_cnt?
                // Or just use AvgCCN. Let's use AvgCCN for now as "Complexity".
                // But "Cyclomatic complexity" usually refers to the code's complexity.
                // Let's try to get the max CCN of any function.

                // Parse function lines
                // Function lines start with a number.
                // 3      1     14      0       3 temp_s1.py@1-3@read_root

                let maxCcn = 0;
                let totalCcn = 0;

                for (const line of lines) {
                    const match = line.match(/^\s*(\d+)\s+(\d+)\s+\d+\s+\d+\s+\d+\s+.*@.*@/);
                    if (match) {
                        const funcCcn = parseInt(match[2], 10);
                        if (funcCcn > maxCcn) maxCcn = funcCcn;
                        totalCcn += funcCcn;
                    }
                }

                ccn = maxCcn > 0 ? maxCcn : 1; // Default to 1 if no functions or flat code
            }

        } catch (e) {
            console.error(`Error processing ${id}:`, e);
        } finally {
            if (fs.existsSync(tempFile)) fs.unlinkSync(tempFile);
        }

        const nesting = calculateNestingDepth(code);

        featuresList.push({
            id,
            lines_of_code: loc,
            num_functions: numFunctions,
            cyclomatic_complexity: ccn,
            max_nesting_depth: nesting
        });

        processed++;
        if (processed % 100 === 0) console.log(`Processed ${processed} rows...`);
    }

    // Write CSV
    const csvStream = csv.format({ headers: true });
    const writeStream = fs.createWriteStream(OUTPUT_CSV);

    csvStream.pipe(writeStream).on('finish', () => console.log('Done writing features.'));

    featuresList.forEach(f => csvStream.write(f));
    csvStream.end();

    console.log(`Saved features to ${OUTPUT_CSV}`);
}

processRows().catch(console.error);
