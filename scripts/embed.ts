import fs from 'fs';
import path from 'path';
import * as csv from 'fast-csv';
import OpenAI from 'openai';

const INPUT_CSV = path.join(__dirname, '../data/codenet_sample.csv');
const OUTPUT_EMBEDDINGS = path.join(__dirname, '../data/embeddings.json'); // Using JSON for simplicity in JS
const OUTPUT_INDEX = path.join(__dirname, '../data/embeddings_index.csv');

const openai = new OpenAI({
    apiKey: process.env.OPENAI_API_KEY,
});

interface CodeNetRow {
    id: string;
    code: string;
}

async function generateEmbeddings() {
    if (!process.env.OPENAI_API_KEY) {
        console.error("Error: OPENAI_API_KEY is not set.");
        process.exit(1);
    }

    const rows: CodeNetRow[] = [];

    await new Promise<void>((resolve, reject) => {
        fs.createReadStream(INPUT_CSV)
            .pipe(csv.parse({ headers: true }))
            .on('error', (error: any) => reject(error))
            .on('data', (row: any) => rows.push(row))
            .on('end', () => resolve());
    });

    console.log(`Loaded ${rows.length} rows for embedding.`);

    const embeddings: number[][] = [];
    const indices: { id: string, index: number }[] = [];

    // Process in batches to avoid rate limits
    const BATCH_SIZE = 20;

    for (let i = 0; i < rows.length; i += BATCH_SIZE) {
        const batch = rows.slice(i, i + BATCH_SIZE);
        const inputs = batch.map(r => r.code.substring(0, 8000)); // Truncate to avoid token limits

        try {
            const response = await openai.embeddings.create({
                model: "text-embedding-3-large",
                input: inputs,
            });

            response.data.forEach((item, idx) => {
                embeddings.push(item.embedding);
                indices.push({ id: batch[idx].id, index: i + idx });
            });

            console.log(`Processed batch ${i / BATCH_SIZE + 1}/${Math.ceil(rows.length / BATCH_SIZE)}`);

        } catch (e) {
            console.error(`Error processing batch starting at ${i}:`, e);
        }
    }

    // Save embeddings as JSON (easier for JS/TS to read than .npy without libraries)
    // If .npy is strictly required, we'd need a buffer writer.
    // For now, JSON is safer for "production-ready" JS without native deps.
    fs.writeFileSync(OUTPUT_EMBEDDINGS, JSON.stringify(embeddings));

    // Save index
    const csvStream = csv.format({ headers: true });
    const writeStream = fs.createWriteStream(OUTPUT_INDEX);
    csvStream.pipe(writeStream);
    indices.forEach(idx => csvStream.write(idx));
    csvStream.end();

    console.log(`Saved embeddings to ${OUTPUT_EMBEDDINGS}`);
    console.log(`Saved index to ${OUTPUT_INDEX}`);
}

generateEmbeddings().catch(console.error);
