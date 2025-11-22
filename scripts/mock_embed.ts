import fs from 'fs';
import path from 'path';
import * as csv from 'fast-csv';

const INPUT_CSV = path.join(__dirname, '../data/molecules_sample.csv');
const OUTPUT_EMBEDDINGS = path.join(__dirname, '../data/embeddings.json');
const OUTPUT_INDEX = path.join(__dirname, '../data/embeddings_index.csv');

interface MoleculeRow {
    id: string;
    smiles: string;
}

async function generateMockEmbeddings() {
    const rows: MoleculeRow[] = [];

    await new Promise<void>((resolve, reject) => {
        fs.createReadStream(INPUT_CSV)
            .pipe(csv.parse({ headers: true }))
            .on('error', (error: any) => reject(error))
            .on('data', (row: any) => rows.push(row))
            .on('end', () => resolve());
    });

    console.log(`Loaded ${rows.length} rows for mock embedding.`);

    const embeddings: number[][] = [];
    const indices: { id: string, index: number }[] = [];

    // Generate random 1536-dim vectors (standard OpenAI size)
    // For VQC we will reduce this later.
    const DIM = 1536;

    for (let i = 0; i < rows.length; i++) {
        const vector = Array.from({ length: DIM }, () => Math.random() * 2 - 1);
        embeddings.push(vector);
        indices.push({ id: rows[i].id, index: i });
    }

    fs.writeFileSync(OUTPUT_EMBEDDINGS, JSON.stringify(embeddings));

    const csvStream = csv.format({ headers: true });
    const writeStream = fs.createWriteStream(OUTPUT_INDEX);
    csvStream.pipe(writeStream);
    indices.forEach(idx => csvStream.write(idx));
    csvStream.end();

    console.log(`Saved MOCK embeddings to ${OUTPUT_EMBEDDINGS}`);
    console.log(`Saved index to ${OUTPUT_INDEX}`);
}

generateMockEmbeddings().catch(console.error);
