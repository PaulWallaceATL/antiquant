import { NextResponse } from 'next/server';
import { predictBoth } from '@/lib/molecular_inference';

export async function POST(request: Request) {
    try {
        const body = await request.json();
        const { smiles } = body;

        if (!smiles) {
            return NextResponse.json({ error: 'SMILES string is required' }, { status: 400 });
        }

        const result = await predictBoth(smiles);

        if (result.error) {
            return NextResponse.json({ error: result.error }, { status: 500 });
        }

        return NextResponse.json(result);
    } catch (error) {
        console.error('API Error:', error);
        return NextResponse.json({ error: 'Internal Server Error' }, { status: 500 });
    }
}
