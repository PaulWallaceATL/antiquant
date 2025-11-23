import { NextResponse } from 'next/server';

// Ketcher info endpoint - required for Ketcher initialization
// Ketcher checks this endpoint to verify the backend is available
export async function GET() {
    return NextResponse.json({ 
        indigoVersion: '1.0.0',
        isAvailable: true,
        version: '1.0.0',
        apiVersion: '1.0',
        capabilities: {
            struct: true,
            convert: true,
            layout: true
        }
    });
}

