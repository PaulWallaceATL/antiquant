import { NextResponse } from 'next/server';

// Ketcher API endpoint - implements required struct service endpoints
// This provides the minimal API that Ketcher needs to function

export async function GET(request: Request) {
    const url = new URL(request.url);
    const pathname = url.pathname;
    const { searchParams } = url;
    const path = searchParams.get('path') || '';
    
    // Handle /api/ketcher/info endpoint (Ketcher checks this on initialization)
    if (pathname.includes('/info') || path.includes('info') || path.includes('version')) {
        return NextResponse.json({ 
            indigoVersion: '1.0.0',
            isAvailable: true,
            version: '1.0.0',
            apiVersion: '1.0'
        });
    }
    
    return NextResponse.json({ 
        message: 'Ketcher API endpoint',
        path,
        pathname 
    });
}

export async function POST(request: Request) {
    try {
        const url = new URL(request.url);
        const pathname = url.pathname;
        const body = await request.json();
        
        // Handle struct service operations that Ketcher requires
        // The key endpoint is for structure conversion/processing
        
        // Check the request path to determine the operation
        if (pathname.includes('struct') || pathname.includes('convert') || pathname.includes('layout') || 
            body.output_format || body.input_format) {
            
            // For structure conversion, we need to return a valid structure
            // Since we don't have Indigo, we'll return the input structure as-is
            // This allows basic editing to work, though advanced features won't
            const inputStruct = body.struct || body.structure || body;
            const outputFormat = body.output_format || 'mol';
            
            return NextResponse.json({ 
                struct: inputStruct,
                format: outputFormat
            });
        }
        
        // Handle info requests
        if (pathname.includes('info') || body.type === 'info') {
            return NextResponse.json({
                indigoVersion: '1.0.0',
                isAvailable: true
            });
        }
        
        // Default response for other operations
        return NextResponse.json({ 
            success: true,
            data: body 
        });
    } catch (error: any) {
        console.error('[Ketcher API] Error:', error);
        return NextResponse.json(
            { error: error.message || 'Internal server error' },
            { status: 500 }
        );
    }
}
