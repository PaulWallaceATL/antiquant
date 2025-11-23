import { NextResponse } from 'next/server';

export async function GET() {
    try {
        // Proxy JSME script from official website (most reliable)
        const jsmeUrl = 'https://peter-ertl.com/jsme/JSME_2024-04-29/jsme.nocache.js';
        
        const response = await fetch(jsmeUrl);
        
        if (!response.ok) {
            return NextResponse.json(
                { error: 'Failed to fetch JSME script' },
                { status: response.status }
            );
        }

        const scriptContent = await response.text();

        // Validate that we got JavaScript, not HTML/JSON error page
        const trimmed = scriptContent.trim();
        if (trimmed.startsWith('<') || trimmed.startsWith('{') || trimmed.startsWith('<!')) {
            console.error('[JSME Proxy] Received non-JavaScript content:', trimmed.substring(0, 200));
            return NextResponse.json(
                { error: 'JSME script source returned HTML/JSON instead of JavaScript' },
                { status: 502 }
            );
        }

        // Return as JavaScript with proper headers
        return new NextResponse(scriptContent, {
            headers: {
                'Content-Type': 'application/javascript; charset=utf-8',
                'Cache-Control': 'public, max-age=86400', // Cache for 1 day
            },
        });
    } catch (error: any) {
        console.error('JSME proxy error:', error);
        return NextResponse.json(
            { error: 'Failed to proxy JSME script' },
            { status: 500 }
        );
    }
}

