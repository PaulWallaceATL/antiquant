import { NextResponse } from 'next/server';
import { processChatMessage } from '@/lib/chat_service';

export async function POST(request: Request) {
    try {
        const body = await request.json();
        const { message, sessionId } = body;

        console.log('[API /chat] Received request:', { message: message?.substring(0, 100), sessionId });

        if (!message) {
            console.error('[API /chat] No message provided');
            return NextResponse.json({ error: 'Message is required' }, { status: 400 });
        }

        console.log('[API /chat] Calling processChatMessage');
        const response = await processChatMessage(message, sessionId);
        console.log('[API /chat] Got response:', { responseLength: response.response?.length, hasError: !!response.error });

        return NextResponse.json(response);
    } catch (error: any) {
        console.error('[API /chat] Error:', error);
        return NextResponse.json(
            { error: error.message || 'Failed to process chat message' },
            { status: 500 }
        );
    }
}

