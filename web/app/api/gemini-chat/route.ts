import { NextResponse } from 'next/server';

export async function POST(request: Request) {
    try {
        const body = await request.json();
        const { message, context, temperature = 0.7 } = body;

        if (!message) {
            return NextResponse.json({ error: 'Message is required' }, { status: 400 });
        }

        const apiKey = process.env.GEMINI_API_KEY;
        if (!apiKey) {
            return NextResponse.json({
                response: "I'm a molecular analysis assistant. I can help you generate SMILES strings, find similar molecules, and analyze molecular properties. However, Gemini API is not configured. Please set GEMINI_API_KEY in .env.local for full functionality."
            });
        }

        // Use Gemini 3 Pro (most advanced) - per https://ai.google.dev/gemini-api/docs/gemini-3
        // Model ID: gemini-3-pro-preview
        const model = 'gemini-3-pro-preview';

        // Enhanced system prompt for molecular analysis
        const systemPrompt = `You are an expert AI assistant specialized in molecular analysis, drug discovery, and pharmaceutical research. Your expertise includes:

- SMILES (Simplified Molecular Input Line Entry System) notation for chemical structures
- Molecular property prediction (solubility, drug-likeness, ADMET properties)
- Drug discovery and lead optimization
- Understanding molecular descriptors (molecular weight, LogP, H-bond donors/acceptors, TPSA, etc.)
- Lipinski's Rule of Five and other drug-likeness rules
- Finding similar molecules using chemical similarity
- Generating novel molecular structures based on therapeutic requirements

When generating SMILES strings:
- Always return valid, syntactically correct SMILES notation
- Consider drug-likeness principles (MW < 500, LogP < 5, etc.)
- Suggest molecules relevant to the therapeutic area mentioned
- Provide 3-5 diverse candidate structures when possible

Be concise, technically accurate, and helpful. When users ask for drug design, provide realistic molecular structures that could be synthesized and tested.`;

        // Build messages array for Gemini
        const messages: Array<{ role: 'user' | 'model'; parts: Array<{ text: string }> }> = [];

        // Add system context as first user message
        messages.push({
            role: 'user',
            parts: [{ text: systemPrompt }]
        });
        messages.push({
            role: 'model',
            parts: [{ text: 'I understand. I\'m ready to help with molecular analysis, SMILES generation, drug discovery, and pharmaceutical research. I can generate valid SMILES strings, explain molecular properties, and assist with drug design questions.' }]
        });

        // Add conversation context if provided
        if (context && Array.isArray(context)) {
            context.forEach((msg: { role: string; content: string }) => {
                if (msg.role === 'user') {
                    messages.push({
                        role: 'user',
                        parts: [{ text: msg.content }]
                    });
                } else if (msg.role === 'assistant') {
                    messages.push({
                        role: 'model',
                        parts: [{ text: msg.content }]
                    });
                }
            });
        }

        // Add current message
        messages.push({
            role: 'user',
            parts: [{ text: message }]
        });

        // Call Gemini API using correct format per https://ai.google.dev/gemini-api/docs/gemini-3
        const response = await fetch(
            `https://generativelanguage.googleapis.com/v1beta/models/${model}:generateContent`,
            {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'x-goog-api-key': apiKey,
                },
                body: JSON.stringify({
                    contents: messages,
                    generationConfig: {
                        temperature,
                        maxOutputTokens: 2048,
                    },
                }),
            }
        );

        if (!response.ok) {
            const errorData = await response.text();
            console.error('Gemini API Error:', errorData);
            let errorMessage = `Gemini API error: ${response.status}`;
            try {
                const errorJson = JSON.parse(errorData);
                errorMessage = errorJson.error?.message || errorMessage;
            } catch (e) {
                // Not JSON, use text as is
            }
            return NextResponse.json({ error: errorMessage }, { status: response.status });
        }

        const data = await response.json();
        const text = data.candidates?.[0]?.content?.parts?.[0]?.text || '';
        
        if (!text) {
            console.error('No text in Gemini response:', JSON.stringify(data, null, 2));
            return NextResponse.json({ 
                error: 'No response from Gemini API. Check API key and quota.' 
            }, { status: 500 });
        }

        return NextResponse.json({ response: text });
    } catch (error: any) {
        console.error('Gemini API Error:', error);
        return NextResponse.json(
            { error: error.message || 'Failed to process chat message' },
            { status: 500 }
        );
    }
}

