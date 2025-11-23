import { NextResponse } from 'next/server';

export async function POST(request: Request) {
    try {
        const body = await request.json();
        const { message, context, temperature = 0.7 } = body;

        if (!message) {
            return NextResponse.json({ error: 'Message is required' }, { status: 400 });
        }

        const apiKey = process.env.OPENAI_API_KEY;
        if (!apiKey) {
            return NextResponse.json({
                response: "I'm a molecular analysis assistant. I can help you generate SMILES strings, find similar molecules, and analyze molecular properties. However, OpenAI API is not configured. Please set OPENAI_API_KEY in .env.local for full functionality."
            });
        }

        // Build conversation history for GPT-5.1 Responses API
        // The Responses API uses a different format - we need to build the full conversation
        let conversationInput = message;
        
        // If we have context, prepend it to the input
        if (context && Array.isArray(context) && context.length > 0) {
            const contextText = context
                .map((msg: { role: string; content: string }) => {
                    if (msg.role === 'user') {
                        return `User: ${msg.content}`;
                    } else if (msg.role === 'assistant') {
                        return `Assistant: ${msg.content}`;
                    }
                    return '';
                })
                .filter(Boolean)
                .join('\n\n');
            conversationInput = `${contextText}\n\nUser: ${message}`;
        }

        // Enhanced system prompt for molecular analysis
        const systemInstruction = `You are an expert AI assistant specialized in molecular analysis, drug discovery, and pharmaceutical research. Your expertise includes:

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

        // Call OpenAI Responses API for GPT-5.1
        const response = await fetch(
            'https://api.openai.com/v1/responses',
            {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'Authorization': `Bearer ${apiKey}`,
                },
                body: JSON.stringify({
                    model: 'gpt-5.1',
                    input: `${systemInstruction}\n\n${conversationInput}`,
                    reasoning: {
                        effort: 'medium', // Use medium reasoning for complex molecular analysis
                    },
                    text: {
                        verbosity: 'medium',
                    },
                }),
            }
        );

        if (!response.ok) {
            const errorData = await response.text();
            console.error('OpenAI API Error:', errorData);
            let errorMessage = `OpenAI API error: ${response.status}`;
            try {
                const errorJson = JSON.parse(errorData);
                errorMessage = errorJson.error?.message || errorMessage;
            } catch (e) {
                // Not JSON, use text as is
            }
            return NextResponse.json({ error: errorMessage }, { status: response.status });
        }

        const data = await response.json();
        const text = data.output_text || data.output?.text || '';
        
        if (!text) {
            console.error('No text in OpenAI response:', JSON.stringify(data, null, 2));
            return NextResponse.json({ 
                error: 'No response from OpenAI API. Check API key and quota.' 
            }, { status: 500 });
        }

        return NextResponse.json({ response: text });
    } catch (error: any) {
        console.error('OpenAI API Error:', error);
        return NextResponse.json(
            { error: error.message || 'Failed to process chat message' },
            { status: 500 }
        );
    }
}
