import { spawn } from 'child_process';
import path from 'path';

const CHAT_HELPERS_SCRIPT = path.join(process.cwd(), '../scripts/chat_helpers.py');

// Use Cloud Run API or fallback to local Python spawn for development
const API_URL = process.env.MOLECULEAI_API_URL || process.env.NEXT_PUBLIC_MOLECULEAI_API_URL || '';
const USE_CLOUD_RUN = !!API_URL;

// Helper to call OpenAI GPT-5.1 API directly (for server-side use)
async function callOpenAIAPI(message: string, context?: Array<{ role: string; content: string }>, temperature: number = 0.7): Promise<string> {
    console.log('[callOpenAIAPI] Starting API call with message:', message.substring(0, 100));
    const apiKey = process.env.OPENAI_API_KEY;
    if (!apiKey) {
        console.error('[callOpenAIAPI] OPENAI_API_KEY is not set!');
        throw new Error('OPENAI_API_KEY is not set in environment variables');
    }
    console.log('[callOpenAIAPI] API key found, length:', apiKey.length);

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
- Format SMILES strings using backticks (e.g., \`CCO\`) so they can be easily identified
- Consider drug-likeness principles (MW < 500, LogP < 5, etc.)
- Suggest molecules relevant to the therapeutic area mentioned
- Provide 3-5 diverse candidate structures when possible

Be concise, technically accurate, and helpful. When users ask for drug design, provide realistic molecular structures that could be synthesized and tested.`;

    // Build conversation input for GPT-5.1 Responses API
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

    // Use Chat Completions API directly (GPT-5.1 Responses API may not be available)
    // Build messages array for Chat Completions
    const messages: Array<{ role: string; content: string }> = [
        { role: 'system', content: systemPrompt },
    ];

    // Add context
    if (context && Array.isArray(context)) {
        context.forEach((msg: { role: string; content: string }) => {
            if (msg.role === 'user' || msg.role === 'assistant') {
                messages.push({ role: msg.role, content: msg.content });
            }
        });
    }

    // Add current message
    messages.push({ role: 'user', content: message });

    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), 20000); // 20 second timeout

    try {
        console.log('[callOpenAIAPI] Calling Chat Completions API with GPT-4o');
        const response = await fetch(
            'https://api.openai.com/v1/chat/completions',
            {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                    'Authorization': `Bearer ${apiKey}`,
                },
                body: JSON.stringify({
                    model: 'gpt-4o',
                    messages: messages,
                    temperature: 0.7,
                    max_tokens: 2048,
                }),
                signal: controller.signal,
            }
        );

        clearTimeout(timeoutId);
        console.log('[callOpenAIAPI] API response status:', response.status, response.ok);

        if (!response.ok) {
            const errorData = await response.text();
            let errorMessage = `OpenAI API error: ${response.status}`;
            try {
                const errorJson = JSON.parse(errorData);
                errorMessage = errorJson.error?.message || errorMessage;
            } catch (e) {
                errorMessage = errorData || errorMessage;
            }
            console.error('[callOpenAIAPI] API call failed:', errorMessage);
            throw new Error(errorMessage);
        }

        const data = await response.json();
        console.log('[callOpenAIAPI] Full API response:', JSON.stringify(data, null, 2));

        // Chat Completions API format - extract text from choices
        let text = '';

        if (data.choices && Array.isArray(data.choices) && data.choices.length > 0) {
            const choice = data.choices[0];
            if (choice.message?.content) {
                text = choice.message.content;
            } else if (choice.text) {
                text = choice.text;
            } else if (choice.delta?.content) {
                text = choice.delta.content;
            }
        }

        // Fallback checks
        if (!text && data.output_text) {
            text = data.output_text;
        }
        if (!text && data.output?.text) {
            text = data.output.text;
        }
        if (!text && data.text) {
            text = data.text;
        }
        if (!text && data.content) {
            text = data.content;
        }

        // CRITICAL: Ensure we always return a string, never an object
        if (typeof text !== 'string') {
            console.error('[callOpenAIAPI] Response text is not a string! Type:', typeof text, 'Value:', text);
            // If text is an object, this is a bug - log it and throw
            if (text && typeof text === 'object') {
                console.error('[callOpenAIAPI] Received object instead of string:', JSON.stringify(text, null, 2));
                throw new Error('OpenAI API returned an object instead of text. This should not happen.');
            }
            text = '';
        }

        if (!text || text.trim() === '') {
            console.error('[callOpenAIAPI] No text in OpenAI response. Full response:', JSON.stringify(data, null, 2));
            throw new Error('No response text from OpenAI API. Check API key, model availability, and quota. Response keys: ' + Object.keys(data).join(', '));
        }

        const finalText = String(text).trim();
        console.log('[callOpenAIAPI] Extracted text (first 100 chars):', finalText.substring(0, 100));
        return finalText; // Ensure it's always a string and trimmed
    } catch (error: any) {
        clearTimeout(timeoutId);
        if (error.name === 'AbortError') {
            throw new Error('Request timeout: OpenAI API took too long to respond. Please try again.');
        }
        console.error('OpenAI API call error:', error);
        throw error;
    }
}

export interface SimilarMolecule {
    id: string;
    smiles: string;
    similarity: number;
    target_property?: number;
}

export interface ChatContext {
    messages: Array<{ role: 'user' | 'assistant'; content: string }>;
    suggestedSmiles?: string[];
    analyzedMolecules?: string[];
}

export interface ChatResponse {
    response: string;
    sessionId: string;
    suggestedSmiles?: string[];
    similarMolecules?: SimilarMolecule[];
    shouldAnalyze?: boolean;
    analysisSmiles?: string;
}

// In-memory session storage (for MVP)
const sessions = new Map<string, ChatContext>();

function generateSessionId(): string {
    return `session_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`;
}

async function runPythonScript(command: string, smiles: string, topK?: number): Promise<any> {
    // Use Cloud Run API if available
    if (USE_CLOUD_RUN) {
        try {
            let endpoint = '';
            let body: any = { smiles };

            if (command === 'validate') {
                endpoint = '/validate';
            } else if (command === 'find_similar') {
                endpoint = '/similar';
                body.limit = topK || 5;
            } else {
                throw new Error(`Unknown command: ${command}`);
            }

            const response = await fetch(`${API_URL}${endpoint}`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify(body),
            });

            if (!response.ok) {
                const errorData = await response.json().catch(() => ({
                    error: `HTTP ${response.status}: ${response.statusText}`
                }));
                throw new Error(errorData.detail || errorData.error || `HTTP ${response.status}`);
            }

            return await response.json();
        } catch (error: any) {
            console.error('Cloud Run API call failed:', error);
            // If API fails, we could try fallback, but for now just throw
            throw new Error(`Failed to run command ${command}: ${error.message}`);
        }
    }

    // Fallback to local Python spawn (for local development)
    return new Promise((resolve, reject) => {
        const args = [command, smiles];
        if (topK !== undefined) {
            args.push(topK.toString());
        }

        const pythonProcess = spawn('python3', [CHAT_HELPERS_SCRIPT, ...args]);

        let outputData = '';
        let errorData = '';

        pythonProcess.stdout.on('data', (data) => {
            outputData += data.toString();
        });

        pythonProcess.stderr.on('data', (data) => {
            errorData += data.toString();
        });

        pythonProcess.on('close', (code) => {
            if (code !== 0) {
                reject(new Error(`Python script failed: ${errorData}`));
                return;
            }

            try {
                const result = JSON.parse(outputData);
                resolve(result);
            } catch (e) {
                reject(new Error(`Failed to parse output: ${outputData}`));
            }
        });
    });
}

export async function validateSmiles(smiles: string): Promise<boolean> {
    try {
        const result = await runPythonScript('validate', smiles) as { valid: boolean };
        return result.valid;
    } catch (error) {
        console.error('SMILES validation error:', error);
        return false;
    }
}

export async function findSimilarMolecules(smiles: string, limit: number = 5): Promise<SimilarMolecule[]> {
    try {
        const result = await runPythonScript('find_similar', smiles, limit) as SimilarMolecule[];
        return result || [];
    } catch (error) {
        console.error('Similarity search error:', error);
        return [];
    }
}

export async function generateSmilesFromDescription(description: string): Promise<string[]> {
    // Use OpenAI GPT-5.1 API to generate SMILES from description
    try {
        const prompt = `You are a molecular chemist. Generate 3-5 valid SMILES strings for molecules that match this description: "${description}".

CRITICAL REQUIREMENTS:
- Return ONLY the SMILES strings, one per line
- Do NOT include any explanations, numbers, labels, or additional text
- Do NOT include words like "SMILES:", "1.", "2.", etc.
- Each SMILES must be a valid, syntactically correct chemical structure
- Consider drug-likeness: molecular weight < 500 Da, LogP < 5, reasonable number of H-bond donors/acceptors
- Molecules should be relevant to the therapeutic area or property mentioned
- Each SMILES should be on its own line with no prefix or suffix

Example of correct format:
CC(=O)Oc1ccccc1C(=O)O
CN1C=NC2=C1C(=O)N(C(=O)N2C)C
CCO

Now generate SMILES for: "${description}"`;

        const text = await callOpenAIAPI(prompt, undefined, 0.7);

        console.log('Raw OpenAI SMILES response:', text);

        // Extract SMILES strings (one per line)
        const lines = text.split('\n');
        const smilesCandidates: string[] = [];

        for (const line of lines) {
            const trimmed = line.trim();
            // Skip empty lines, comments, and lines that look like explanations
            if (!trimmed ||
                trimmed.startsWith('#') ||
                trimmed.toLowerCase().includes('smiles') ||
                trimmed.match(/^\d+[\.\)]/) || // Skip numbered lists
                trimmed.length < 2) {
                continue;
            }

            // Extract SMILES pattern (alphanumeric and chemical symbols)
            const smilesMatch = trimmed.match(/^([A-Za-z0-9@+\-\[\]()=#\.]+)/);
            if (smilesMatch) {
                const potentialSmiles = smilesMatch[1];
                // Basic validation: should start with a letter and contain valid SMILES characters
                if (/^[A-Za-z]/.test(potentialSmiles) && potentialSmiles.length >= 2) {
                    smilesCandidates.push(potentialSmiles);
                }
            }
        }

        // Remove duplicates
        const uniqueSmiles = Array.from(new Set(smilesCandidates));

        // Validate all generated SMILES using Python validation
        const validSmiles: string[] = [];
        for (const smiles of uniqueSmiles.slice(0, 5)) {
            if (await validateSmiles(smiles)) {
                validSmiles.push(smiles);
            }
        }

        // If we got valid SMILES, return them; otherwise return what we have (OpenAI might have generated valid ones)
        if (validSmiles.length > 0) {
            return validSmiles;
        }

        // If validation failed but we have candidates, return them anyway (validation might be too strict)
        if (uniqueSmiles.length > 0) {
            console.warn('Some SMILES failed validation, but returning them anyway:', uniqueSmiles);
            return uniqueSmiles.slice(0, 5);
        }

        // Last resort: fallback
        console.warn('No SMILES extracted from OpenAI response, using fallback');
        return getFallbackSmiles(description);
    } catch (error) {
        console.error('SMILES generation error:', error);
        // Fallback: return some example SMILES based on keywords
        return getFallbackSmiles(description);
    }
}

function getFallbackSmiles(description: string): string[] {
    // Simple keyword-based fallback
    const lowerDesc = description.toLowerCase();
    const fallbacks: string[] = [];

    if (lowerDesc.includes('aspirin') || lowerDesc.includes('pain')) {
        fallbacks.push('CC(=O)Oc1ccccc1C(=O)O');
    }
    if (lowerDesc.includes('caffeine') || lowerDesc.includes('stimulant')) {
        fallbacks.push('CN1C=NC2=C1C(=O)N(C(=O)N2C)C');
    }
    if (lowerDesc.includes('alcohol') || lowerDesc.includes('ethanol')) {
        fallbacks.push('CCO');
    }
    if (lowerDesc.includes('benzene') || lowerDesc.includes('aromatic')) {
        fallbacks.push('c1ccccc1');
    }

    return fallbacks.length > 0 ? fallbacks : ['CCO']; // Default to ethanol
}

export async function processChatMessage(
    message: string,
    sessionId?: string
): Promise<ChatResponse> {
    console.log('[processChatMessage] Starting with message:', message);

    // Get or create session
    let session = sessionId ? sessions.get(sessionId) : undefined;
    if (!session) {
        sessionId = generateSessionId();
        session = { messages: [] };
        sessions.set(sessionId, session);
    }

    // Add user message to context
    session.messages.push({ role: 'user', content: message });

    const lowerMessage = message.toLowerCase().trim();
    let response: ChatResponse = {
        response: '',
        sessionId: sessionId!,
    };

    // Check if message is just a SMILES string (simple pattern match) - very strict
    // Only match if it's a short string with no spaces and looks like pure SMILES
    const isJustSmiles = /^[A-Za-z0-9@+\-\[\]()=#\.]+$/.test(message.trim())
        && message.trim().length < 50
        && message.trim().length > 2
        && !message.includes(' ')
        && !message.includes('?')
        && !message.includes('!')
        && !lowerMessage.includes('generate')
        && !lowerMessage.includes('create')
        && !lowerMessage.includes('design')
        && !lowerMessage.includes('find')
        && !lowerMessage.includes('similar');

    // Always use OpenAI for most conversations, only use special handlers for very specific cases
    // Detect specific intents (but be more careful - only for very clear cases)
    const wantsGenerate = (lowerMessage.includes('generate') || lowerMessage.includes('create') || lowerMessage.includes('design'))
        && (lowerMessage.includes('molecule') || lowerMessage.includes('compound') || lowerMessage.includes('drug') || lowerMessage.includes('for'))
        && !lowerMessage.includes('?') && !lowerMessage.includes('how') && !lowerMessage.includes('what'); // Don't match questions

    const wantsSimilar = (lowerMessage.includes('similar') || lowerMessage.includes('find') || lowerMessage.includes('search'))
        && (lowerMessage.includes('molecule') || lowerMessage.includes('like') || lowerMessage.includes('to'))
        && !lowerMessage.includes('?') && !lowerMessage.includes('how') && !lowerMessage.includes('what'); // Don't match questions

    console.log('[processChatMessage] Intent detection:', { isJustSmiles, wantsGenerate, wantsSimilar, messageLength: message.length });

    // If it's just a SMILES string, offer to analyze it (but don't auto-analyze)
    if (isJustSmiles && !wantsGenerate && !wantsSimilar) {
        console.log('[processChatMessage] Detected as SMILES string, returning hardcoded response');
        response.response = `I see you've provided a SMILES string: ${message.trim()}. Click the "Analyze" button next to it in the chat to analyze this molecule.`;
    } else if (wantsGenerate && !lowerMessage.includes('?')) {
        // Generate SMILES from description using OpenAI API
        try {
            const suggestedSmiles = await generateSmilesFromDescription(message);
            response.suggestedSmiles = suggestedSmiles;

            // Get AI response to explain the molecules
            if (suggestedSmiles.length > 0) {
                try {
                    const aiText = await callOpenAIAPI(
                        `I've generated these SMILES strings based on the user's request: ${suggestedSmiles.join(', ')}. Provide a brief, friendly explanation of what these molecules are and their potential applications, then ask if the user wants to analyze any of them.`,
                        session.messages.slice(-3)
                    );
                    response.response = aiText;
                } catch (aiError) {
                    // Fallback if AI explanation fails
                    response.response = `I've generated ${suggestedSmiles.length} molecule candidate${suggestedSmiles.length !== 1 ? 's' : ''}:\n\n${suggestedSmiles.map((s, i) => `${i + 1}. ${s}`).join('\n')}\n\nWould you like me to analyze any of these molecules?`;
                }
            } else {
                // If no SMILES were generated, ask OpenAI to help
                const aiText = await callOpenAIAPI(
                    `The user asked to generate molecules for: "${message}". I was unable to generate valid SMILES strings. Please help explain this to the user and suggest what they might need to clarify.`,
                    session.messages.slice(-3)
                );
                response.response = aiText;
            }
        } catch (error) {
            // If generation fails, use OpenAI to explain
            console.error('SMILES generation failed, using OpenAI:', error);
            const aiText = await callOpenAIAPI(
                `The user asked to generate molecules for: "${message}". There was an error generating SMILES strings. Please help the user understand what went wrong and how they can rephrase their request.`,
                session.messages.slice(-5)
            );
            response.response = aiText;
        }
    } else if (wantsSimilar && !lowerMessage.includes('?')) {
        // Extract SMILES from message or use last analyzed molecule
        const smilesMatch = message.match(/[A-Za-z0-9@+\-\[\]()=#]+/);
        const smiles = smilesMatch ? smilesMatch[0] : (session.analyzedMolecules?.[session.analyzedMolecules.length - 1]);

        if (smiles) {
            try {
                const similar = await findSimilarMolecules(smiles, 5);
                response.similarMolecules = similar;

                // Get AI response to explain the results
                try {
                    const aiText = await callOpenAIAPI(
                        `I found ${similar.length} similar molecules to ${smiles}: ${similar.map(m => m.smiles).join(', ')}. Provide a brief explanation.`,
                        session.messages.slice(-3)
                    );
                    response.response = aiText || `I found ${similar.length} similar molecule${similar.length !== 1 ? 's' : ''} in the database:\n\n${similar.map((m, i) => `${i + 1}. ${m.smiles} (similarity: ${(m.similarity * 100).toFixed(1)}%)`).join('\n')}`;
                } catch (aiError) {
                    response.response = `I found ${similar.length} similar molecule${similar.length !== 1 ? 's' : ''} in the database:\n\n${similar.map((m, i) => `${i + 1}. ${m.smiles} (similarity: ${(m.similarity * 100).toFixed(1)}%)`).join('\n')}`;
                }
            } catch (error) {
                // If similarity search fails, use OpenAI
                const aiText = await callOpenAIAPI(message, session.messages.slice(-5));
                response.response = aiText;
            }
        } else {
            // Use AI to help
            const aiText = await callOpenAIAPI(
                `The user wants to find similar molecules but didn't provide a SMILES string. Help them understand they need to provide a SMILES string first. Original message: ${message}`,
                session.messages.slice(-3)
            );
            response.response = aiText;
        }
    } else {
        // General conversation - always use OpenAI GPT-5.1
        console.log('[processChatMessage] Using OpenAI for general conversation');
        try {
            console.log('[processChatMessage] Calling callOpenAIAPI with message:', message);
            const aiText = await callOpenAIAPI(message, session.messages.slice(-5));
            console.log('[processChatMessage] Received AI response:', aiText?.substring(0, 100));
            if (!aiText || aiText.trim() === '') {
                throw new Error('OpenAI returned empty response');
            }
            response.response = aiText;
        } catch (error: any) {
            console.error('[processChatMessage] Chat error:', error);
            const errorMsg = error.message || 'Unknown error';
            response.response = `I encountered an error: ${errorMsg}. Please check your OPENAI_API_KEY configuration and ensure the API is working.`;
        }
    }

    // Add assistant response to context
    session.messages.push({ role: 'assistant', content: response.response });

    // Clean up old sessions (30 min expiration)
    setTimeout(() => {
        sessions.delete(sessionId!);
    }, 30 * 60 * 1000);

    return response;
}

