'use client';

import React, { useState, useRef, useEffect } from 'react';
import { Send, Loader2, Beaker, Search } from 'lucide-react';
import { ChatResponse, SimilarMolecule } from '../lib/chat_service';

interface ChatMessage {
    role: 'user' | 'assistant';
    content: string;
    suggestedSmiles?: string[];
    similarMolecules?: SimilarMolecule[];
}

interface ChatInterfaceProps {
    onAnalyzeSmiles?: (smiles: string) => void;
}

const STARTER_PROMPTS = [
    "Generate a molecule for treating diabetes",
    "Find molecules similar to aspirin",
    "Create a pain relief compound",
    "Design an anti-inflammatory drug",
    "Generate a molecule with high solubility",
    "Find similar molecules to caffeine",
];

export default function ChatInterface({ onAnalyzeSmiles }: ChatInterfaceProps) {
    const [messages, setMessages] = useState<ChatMessage[]>([
        {
            role: 'assistant',
            content: "Hello! I'm your molecular analysis assistant. I can help you:\n\n• Generate SMILES strings from descriptions\n• Find similar molecules in our database\n• Analyze molecular properties\n\nTry one of the starter prompts below or ask me anything!"
        }
    ]);
    const [input, setInput] = useState('');
    const [loading, setLoading] = useState(false);
    const [sessionId, setSessionId] = useState<string | undefined>();
    const messagesEndRef = useRef<HTMLDivElement>(null);

    const scrollToBottom = () => {
        messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
    };

    useEffect(() => {
        scrollToBottom();
    }, [messages]);

    const sendMessage = async (messageText: string) => {
        if (!messageText.trim() || loading) return;

        const userMessage = messageText.trim();
        setLoading(true);

        // Add user message to UI
        setMessages(prev => [...prev, { role: 'user', content: userMessage }]);

        try {
            const response = await fetch('/api/chat', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ message: userMessage, sessionId }),
            });

            if (!response.ok) {
                const errorData = await response.json().catch(() => ({ error: 'Unknown error' }));
                throw new Error(errorData.error || `HTTP ${response.status}`);
            }

            const data: ChatResponse = await response.json();

            if (data.sessionId) {
                setSessionId(data.sessionId);
            }

            // Add assistant response - ensure content is always a string
            let responseContent = data.response;
            if (typeof responseContent !== 'string') {
                console.error('Response is not a string:', typeof responseContent, responseContent);
                // If it's an object, try to extract a meaningful string
                if (responseContent && typeof responseContent === 'object') {
                    responseContent = JSON.stringify(responseContent);
                } else {
                    responseContent = 'I apologize, but I did not receive a valid response.';
                }
            }
            
            const assistantMessage: ChatMessage = {
                role: 'assistant',
                content: String(responseContent || 'I apologize, but I did not receive a valid response.'),
                suggestedSmiles: data.suggestedSmiles,
                similarMolecules: data.similarMolecules,
            };

            setMessages(prev => [...prev, assistantMessage]);

            // DO NOT auto-analyze - only analyze when user clicks the "Analyze" button
            // Removed auto-analysis to prevent SMILES input from being updated by chat text
        } catch (error: any) {
            console.error('Chat error:', error);
            const errorMessage = error.message || 'Sorry, I encountered an error. Please try again.';
            setMessages(prev => [...prev, {
                role: 'assistant',
                content: `Error: ${errorMessage}`
            }]);
        } finally {
            setLoading(false);
        }
    };

    const handleSend = async () => {
        if (!input.trim() || loading) return;
        const message = input.trim();
        setInput('');
        await sendMessage(message);
    };

    const handleAnalyze = (smiles: string) => {
        // Only analyze when user explicitly clicks the Analyze button
        // Validate it's a proper SMILES string
        const trimmedSmiles = smiles.trim();
        if (!onAnalyzeSmiles || !trimmedSmiles) {
            return;
        }
        
        // Strict SMILES validation - must contain only valid SMILES characters
        // SMILES can contain: letters, numbers, @, +, -, [, ], (, ), =, #, and common elements
        const smilesPattern = /^[A-Za-z0-9@+\-\[\]()=#\.]+$/;
        if (!smilesPattern.test(trimmedSmiles) || trimmedSmiles.length < 2) {
            console.warn('Invalid SMILES string:', trimmedSmiles);
            return;
        }
        
        // Only proceed if it looks like a valid SMILES (not chat text)
        // SMILES typically start with a capital letter (element) or lowercase (aromatic)
        const firstChar = trimmedSmiles[0];
        if (!/[A-Za-z]/.test(firstChar)) {
            console.warn('SMILES must start with a letter:', trimmedSmiles);
            return;
        }
        
        onAnalyzeSmiles(trimmedSmiles);
    };

    return (
        <div className="flex flex-col h-full bg-white dark:bg-gray-900">
            {/* Messages */}
            <div className="flex-1 overflow-y-auto p-4 space-y-4">
                {messages.map((msg, idx) => (
                    <div
                        key={idx}
                        className={`flex ${msg.role === 'user' ? 'justify-end' : 'justify-start'}`}
                    >
                        <div
                            className={`max-w-[80%] rounded-lg px-4 py-2 ${
                                msg.role === 'user'
                                    ? 'bg-blue-600 text-white'
                                    : 'bg-gray-100 dark:bg-gray-800 text-gray-900 dark:text-gray-100'
                            }`}
                        >
                            <p className="whitespace-pre-wrap text-sm">{msg.content}</p>

                            {/* Suggested SMILES */}
                            {msg.suggestedSmiles && msg.suggestedSmiles.length > 0 && (
                                <div className="mt-3 space-y-2">
                                    <p className="text-xs font-semibold mb-2">Suggested Molecules:</p>
                                    {msg.suggestedSmiles.map((smiles, i) => (
                                        <div
                                            key={i}
                                            className="flex items-center justify-between bg-white dark:bg-gray-700 rounded p-2"
                                        >
                                            <code className="text-xs font-mono text-gray-800 dark:text-gray-200">
                                                {smiles}
                                            </code>
                                            <button
                                                onClick={() => handleAnalyze(smiles)}
                                                className="ml-2 px-2 py-1 bg-blue-600 text-white text-xs rounded hover:bg-blue-700 flex items-center gap-1"
                                            >
                                                <Beaker className="w-3 h-3" />
                                                Analyze
                                            </button>
                                        </div>
                                    ))}
                                </div>
                            )}

                            {/* Similar Molecules */}
                            {msg.similarMolecules && msg.similarMolecules.length > 0 && (
                                <div className="mt-3 space-y-2">
                                    <p className="text-xs font-semibold mb-2">Similar Molecules Found:</p>
                                    {msg.similarMolecules.map((mol, i) => (
                                        <div
                                            key={i}
                                            className="flex items-center justify-between bg-white dark:bg-gray-700 rounded p-2"
                                        >
                                            <div className="flex-1">
                                                <code className="text-xs font-mono text-gray-800 dark:text-gray-200 block">
                                                    {mol.smiles}
                                                </code>
                                                <span className="text-xs text-gray-500 dark:text-gray-400">
                                                    Similarity: {(mol.similarity * 100).toFixed(1)}%
                                                </span>
                                            </div>
                                            <button
                                                onClick={() => handleAnalyze(mol.smiles)}
                                                className="ml-2 px-2 py-1 bg-blue-600 text-white text-xs rounded hover:bg-blue-700 flex items-center gap-1"
                                            >
                                                <Search className="w-3 h-3" />
                                                Analyze
                                            </button>
                                        </div>
                                    ))}
                                </div>
                            )}
                        </div>
                    </div>
                ))}

                {loading && (
                    <div className="flex justify-start">
                        <div className="bg-gray-100 dark:bg-gray-800 rounded-lg px-4 py-2">
                            <Loader2 className="w-4 h-4 animate-spin text-gray-600 dark:text-gray-400" />
                        </div>
                    </div>
                )}

                <div ref={messagesEndRef} />
            </div>

            {/* Starter Prompts */}
            {messages.length === 1 && (
                <div className="border-t border-gray-200 dark:border-gray-800 p-4">
                    <p className="text-xs text-gray-500 dark:text-gray-400 mb-2">Try these prompts:</p>
                    <div className="flex flex-wrap gap-2">
                        {STARTER_PROMPTS.map((prompt, idx) => (
                            <button
                                key={idx}
                                onClick={() => sendMessage(prompt)}
                                disabled={loading}
                                className="px-3 py-1.5 text-xs bg-gray-100 dark:bg-gray-800 hover:bg-gray-200 dark:hover:bg-gray-700 rounded-lg transition-colors text-left disabled:opacity-50 disabled:cursor-not-allowed"
                            >
                                {prompt}
                            </button>
                        ))}
                    </div>
                </div>
            )}

            {/* Input */}
            <div className="border-t border-gray-200 dark:border-gray-800 p-4">
                <div className="flex gap-2">
                    <input
                        type="text"
                        value={input}
                        onChange={(e) => setInput(e.target.value)}
                        onKeyPress={(e) => {
                            if (e.key === 'Enter' && !e.shiftKey) {
                                e.preventDefault();
                                handleSend();
                            }
                        }}
                        placeholder="Ask me to generate molecules, find similar ones, or analyze properties..."
                        className="flex-1 px-4 py-2 border border-gray-300 dark:border-gray-700 rounded-lg focus:ring-2 focus:ring-blue-500 outline-none bg-white dark:bg-gray-800 text-gray-900 dark:text-gray-100"
                        disabled={loading}
                    />
                    <button
                        onClick={handleSend}
                        disabled={loading || !input.trim()}
                        className="px-4 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700 disabled:opacity-50 disabled:cursor-not-allowed flex items-center gap-2"
                    >
                        {loading ? (
                            <Loader2 className="w-4 h-4 animate-spin" />
                        ) : (
                            <Send className="w-4 h-4" />
                        )}
                    </button>
                </div>
            </div>
        </div>
    );
}

