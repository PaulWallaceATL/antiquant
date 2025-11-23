'use client';

import React, { useEffect, useRef, useState, useMemo } from 'react';
import { Ketcher, RemoteStructServiceProvider } from 'ketcher-core';
import { Editor } from 'ketcher-react';
import 'ketcher-react/dist/index.css';

interface KetcherEditorProps {
    initialSmiles?: string;
    onSmilesChange: (smiles: string) => void;
    width?: number;
    height?: number;
}

export default function KetcherEditor({ 
    initialSmiles = '', 
    onSmilesChange, 
    width = 600, 
    height = 400 
}: KetcherEditorProps) {
    const ketcherRef = useRef<Ketcher | null>(null);
    const [error, setError] = useState<string | null>(null);
    const containerRef = useRef<HTMLDivElement>(null);
    const isInitializedRef = useRef(false);
    
    // Memoize the struct service provider to avoid recreating it on every render
    // This prevents React warnings about state updates during render
    const structServiceProvider = useMemo(() => {
        return new RemoteStructServiceProvider('/api/ketcher');
    }, []);

    const handleOnInit = async (ketcher: Ketcher) => {
        // Prevent multiple initializations
        if (isInitializedRef.current) {
            return;
        }
        isInitializedRef.current = true;

        try {
            ketcherRef.current = ketcher;
            
            // Wait for Ketcher to fully initialize - use a longer delay
            // Ketcher needs time to set up its logger, editor, Redux store, and all internal components
            await new Promise(resolve => setTimeout(resolve, 1000));
            
            // Set initial SMILES if provided - use a more robust approach
            if (initialSmiles && initialSmiles.trim()) {
                // Wait for editor to be fully ready, then set SMILES with retries
                const setSmilesWithRetry = async (retries = 5) => {
                    for (let i = 0; i < retries; i++) {
                        try {
                            // Increasing delay with each retry
                            await new Promise(resolve => setTimeout(resolve, 500 + (i * 200)));
                            
                            // Check if ketcher is ready by trying to get current structure
                            try {
                                await ketcher.getSmiles();
                            } catch (checkError) {
                                // If we can't even get SMILES, ketcher isn't ready yet
                                if (i < retries - 1) continue;
                            }
                            
                            // Try to set the molecule
                            await ketcher.setMolecule(initialSmiles);
                            console.log('[KetcherEditor] ✓ Successfully set initial SMILES:', initialSmiles);
                            return; // Success, exit retry loop
                        } catch (e: any) {
                            const errorMsg = e?.message || String(e);
                            if (errorMsg.includes('KetcherLogger') || errorMsg.includes('initialized')) {
                                // Still initializing, continue retrying
                                if (i < retries - 1) {
                                    console.log(`[KetcherEditor] Ketcher still initializing, retry ${i + 1}/${retries}...`);
                                    continue;
                                }
                            }
                            if (i === retries - 1) {
                                console.warn('[KetcherEditor] Could not set initial SMILES after retries:', e);
                            }
                        }
                    }
                };
                
                // Start the retry process
                setSmilesWithRetry();
            }

            // Listen for structure changes - wait for editor to be fully ready
            setTimeout(() => {
                try {
                    ketcher.editor.subscribe('change', async () => {
                        requestAnimationFrame(async () => {
                            try {
                                const smiles = await ketcher.getSmiles();
                                if (smiles && smiles.trim()) {
                                    onSmilesChange(smiles);
                                }
                            } catch (e) {
                                // Ignore errors during structure changes
                            }
                        });
                    });
                    console.log('[KetcherEditor] ✓ Subscribed to structure changes');
                } catch (e) {
                    console.warn('[KetcherEditor] Could not subscribe to changes:', e);
                }
            }, 1500);
        } catch (e: any) {
            console.error('[KetcherEditor] Initialization error:', e);
            // Defer error state update to avoid render-time updates
            setTimeout(() => {
                setError('Failed to initialize molecular editor. Please use Text Mode.');
            }, 0);
        }
    };

    // Update SMILES when initialSmiles prop changes
    useEffect(() => {
        if (ketcherRef.current && initialSmiles && initialSmiles.trim() && isInitializedRef.current) {
            // Use setTimeout to avoid render-time updates
            const timeoutId = setTimeout(async () => {
                try {
                    const currentSmiles = await ketcherRef.current!.getSmiles();
                    if (initialSmiles !== currentSmiles) {
                        await ketcherRef.current!.setMolecule(initialSmiles);
                    }
                } catch (e) {
                    // Silently fail - user might be editing
                }
            }, 0);
            
            return () => clearTimeout(timeoutId);
        }
    }, [initialSmiles]);

    if (error) {
        return (
            <div className="p-4 border rounded bg-red-50 dark:bg-red-900/20">
                <p className="text-sm text-red-600 dark:text-red-400 font-medium mb-2">{error}</p>
                <p className="text-xs text-gray-600 dark:text-gray-400">
                    Please use Text Mode to enter SMILES strings directly.
                </p>
            </div>
        );
    }

    return (
        <div 
            ref={containerRef}
            className="border rounded-lg overflow-hidden bg-white dark:bg-gray-900" 
            style={{ width, height }}
        >
            <Editor 
                onInit={handleOnInit}
                staticResourcesUrl="/"
                structServiceProvider={structServiceProvider}
                errorHandler={(message: any) => {
                    // Convert message to string if it's not already
                    const messageStr = typeof message === 'string' ? message : String(message || 'Unknown error');
                    console.error('[KetcherEditor] Error:', messageStr);
                    // Defer error state updates to avoid render-time updates
                    setTimeout(() => {
                        // Only show user-facing errors, not internal API errors
                        if (messageStr && typeof messageStr === 'string' && 
                            !messageStr.includes('createStructService') && 
                            !messageStr.includes('api') && 
                            !messageStr.includes('structService') && 
                            !messageStr.includes('Indigo') &&
                            !messageStr.includes('KetcherLogger')) {
                            setError(messageStr);
                        }
                    }, 0);
                }}
            />
        </div>
    );
}
