'use client';

import React, { useState, useEffect, useCallback } from 'react';
import dynamic from 'next/dynamic';

// Dynamically import Ketcher editor to avoid SSR issues
const KetcherEditor = dynamic(() => import('./KetcherEditor'), {
    ssr: false,
    loading: () => (
        <div className="h-[400px] flex items-center justify-center border rounded-lg bg-gray-50 dark:bg-gray-800">
            <p className="text-sm text-gray-600 dark:text-gray-400">Loading molecular editor...</p>
        </div>
    ),
});

interface MolecularEditorProps {
    initialSmiles: string;
    onChange: (smiles: string) => void;
}

export default function MolecularEditor({ initialSmiles, onChange }: MolecularEditorProps) {
    // Default to text mode - it's more reliable, works offline, and isn't blocked by service workers
    const [showVisualEditor, setShowVisualEditor] = useState(false);
    const [smiles, setSmiles] = useState(initialSmiles);
    const [hasServiceWorker, setHasServiceWorker] = useState(false);

    // Check for service workers on mount
    useEffect(() => {
        if ('serviceWorker' in navigator) {
            navigator.serviceWorker.getRegistrations().then(registrations => {
                if (registrations.length > 0) {
                    setHasServiceWorker(true);
                }
            });
        }
    }, []);

    // Sync with external changes
    useEffect(() => {
        setSmiles(initialSmiles);
    }, [initialSmiles]);

    // Debounced onChange to avoid too many updates
    const debouncedOnChange = useCallback(
        (() => {
            let timeout: NodeJS.Timeout;
            return (newSmiles: string) => {
                clearTimeout(timeout);
                timeout = setTimeout(() => {
                    onChange(newSmiles);
                }, 300);
            };
        })(),
        [onChange]
    );

    const handleTextChange = (e: React.ChangeEvent<HTMLTextAreaElement>) => {
        const newSmiles = e.target.value.trim();
        setSmiles(newSmiles);
        debouncedOnChange(newSmiles);
    };

    const handleVisualChange = (newSmiles: string) => {
        if (newSmiles && newSmiles.trim()) {
            setSmiles(newSmiles.trim());
            onChange(newSmiles.trim());
        }
    };

    return (
        <div className="space-y-4">
            <div className="flex justify-between items-center">
                <h3 className="font-semibold text-sm">Molecular Editor</h3>
                <button
                    onClick={() => setShowVisualEditor(!showVisualEditor)}
                    className="text-xs text-blue-600 dark:text-blue-400 hover:underline font-medium transition-colors"
                    disabled={hasServiceWorker}
                    title={hasServiceWorker ? 'Visual Mode disabled due to active service worker' : ''}
                >
                    {showVisualEditor ? 'Switch to Text Mode' : 'Switch to Visual Mode'}
                </button>
            </div>
            
            {!showVisualEditor && (
                <div className="text-xs text-green-600 dark:text-green-400 bg-green-50 dark:bg-green-900/20 p-2 rounded border border-green-200 dark:border-green-800">
                    <strong>✓ Text Mode:</strong> The most reliable way to enter SMILES. Paste SMILES strings here or type them directly. This works perfectly without any backend requirements.
                </div>
            )}

            {showVisualEditor ? (
                <div className="space-y-2">
                    <div className="p-4 border rounded-lg bg-amber-50 dark:bg-amber-900/20 border-amber-200 dark:border-amber-800">
                        <p className="text-sm text-amber-800 dark:text-amber-200 font-medium mb-2">
                            ⚠️ Visual Editor Limitations
                        </p>
                        <p className="text-xs text-amber-700 dark:text-amber-300 mb-3">
                            The visual editor requires a full Indigo backend service to render structures from SMILES. 
                            Without it, structures may not display correctly. <strong>Text Mode below works perfectly</strong> and is the recommended way to enter SMILES strings.
                        </p>
                        <p className="text-xs text-amber-600 dark:text-amber-400">
                            You can still use the visual editor to draw new structures, but SMILES input may not render.
                        </p>
                    </div>
                    <KetcherEditor
                        initialSmiles={smiles}
                        onSmilesChange={handleVisualChange}
                        width={600}
                        height={400}
                    />
                </div>
            ) : (
                <div className="space-y-2">
                    <textarea
                        value={smiles}
                        onChange={handleTextChange}
                        className="w-full h-32 p-4 font-mono text-sm border rounded-lg focus:ring-2 focus:ring-blue-500 outline-none resize-none bg-white dark:bg-gray-800 text-gray-900 dark:text-gray-100"
                        placeholder="Enter SMILES notation (e.g., CCO for ethanol)..."
                    />
                    <p className="text-xs text-gray-500 dark:text-gray-400">
                        Enter a SMILES string or switch to visual mode to draw the molecule
                    </p>
                </div>
            )}
        </div>
    );
}
