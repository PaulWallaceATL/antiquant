'use client';

import React, { useEffect, useRef, useState } from 'react';

declare global {
    interface Window {
        JSME: any;
    }
}

interface JsmeEditorProps {
    initialSmiles?: string;
    onSmilesChange: (smiles: string) => void;
    width?: number;
    height?: number;
}

export default function JsmeEditor({ 
    initialSmiles = '', 
    onSmilesChange, 
    width = 400, 
    height = 300 
}: JsmeEditorProps) {
    const containerRef = useRef<HTMLDivElement>(null);
    const editorRef = useRef<any>(null);
    const [loaded, setLoaded] = useState(false);
    const [error, setError] = useState<string | null>(null);
    const editorIdRef = useRef<string | null>(null);
    const isInitializingRef = useRef(false);

    // Load JSME script dynamically (client-side only)
    useEffect(() => {
        if (typeof window === 'undefined') return;

        // Check for service workers that might interfere - if present, disable visual editor immediately
        const checkServiceWorkers = async () => {
            if ('serviceWorker' in navigator) {
                try {
                    const registrations = await navigator.serviceWorker.getRegistrations();
                    if (registrations.length > 0) {
                        console.warn('[JsmeEditor] Service workers detected - blocking visual editor');
                        setError('⚠️ Service Worker Detected: An active service worker is blocking external scripts. The visual editor cannot load. Please unregister service workers in DevTools (Application → Service Workers → Unregister) or use Text Mode below, which works offline.');
                        return true; // Service worker found, abort
                    }
                } catch (e) {
                    console.warn('[JsmeEditor] Error checking service workers:', e);
                }
            }
            return false; // No service worker, continue
        };

        // Check service workers first, then proceed if none found
        (async () => {
            const hasServiceWorker = await checkServiceWorkers();
            if (hasServiceWorker) {
                return; // Don't attempt to load JSME if service worker is active
            }

            // Continue with normal loading...
            // Check if JSME is already loaded
            if (window.JSME) {
                console.log('[JsmeEditor] JSME already loaded');
                setLoaded(true);
                return;
            }

        // Check if script is already being loaded
        const existingScript = document.querySelector('script[src*="jsme"], script[id*="jsme"]') as HTMLScriptElement | null;
        if (existingScript) {
            console.log('[JsmeEditor] JSME script already exists, waiting for load...');
            if (window.JSME && typeof window.JSME === 'function') {
                setLoaded(true);
                return;
            }
            const handleLoad = () => {
                if (window.JSME) {
                    setLoaded(true);
                }
            };
            existingScript.addEventListener('load', handleLoad);
            existingScript.addEventListener('error', () => {
                setError('Failed to load JSME editor script.');
            });
            return () => {
                existingScript.removeEventListener('load', handleLoad);
            };
        }

        console.log('[JsmeEditor] Loading JSME script...');
        
        // Strategy: Try proxy first (bypasses service worker), then direct URLs
        const loadViaProxy = async () => {
            try {
                console.log('[JsmeEditor] Attempting to load via proxy API (bypasses service worker)...');
                const response = await fetch('/api/jsme-proxy', {
                    cache: 'no-store', // Don't cache to avoid stale responses
                });
                
                if (!response.ok) {
                    throw new Error(`Proxy returned ${response.status}`);
                }
                
                const scriptContent = await response.text();
                
                if (!scriptContent || scriptContent.length < 100) {
                    throw new Error('Proxy returned empty or invalid script');
                }
                
                // Check if response is actually JavaScript (not HTML error page or JSON)
                const trimmed = scriptContent.trim();
                if (trimmed.startsWith('<') || trimmed.startsWith('{') || trimmed.startsWith('<!')) {
                    console.error('[JsmeEditor] Proxy returned non-JavaScript content:', trimmed.substring(0, 200));
                    throw new Error('Proxy returned HTML/JSON instead of JavaScript - likely an error page');
                }
                
                // Inject script content directly (bypasses service worker and CORS)
                const script = document.createElement('script');
                script.textContent = scriptContent;
                script.type = 'text/javascript';
                script.id = 'jsme-proxy-loaded';
                document.head.appendChild(script);
                
                console.log('[JsmeEditor] Script injected via proxy, waiting for JSME to initialize...');
                
                // Wait for JSME to initialize
                const checkInterval = setInterval(() => {
                    if (window.JSME && typeof window.JSME === 'function') {
                        console.log('[JsmeEditor] ✓ JSME loaded successfully via proxy!');
                        clearInterval(checkInterval);
                        setLoaded(true);
                        setError(null);
                    }
                }, 200);
                
                setTimeout(() => {
                    clearInterval(checkInterval);
                    if (!window.JSME) {
                        console.warn('[JsmeEditor] JSME not available after proxy load, trying direct URLs...');
                        loadViaDirectURLs();
                    }
                }, 10000);
            } catch (proxyError) {
                console.warn('[JsmeEditor] Proxy load failed:', proxyError);
                loadViaDirectURLs();
            }
        };
        
        const loadViaDirectURLs = () => {
            // Use official JSME website - these URLs work but may be blocked by service workers
            const cdnSources = [
                'https://peter-ertl.com/jsme/JSME_2024-04-29/jsme.nocache.js',
                'https://peter-ertl.com/jsme/JSME_2023-07-31/jsme.nocache.js',
                'https://peter-ertl.com/jsme/JSME_2022-09-26/jsme.nocache.js',
            ];

            let attemptCount = 0;
            const maxAttempts = cdnSources.length * 2; // Try each source twice

            const tryLoadScript = (sourceIndex: number) => {
                attemptCount++;
                
                if (sourceIndex >= cdnSources.length) {
                if (attemptCount >= maxAttempts) {
                    setError('JSME editor failed to load after multiple attempts. The visual editor requires an internet connection. Please use Text Mode to enter SMILES strings directly, or check your network connection.');
                    return;
                }
                // Retry from beginning
                setTimeout(() => tryLoadScript(0), 1000);
                return;
                }

                console.log(`[JsmeEditor] Attempting to load from: ${cdnSources[sourceIndex]} (attempt ${attemptCount})`);

                // Check if this script is already being loaded
                const existingScriptId = `jsme-script-${sourceIndex}`;
                const existing = document.getElementById(existingScriptId);
                if (existing) {
                    console.log('[JsmeEditor] Script already exists, skipping duplicate load');
                    if (sourceIndex < cdnSources.length - 1) {
                        tryLoadScript(sourceIndex + 1);
                    }
                    return;
                }

                const script = document.createElement('script');
                script.src = cdnSources[sourceIndex];
                script.async = true;
                script.id = existingScriptId;
                // Don't set crossOrigin - script tags don't need it and it can cause MIME type issues
                
                let checkInterval: NodeJS.Timeout;
                let timeoutId: NodeJS.Timeout;
                let loaded = false;
                let errorHandled = false;
                
                const handleSuccess = () => {
                    if (loaded) return;
                    loaded = true;
                    console.log('[JsmeEditor] ✓ JSME successfully loaded and available!');
                    if (checkInterval) clearInterval(checkInterval);
                    if (timeoutId) clearTimeout(timeoutId);
                    setLoaded(true);
                    setError(null);
                };
                
                const handleError = () => {
                    if (errorHandled) return;
                    errorHandled = true;
                    console.warn('[JsmeEditor] Script error from:', cdnSources[sourceIndex]);
                    if (checkInterval) clearInterval(checkInterval);
                    if (timeoutId) clearTimeout(timeoutId);
                    
                    // Remove failed script
                    if (script.parentNode) {
                        try {
                            script.parentNode.removeChild(script);
                        } catch (removeError) {
                            // Ignore removal errors
                        }
                    }
                    
                    // Try next source
                    if (sourceIndex < cdnSources.length - 1) {
                        setTimeout(() => tryLoadScript(sourceIndex + 1), 1000);
                    } else if (attemptCount < maxAttempts) {
                        // Retry from beginning
                        setTimeout(() => tryLoadScript(0), 2000);
                    } else {
                        setError('JSME editor could not be loaded after multiple attempts. This may be due to network restrictions, ad blockers, or firewall settings. Please use Text Mode to enter SMILES strings directly - it works offline and is more reliable.');
                    }
                };
                
                script.onload = () => {
                    console.log('[JsmeEditor] Script load event fired from:', cdnSources[sourceIndex]);
                    
                    // Wait for JSME to actually be available (might take a moment)
                    checkInterval = setInterval(() => {
                        if (window.JSME && typeof window.JSME === 'function') {
                            handleSuccess();
                        }
                    }, 200);

                    // Timeout after 10 seconds
                    timeoutId = setTimeout(() => {
                        if (!window.JSME && !loaded) {
                            console.warn('[JsmeEditor] JSME not available after 10s, trying next source...');
                            handleError();
                        }
                    }, 10000);
                };

                script.onerror = handleError;

                // Add to document
                try {
                    document.head.appendChild(script);
                } catch (appendError) {
                    console.error('[JsmeEditor] Failed to append script:', appendError);
                    if (sourceIndex < cdnSources.length - 1) {
                        tryLoadScript(sourceIndex + 1);
                    } else {
                        setError('Failed to load JSME editor. Please use Text Mode.');
                    }
                }
            };

            // Start loading direct URLs
            tryLoadScript(0);
        };
        
            // Try proxy first (bypasses service worker), then fallback to direct URLs
            loadViaProxy();
        })();

        // Cleanup function
        return () => {
            // Don't remove the script as it might be used by other instances
        };
    }, []);

    // Initialize editor when script is loaded
    useEffect(() => {
        if (!loaded || !containerRef.current || isInitializingRef.current) return;
        
        // If editor already exists, just update SMILES
        if (editorRef.current) {
            try {
                if (initialSmiles && initialSmiles !== editorRef.current.getSmiles()) {
                    editorRef.current.readMolFile(initialSmiles);
                }
            } catch (e) {
                console.warn('Failed to update JSME with initial SMILES:', e);
            }
            return;
        }

        // Create unique ID for this editor instance
        const id = `jsme-editor-${Date.now()}-${Math.random().toString(36).substr(2, 9)}`;
        editorIdRef.current = id;

        // Create container div
        const div = document.createElement('div');
        div.id = id;
        div.style.width = `${width}px`;
        div.style.height = `${height}px`;
        containerRef.current.innerHTML = ''; // Clear any existing content
        containerRef.current.appendChild(div);

        isInitializingRef.current = true;

        try {
            // Verify JSME is actually a constructor function
            if (typeof window.JSME !== 'function') {
                throw new Error('JSME is not a function. It may not have loaded correctly.');
            }

            // Initialize JSME editor
            const editor = new window.JSME(id, width, height, {
                options: 'oldlook,star'
            });
            
            if (!editor) {
                throw new Error('Failed to create JSME editor instance');
            }
            
            editorRef.current = editor;

            // Set initial SMILES if provided (with delay to ensure editor is ready)
            if (initialSmiles && initialSmiles.trim()) {
                setTimeout(() => {
                    try {
                        // Try different methods to set SMILES
                        if (typeof editor.readMolFile === 'function') {
                            editor.readMolFile(initialSmiles);
                        } else if (typeof editor.setMolFile === 'function') {
                            editor.setMolFile(initialSmiles);
                        } else if (typeof editor.readSmiles === 'function') {
                            editor.readSmiles(initialSmiles);
                        } else if (typeof editor.setSmiles === 'function') {
                            editor.setSmiles(initialSmiles);
                        }
                    } catch (e) {
                        console.warn('[JsmeEditor] Could not set initial SMILES:', e);
                        // Don't fail if we can't set initial SMILES - user can draw it
                    }
                }, 200);
            }

            // Set up callback for structure changes
            const handleChange = () => {
                try {
                    const smiles = editor.getSmiles();
                    if (smiles && smiles.trim()) {
                        onSmilesChange(smiles);
                    }
                } catch (e) {
                    console.warn('Error getting SMILES from JSME:', e);
                }
            };

            editor.setCallBack('AfterStructureModified', handleChange);
            
            isInitializingRef.current = false;
            setError(null);
        } catch (e) {
            console.error('JSME initialization error:', e);
            setError(`Failed to initialize editor: ${e instanceof Error ? e.message : 'Unknown error'}`);
            isInitializingRef.current = false;
        }
    }, [loaded, width, height, initialSmiles]); // Removed onSmilesChange from deps

    // Cleanup on unmount
    useEffect(() => {
        return () => {
            if (editorRef.current) {
                try {
                    // JSME doesn't have a formal destroy method, but we can clear the container
                    if (containerRef.current) {
                        containerRef.current.innerHTML = '';
                    }
                } catch (e) {
                    console.warn('Error cleaning up JSME:', e);
                }
                editorRef.current = null;
            }
            isInitializingRef.current = false;
        };
    }, []);

    // Update SMILES when initialSmiles prop changes (but not during initialization)
    useEffect(() => {
        if (!editorRef.current || isInitializingRef.current || !initialSmiles) return;
        
        try {
            const currentSmiles = editorRef.current.getSmiles();
            if (initialSmiles !== currentSmiles) {
                editorRef.current.readMolFile(initialSmiles);
            }
        } catch (e) {
            // Silently fail - user might be editing
        }
    }, [initialSmiles]);

    if (error) {
        return (
            <div className="p-4 border rounded bg-red-50 dark:bg-red-900/20">
                <p className="text-sm text-red-600 dark:text-red-400 font-medium mb-2">{error}</p>
                <div className="space-y-2">
                    <button
                        onClick={async () => {
                            setError(null);
                            setLoaded(false);
                            
                            // Unregister service workers that might be interfering
                            if ('serviceWorker' in navigator) {
                                try {
                                    const registrations = await navigator.serviceWorker.getRegistrations();
                                    for (const registration of registrations) {
                                        await registration.unregister();
                                        console.log('[JsmeEditor] Unregistered service worker');
                                    }
                                } catch (e) {
                                    console.warn('[JsmeEditor] Could not unregister service workers:', e);
                                }
                            }
                            
                            // Clear JSME scripts
                            const scripts = document.querySelectorAll('script[src*="jsme"]');
                            scripts.forEach(s => s.remove());
                            delete (window as any).JSME;
                            
                            // Hard reload to clear cache
                            window.location.reload();
                        }}
                        className="block text-xs text-red-600 dark:text-red-400 underline hover:no-underline"
                    >
                        Reload Page (Clears Cache & Service Workers)
                    </button>
                    <p className="text-xs text-gray-600 dark:text-gray-400">
                        <strong>Service Worker Issue:</strong> If you see CORS errors, a service worker may be blocking requests. 
                        Open DevTools → Application → Service Workers and click "Unregister", then reload. 
                        Or use <strong>Text Mode</strong> which works offline and doesn't require external scripts.
                    </p>
                </div>
            </div>
        );
    }

    return (
        <div className="w-full h-full flex flex-col">
            <div 
                ref={containerRef} 
                className="flex-1 border rounded bg-white"
                style={{ minHeight: `${height}px` }}
            />
            {!loaded && (
                <div className="absolute inset-0 flex items-center justify-center bg-gray-50 dark:bg-gray-800">
                    <p className="text-sm text-gray-600 dark:text-gray-400">Loading molecular editor...</p>
                </div>
            )}
        </div>
    );
}
