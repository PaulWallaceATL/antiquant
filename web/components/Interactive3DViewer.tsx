'use client';

import React, { useEffect, useRef } from 'react';

interface Interactive3DViewerProps {
    sdf: string;
    title?: string;
    viewerKey?: string; // Force re-render when this changes
}

declare global {
    interface Window {
        $3Dmol: any;
    }
}

export default function Interactive3DViewer({ sdf, title = "Molecule", viewerKey = "default" }: Interactive3DViewerProps) {
    const containerRef = useRef<HTMLDivElement>(null);
    const viewerRef = useRef<any>(null);

    useEffect(() => {
        // Load 3Dmol if not already loaded
        if (typeof window !== 'undefined' && !window.$3Dmol) {
            const script = document.createElement('script');
            script.src = 'https://3Dmol.csb.pitt.edu/build/3Dmol-min.js';
            script.async = true;
            script.onload = () => {
                console.log('3Dmol.js loaded, initializing...');
                setTimeout(initViewer, 100);
            };
            document.body.appendChild(script);
        } else if (window.$3Dmol) {
            setTimeout(initViewer, 100);
        }

        return () => {
            if (viewerRef.current) {
                viewerRef.current = null;
            }
        };
    }, [sdf, viewerKey]); // Re-run when sdf or key changes

    const initViewer = () => {
        if (!containerRef.current || !window.$3Dmol || !sdf) {
            console.warn('Cannot initialize: missing requirements');
            return;
        }

        try {
            // Clear container
            containerRef.current.innerHTML = '';

            // Create viewer
            const viewer = window.$3Dmol.createViewer(containerRef.current, {
                backgroundColor: 'white',
            });

            // Add model from SDF
            viewer.addModel(sdf, 'sdf');

            // Set style
            viewer.setStyle({}, { stick: { radius: 0.15 }, sphere: { scale: 0.25 } });

            // Center and render
            viewer.zoomTo();
            viewer.render();

            // Enable rotation
            viewer.spin(true);

            viewerRef.current = viewer;
            console.log('3D viewer initialized successfully!');
        } catch (error) {
            console.error('3D Viewer initialization error:', error);
        }
    };

    return (
        <div className="relative w-full h-full min-h-[400px] bg-gradient-to-br from-gray-50 to-gray-100 rounded-lg border border-gray-300">
            {title && (
                <div className="absolute top-3 left-3 bg-blue-600 text-white px-3 py-1.5 rounded-md text-sm font-semibold shadow-lg z-10">
                    {title}
                </div>
            )}
            <div
                ref={containerRef}
                className="w-full h-full rounded-lg"
                style={{ minHeight: '400px' }}
            />
            <div className="absolute bottom-3 right-3 text-xs text-gray-500 bg-white/80 px-2 py-1 rounded">
                Drag to rotate • Scroll to zoom
            </div>
        </div>
    );
}
