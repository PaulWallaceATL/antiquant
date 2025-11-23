'use client';

import React from 'react';

interface Molecule3DViewerProps {
    structure?: {
        sdf: string;
        image_base64?: string;
        atoms?: any[];
        bonds?: any[];
    };
    title?: string;
}

export default function Molecule3DViewer({ structure, title = "Molecule" }: Molecule3DViewerProps) {
    if (!structure) {
        return <div className="w-full h-full min-h-[300px] flex items-center justify-center bg-gray-50 rounded border">
            <span className="text-gray-500">No structure data available</span>
        </div>;
    }

    return (
        <div className="relative w-full h-full min-h-[300px] bg-white rounded border border-gray-200 flex items-center justify-center p-4">
            {title && (
                <div className="absolute top-2 left-2 bg-blue-600 text-white px-3 py-1 rounded text-xs font-semibold z-10">
                    {title}
                </div>
            )}

            {structure.image_base64 ? (
                <div className="flex flex-col items-center">
                    <img
                        src={`data:image/png;base64,${structure.image_base64}`}
                        alt="Molecular Structure"
                        className="max-w-full h-auto border border-gray-200 rounded shadow-sm"
                        style={{ maxHeight: '400px' }}
                    />
                    <div className="mt-2 text-xs text-gray-500">
                        {structure.atoms ? `${structure.atoms.length} atoms` : ''}
                        {structure.bonds ? `, ${structure.bonds.length} bonds` : ''}
                    </div>
                </div>
            ) : (
                <div className="text-center text-gray-500">
                    <p className="mb-2">No visualization available</p>
                    <p className="text-xs">{structure.atoms?.length || 0} atoms detected</p>
                </div>
            )}
        </div>
    );
}
