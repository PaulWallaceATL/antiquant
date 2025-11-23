'use client';

import React, { useState } from 'react';
import { Download, Loader2 } from 'lucide-react';
import { MolecularAnalysisResult } from '../lib/molecular_inference';
import { downloadReport } from '../lib/report_generator';

interface ExportButtonProps {
    analysis: MolecularAnalysisResult;
    className?: string;
}

export default function ExportButton({ analysis, className = '' }: ExportButtonProps) {
    const [isGenerating, setIsGenerating] = useState(false);

    const handleExport = async () => {
        if (!analysis || analysis.error) {
            alert('No analysis data available to export');
            return;
        }

        setIsGenerating(true);
        try {
            await downloadReport(analysis);
        } catch (error) {
            console.error('Export failed:', error);
            alert('Failed to generate PDF report. Please try again.');
        } finally {
            setIsGenerating(false);
        }
    };

    return (
        <button
            onClick={handleExport}
            disabled={isGenerating || !analysis || !!analysis.error}
            className={`px-4 py-2 rounded-md font-medium transition-all flex items-center gap-2 ${
                isGenerating || !analysis || !!analysis.error
                    ? 'opacity-50 cursor-not-allowed bg-gray-300 text-gray-600'
                    : 'bg-blue-600 text-white hover:bg-blue-700 shadow-lg shadow-blue-500/50 hover:scale-105 active:scale-95'
            } ${className}`}
        >
            {isGenerating ? (
                <>
                    <Loader2 className="w-4 h-4 animate-spin" />
                    Generating...
                </>
            ) : (
                <>
                    <Download className="w-4 h-4" />
                    Export PDF
                </>
            )}
        </button>
    );
}

