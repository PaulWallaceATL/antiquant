'use client';

import React from 'react';
import { CheckCircle, XCircle, AlertTriangle, Activity, Shield, Zap } from 'lucide-react';
import { DrugLikeness } from '../lib/molecular_inference';
import InfoTooltip from './InfoTooltip';
import { TOOLTIPS } from '../lib/tooltips';

interface DrugLikenessDashboardProps {
    data: DrugLikeness;
}

const StatusIcon = ({ status }: { status: boolean | string }) => {
    if (status === true || status === 'High' || status === 'Low Risk') {
        return <CheckCircle className="w-5 h-5 text-green-500" />;
    } else if (status === false || status === 'Low' || status === 'High Risk') {
        return <XCircle className="w-5 h-5 text-red-500" />;
    } else {
        return <AlertTriangle className="w-5 h-5 text-yellow-500" />;
    }
};

const ScoreGauge = ({ score }: { score: number }) => {
    const percentage = (score / 10) * 100;
    const color = score >= 7 ? 'text-green-500' : score >= 4 ? 'text-yellow-500' : 'text-red-500';
    const strokeColor = score >= 7 ? '#22c55e' : score >= 4 ? '#eab308' : '#ef4444';

    return (
        <div className="relative w-32 h-32 flex items-center justify-center">
            <svg className="w-full h-full transform -rotate-90">
                <circle
                    cx="64"
                    cy="64"
                    r="56"
                    stroke="#e5e7eb"
                    strokeWidth="12"
                    fill="none"
                />
                <circle
                    cx="64"
                    cy="64"
                    r="56"
                    stroke={strokeColor}
                    strokeWidth="12"
                    fill="none"
                    strokeDasharray={351.86}
                    strokeDashoffset={351.86 - (351.86 * percentage) / 100}
                    className="transition-all duration-1000 ease-out"
                />
            </svg>
            <div className="absolute flex flex-col items-center">
                <span className={`text-3xl font-bold ${color}`}>{score}</span>
                <span className="text-xs text-gray-500 uppercase">Score</span>
            </div>
        </div>
    );
};

export default function DrugLikenessDashboard({ data }: DrugLikenessDashboardProps) {
    return (
        <div className="space-y-6">
            {/* Top Section: Score & Quick Stats */}
            <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
                <div className="bg-white dark:bg-gray-900 p-6 rounded-xl border border-gray-200 dark:border-gray-800 flex flex-col items-center justify-center">
                    <div className="flex items-center gap-2 mb-4">
                        <h3 className="text-sm font-medium text-gray-500">Overall Drug Score</h3>
                        <InfoTooltip title={TOOLTIPS.drugScore.title} content={TOOLTIPS.drugScore.content} />
                    </div>
                    <ScoreGauge score={data.drug_score} />
                </div>

                <div className="col-span-2 grid grid-cols-2 gap-4">
                    <div className="bg-blue-50 dark:bg-blue-900/20 p-4 rounded-xl border border-blue-100 dark:border-blue-800">
                        <div className="flex items-center gap-3 mb-2">
                            <Shield className="w-5 h-5 text-blue-600" />
                            <h3 className="font-semibold text-blue-900 dark:text-blue-100">Lipinski's Rule of 5</h3>
                            <InfoTooltip title={TOOLTIPS.lipinski.title} content={TOOLTIPS.lipinski.content} />
                        </div>
                        <div className="flex items-center gap-2 mb-1">
                            <StatusIcon status={data.lipinski.pass} />
                            <span className="font-medium">{data.lipinski.pass ? 'Passed' : 'Failed'}</span>
                        </div>
                        <p className="text-sm text-gray-600 dark:text-gray-400">
                            {data.lipinski.violations} violation(s) found
                        </p>
                    </div>

                    <div className="bg-purple-50 dark:bg-purple-900/20 p-4 rounded-xl border border-purple-100 dark:border-purple-800">
                        <div className="flex items-center gap-3 mb-2">
                            <Zap className="w-5 h-5 text-purple-600" />
                            <h3 className="font-semibold text-purple-900 dark:text-purple-100">Synthesizability</h3>
                            <InfoTooltip title={TOOLTIPS.syntheticAccessibility.title} content={TOOLTIPS.syntheticAccessibility.content} />
                        </div>
                        <div className="text-2xl font-bold text-purple-700 dark:text-purple-300">
                            {data.synthetic_accessibility}/10
                        </div>
                        <p className="text-sm text-gray-600 dark:text-gray-400">
                            {data.synthetic_accessibility <= 3 ? 'Easy' : data.synthetic_accessibility <= 7 ? 'Moderate' : 'Difficult'} to synthesize
                        </p>
                    </div>
                </div>
            </div>

            {/* Detailed Rules Grid */}
            <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                {/* Lipinski Details */}
                <div className="bg-white dark:bg-gray-900 rounded-xl border border-gray-200 dark:border-gray-800 overflow-hidden">
                    <div className="px-6 py-4 border-b border-gray-200 dark:border-gray-800 bg-gray-50 dark:bg-gray-800/50">
                        <h3 className="font-semibold">Lipinski Rules Breakdown</h3>
                    </div>
                    <div className="p-6 space-y-4">
                        {[
                            { label: 'Molecular Weight', val: data.lipinski.details.molecular_weight, unit: 'Da', tooltip: TOOLTIPS.molecularWeight },
                            { label: 'LogP', val: data.lipinski.details.logp, unit: '', tooltip: TOOLTIPS.logp },
                            { label: 'H-Bond Donors', val: data.lipinski.details.h_donors, unit: '', tooltip: TOOLTIPS.hDonors },
                            { label: 'H-Bond Acceptors', val: data.lipinski.details.h_acceptors, unit: '', tooltip: TOOLTIPS.hAcceptors },
                        ].map((item) => (
                            <div key={item.label} className="flex justify-between items-center">
                                <div className="flex items-center gap-2">
                                    <span className="text-sm text-gray-600 dark:text-gray-400">{item.label}</span>
                                    <InfoTooltip title={item.tooltip.title} content={item.tooltip.content} />
                                </div>
                                <div className="flex items-center gap-3">
                                    <span className="font-mono text-sm">
                                        {item.val.value.toFixed(1)} {item.unit}
                                    </span>
                                    <StatusIcon status={item.val.pass} />
                                </div>
                            </div>
                        ))}
                    </div>
                </div>

                {/* ADMET Prediction */}
                <div className="bg-white dark:bg-gray-900 rounded-xl border border-gray-200 dark:border-gray-800 overflow-hidden">
                    <div className="px-6 py-4 border-b border-gray-200 dark:border-gray-800 bg-gray-50 dark:bg-gray-800/50 flex items-center gap-2">
                        <h3 className="font-semibold">ADMET Profile</h3>
                        <InfoTooltip title={TOOLTIPS.admet.title} content={TOOLTIPS.admet.content} />
                    </div>
                    <div className="p-6 space-y-4">
                        <div className="flex justify-between items-center">
                            <div className="flex items-center gap-2">
                                <span className="text-sm text-gray-600 dark:text-gray-400">BBB Permeability</span>
                                <InfoTooltip title={TOOLTIPS.bbbPermeability.title} content={TOOLTIPS.bbbPermeability.content} />
                            </div>
                            <span className={`text-sm font-medium px-2 py-1 rounded ${data.bbb_permeability === 'High' ? 'bg-green-100 text-green-800' : 'bg-yellow-100 text-yellow-800'
                                }`}>
                                {data.bbb_permeability}
                            </span>
                        </div>
                        <div className="flex justify-between items-center">
                            <div className="flex items-center gap-2">
                                <span className="text-sm text-gray-600 dark:text-gray-400">Hepatotoxicity</span>
                                <InfoTooltip title={TOOLTIPS.hepatotoxicity.title} content={TOOLTIPS.hepatotoxicity.content} />
                            </div>
                            <span className={`text-sm font-medium px-2 py-1 rounded ${data.admet.toxicity.hepatotoxicity === 'Low Risk' ? 'bg-green-100 text-green-800' : 'bg-red-100 text-red-800'
                                }`}>
                                {data.admet.toxicity.hepatotoxicity}
                            </span>
                        </div>
                        <div className="flex justify-between items-center">
                            <div className="flex items-center gap-2">
                                <span className="text-sm text-gray-600 dark:text-gray-400">CYP450 Inhibitor</span>
                                <InfoTooltip title={TOOLTIPS.cyp450.title} content={TOOLTIPS.cyp450.content} />
                            </div>
                            <span className={`text-sm font-medium px-2 py-1 rounded ${data.admet.metabolism.cyp450_inhibitor === 'Unlikely' ? 'bg-green-100 text-green-800' : 'bg-yellow-100 text-yellow-800'
                                }`}>
                                {data.admet.metabolism.cyp450_inhibitor}
                            </span>
                        </div>
                        <div className="flex justify-between items-center">
                            <div className="flex items-center gap-2">
                                <span className="text-sm text-gray-600 dark:text-gray-400">Mutagenicity (Ames)</span>
                                <InfoTooltip title={TOOLTIPS.amesMutagenicity.title} content={TOOLTIPS.amesMutagenicity.content} />
                            </div>
                            <span className={`text-sm font-medium px-2 py-1 rounded ${data.admet.toxicity.ames_mutagenicity === 'Low Risk' ? 'bg-green-100 text-green-800' : 'bg-red-100 text-red-800'
                                }`}>
                                {data.admet.toxicity.ames_mutagenicity}
                            </span>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    );
}
