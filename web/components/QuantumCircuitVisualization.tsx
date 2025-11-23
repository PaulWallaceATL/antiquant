'use client';

import React from 'react';
import { Zap } from 'lucide-react';
import InfoTooltip from './InfoTooltip';
import { TOOLTIPS } from '../lib/tooltips';

interface QuantumCircuitProps {
    circuit: {
        n_qubits: number;
        n_layers: number;
        operations: Array<{
            gate: string;
            qubit?: number;
            parameter?: number;
            parameters?: number[];
            control?: number;
            target?: number;
            layer: string;
        }>;
        expectation_value: number;
        final_state_probabilities: number[];
        input_encoding: number[];
    };
}

export default function QuantumCircuitVisualization({ circuit }: QuantumCircuitProps) {
    const { n_qubits, operations, expectation_value, final_state_probabilities, input_encoding } = circuit;

    // Group operations by layer
    const layers: Map<string, typeof operations> = new Map();
    operations.forEach(op => {
        if (!layers.has(op.layer)) {
            layers.set(op.layer, []);
        }
        layers.get(op.layer)!.push(op);
    });

    const layerNames = Array.from(layers.keys());

    // Calculate grid dimensions
    const wireSpacing = 60;
    const gateWidth = 80;
    const gateSpacing = 40;
    const leftMargin = 100;
    const topMargin = 50;

    const svgWidth = leftMargin + (layerNames.length * (gateWidth + gateSpacing)) + 100;
    const svgHeight = topMargin + (n_qubits * wireSpacing) + 100;

    return (
        <div className="space-y-4">
            {/* Circuit SVG */}
            <div className="bg-white dark:bg-gray-900 border-2 border-purple-200 dark:border-purple-800 rounded-lg p-6 overflow-x-auto">
                <div className="flex items-center gap-2 mb-4">
                    <Zap className="w-5 h-5 text-purple-600" />
                    <h3 className="font-semibold text-lg">Quantum Circuit Execution</h3>
                </div>

                <svg width={svgWidth} height={svgHeight} className="mx-auto">
                    {/* Draw qubit wires */}
                    {Array.from({ length: n_qubits }).map((_, i) => (
                        <g key={`qubit-${i}`}>
                            {/* Wire line */}
                            <line
                                x1={leftMargin}
                                y1={topMargin + i * wireSpacing}
                                x2={svgWidth - 50}
                                y2={topMargin + i * wireSpacing}
                                stroke="#666"
                                strokeWidth="2"
                            />
                            {/* Qubit label */}
                            <text
                                x={50}
                                y={topMargin + i * wireSpacing + 5}
                                fill="#888"
                                fontSize="14"
                                fontFamily="monospace"
                            >
                                q[{i}]
                            </text>
                            {/* Initial state */}
                            <circle
                                cx={leftMargin - 20}
                                cy={topMargin + i * wireSpacing}
                                r="8"
                                fill="#3b82f6"
                                stroke="#1e40af"
                                strokeWidth="2"
                            />
                        </g>
                    ))}

                    {/* Draw gates */}
                    {layerNames.map((layerName, layerIdx) => {
                        const layerOps = layers.get(layerName)!;
                        const x = leftMargin + layerIdx * (gateWidth + gateSpacing);

                        return (
                            <g key={`layer-${layerName}`}>
                                {layerOps.map((op, opIdx) => {
                                    if (op.gate === 'RX' && op.qubit !== undefined) {
                                        const y = topMargin + op.qubit * wireSpacing;
                                        return (
                                            <g key={`${layerName}-${opIdx}`}>
                                                <rect
                                                    x={x}
                                                    y={y - 20}
                                                    width={60}
                                                    height={40}
                                                    fill="#3b82f6"
                                                    stroke="#1e40af"
                                                    strokeWidth="2"
                                                    rx="4"
                                                />
                                                <text
                                                    x={x + 30}
                                                    y={y + 5}
                                                    fill="white"
                                                    fontSize="12"
                                                    fontWeight="bold"
                                                    textAnchor="middle"
                                                >
                                                    RX
                                                </text>
                                                <text
                                                    x={x + 30}
                                                    y={y - 25}
                                                    fill="#666"
                                                    fontSize="9"
                                                    textAnchor="middle"
                                                >
                                                    θ={op.parameter?.toFixed(2)}
                                                </text>
                                            </g>
                                        );
                                    } else if (op.gate === 'Rot' && op.qubit !== undefined) {
                                        const y = topMargin + op.qubit * wireSpacing;
                                        return (
                                            <g key={`${layerName}-${opIdx}`}>
                                                <rect
                                                    x={x}
                                                    y={y - 20}
                                                    width={60}
                                                    height={40}
                                                    fill="#8b5cf6"
                                                    stroke="#6d28d9"
                                                    strokeWidth="2"
                                                    rx="4"
                                                />
                                                <text
                                                    x={x + 30}
                                                    y={y + 5}
                                                    fill="white"
                                                    fontSize="11"
                                                    fontWeight="bold"
                                                    textAnchor="middle"
                                                >
                                                    Rot
                                                </text>
                                            </g>
                                        );
                                    } else if (op.gate === 'CNOT' && op.control !== undefined && op.target !== undefined) {
                                        const y1 = topMargin + op.control * wireSpacing;
                                        const y2 = topMargin + op.target * wireSpacing;
                                        return (
                                            <g key={`${layerName}-${opIdx}`}>
                                                {/* Control qubit (dot) */}
                                                <circle
                                                    cx={x + 30}
                                                    cy={y1}
                                                    r="6"
                                                    fill="#ec4899"
                                                />
                                                {/* Vertical line */}
                                                <line
                                                    x1={x + 30}
                                                    y1={y1}
                                                    x2={x + 30}
                                                    y2={y2}
                                                    stroke="#ec4899"
                                                    strokeWidth="2"
                                                />
                                                {/* Target qubit (⊕) */}
                                                <circle
                                                    cx={x + 30}
                                                    cy={y2}
                                                    r="12"
                                                    fill="none"
                                                    stroke="#ec4899"
                                                    strokeWidth="2"
                                                />
                                                <line
                                                    x1={x + 18}
                                                    y1={y2}
                                                    x2={x + 42}
                                                    y2={y2}
                                                    stroke="#ec4899"
                                                    strokeWidth="2"
                                                />
                                                <line
                                                    x1={x + 30}
                                                    y1={y2 - 12}
                                                    x2={x + 30}
                                                    y2={y2 + 12}
                                                    stroke="#ec4899"
                                                    strokeWidth="2"
                                                />
                                            </g>
                                        );
                                    }
                                    return null;
                                })}
                            </g>
                        );
                    })}

                    {/* Measurement at the end */}
                    {Array.from({ length: n_qubits }).map((_, i) => {
                        const x = svgWidth - 80;
                        const y = topMargin + i * wireSpacing;
                        return (
                            <g key={`measure-${i}`}>
                                <rect
                                    x={x}
                                    y={y - 20}
                                    width={50}
                                    height={40}
                                    fill="#10b981"
                                    stroke="#059669"
                                    strokeWidth="2"
                                    rx="4"
                                />
                                <text
                                    x={x + 25}
                                    y={y + 5}
                                    fill="white"
                                    fontSize="10"
                                    fontWeight="bold"
                                    textAnchor="middle"
                                >
                                    M
                                </text>
                            </g>
                        );
                    })}
                </svg>
            </div>

            {/* Legend */}
            <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-sm">
                <div className="flex items-center gap-2">
                    <div className="w-8 h-8 bg-blue-500 rounded border-2 border-blue-700"></div>
                    <span>RX Gate (Encoding)</span>
                    <InfoTooltip title={TOOLTIPS.rxGate.title} content={TOOLTIPS.rxGate.content} />
                </div>
                <div className="flex items-center gap-2">
                    <div className="w-8 h-8 bg-purple-500 rounded border-2 border-purple-700"></div>
                    <span>Rotation Gates</span>
                    <InfoTooltip title={TOOLTIPS.rotGate.title} content={TOOLTIPS.rotGate.content} />
                </div>
                <div className="flex items-center gap-2">
                    <div className="w-8 h-8 bg-pink-500 rounded-full border-2 border-pink-700"></div>
                    <span>CNOT (Entangle)</span>
                    <InfoTooltip title={TOOLTIPS.cnotGate.title} content={TOOLTIPS.cnotGate.content} />
                </div>
                <div className="flex items-center gap-2">
                    <div className="w-8 h-8 bg-green-500 rounded border-2 border-green-700"></div>
                    <span>Measurement</span>
                    <InfoTooltip title={TOOLTIPS.measurement.title} content={TOOLTIPS.measurement.content} />
                </div>
            </div>

            {/* Quantum State Info */}
            <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                <div className="bg-gradient-to-br from-purple-50 to-pink-50 dark:from-purple-950/30 dark:to-pink-950/30 p-4 rounded-lg border border-purple-200 dark:border-purple-800">
                    <div className="flex items-center gap-2 mb-2">
                        <h4 className="font-semibold text-sm">Expectation Value</h4>
                        <InfoTooltip title={TOOLTIPS.expectationValue.title} content={TOOLTIPS.expectationValue.content} />
                    </div>
                    <div className="text-3xl font-bold text-purple-600 dark:text-purple-400">
                        {expectation_value.toFixed(4)}
                    </div>
                    <p className="text-xs text-gray-600 dark:text-gray-400 mt-1">⟨Z⟩ measurement on qubit 0</p>
                </div>

                <div className="bg-gradient-to-br from-blue-50 to-purple-50 dark:from-blue-950/30 dark:to-purple-950/30 p-4 rounded-lg border border-blue-200 dark:border-blue-800">
                    <div className="flex items-center gap-2 mb-2">
                        <h4 className="font-semibold text-sm">Input Encoding</h4>
                        <InfoTooltip title={TOOLTIPS.inputEncoding.title} content={TOOLTIPS.inputEncoding.content} />
                    </div>
                    <div className="grid grid-cols-4 gap-1 text-xs font-mono">
                        {input_encoding.map((val, idx) => (
                            <div key={idx} className="bg-white dark:bg-gray-800 px-2 py-1 rounded text-center">
                                {val.toFixed(2)}
                            </div>
                        ))}
                    </div>
                    <p className="text-xs text-gray-600 dark:text-gray-400 mt-2">PCA-reduced embedding (4D)</p>
                </div>
            </div>

            {/* State Probabilities */}
            <div className="bg-white dark:bg-gray-900 border border-gray-200 dark:border-gray-700 rounded-lg p-4">
                <h4 className="font-semibold mb-3 text-sm">Final Quantum State Probabilities</h4>
                <div className="grid grid-cols-8 gap-2">
                    {final_state_probabilities.map((prob, idx) => (
                        <div key={idx} className="text-center">
                            <div
                                className="bg-gradient-to-t from-purple-500 to-pink-500 rounded-t"
                                style={{ height: `${Math.max(prob * 200, 4)}px` }}
                            ></div>
                            <div className="text-xs text-gray-600 dark:text-gray-400 mt-1">
                                |{idx.toString(2).padStart(4, '0')}⟩
                            </div>
                            <div className="text-xs font-mono">{(prob * 100).toFixed(1)}%</div>
                        </div>
                    ))}
                </div>
            </div>
        </div>
    );
}
