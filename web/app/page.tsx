'use client';

import React, { useState } from 'react';
import { Play, Pill, Beaker, Layers, Zap, Atom } from 'lucide-react';
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from 'recharts';

// Lightweight UI components
const Button = ({ children, onClick, disabled, className, variant = 'primary' }: any) => (
  <button
    onClick={onClick}
    disabled={disabled}
    className={`px-4 py-2 rounded-md font-medium transition-colors flex items-center gap-2 ${disabled ? 'opacity-50 cursor-not-allowed' : ''
      } ${variant === 'primary'
        ? 'bg-blue-600 text-white hover:bg-blue-700'
        : 'bg-gray-100 text-gray-900 hover:bg-gray-200 dark:bg-gray-800 dark:text-gray-100'
      } ${className}`}
  >
    {children}
  </button>
);

const Card = ({ children, className }: any) => (
  <div className={`bg-white dark:bg-gray-900 border border-gray-200 dark:border-gray-800 rounded-lg shadow-sm ${className}`}>
    {children}
  </div>
);

const Badge = ({ children, variant = 'default' }: any) => {
  const colors = {
    default: 'bg-gray-100 text-gray-800 dark:bg-gray-700 dark:text-gray-200',
    success: 'bg-green-100 text-green-800 dark:bg-green-900 dark:text-green-200',
    warning: 'bg-yellow-100 text-yellow-800 dark:bg-yellow-900 dark:text-yellow-200',
    info: 'bg-blue-100 text-blue-800 dark:bg-blue-900 dark:text-blue-200',
  };
  return (
    <span className={`px-2.5 py-0.5 rounded-full text-xs font-medium ${colors[variant as keyof typeof colors]}`}>
      {children}
    </span>
  );
};

export default function Home() {
  const [smiles, setSmiles] = useState<string>('CCO');  // Ethanol
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState<any>(null);
  const [mode, setMode] = useState<'classical' | 'quantum'>('classical');

  const analyze = async () => {
    setLoading(true);
    try {
      const endpoint = mode === 'classical' ? '/api/analyze' : '/api/analyze-quantum';
      const res = await fetch(endpoint, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ smiles }),
      });
      const data = await res.json();
      setResult(data);
    } catch (error) {
      console.error(error);
      alert('Analysis failed');
    } finally {
      setLoading(false);
    }
  };

  const exampleMolecules = [
    { name: 'Ethanol', smiles: 'CCO' },
    { name: 'Aspirin', smiles: 'CC(=O)Oc1ccccc1C(=O)O' },
    { name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' },
    { name: 'Benzene', smiles: 'c1ccccc1' },
  ];

  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-50 via-white to-purple-50 dark:from-gray-950 dark:via-gray-900 dark:to-blue-950 text-gray-900 dark:text-gray-100 p-8 font-sans">
      <header className="mb-8 flex justify-between items-start">
        <div>
          <h1 className="text-4xl font-bold tracking-tight flex items-center gap-3 mb-2">
            <Pill className="w-10 h-10 text-blue-600 dark:text-blue-400" />
            MoleculeAI
          </h1>
          <p className="text-gray-600 dark:text-gray-400 text-lg">
            Hybrid Classical-Quantum Molecular Property Predictor
          </p>
          <p className="text-sm text-gray-500 dark:text-gray-500 mt-1">
            Predicting Log Solubility (mol/L) using XGBoost + Quantum VQC
          </p>
        </div>
        <div className="flex gap-2">
          <Button
            variant={mode === 'classical' ? 'primary' : 'secondary'}
            onClick={() => setMode('classical')}
          >
            <Layers className="w-4 h-4" /> Classical
          </Button>
          <Button
            variant={mode === 'quantum' ? 'primary' : 'secondary'}
            onClick={() => setMode('quantum')}
            className={mode === 'quantum' ? 'bg-purple-600 hover:bg-purple-700 text-white' : ''}
          >
            <Atom className="w-4 h-4" /> Quantum
          </Button>
        </div>
      </header>

      <main className="grid grid-cols-1 lg:grid-cols-2 gap-8">
        {/* Left Column: Input */}
        <div className="space-y-4">
          <Card className="overflow-hidden border-2 border-blue-100 dark:border-blue-900/50">
            <div className="bg-gradient-to-r from-blue-500 to-purple-600 px-4 py-3 flex justify-between items-center">
              <span className="text-sm font-semibold flex items-center gap-2 text-white">
                <Beaker className="w-4 h-4" /> Enter SMILES String
              </span>
              <Badge variant="info">SMILES Format</Badge>
            </div>
            <div className="p-4">
              <textarea
                value={smiles}
                onChange={(e) => setSmiles(e.target.value)}
                className="w-full h-32 bg-gray-50 dark:bg-gray-800 border border-gray-300 dark:border-gray-700 rounded-lg p-4 font-mono text-lg focus:ring-2 focus:ring-blue-500 focus:border-transparent outline-none"
                placeholder="Enter SMILES notation (e.g., CCO for ethanol)"
              />
            </div>
          </Card>

          <Card className="p-4">
            <h3 className="text-sm font-semibold mb-3 text-gray-700 dark:text-gray-300">Example Molecules</h3>
            <div className="grid grid-cols-2 gap-2">
              {exampleMolecules.map((mol) => (
                <button
                  key={mol.smiles}
                  onClick={() => setSmiles(mol.smiles)}
                  className="px-3 py-2 bg-gray-100 dark:bg-gray-800 hover:bg-blue-100 dark:hover:bg-blue-900/30 rounded-md text-sm text-left transition-colors"
                >
                  <div className="font-medium">{mol.name}</div>
                  <div className="text-xs text-gray-500 dark:text-gray-400 font-mono">{mol.smiles}</div>
                </button>
              ))}
            </div>
          </Card>

          <Button
            onClick={analyze}
            disabled={loading || !smiles}
            className={`w-full py-4 text-lg shadow-lg ${mode === 'quantum' ? 'bg-purple-600 hover:bg-purple-700' : ''}`}
          >
            {loading ? 'Analyzing...' : (
              <>
                <Play className="w-5 h-5" /> Predict Solubility
              </>
            )}
          </Button>
        </div>

        {/* Right Column: Results */}
        <div className="space-y-6">
          {result ? (
            result.error ? (
              <Card className="p-6 border-2 border-red-200 dark:border-red-900">
                <div className="text-red-600 dark:text-red-400 text-center">
                  <p className="font-semibold">Error</p>
                  <p className="text-sm mt-2">{result.error}</p>
                </div>
              </Card>
            ) : (
              <>
                {/* Main Prediction Card */}
                <Card className="p-6 bg-gradient-to-br from-blue-50 to-purple-50 dark:from-blue-950/50 dark:to-purple-950/50 border-2 border-blue-200 dark:border-blue-800">
                  <h3 className="text-sm font-medium text-gray-600 dark:text-gray-400 uppercase tracking-wider mb-4">Prediction Result</h3>
                  <div className="flex items-end gap-4 mb-6">
                    <div className="text-6xl font-bold text-blue-600 dark:text-blue-400">
                      {result.prediction}
                    </div>
                    <div className="text-lg text-gray-600 dark:text-gray-400 mb-2">
                      Log(mol/L)
                    </div>
                  </div>

                  <div className="space-y-2 text-sm">
                    <div className="flex justify-between">
                      <span className="text-gray-600 dark:text-gray-400">Property:</span>
                      <span className="font-medium">{result.property}</span>
                    </div>
                    <div className="flex justify-between">
                      <span className="text-gray-600 dark:text-gray-400">Model:</span>
                      <Badge variant={mode === 'quantum' ? 'info' : 'success'}>{result.model_used}</Badge>
                    </div>
                    <div className="flex justify-between">
                      <span className="text-gray-600 dark:text-gray-400">SMILES:</span>
                      <span className="font-mono text-xs">{result.smiles}</span>
                    </div>
                  </div>
                </Card>

                {/* Molecular Features Table */}
                <Card className="overflow-hidden">
                  <div className="px-6 py-4 bg-gray-50 dark:bg-gray-800/50 border-b border-gray-200 dark:border-gray-700">
                    <h3 className="font-semibold flex items-center gap-2">
                      <Beaker className="w-4 h-4" /> Molecular Descriptors
                    </h3>
                  </div>
                  <div className="p-0">
                    <table className="w-full text-sm text-left">
                      <thead className="bg-gray-50 dark:bg-gray-800/30 text-gray-600 dark:text-gray-400">
                        <tr>
                          <th className="px-6 py-3 font-medium">Property</th>
                          <th className="px-6 py-3 font-medium text-right">Value</th>
                        </tr>
                      </thead>
                      <tbody className="divide-y divide-gray-200 dark:divide-gray-700">
                        <tr>
                          <td className="px-6 py-3">Molecular Weight</td>
                          <td className="px-6 py-3 text-right font-mono">{result.features.molecular_weight.toFixed(2)}</td>
                        </tr>
                        <tr>
                          <td className="px-6 py-3">LogP</td>
                          <td className="px-6 py-3 text-right font-mono">{result.features.logp.toFixed(2)}</td>
                        </tr>
                        <tr>
                          <td className="px-6 py-3">H-Bond Donors</td>
                          <td className="px-6 py-3 text-right font-mono">{result.features.num_h_donors}</td>
                        </tr>
                        <tr>
                          <td className="px-6 py-3">H-Bond Acceptors</td>
                          <td className="px-6 py-3 text-right font-mono">{result.features.num_h_acceptors}</td>
                        </tr>
                        <tr>
                          <td className="px-6 py-3">Rotatable Bonds</td>
                          <td className="px-6 py-3 text-right font-mono">{result.features.num_rotatable_bonds}</td>
                        </tr>
                        <tr>
                          <td className="px-6 py-3">Aromatic Rings</td>
                          <td className="px-6 py-3 text-right font-mono">{result.features.num_aromatic_rings}</td>
                        </tr>
                        <tr>
                          <td className="px-6 py-3">TPSA</td>
                          <td className="px-6 py-3 text-right font-mono">{result.features.tpsa.toFixed(2)}</td>
                        </tr>
                        <tr>
                          <td className="px-6 py-3">Total Atoms</td>
                          <td className="px-6 py-3 text-right font-mono">{result.features.num_atoms}</td>
                        </tr>
                      </tbody>
                    </table>
                  </div>
                </Card>

                {/* Chart */}
                <Card className="p-6">
                  <h3 className="font-semibold mb-6">Descriptor Visualization</h3>
                  <div className="h-64 w-full">
                    <ResponsiveContainer width="100%" height="100%">
                      <BarChart data={[
                        { name: 'MW/10', value: result.features.molecular_weight / 10 },
                        { name: 'LogP', value: result.features.logp },
                        { name: 'HBD', value: result.features.num_h_donors },
                        { name: 'HBA', value: result.features.num_h_acceptors },
                        { name: 'RotB', value: result.features.num_rotatable_bonds },
                        { name: 'AromR', value: result.features.num_aromatic_rings },
                      ]}>
                        <CartesianGrid strokeDasharray="3 3" stroke="#333" vertical={false} />
                        <XAxis dataKey="name" stroke="#888" />
                        <YAxis stroke="#888" />
                        <Tooltip
                          contentStyle={{ backgroundColor: '#1f2937', border: 'none', borderRadius: '8px' }}
                          itemStyle={{ color: '#fff' }}
                        />
                        <Bar dataKey="value" fill={mode === 'quantum' ? '#9333ea' : '#2563eb'} radius={[4, 4, 0, 0]} />
                      </BarChart>
                    </ResponsiveContainer>
                  </div>
                </Card>
              </>
            )
          ) : (
            <div className="h-full flex flex-col items-center justify-center text-gray-400 border-2 border-dashed border-gray-300 dark:border-gray-700 rounded-lg p-12 min-h-[500px]">
              <Zap className="w-16 h-16 mb-4 opacity-20" />
              <p className="text-lg">Enter a SMILES string and run analysis</p>
              <p className="text-sm mt-2">Try examples like CCO (ethanol) or c1ccccc1 (benzene)</p>
            </div>
          )}
        </div>
      </main>
    </div>
  );
}
