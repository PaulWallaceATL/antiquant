'use client';

import React, { useState } from 'react';
import { Play, Pill, Beaker, Layers, Atom, Activity, GitCompare, Shield, MessageSquare, TrendingUp, Download, X } from 'lucide-react';
import { BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, RadarChart, PolarGrid, PolarAngleAxis, PolarRadiusAxis, Radar } from 'recharts';
import QuantumCircuitVisualization from '../components/QuantumCircuitVisualization';
import DrugLikenessDashboard from '../components/DrugLikenessDashboard';
import Interactive3DViewer from '../components/Interactive3DViewer';
import MolecularEditor from '../components/MolecularEditor';
import ChatInterface from '../components/ChatInterface';
import MarketAnalysisDashboard from '../components/MarketAnalysisDashboard';
import ExportButton from '../components/ExportButton';
import InfoTooltip from '../components/InfoTooltip';
import { TOOLTIPS } from '../lib/tooltips';

// UI components - Mobile Optimized
const Button = ({ children, onClick, disabled, className, variant = 'primary' }: any) => (
  <button
    onClick={onClick}
    disabled={disabled}
    className={`px-3 sm:px-4 py-2.5 sm:py-2 rounded-md font-medium transition-all flex items-center justify-center gap-2 min-h-[44px] touch-manipulation ${disabled ? 'opacity-50 cursor-not-allowed' : 'hover:scale-105 active:scale-95'
      } ${variant === 'primary'
        ? 'bg-blue-600 text-white hover:bg-blue-700 shadow-lg shadow-blue-500/50'
        : variant === 'compare'
          ? 'bg-gradient-to-r from-blue-600 to-purple-600 text-white hover:from-blue-700 hover:to-purple-700 shadow-lg shadow-purple-500/50'
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
    info: 'bg-blue-100 text-blue-800 dark:bg-blue-900 dark:text-blue-200',
    purple: 'bg-purple-100 text-purple-800 dark:bg-purple-900 dark:text-purple-200',
  };
  return (
    <span className={`px-2.5 py-0.5 rounded-full text-xs font-medium ${colors[variant as keyof typeof colors]}`}>
      {children}
    </span>
  );
};

export default function Home() {
  const [smiles, setSmiles] = useState<string>('CCO');
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState<any>(null);
  const [mode, setMode] = useState<'classical' | 'quantum' | 'compare'>('quantum');


  const [activeTab, setActiveTab] = useState<'druglikeness' | '3d' | 'overview' | 'circuit' | 'details' | 'editor' | 'chat' | 'market'>('druglikeness');
  const [showChat, setShowChat] = useState(false);

  const analyze = async () => {
    setLoading(true);
    setResult(null);
    try {
      const endpoint = mode === 'compare' ? '/api/compare' : mode === 'quantum' ? '/api/analyze-quantum' : '/api/analyze';
      const res = await fetch(endpoint, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ smiles }),
      });
      const data = await res.json();
      console.log("API Response:", data);
      console.log("structure_3d exists?", !!data.structure_3d);
      console.log("structure_3d content:", data.structure_3d);
      setResult(data);
      setActiveTab('druglikeness');
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

  // Helper to get result data based on mode
  const getDisplayData = () => {
    if (!result || result.error) return null;

    if (mode === 'compare') {
      return { classical: result.classical, quantum: result.quantum, features: result.features };
    } else {
      return { single: result, features: result.features };
    }
  };

  const displayData = getDisplayData();

  const handleChatAnalyze = (smilesToAnalyze: string) => {
    // Strict validation - only accept valid SMILES strings
    const trimmed = smilesToAnalyze.trim();
    if (!trimmed || trimmed.length < 2) {
      return;
    }
    
    // Must be valid SMILES pattern
    const smilesPattern = /^[A-Za-z0-9@+\-\[\]()=#\.]+$/;
    if (!smilesPattern.test(trimmed)) {
      console.warn('Invalid SMILES format:', trimmed);
      return;
    }
    
    // Must start with a letter (element symbol)
    if (!/[A-Za-z]/.test(trimmed[0])) {
      console.warn('SMILES must start with a letter:', trimmed);
      return;
    }
    
    // Only update if it passes all validation
    setSmiles(trimmed);
    setShowChat(false); // Close chat when analyzing
    setActiveTab('druglikeness');
    // Trigger analysis
    setTimeout(() => {
      analyze();
    }, 100);
  };

  const [mobileMenuOpen, setMobileMenuOpen] = useState(false);

  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-50 via-white to-purple-50 dark:from-gray-950 dark:via-gray-900 dark:to-blue-950 text-gray-900 dark:text-gray-100 overflow-x-hidden">
      {/* Header - Mobile Responsive */}
      <header className="border-b border-gray-200 dark:border-gray-800 bg-white/95 dark:bg-gray-900/95 backdrop-blur-lg sticky top-0 z-40 shadow-sm">
        <div className="px-4 sm:px-6 py-3 sm:py-4">
          <div className="flex justify-between items-center">
            <div className="flex items-center gap-2 sm:gap-3">
              <Pill className="w-7 h-7 sm:w-8 sm:h-8 text-blue-600 dark:text-blue-400 flex-shrink-0" />
              <div className="min-w-0">
                <h1 className="text-xl sm:text-2xl font-bold tracking-tight truncate">MoleculeAI</h1>
                <p className="text-[10px] sm:text-xs text-gray-600 dark:text-gray-400 truncate">Quantum-Enhanced Drug Discovery</p>
              </div>
            </div>
            {/* Mobile Menu Toggle */}
            <button
              onClick={() => setMobileMenuOpen(!mobileMenuOpen)}
              className="lg:hidden p-2 rounded-lg hover:bg-gray-100 dark:hover:bg-gray-800 transition-colors"
              aria-label="Toggle menu"
            >
              {mobileMenuOpen ? (
                <X className="w-6 h-6" />
              ) : (
                <Activity className="w-6 h-6" />
              )}
            </button>
            {/* Desktop Buttons */}
            <div className="hidden lg:flex gap-2 flex-wrap">
              <Button
                variant={mode === 'quantum' ? 'primary' : 'secondary'}
                onClick={() => setMode('quantum')}
                className={`text-xs sm:text-sm ${mode === 'quantum' ? 'bg-purple-600 hover:bg-purple-700 text-white shadow-lg shadow-purple-500/50' : ''}`}
              >
                <Atom className="w-4 h-4" /> <span className="hidden sm:inline">Quantum</span>
              </Button>
              {result && !result.error && (
                <ExportButton analysis={result} />
              )}
              <Button
                variant={showChat ? 'primary' : 'secondary'}
                onClick={() => setShowChat(!showChat)}
                className="text-xs sm:text-sm"
              >
                <MessageSquare className="w-4 h-4" /> <span className="hidden xl:inline">Chat</span>
              </Button>
            </div>
          </div>
          
          {/* Mobile Menu */}
          {mobileMenuOpen && (
            <div className="lg:hidden mt-4 pt-4 border-t border-gray-200 dark:border-gray-800 space-y-2">
              <div className="grid grid-cols-3 gap-2">
                <Button
                  variant={mode === 'quantum' ? 'primary' : 'secondary'}
                  onClick={() => { setMode('quantum'); setMobileMenuOpen(false); }}
                  className={`text-xs py-2 ${mode === 'quantum' ? 'bg-purple-600 hover:bg-purple-700 text-white' : ''}`}
                >
                  <Atom className="w-4 h-4" /> Quantum
                </Button>
                {result && !result.error && (
                  <div className="col-span-1">
                    <ExportButton analysis={result} />
                  </div>
                )}
                <Button
                  variant={showChat ? 'primary' : 'secondary'}
                  onClick={() => { setShowChat(!showChat); setMobileMenuOpen(false); }}
                  className="text-xs py-2"
                >
                  <MessageSquare className="w-4 h-4" /> Chat
                </Button>
              </div>
            </div>
          )}
        </div>
      </header>

      {/* Main Content - Mobile Responsive */}
      <div className="min-h-[calc(100vh-73px)] lg:h-[calc(100vh-73px)] flex flex-col lg:grid lg:grid-cols-5 gap-2 sm:gap-4 p-2 sm:p-4">
        {/* Left Panel - Input - Mobile Optimized */}
        <div className="lg:col-span-1 flex flex-col lg:block lg:space-y-3 sm:lg:space-y-4 lg:overflow-y-auto lg:max-h-none">
          {/* Mobile: Sticky Input Section */}
          <div className="lg:static sticky top-[73px] z-30 bg-white dark:bg-gray-900 lg:bg-transparent lg:dark:bg-transparent pb-2 lg:pb-0 border-b border-gray-200 dark:border-gray-800 lg:border-0 mb-2 lg:mb-0 space-y-2 lg:space-y-3 sm:lg:space-y-4">
            <Card className="border-2 border-blue-100 dark:border-blue-900/50 lg:shadow-sm">
              <div className="bg-gradient-to-r from-blue-500 to-purple-600 px-2 sm:px-3 py-1.5 sm:py-2">
                <span className="text-[10px] sm:text-xs font-semibold flex items-center gap-1.5 sm:gap-2 text-white">
                  <Beaker className="w-3.5 h-3.5 sm:w-4 sm:h-4" /> SMILES Input
                  <InfoTooltip title={TOOLTIPS.smiles.title} content={TOOLTIPS.smiles.content} className="text-white/80 hover:text-white" />
                </span>
              </div>
              <div className="p-2 sm:p-3">
                <textarea
                  value={smiles}
                  onChange={(e) => setSmiles(e.target.value)}
                  className="w-full h-20 sm:h-24 lg:h-32 bg-gray-50 dark:bg-gray-800 border border-gray-300 dark:border-gray-700 rounded-lg p-2 sm:p-3 font-mono text-[11px] sm:text-xs lg:text-sm focus:ring-2 focus:ring-blue-500 outline-none resize-none touch-manipulation"
                  placeholder="Enter SMILES..."
                  aria-label="SMILES input"
                />
              </div>
            </Card>

            {/* Mobile: Horizontal Scrollable Examples */}
            <Card className="p-2 sm:p-3 lg:shadow-sm">
              <h3 className="text-[10px] sm:text-xs font-semibold mb-1.5 sm:mb-2 text-gray-700 dark:text-gray-300">Examples</h3>
              {/* Mobile: Horizontal scroll, Desktop: Vertical stack */}
              <div className="flex lg:flex-col gap-2 lg:space-y-2 overflow-x-auto lg:overflow-x-visible pb-1 lg:pb-0 -mx-2 lg:mx-0 px-2 lg:px-0 snap-x lg:snap-none scrollbar-hide">
                {exampleMolecules.map((mol) => (
                  <button
                    key={mol.smiles}
                    onClick={() => setSmiles(mol.smiles)}
                    className="flex-shrink-0 lg:w-full px-2.5 sm:px-3 py-2 sm:py-2.5 lg:py-2 bg-gray-100 dark:bg-gray-800 hover:bg-blue-100 dark:hover:bg-blue-900/30 active:bg-blue-200 dark:active:bg-blue-900/50 rounded text-[10px] sm:text-xs text-left transition-colors touch-manipulation min-h-[40px] sm:min-h-[44px] snap-start lg:snap-none"
                    aria-label={`Load ${mol.name} example`}
                  >
                    <div className="font-medium truncate">{mol.name}</div>
                    <div className="text-[9px] sm:text-[10px] text-gray-500 dark:text-gray-400 font-mono truncate max-w-[140px] sm:max-w-none">{mol.smiles}</div>
                  </button>
                ))}
              </div>
            </Card>

            <Button
              onClick={analyze}
              disabled={loading || !smiles}
              className={`w-full py-2.5 sm:py-3 lg:py-3 text-xs sm:text-sm min-h-[44px] sm:min-h-[48px] ${mode === 'compare' ? 'bg-gradient-to-r from-blue-600 to-purple-600' : mode === 'quantum' ? 'bg-purple-600 hover:bg-purple-700' : ''}`}
              aria-label={loading ? 'Analyzing molecule...' : 'Analyze molecule'}
            >
              {loading ? (
                <>
                  <Activity className="w-3.5 h-3.5 sm:w-4 sm:h-4 animate-spin" /> <span className="text-xs sm:text-sm">Analyzing...</span>
                </>
              ) : (
                <>
                  <Play className="w-3.5 h-3.5 sm:w-4 sm:h-4" /> <span className="text-xs sm:text-sm">{mode === 'compare' ? 'Compare' : 'Predict'}</span>
                </>
              )}
            </Button>
          </div>
        </div>

        {/* Right Panel - Results */}
        <div className="lg:col-span-4 flex flex-col overflow-hidden min-h-0 flex-1">
          {result && !result.error && displayData ? (
            <>
              {/* Top: Result Summary */}
              <div className="mb-3 sm:mb-4">
                {mode === 'compare' ? (
                  <div className="grid grid-cols-1 sm:grid-cols-2 gap-3 sm:gap-4">
                    <Card className="p-4 bg-gradient-to-br from-blue-50 to-cyan-50 dark:from-blue-950/50 dark:to-cyan-950/50 border-2 border-blue-200 dark:border-blue-800">
                      <div className="flex justify-between items-start mb-3">
                        <div className="flex items-center gap-2">
                          <h3 className="text-xs font-medium text-gray-600 dark:text-gray-400 uppercase">Classical XGBoost</h3>
                          <InfoTooltip title={TOOLTIPS.xgboost.title} content={TOOLTIPS.xgboost.content} />
                        </div>
                        <Badge variant="success">✓ {(displayData.classical.confidence * 100).toFixed(0)}%</Badge>
                      </div>
                      <div className="text-4xl sm:text-5xl font-bold text-blue-600 dark:text-blue-400 mb-1">
                        {displayData.classical.prediction}
                      </div>
                      <div className="text-sm text-gray-600 dark:text-gray-400">Log Solubility (mol/L)</div>
                    </Card>

                    <Card className="p-4 bg-gradient-to-br from-purple-50 to-pink-50 dark:from-purple-950/50 dark:to-pink-950/50 border-2 border-purple-200 dark:border-purple-800">
                      <div className="flex justify-between items-start mb-3">
                        <div className="flex items-center gap-2">
                          <h3 className="text-xs font-medium text-gray-600 dark:text-gray-400 uppercase">Quantum VQC</h3>
                          <InfoTooltip title={TOOLTIPS.vqc.title} content={TOOLTIPS.vqc.content} />
                        </div>
                        <Badge variant="purple">⚛ {(displayData.quantum.confidence * 100).toFixed(0)}%</Badge>
                      </div>
                      <div className="text-4xl sm:text-5xl font-bold text-purple-600 dark:text-purple-400 mb-1">
                        {displayData.quantum.prediction}
                      </div>
                      <div className="text-sm text-gray-600 dark:text-gray-400">Log Solubility (mol/L)</div>
                    </Card>
                  </div>
                ) : (
                  <Card className="bg-gradient-to-br from-blue-50 to-purple-50 dark:from-blue-950/50 dark:to-purple-950/50 border-2 border-blue-200 dark:border-blue-800 p-4">
                    <div className="flex justify-between items-center">
                      <div>
                        <div className="text-4xl sm:text-5xl font-bold text-blue-600 dark:text-blue-400 mb-1">
                          {displayData.single.prediction}
                        </div>
                        <div className="text-sm text-gray-600 dark:text-gray-400">Log Solubility (mol/L)</div>
                      </div>
                      <Badge variant={mode === 'quantum' ? 'purple' : 'success'}>{displayData.single.model_used}</Badge>
                    </div>
                  </Card>
                )}
              </div>

              {/* Tabs - Mobile Scrollable */}
              <div className="flex gap-1 sm:gap-2 mb-3 overflow-x-auto pb-1 scrollbar-hide -mx-3 sm:mx-0 px-3 sm:px-0 snap-x">
                <button
                  onClick={() => setActiveTab('druglikeness')}
                  className={`px-3 sm:px-4 py-2 rounded-t-lg text-xs sm:text-sm font-medium transition-colors whitespace-nowrap snap-start touch-manipulation ${activeTab === 'druglikeness'
                    ? 'bg-white dark:bg-gray-900 border-t border-x border-gray-200 dark:border-gray-800'
                    : 'bg-gray-100 dark:bg-gray-800 text-gray-600 hover:bg-gray-200'
                    }`}
                >
                  <Shield className="w-3 h-3 sm:w-4 sm:h-4 inline mr-1" /> <span className="hidden sm:inline">Drug </span>Likeness
                </button>
                <button
                  onClick={() => setActiveTab('3d')}
                  className={`px-3 sm:px-4 py-2 rounded-t-lg text-xs sm:text-sm font-medium transition-colors whitespace-nowrap snap-start touch-manipulation ${activeTab === '3d'
                    ? 'bg-white dark:bg-gray-900 border-t border-x border-gray-200 dark:border-gray-800'
                    : 'bg-gray-100 dark:bg-gray-800 text-gray-600 hover:bg-gray-200'
                    }`}
                >
                  <Pill className="w-3 h-3 sm:w-4 sm:h-4 inline mr-1" /> 3D
                </button>
                <button
                  onClick={() => setActiveTab('overview')}
                  className={`px-3 sm:px-4 py-2 rounded-t-lg text-xs sm:text-sm font-medium transition-colors whitespace-nowrap snap-start touch-manipulation ${activeTab === 'overview'
                    ? 'bg-white dark:bg-gray-900 border-t border-x border-gray-200 dark:border-gray-800'
                    : 'bg-gray-100 dark:bg-gray-800 text-gray-600 hover:bg-gray-200'
                    }`}
                >
                  <Activity className="w-3 h-3 sm:w-4 sm:h-4 inline mr-1" /> Overview
                </button>
                {mode === 'quantum' && displayData?.single?.quantum_circuit && (
                  <button
                    onClick={() => setActiveTab('circuit')}
                    className={`px-3 sm:px-4 py-2 rounded-t-lg text-xs sm:text-sm font-medium transition-colors whitespace-nowrap snap-start touch-manipulation ${activeTab === 'circuit'
                      ? 'bg-white dark:bg-gray-900 border-t border-x border-gray-200 dark:border-gray-800'
                      : 'bg-gray-100 dark:bg-gray-800 text-gray-600 hover:bg-gray-200'
                      }`}
                  >
                    <Atom className="w-3 h-3 sm:w-4 sm:h-4 inline mr-1" /> <span className="hidden sm:inline">Quantum </span>Circuit
                    <InfoTooltip title={TOOLTIPS.quantumCircuit.title} content={TOOLTIPS.quantumCircuit.content} className="ml-1 hidden sm:inline" />
                  </button>
                )}
                <button
                  onClick={() => setActiveTab('details')}
                  className={`px-3 sm:px-4 py-2 rounded-t-lg text-xs sm:text-sm font-medium transition-colors whitespace-nowrap snap-start touch-manipulation ${activeTab === 'details'
                    ? 'bg-white dark:bg-gray-900 border-t border-x border-gray-200 dark:border-gray-800'
                    : 'bg-gray-100 dark:bg-gray-800 text-gray-600 hover:bg-gray-200'
                    }`}
                >
                  <Beaker className="w-3 h-3 sm:w-4 sm:h-4 inline mr-1" /> Details
                </button>
                {result?.marketAnalysis && (
                  <button
                    onClick={() => setActiveTab('market')}
                    className={`px-3 sm:px-4 py-2 rounded-t-lg text-xs sm:text-sm font-medium transition-colors whitespace-nowrap snap-start touch-manipulation ${activeTab === 'market'
                      ? 'bg-white dark:bg-gray-900 border-t border-x border-gray-200 dark:border-gray-800'
                      : 'bg-gray-100 dark:bg-gray-800 text-gray-600 hover:bg-gray-200'
                      }`}
                  >
                    <TrendingUp className="w-3 h-3 sm:w-4 sm:h-4 inline mr-1" /> Market
                  </button>
                )}
              </div>

              {/* Tab Content */}
              <div className="flex-1 overflow-y-auto">
                {activeTab === 'druglikeness' && result.drug_likeness && (
                  <DrugLikenessDashboard data={result.drug_likeness} />
                )}

                {activeTab === '3d' && (
                  result.structure_3d ? (
                    result.structure_3d.sdf ? (
                      <Card className="p-4">
                        <h3 className="font-semibold mb-3 text-sm">Interactive 3D Molecular Structure</h3>
                        {mode === 'compare' ? (
                          <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
                            <Interactive3DViewer
                              sdf={result.structure_3d.sdf}
                              title="Classical"
                              viewerKey={`classical-${result.smiles}`}
                            />
                            <Interactive3DViewer
                              sdf={result.structure_3d.sdf}
                              title="Quantum"
                              viewerKey={`quantum-${result.smiles}`}
                            />
                          </div>
                        ) : (
                          <Interactive3DViewer
                            sdf={result.structure_3d.sdf}
                            viewerKey={`single-${result.smiles}`}
                          />
                        )}
                      </Card>
                    ) : (
                      <Card className="p-8 h-96 flex items-center justify-center">
                        <div className="text-center text-gray-500">
                          <p>No SDF data available</p>
                        </div>
                      </Card>
                    )
                  ) : (
                    <Card className="p-8 h-96 flex items-center justify-center">
                      <div className="text-center text-gray-500">
                        <p className="font-semibold mb-2">No 3D Structure Data</p>
                        <p className="text-sm">Check console for API response</p>
                      </div>
                    </Card>
                  )
                )}

                {activeTab === 'overview' && (
                  <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
                    <Card className="p-4">
                      <h3 className="font-semibold mb-4 text-sm">Molecular Profile</h3>
                      <ResponsiveContainer width="100%" height={250}>
                        <RadarChart data={[
                          { property: 'MW/100', value: displayData.features.molecular_weight / 100 },
                          { property: 'LogP', value: Math.max(displayData.features.logp, 0) },
                          { property: 'HBD', value: displayData.features.num_h_donors },
                          { property: 'HBA', value: displayData.features.num_h_acceptors },
                          { property: 'RotB', value: displayData.features.num_rotatable_bonds },
                          { property: 'TPSA/50', value: displayData.features.tpsa / 50 },
                        ]}>
                          <PolarGrid stroke="#333" />
                          <PolarAngleAxis dataKey="property" tick={{ fill: '#888', fontSize: 11 }} />
                          <PolarRadiusAxis angle={90} domain={[0, 5]} tick={{ fill: '#888' }} />
                          <Radar name="Value" dataKey="value" stroke={mode === 'quantum' ? '#9333ea' : '#2563eb'} fill={mode === 'quantum' ? '#9333ea' : '#2563eb'} fillOpacity={0.6} />
                        </RadarChart>
                      </ResponsiveContainer>
                    </Card>

                    <Card className="p-4">
                      <h3 className="font-semibold mb-4 text-sm">Descriptor Values</h3>
                      <ResponsiveContainer width="100%" height={250}>
                        <BarChart data={[
                          { name: 'MW/10', value: displayData.features.molecular_weight / 10 },
                          { name: 'LogP', value: displayData.features.logp },
                          { name: 'HBD', value: displayData.features.num_h_donors },
                          { name: 'HBA', value: displayData.features.num_h_acceptors },
                          { name: 'RotB', value: displayData.features.num_rotatable_bonds },
                        ]}>
                          <CartesianGrid strokeDasharray="3 3" stroke="#333" vertical={false} />
                          <XAxis dataKey="name" stroke="#888" tick={{ fontSize: 11 }} />
                          <YAxis stroke="#888" tick={{ fontSize: 11 }} />
                          <Tooltip contentStyle={{ backgroundColor: '#1f2937', border: 'none', borderRadius: '8px' }} />
                          <Bar dataKey="value" fill={mode === 'quantum' ? '#9333ea' : '#2563eb'} radius={[4, 4, 0, 0]} />
                        </BarChart>
                      </ResponsiveContainer>
                    </Card>
                  </div>
                )}

                {activeTab === 'circuit' && mode === 'quantum' && displayData.single?.quantum_circuit && (
                  <QuantumCircuitVisualization circuit={displayData.single.quantum_circuit} />
                )}

                {activeTab === 'details' && (
                  <Card className="overflow-hidden">
                    <div className="overflow-x-auto -mx-4 sm:mx-0">
                      <table className="w-full text-sm min-w-full">
                        <thead className="bg-gray-50 dark:bg-gray-800/30 text-gray-600 dark:text-gray-400">
                          <tr>
                            <th className="px-4 sm:px-6 py-3 font-medium text-left">Property</th>
                            <th className="px-4 sm:px-6 py-3 font-medium text-right">Value</th>
                          </tr>
                        </thead>
                        <tbody className="divide-y divide-gray-200 dark:divide-gray-700">
                          {displayData.features.pubchem_cid && (
                            <tr className="bg-blue-50 dark:bg-blue-900/20">
                              <td className="px-4 sm:px-6 py-3">
                                <div className="flex items-center gap-2">
                                  PubChem CID
                                </div>
                              </td>
                              <td className="px-4 sm:px-6 py-3 text-right">
                                {displayData.features.pubchem_url ? (
                                  <a href={displayData.features.pubchem_url} target="_blank" rel="noopener noreferrer" className="font-mono text-blue-600 dark:text-blue-400 hover:underline">
                                    {displayData.features.pubchem_cid}
                                  </a>
                                ) : (
                                  <span className="font-mono">{displayData.features.pubchem_cid}</span>
                                )}
                              </td>
                            </tr>
                          )}
                          {displayData.features.iupac_name && (
                            <tr>
                              <td className="px-4 sm:px-6 py-3">
                                <div className="flex items-center gap-2">
                                  IUPAC Name
                                </div>
                              </td>
                              <td className="px-4 sm:px-6 py-3 text-right text-xs sm:text-sm break-words">{displayData.features.iupac_name}</td>
                            </tr>
                          )}
                          {displayData.features.molecule_name && (
                            <tr>
                              <td className="px-4 sm:px-6 py-3">
                                <div className="flex items-center gap-2">
                                  Common Name
                                </div>
                              </td>
                              <td className="px-4 sm:px-6 py-3 text-right font-medium">{displayData.features.molecule_name}</td>
                            </tr>
                          )}
                          <tr>
                            <td className="px-4 sm:px-6 py-3">
                            <div className="flex items-center gap-2">
                              Molecular Weight
                              <InfoTooltip title={TOOLTIPS.molecularWeight.title} content={TOOLTIPS.molecularWeight.content} />
                            </div>
                          </td>
                          <td className="px-4 sm:px-6 py-3 text-right font-mono">{displayData.features.molecular_weight.toFixed(2)}</td>
                        </tr>
                        <tr>
                          <td className="px-4 sm:px-6 py-3">
                            <div className="flex items-center gap-2">
                              LogP
                              <InfoTooltip title={TOOLTIPS.logp.title} content={TOOLTIPS.logp.content} />
                            </div>
                          </td>
                          <td className="px-4 sm:px-6 py-3 text-right font-mono">{displayData.features.logp.toFixed(2)}</td>
                        </tr>
                        <tr>
                          <td className="px-4 sm:px-6 py-3">
                            <div className="flex items-center gap-2">
                              H-Bond Donors
                              <InfoTooltip title={TOOLTIPS.hDonors.title} content={TOOLTIPS.hDonors.content} />
                            </div>
                          </td>
                          <td className="px-4 sm:px-6 py-3 text-right font-mono">{displayData.features.num_h_donors}</td>
                        </tr>
                        <tr>
                          <td className="px-4 sm:px-6 py-3">
                            <div className="flex items-center gap-2">
                              H-Bond Acceptors
                              <InfoTooltip title={TOOLTIPS.hAcceptors.title} content={TOOLTIPS.hAcceptors.content} />
                            </div>
                          </td>
                          <td className="px-4 sm:px-6 py-3 text-right font-mono">{displayData.features.num_h_acceptors}</td>
                        </tr>
                        <tr>
                          <td className="px-4 sm:px-6 py-3">
                            <div className="flex items-center gap-2">
                              Rotatable Bonds
                              <InfoTooltip title={TOOLTIPS.rotatableBonds.title} content={TOOLTIPS.rotatableBonds.content} />
                            </div>
                          </td>
                          <td className="px-4 sm:px-6 py-3 text-right font-mono">{displayData.features.num_rotatable_bonds}</td>
                        </tr>
                        <tr>
                          <td className="px-4 sm:px-6 py-3">
                            <div className="flex items-center gap-2">
                              Aromatic Rings
                              <InfoTooltip title={TOOLTIPS.aromaticRings.title} content={TOOLTIPS.aromaticRings.content} />
                            </div>
                          </td>
                          <td className="px-4 sm:px-6 py-3 text-right font-mono">{displayData.features.num_aromatic_rings}</td>
                        </tr>
                        <tr>
                          <td className="px-4 sm:px-6 py-3">
                            <div className="flex items-center gap-2">
                              TPSA
                              <InfoTooltip title={TOOLTIPS.tpsa.title} content={TOOLTIPS.tpsa.content} />
                            </div>
                          </td>
                          <td className="px-4 sm:px-6 py-3 text-right font-mono">{displayData.features.tpsa.toFixed(2)}</td>
                        </tr>
                        <tr>
                          <td className="px-4 sm:px-6 py-3">
                            <div className="flex items-center gap-2">
                              Total Atoms
                              <InfoTooltip title={TOOLTIPS.totalAtoms.title} content={TOOLTIPS.totalAtoms.content} />
                            </div>
                          </td>
                          <td className="px-4 sm:px-6 py-3 text-right font-mono">{displayData.features.num_atoms}</td>
                        </tr>
                        {displayData.features.synonyms && displayData.features.synonyms.length > 0 && (
                          <tr>
                            <td className="px-4 sm:px-6 py-3">
                              <div className="flex items-center gap-2">
                                Synonyms
                              </div>
                            </td>
                            <td className="px-4 sm:px-6 py-3 text-right text-xs sm:text-sm">
                              <div className="flex flex-wrap gap-1 justify-end">
                                {displayData.features.synonyms.slice(0, 5).map((syn: string, idx: number) => (
                                  <span key={idx} className="px-2 py-0.5 bg-gray-100 dark:bg-gray-800 rounded text-[10px] sm:text-xs">
                                    {syn}
                                  </span>
                                ))}
                              </div>
                            </td>
                          </tr>
                        )}
                        </tbody>
                      </table>
                    </div>
                  </Card>
                )}

                {activeTab === 'market' && result.marketAnalysis && (
                  <MarketAnalysisDashboard data={result.marketAnalysis} />
                )}

                {activeTab === 'editor' && (
                  <Card className="p-6">
                    <MolecularEditor
                      initialSmiles={smiles}
                      onChange={(newSmiles) => {
                        setSmiles(newSmiles);
                        // Optionally auto-analyze when SMILES changes from editor
                        // Uncomment if you want auto-analysis:
                        // if (newSmiles && newSmiles.trim()) {
                        //   setTimeout(() => analyze(), 500);
                        // }
                      }}
                    />
                    <div className="mt-4 flex gap-2">
                      <Button
                        onClick={() => {
                          if (smiles && smiles.trim()) {
                            analyze();
                          }
                        }}
                        disabled={!smiles || !smiles.trim() || loading}
                      >
                        <Play className="w-4 h-4" />
                        Analyze Current SMILES
                      </Button>
                      <Button
                        variant="secondary"
                        onClick={() => setSmiles('CCO')}
                      >
                        Reset to Example
                      </Button>
                    </div>
                  </Card>
                )}
              </div>
            </>
          ) : result?.error ? (
            <Card className="p-6 sm:p-8 border-2 border-red-200 dark:border-red-900 flex items-center justify-center h-full min-h-[300px]">
              <div className="text-red-600 dark:text-red-400 text-center px-4">
                <p className="font-semibold text-base sm:text-lg mb-2">Error</p>
                <p className="text-xs sm:text-sm break-words">{result.error}</p>
              </div>
            </Card>
          ) : (
            <div className="h-full flex flex-col items-center justify-center text-gray-400 border-2 border-dashed border-gray-300 dark:border-gray-700 rounded-lg p-6 sm:p-8">
              <GitCompare className="w-16 h-16 sm:w-20 sm:h-20 mb-4 opacity-20" />
              <p className="text-lg sm:text-xl font-medium text-center">Ready to Analyze</p>
              <p className="text-xs sm:text-sm mt-2 text-center max-w-md">Enter a SMILES string and click {mode === 'compare' ? 'Compare Models' : 'Predict'}</p>
            </div>
          )}
        </div>
      </div>

      {/* Chat Panel - Floating Overlay - Mobile Optimized */}
      {showChat && (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/50 backdrop-blur-sm p-0 sm:p-4" onClick={() => setShowChat(false)}>
          <div 
            className="w-full h-full sm:w-full sm:max-w-2xl sm:h-[80vh] bg-white dark:bg-gray-900 sm:rounded-lg shadow-2xl border-0 sm:border border-gray-200 dark:border-gray-800 flex flex-col"
            onClick={(e) => e.stopPropagation()}
          >
            <div className="flex items-center justify-between p-4 border-b border-gray-200 dark:border-gray-800">
              <div className="flex items-center gap-2">
                <MessageSquare className="w-5 h-5 text-blue-600 dark:text-blue-400" />
                <h2 className="text-lg font-semibold">AI Assistant</h2>
              </div>
              <button
                onClick={() => setShowChat(false)}
                className="p-2 hover:bg-gray-100 dark:hover:bg-gray-800 rounded-lg transition-colors"
              >
                <X className="w-5 h-5" />
              </button>
            </div>
            <div className="flex-1 overflow-hidden">
              <ChatInterface onAnalyzeSmiles={handleChatAnalyze} />
            </div>
          </div>
        </div>
      )}
    </div>
  );
}
