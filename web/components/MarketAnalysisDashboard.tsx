'use client';

import React from 'react';
import { TrendingUp, FileText, Building2, AlertCircle, CheckCircle, XCircle } from 'lucide-react';
import { MarketAnalysis } from '../lib/molecular_inference';

interface MarketAnalysisDashboardProps {
    data: MarketAnalysis;
}

const StatusBadge = ({ status }: { status: string }) => {
    const isPositive = status.toLowerCase().includes('favorable') || status.toLowerCase().includes('approved');
    return (
        <span className={`px-3 py-1 rounded-full text-xs font-medium ${
            isPositive 
                ? 'bg-green-100 text-green-800 dark:bg-green-900 dark:text-green-200'
                : status.toLowerCase().includes('moderate')
                ? 'bg-yellow-100 text-yellow-800 dark:bg-yellow-900 dark:text-yellow-200'
                : 'bg-red-100 text-red-800 dark:bg-red-900 dark:text-red-200'
        }`}>
            {status}
        </span>
    );
};

export default function MarketAnalysisDashboard({ data }: MarketAnalysisDashboardProps) {
    return (
        <div className="space-y-6">
            {/* Data Source Info */}
            {(data as any).dataSource && (
                <div className="bg-blue-50 dark:bg-blue-900/20 border border-blue-200 dark:border-blue-800 rounded-lg p-3 sm:p-4">
                    <p className="text-xs sm:text-sm text-blue-800 dark:text-blue-200">
                        <span className="font-semibold">Data Sources:</span> {(data as any).dataSource}
                    </p>
                </div>
            )}
            {/* Market Size Overview */}
            <div className="grid grid-cols-1 lg:grid-cols-2 gap-4 sm:gap-6">
                <div className="bg-white dark:bg-gray-900 p-4 sm:p-6 rounded-xl border border-gray-200 dark:border-gray-800">
                    <div className="flex items-center gap-2 sm:gap-3 mb-3 sm:mb-4">
                        <TrendingUp className="w-5 h-5 sm:w-6 sm:h-6 text-blue-600 flex-shrink-0" />
                        <h3 className="text-sm sm:text-lg font-semibold">Market Size</h3>
                    </div>
                    <div className="space-y-3">
                        <div>
                            <p className="text-xs sm:text-sm text-gray-600 dark:text-gray-400">Therapeutic Area</p>
                            <p className="text-lg sm:text-xl font-bold text-gray-900 dark:text-gray-100">{data.marketSize.therapeuticArea}</p>
                        </div>
                        <div>
                            <p className="text-xs sm:text-sm text-gray-600 dark:text-gray-400">Estimated Market Size</p>
                            <p className="text-2xl sm:text-3xl font-bold text-blue-600 dark:text-blue-400">
                                ${data.marketSize.estimatedMarketSize}B
                            </p>
                            <p className="text-[10px] sm:text-xs text-gray-500 dark:text-gray-400">{data.marketSize.currency} ({data.marketSize.year})</p>
                        </div>
                        <div>
                            <p className="text-xs sm:text-sm text-gray-600 dark:text-gray-400">Annual Growth Rate</p>
                            <p className="text-base sm:text-lg font-semibold text-green-600 dark:text-green-400">{data.marketSize.growthRate}</p>
                        </div>
                    </div>
                </div>

                <div className="bg-white dark:bg-gray-900 p-4 sm:p-6 rounded-xl border border-gray-200 dark:border-gray-800">
                    <div className="flex items-center gap-2 sm:gap-3 mb-3 sm:mb-4">
                        <AlertCircle className="w-5 h-5 sm:w-6 sm:h-6 text-purple-600 flex-shrink-0" />
                        <h3 className="text-sm sm:text-lg font-semibold">Regulatory Status</h3>
                    </div>
                    <div className="flex flex-col items-center justify-center h-full">
                        <StatusBadge status={data.regulatoryStatus} />
                        <p className="text-sm text-gray-600 dark:text-gray-400 mt-4 text-center">
                            Based on drug-likeness criteria and molecular properties
                        </p>
                    </div>
                </div>
            </div>

            {/* Competitors */}
            <div className="bg-white dark:bg-gray-900 rounded-xl border border-gray-200 dark:border-gray-800 overflow-hidden">
                <div className="px-4 sm:px-6 py-3 sm:py-4 border-b border-gray-200 dark:border-gray-800 bg-gray-50 dark:bg-gray-800/50">
                    <div className="flex items-center gap-2 sm:gap-3 flex-wrap">
                        <Building2 className="w-4 h-4 sm:w-5 sm:h-5 text-gray-600 dark:text-gray-400 flex-shrink-0" />
                        <h3 className="text-sm sm:text-base font-semibold">Competitor Drugs</h3>
                        <span className="ml-auto text-xs sm:text-sm text-gray-500 dark:text-gray-400">
                            {data.competitors.length} similar molecule{data.competitors.length !== 1 ? 's' : ''} found
                        </span>
                    </div>
                </div>
                <div className="overflow-x-auto -mx-4 sm:mx-0 scrollbar-hide">
                    {data.competitors.length > 0 ? (
                        <table className="w-full text-xs sm:text-sm min-w-full">
                            <thead className="bg-gray-50 dark:bg-gray-800/30 text-gray-600 dark:text-gray-400">
                                <tr>
                                    <th className="px-4 sm:px-6 py-2 sm:py-3 font-medium text-left">Drug Name</th>
                                    <th className="px-4 sm:px-6 py-2 sm:py-3 font-medium text-left">Indication</th>
                                    <th className="px-4 sm:px-6 py-2 sm:py-3 font-medium text-center">Similarity</th>
                                    <th className="px-4 sm:px-6 py-2 sm:py-3 font-medium text-left">Status</th>
                                    <th className="px-4 sm:px-6 py-2 sm:py-3 font-medium text-left">Company</th>
                                </tr>
                            </thead>
                            <tbody className="divide-y divide-gray-200 dark:divide-gray-700">
                                {data.competitors.map((competitor, idx) => (
                                    <tr key={idx} className="hover:bg-gray-50 dark:hover:bg-gray-800/50">
                                        <td className="px-4 sm:px-6 py-3 sm:py-4 font-medium text-xs sm:text-sm">{competitor.name}</td>
                                        <td className="px-4 sm:px-6 py-3 sm:py-4 text-gray-600 dark:text-gray-400 text-xs sm:text-sm">{competitor.indication}</td>
                                        <td className="px-4 sm:px-6 py-3 sm:py-4 text-center">
                                            <span className="px-2 py-1 bg-blue-100 dark:bg-blue-900/30 text-blue-800 dark:text-blue-200 rounded text-xs sm:text-sm whitespace-nowrap">
                                                {(competitor.similarity * 100).toFixed(1)}%
                                            </span>
                                        </td>
                                        <td className="px-4 sm:px-6 py-3 sm:py-4">
                                            {competitor.status === 'Approved' ? (
                                                <span className="flex items-center gap-1 text-green-600 dark:text-green-400 text-xs sm:text-sm">
                                                    <CheckCircle className="w-3 h-3 sm:w-4 sm:h-4" />
                                                    {competitor.status}
                                                </span>
                                            ) : (
                                                <span className="text-gray-600 dark:text-gray-400 text-xs sm:text-sm">{competitor.status}</span>
                                            )}
                                        </td>
                                        <td className="px-4 sm:px-6 py-3 sm:py-4 text-gray-600 dark:text-gray-400 text-xs sm:text-sm">{competitor.company}</td>
                                    </tr>
                                ))}
                            </tbody>
                        </table>
                    ) : (
                        <div className="p-8 text-center text-gray-500 dark:text-gray-400">
                            <p>No similar competitor drugs found in database</p>
                            <p className="text-xs mt-2">This may indicate a novel compound</p>
                        </div>
                    )}
                </div>
            </div>

            {/* Patents */}
            <div className="bg-white dark:bg-gray-900 rounded-xl border border-gray-200 dark:border-gray-800 overflow-hidden">
                <div className="px-4 sm:px-6 py-3 sm:py-4 border-b border-gray-200 dark:border-gray-800 bg-gray-50 dark:bg-gray-800/50">
                    <div className="flex items-center gap-2 sm:gap-3 flex-wrap">
                        <FileText className="w-4 h-4 sm:w-5 sm:h-5 text-gray-600 dark:text-gray-400 flex-shrink-0" />
                        <h3 className="text-sm sm:text-base font-semibold">Patent Landscape</h3>
                        <span className="ml-auto text-xs sm:text-sm text-gray-500 dark:text-gray-400">
                            {data.patents.length} relevant patent{data.patents.length !== 1 ? 's' : ''}
                        </span>
                    </div>
                </div>
                <div className="p-4 sm:p-6">
                    {data.patents.length > 0 ? (
                        <div className="space-y-3 sm:space-y-4">
                            {data.patents.map((patent, idx) => (
                                <div key={idx} className="border border-gray-200 dark:border-gray-700 rounded-lg p-3 sm:p-4 hover:bg-gray-50 dark:hover:bg-gray-800/50">
                                    <div className="flex items-start justify-between mb-2 gap-2 sm:gap-4">
                                        <div className="min-w-0 flex-1">
                                            <p className="font-medium text-xs sm:text-sm text-gray-900 dark:text-gray-100 break-words">{patent.patentNumber}</p>
                                            <p className="text-xs sm:text-sm text-gray-600 dark:text-gray-400 mt-1 break-words">{patent.title}</p>
                                        </div>
                                        <span className={`px-2 py-1 rounded text-[10px] sm:text-xs font-medium whitespace-nowrap flex-shrink-0 ${
                                            patent.relevance === 'High'
                                                ? 'bg-red-100 text-red-800 dark:bg-red-900 dark:text-red-200'
                                                : 'bg-yellow-100 text-yellow-800 dark:bg-yellow-900 dark:text-yellow-200'
                                        }`}>
                                            {patent.relevance}
                                        </span>
                                    </div>
                                    <div className="grid grid-cols-2 sm:grid-cols-4 gap-3 sm:gap-4 mt-3 text-xs sm:text-sm">
                                        <div>
                                            <p className="text-[10px] sm:text-xs text-gray-500 dark:text-gray-400">Filing Date</p>
                                            <p className="font-medium text-xs sm:text-sm">{patent.filingDate}</p>
                                        </div>
                                        <div>
                                            <p className="text-[10px] sm:text-xs text-gray-500 dark:text-gray-400">Status</p>
                                            <p className="font-medium text-xs sm:text-sm">{patent.status}</p>
                                        </div>
                                        <div>
                                            <p className="text-[10px] sm:text-xs text-gray-500 dark:text-gray-400">Assignee</p>
                                            <p className="font-medium text-xs sm:text-sm truncate">{patent.assignee}</p>
                                        </div>
                                        <div>
                                            <p className="text-[10px] sm:text-xs text-gray-500 dark:text-gray-400">Relevance</p>
                                            <p className="font-medium text-xs sm:text-sm">{patent.relevance}</p>
                                        </div>
                                    </div>
                                </div>
                            ))}
                        </div>
                    ) : (
                        <div className="p-6 sm:p-8 text-center text-gray-500 dark:text-gray-400">
                            <p className="text-sm sm:text-base">No relevant patents found</p>
                            <p className="text-xs mt-2">This may indicate a clear IP landscape</p>
                        </div>
                    )}
                </div>
            </div>
        </div>
    );
}

