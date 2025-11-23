import { MolecularAnalysisResult } from './molecular_inference';

// Type definitions for jsPDF (avoiding direct import for SSR)
type jsPDFType = any;
type AutoTableType = any;

export interface ReportSection {
    title: string;
    content: (doc: jsPDFType, data: MolecularAnalysisResult, autoTable: AutoTableType) => void;
}

export function getReportSections(jsPDF: jsPDFType, autoTable: AutoTableType): ReportSection[] {
    return [
    {
        title: 'Cover Page',
        content: (doc, data, autoTable) => {
            doc.setFontSize(24);
            doc.setTextColor(37, 99, 235); // Blue
            doc.text('MoleculeAI', 105, 50, { align: 'center' });
            
            doc.setFontSize(16);
            doc.setTextColor(100, 100, 100);
            doc.text('Molecular Property Analysis Report', 105, 65, { align: 'center' });
            
            doc.setFontSize(12);
            doc.setTextColor(0, 0, 0);
            doc.text(`SMILES: ${data.smiles}`, 105, 100, { align: 'center' });
            
            doc.setFontSize(10);
            doc.setTextColor(150, 150, 150);
            const date = new Date().toLocaleDateString('en-US', { 
                year: 'numeric', 
                month: 'long', 
                day: 'numeric' 
            });
            doc.text(`Generated: ${date}`, 105, 120, { align: 'center' });
        }
    },
    {
        title: 'Executive Summary',
        content: (doc, data, autoTable) => {
            doc.addPage();
            doc.setFontSize(18);
            doc.setTextColor(37, 99, 235);
            doc.text('Executive Summary', 20, 30);
            
            doc.setFontSize(11);
            doc.setTextColor(0, 0, 0);
            
            let yPos = 50;
            
            // Key Prediction
            doc.setFontSize(12);
            doc.setFont(undefined, 'bold');
            doc.text('Key Prediction:', 20, yPos);
            yPos += 10;
            
            doc.setFont(undefined, 'normal');
            if (data.classical && data.quantum) {
                doc.text(`Solubility (Classical): ${data.classical.prediction} log(mol/L)`, 25, yPos);
                yPos += 7;
                doc.text(`Solubility (Quantum): ${data.quantum.prediction} log(mol/L)`, 25, yPos);
            } else if (data.prediction) {
                doc.text(`Solubility: ${data.prediction} log(mol/L)`, 25, yPos);
            }
            yPos += 15;
            
            // Drug Score
            if (data.drug_likeness) {
                doc.setFont(undefined, 'bold');
                doc.text('Drug-Likeness Score:', 20, yPos);
                yPos += 10;
                doc.setFont(undefined, 'normal');
                doc.text(`${data.drug_likeness.drug_score}/10`, 25, yPos);
                yPos += 7;
                doc.text(`Lipinski's Rule: ${data.drug_likeness.lipinski.pass ? 'PASS' : 'FAIL'} (${data.drug_likeness.lipinski.violations} violations)`, 25, yPos);
                yPos += 15;
            }
            
            // Market Analysis Summary
            if (data.marketAnalysis) {
                doc.setFont(undefined, 'bold');
                doc.text('Market Overview:', 20, yPos);
                yPos += 10;
                doc.setFont(undefined, 'normal');
                doc.text(`Therapeutic Area: ${data.marketAnalysis.marketSize.therapeuticArea}`, 25, yPos);
                yPos += 7;
                doc.text(`Market Size: $${data.marketAnalysis.marketSize.estimatedMarketSize}B ${data.marketAnalysis.marketSize.currency}`, 25, yPos);
                yPos += 7;
                doc.text(`Growth Rate: ${data.marketAnalysis.marketSize.growthRate}`, 25, yPos);
                yPos += 7;
                doc.text(`Competitors Found: ${data.marketAnalysis.competitors.length}`, 25, yPos);
                yPos += 7;
                doc.text(`Relevant Patents: ${data.marketAnalysis.patents.length}`, 25, yPos);
            }
        }
    },
    {
        title: 'Molecular Descriptors',
        content: (doc, data, autoTable) => {
            doc.addPage();
            doc.setFontSize(18);
            doc.setTextColor(37, 99, 235);
            doc.text('Molecular Descriptors', 20, 30);
            
            if (data.features) {
                autoTable(doc, {
                    startY: 40,
                    head: [['Property', 'Value']],
                    body: [
                        ['Molecular Weight', `${data.features.molecular_weight.toFixed(2)} Da`],
                        ['LogP', data.features.logp.toFixed(2)],
                        ['H-Bond Donors', data.features.num_h_donors.toString()],
                        ['H-Bond Acceptors', data.features.num_h_acceptors.toString()],
                        ['Rotatable Bonds', data.features.num_rotatable_bonds.toString()],
                        ['Aromatic Rings', data.features.num_aromatic_rings.toString()],
                        ['TPSA', `${data.features.tpsa.toFixed(2)}`],
                        ['Total Atoms', data.features.num_atoms.toString()],
                    ],
                    theme: 'striped',
                    headStyles: { fillColor: [37, 99, 235] },
                });
            }
        }
    },
    {
        title: 'Drug-Likeness Analysis',
        content: (doc, data, autoTable) => {
            if (!data.drug_likeness) return;
            
            doc.addPage();
            doc.setFontSize(18);
            doc.setTextColor(37, 99, 235);
            doc.text('Drug-Likeness Analysis', 20, 30);
            
            let yPos = 50;
            doc.setFontSize(11);
            doc.setTextColor(0, 0, 0);
            
            // Lipinski's Rule
            doc.setFont(undefined, 'bold');
            doc.text('Lipinski\'s Rule of Five:', 20, yPos);
            yPos += 10;
            doc.setFont(undefined, 'normal');
            
            const lipinski = data.drug_likeness.lipinski.details;
            doc.text(`Molecular Weight: ${lipinski.molecular_weight.value.toFixed(1)} / ${lipinski.molecular_weight.limit} Da ${lipinski.molecular_weight.pass ? '✓' : '✗'}`, 25, yPos);
            yPos += 7;
            doc.text(`LogP: ${lipinski.logp.value.toFixed(2)} / ${lipinski.logp.limit} ${lipinski.logp.pass ? '✓' : '✗'}`, 25, yPos);
            yPos += 7;
            doc.text(`H-Bond Donors: ${lipinski.h_donors.value} / ${lipinski.h_donors.limit} ${lipinski.h_donors.pass ? '✓' : '✗'}`, 25, yPos);
            yPos += 7;
            doc.text(`H-Bond Acceptors: ${lipinski.h_acceptors.value} / ${lipinski.h_acceptors.limit} ${lipinski.h_acceptors.pass ? '✓' : '✗'}`, 25, yPos);
            yPos += 15;
            
            // Veber's Rule
            doc.setFont(undefined, 'bold');
            doc.text('Veber\'s Rule:', 20, yPos);
            yPos += 10;
            doc.setFont(undefined, 'normal');
            doc.text(`Rotatable Bonds: ${data.drug_likeness.veber.rotatable_bonds} / 10 ${data.drug_likeness.veber.rotatable_bonds <= 10 ? '✓' : '✗'}`, 25, yPos);
            yPos += 7;
            doc.text(`TPSA: ${data.drug_likeness.veber.tpsa.toFixed(1)} / 140 ${data.drug_likeness.veber.tpsa <= 140 ? '✓' : '✗'}`, 25, yPos);
            yPos += 15;
            
            // ADMET Summary
            if (data.drug_likeness.admet) {
                doc.setFont(undefined, 'bold');
                doc.text('ADMET Profile:', 20, yPos);
                yPos += 10;
                doc.setFont(undefined, 'normal');
                doc.text(`BBB Permeability: ${data.drug_likeness.bbb_permeability}`, 25, yPos);
                yPos += 7;
                doc.text(`Hepatotoxicity: ${data.drug_likeness.admet.toxicity.hepatotoxicity}`, 25, yPos);
                yPos += 7;
                doc.text(`CYP450 Inhibitor: ${data.drug_likeness.admet.metabolism.cyp450_inhibitor}`, 25, yPos);
            }
        }
    },
    {
        title: 'Market Analysis',
        content: (doc, data, autoTable) => {
            if (!data.marketAnalysis) return;
            
            doc.addPage();
            doc.setFontSize(18);
            doc.setTextColor(37, 99, 235);
            doc.text('Market Analysis', 20, 30);
            
            let yPos = 50;
            doc.setFontSize(11);
            doc.setTextColor(0, 0, 0);
            
            // Market Size
            doc.setFont(undefined, 'bold');
            doc.text('Market Size:', 20, yPos);
            yPos += 10;
            doc.setFont(undefined, 'normal');
            doc.text(`Therapeutic Area: ${data.marketAnalysis.marketSize.therapeuticArea}`, 25, yPos);
            yPos += 7;
            doc.text(`Estimated Market: $${data.marketAnalysis.marketSize.estimatedMarketSize}B ${data.marketAnalysis.marketSize.currency}`, 25, yPos);
            yPos += 7;
            doc.text(`Growth Rate: ${data.marketAnalysis.marketSize.growthRate}`, 25, yPos);
            yPos += 15;
            
            // Competitors
            if (data.marketAnalysis.competitors.length > 0) {
                doc.setFont(undefined, 'bold');
                doc.text('Competitor Drugs:', 20, yPos);
                yPos += 10;
                
                autoTable(doc, {
                    startY: yPos,
                    head: [['Drug Name', 'Indication', 'Similarity', 'Status', 'Company']],
                    body: data.marketAnalysis.competitors.map(c => [
                        c.name,
                        c.indication,
                        `${(c.similarity * 100).toFixed(1)}%`,
                        c.status,
                        c.company
                    ]),
                    theme: 'striped',
                    headStyles: { fillColor: [37, 99, 235] },
                });
            }
            
            // Patents
            if (data.marketAnalysis.patents.length > 0) {
                const lastY = (doc as any).lastAutoTable.finalY || yPos;
                doc.setFont(undefined, 'bold');
                doc.text('Relevant Patents:', 20, lastY + 20);
                
                autoTable(doc, {
                    startY: lastY + 30,
                    head: [['Patent Number', 'Title', 'Status', 'Relevance']],
                    body: data.marketAnalysis.patents.map(p => [
                        p.patentNumber,
                        p.title.substring(0, 40) + (p.title.length > 40 ? '...' : ''),
                        p.status,
                        p.relevance
                    ]),
                    theme: 'striped',
                    headStyles: { fillColor: [37, 99, 235] },
                });
            }
            
            // Regulatory Status
            const lastY = (doc as any).lastAutoTable?.finalY || yPos + 50;
            doc.setFont(undefined, 'bold');
            doc.text('Regulatory Status:', 20, lastY + 20);
            doc.setFont(undefined, 'normal');
            doc.text(data.marketAnalysis.regulatoryStatus, 25, lastY + 30);
        }
    }
    ];
}

