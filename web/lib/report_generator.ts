import { MolecularAnalysisResult } from './molecular_inference';
import { getReportSections } from './report_template';

export async function generateExecutiveReport(analysis: MolecularAnalysisResult): Promise<Blob> {
    // Dynamic import to avoid SSR issues with jspdf
    const [{ default: jsPDF }, { default: autoTable }] = await Promise.all([
        import('jspdf'),
        import('jspdf-autotable')
    ]);
    
    return new Promise((resolve, reject) => {
        try {
            const doc = new jsPDF({
                orientation: 'portrait',
                unit: 'mm',
                format: 'a4'
            });

            // Set default font
            doc.setFont('helvetica');

            // Get report sections with the imported modules
            const reportSections = getReportSections(jsPDF, autoTable);

            // Generate each section
            reportSections.forEach((section, index) => {
                try {
                    section.content(doc, analysis, autoTable);
                } catch (error) {
                    console.error(`Error generating section "${section.title}":`, error);
                }
            });

            // Generate PDF blob
            const pdfBlob = doc.output('blob');
            resolve(pdfBlob);
        } catch (error) {
            reject(error);
        }
    });
}

export function downloadReport(analysis: MolecularAnalysisResult, filename?: string): Promise<void> {
    return new Promise((resolve, reject) => {
        generateExecutiveReport(analysis)
            .then(blob => {
                const url = URL.createObjectURL(blob);
                const link = document.createElement('a');
                link.href = url;
                
                // Generate filename
                const defaultFilename = filename || `MoleculeAI_Report_${analysis.smiles.replace(/[^a-zA-Z0-9]/g, '_')}_${Date.now()}.pdf`;
                link.download = defaultFilename;
                
                document.body.appendChild(link);
                link.click();
                document.body.removeChild(link);
                
                // Clean up
                setTimeout(() => URL.revokeObjectURL(url), 100);
                
                resolve();
            })
            .catch(reject);
    });
}

