import myvariant
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sys
import warnings
import os

# Suppress warnings
warnings.filterwarnings('ignore')

class ClinVarTool:
    def __init__(self):
        self.mv = myvariant.MyVariantInfo()
        self.current_df = None

    def fetch_data(self, gene_symbol, size=400):
        """Fetch real-time variant data from MyVariant.info (ClinVar)."""
        print(f"\n[Fetching] Retrieving up to {size} variants for gene: {gene_symbol.upper()}...")
        try:
            # Query MyVariant.info for ClinVar data
            results = self.mv.query(f'clinvar.gene.symbol:{gene_symbol}', fields='clinvar', size=size)
            
            if 'hits' not in results or not results['hits']:
                print(f"[Error] No data found for gene: {gene_symbol}")
                return None
            
            data = []
            for hit in results['hits']:
                clinvar = hit.get('clinvar', {})
                # Extract clinical significance - can be complex in real data
                sig = "Unknown"
                rcv = clinvar.get('rcv')
                if isinstance(rcv, list) and rcv:
                    sig = rcv[0].get('clinical_significance', 'Unknown')
                elif isinstance(rcv, dict):
                    sig = rcv.get('clinical_significance', 'Unknown')
                
                # Extract position and chromosome
                hg19 = clinvar.get('hg19', {})
                pos = hg19.get('start', 'Unknown')
                chrom = clinvar.get('chrom', 'Unknown')
                
                data.append({
                    'VariantID': hit.get('_id'),
                    'Gene': gene_symbol.upper(),
                    'Chromosome': chrom,
                    'Position': pos,
                    'ClinicalSignificance': sig
                })
            
            self.current_df = pd.DataFrame(data)
            print(f"[Success] Loaded {len(self.current_df)} real variants.")
            return self.current_df
            
        except Exception as e:
            print(f"[Error] Failed to fetch data: {e}")
            return None

    def plot_interactive(self, df, gene):
        """Create an interactive plot using Plotly."""
        if df is None or df.empty:
            print("[Error] No data to plot.")
            return

        print(f"[Plotting] Generating interactive visualization for {gene}...")
        
        # Create a subplot with two charts
        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=("Clinical Significance Distribution", "Variant Positions on Chromosome"),
            specs=[[{"type": "bar"}, {"type": "scatter"}]]
        )

        # 1. Bar chart for Clinical Significance
        sig_counts = df['ClinicalSignificance'].value_counts().reset_index()
        sig_counts.columns = ['Significance', 'Count']
        
        fig.add_trace(
            go.Bar(
                x=sig_counts['Significance'],
                y=sig_counts['Count'],
                name="Significance",
                marker_color='indianred'
            ),
            row=1, col=1
        )

        # 2. Scatter plot for positions
        # Filter out 'Unknown' positions for plotting
        plot_df = df[df['Position'] != 'Unknown'].copy()
        plot_df['Position'] = pd.to_numeric(plot_df['Position'])
        
        fig.add_trace(
            go.Scatter(
                x=plot_df['Position'],
                y=[1] * len(plot_df),
                mode='markers',
                marker=dict(size=10, color=plot_df.index, colorscale='Viridis'),
                text=plot_df['VariantID'],
                hovertemplate="<b>Variant:</b> %{text}<br><b>Position:</b> %{x}<br><b>Significance:</b> " + plot_df['ClinicalSignificance'],
                name="Positions"
            ),
            row=1, col=2
        )

        fig.update_layout(
            title_text=f"Interactive ClinVar Analysis: {gene}",
            showlegend=False,
            height=500,
            template="plotly_white"
        )
        
        # In a sandbox environment, we save to HTML and provide the file/URL
        output_file = f"interactive_analysis_{gene}.html"
        fig.write_html(output_file)
        print(f"[Success] Interactive graph saved as: {output_file}")
        return output_file

def main():
    tool = ClinVarTool()
    
    print("\n" + "="*50)
    print("   REAL-TIME INTERACTIVE CLINVAR STUDY TOOL")
    print("="*50)
    print("This tool connects to MyVariant.info (NCBI/ClinVar) for real data.")
    
    while True:
        try:
            gene = input("\nEnter Gene Symbol (e.g., TP53, BRCA1, CFTR) or 'exit': ").strip()
            if gene.lower() == 'exit':
                print("Exiting. Goodbye!")
                break
            
            if not gene:
                continue
                
            df = tool.fetch_data(gene)
            if df is not None:
                print("\nData Preview (First 5 rows):")
                print(df.head())
                
                html_file = tool.plot_interactive(df, gene.upper())
                print(f"\n[Action] Please open '{html_file}' in your browser to interact with the graph.")
                
        except KeyboardInterrupt:
            break
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()
