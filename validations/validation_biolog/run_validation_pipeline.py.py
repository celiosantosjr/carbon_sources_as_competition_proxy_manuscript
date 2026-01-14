"""
Simple runner script for the CaCo-BIOLOG validation pipeline.
"""

import subprocess
import sys
import os

def main():
    # Configuration
    genome_list = "data/genome_list.txt"  # Isolates genome list
    biolog_data = "data/csource.tsv"  # Isolates BIOLOG result
    output_dir = "validation_pipeline_results"
    
    print("=" * 60)
    print("CaCo-BIOLOG Validation Pipeline")
    print("=" * 60)
    print()
    
    # Check if files exist
    if not os.path.exists(genome_list):
        print(f"Error: {genome_list} not found.")
        print("Please create a genome_list.txt file with paths to your genome files.")
        sys.exit(1)
    
    if not os.path.exists(biolog_data):
        print(f"Error: {biolog_data} not found.")
        print("Please ensure your BIOLOG data file is named csource.tsv")
        sys.exit(1)
    
    # Run the pipeline
    cmd = [
        sys.executable, "caco_biolog_validation_pipeline.py",
        "-g", genome_list,
        "-b", biolog_data,
        "-o", output_dir
    ]
    
    print(f"Running: {' '.join(cmd)}")
    print()
    
    result = subprocess.run(cmd, capture_output=False, text=True)
    
    if result.returncode != 0:
        print("\nPipeline failed with error code:", result.returncode)
        sys.exit(1)
    
    print("\n" + "=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)
    print(f"All results saved to: {os.path.abspath(output_dir)}")
    print()
    print("Key files generated:")
    print(f"  1. {output_dir}/caco_biolog_validation_figure.png (main figure)")
    print(f"  2. {output_dir}/validation_summary_report.txt (detailed report)")
    print(f"  3. {output_dir}/genome_predictions/ (CaCo results from genomes)")
    print(f"  4. {output_dir}/biolog_predictions/ (RPS from BIOLOG data)")
    print(f"  5. {output_dir}/comparison_results/ (comparison data)")
    
    # Try to show the figure
    fig_path = os.path.join(output_dir, "caco_biolog_validation_figure.png")
    if os.path.exists(fig_path):
        print(f"\nFigure saved: {fig_path}")
        
        # Try to open on macOS
        if sys.platform == "darwin":
            try:
                subprocess.run(["open", fig_path])
            except:
                pass
        # Try to open on Windows
        elif sys.platform == "win32":
            try:
                os.startfile(fig_path)
            except:
                pass
        # Try to open on Linux
        else:
            try:
                subprocess.run(["xdg-open", fig_path])
            except:
                pass

if __name__ == "__main__":
    main()