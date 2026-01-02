"""
AMR Thesis Project - Central Pipeline Orchestrator
===================================================

This is the SINGLE SOURCE OF TRUTH for running all pipeline operations.
All commands go through main.py - never run scripts directly.

Usage:
    python main.py --pipeline          # Run core data pipeline
    python main.py --validate          # Run component validation scripts
    python main.py --sensitivity       # Run methodology sensitivity analysis
    python main.py --analyze           # Run analysis modules
    python main.py --viz               # Regenerate visualizations
    python main.py --app               # Launch Streamlit dashboard
    python main.py --all               # Run everything in sequence

Developed by: Queshue
Date: 2024-06-15 -> 2025-12-28
"""

import os
import sys
import argparse
import subprocess
import logging
from pathlib import Path
from datetime import datetime
import time
from contextlib import redirect_stdout
import io

# =============================================================================
# PROJECT SETUP
# =============================================================================
PROJECT_ROOT = Path(__file__).parent.resolve()
sys.path.insert(0, str(PROJECT_ROOT / 'src'))
sys.path.insert(0, str(PROJECT_ROOT / 'scripts'))

# Import config for paths
from config import (
    PROCESSED_DATA_DIR, FIGURES_DIR, MODELS_DIR, ARTIFACTS_DIR, RAW_DATA_DIR,
    RANDOM_STATE, MIN_ANTIBIOTIC_COVERAGE, MAX_ISOLATE_MISSING
)


# =============================================================================
# TERMINAL COLOR & STYLING
# =============================================================================
class Colors:
    """ANSI color codes for terminal output with Windows fallback."""
    
    # Check if terminal supports colors
    ENABLED = True
    
    @classmethod
    def init(cls):
        """Initialize color support, especially for Windows."""
        if sys.platform == 'win32':
            try:
                # Enable ANSI escape codes on Windows 10+
                import ctypes
                kernel32 = ctypes.windll.kernel32
                kernel32.SetConsoleMode(kernel32.GetStdHandle(-11), 7)
            except Exception:
                cls.ENABLED = False
    
    @classmethod
    def colorize(cls, text: str, *codes) -> str:
        """Apply color codes to text if colors are enabled."""
        if not cls.ENABLED:
            return text
        code_str = ';'.join(str(c) for c in codes)
        return f"\033[{code_str}m{text}\033[0m"
    
    # Colors
    BLACK = 30
    RED = 31
    GREEN = 32
    YELLOW = 33
    BLUE = 34
    MAGENTA = 35
    CYAN = 36
    WHITE = 37
    
    # Bright colors
    BRIGHT_BLACK = 90
    BRIGHT_RED = 91
    BRIGHT_GREEN = 92
    BRIGHT_YELLOW = 93
    BRIGHT_BLUE = 94
    BRIGHT_MAGENTA = 95
    BRIGHT_CYAN = 96
    BRIGHT_WHITE = 97
    
    # Styles
    BOLD = 1
    DIM = 2
    ITALIC = 3
    UNDERLINE = 4


# Initialize colors
Colors.init()


# =============================================================================
# ASCII ART & VISUAL ELEMENTS
# =============================================================================
BANNER = """
╔═══════════════════════════════════════════════════════════════════════╗
║                                                                       ║
║      █████╗ ███╗   ███╗██████╗     ██████╗ ██╗██████╗ ███████╗        ║
║     ██╔══██╗████╗ ████║██╔══██╗    ██╔══██╗██║██╔══██╗██╔════╝        ║
║     ███████║██╔████╔██║██████╔╝    ██████╔╝██║██████╔╝█████╗          ║
║     ██╔══██║██║╚██╔╝██║██╔══██╗    ██╔═══╝ ██║██╔═══╝ ██╔══╝          ║
║     ██║  ██║██║ ╚═╝ ██║██║  ██║    ██║     ██║██║     ███████╗        ║
║     ╚═╝  ╚═╝╚═╝     ╚═╝╚═╝  ╚═╝    ╚═╝     ╚═╝╚═╝     ╚══════╝        ║
║                                                                       ║
║     Antimicrobial Resistance Pattern Recognition Pipeline v2.0        ║
║     ─────────────────────────────────────────────────────────         ║
║     Thesis Project | Data-Driven AMR Surveillance                     ║
║                                                                       ║
╚═══════════════════════════════════════════════════════════════════════╝
"""

MINI_BANNER = """
┌───────────────────────────────────────────────────────────────────────┐
│  AMR PIPELINE  ║  Antimicrobial Resistance Pattern Recognition        │
└───────────────────────────────────────────────────────────────────────┘
"""


# =============================================================================
# ENHANCED LOGGING SYSTEM
# =============================================================================
class Console:
    """Enhanced console output with colors, formatting, and timing."""
    
    # Unicode symbols (with fallbacks)
    ICONS = {
        'success': '✓',
        'error': '✗',
        'warning': '⚠',
        'info': '→',
        'bullet': '•',
        'arrow': '▶',
        'check': '✔',
        'cross': '✘',
        'star': '★',
        'circle': '●',
        'square': '■',
        'diamond': '◆',
        'progress': '░▒▓█',
    }
    
    def __init__(self):
        self.phase_start_time = None
        self.total_start_time = None
        self.current_step = 0
        self.total_steps = 0
        self.phase_name = ""
        self.step_times = []
        self.phase_results = []
        
        # Configure base logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(message)s',
            handlers=[logging.StreamHandler(sys.stdout)]
        )
        self.logger = logging.getLogger('amr_pipeline')
    
    def _format_time(self, seconds: float) -> str:
        """Format time duration to human-readable string."""
        if seconds < 60:
            return f"{seconds:.1f}s"
        elif seconds < 3600:
            mins = int(seconds // 60)
            secs = seconds % 60
            return f"{mins}m {secs:.1f}s"
        else:
            hours = int(seconds // 3600)
            mins = int((seconds % 3600) // 60)
            return f"{hours}h {mins}m"
    
    def _create_box(self, text: str, width: int = 75, style: str = 'double') -> str:
        """Create a boxed text element."""
        chars = {
            'double': {'tl': '╔', 'tr': '╗', 'bl': '╚', 'br': '╝', 'h': '═', 'v': '║'},
            'single': {'tl': '┌', 'tr': '┐', 'bl': '└', 'br': '┘', 'h': '─', 'v': '│'},
            'heavy':  {'tl': '┏', 'tr': '┓', 'bl': '┗', 'br': '┛', 'h': '━', 'v': '┃'},
            'round':  {'tl': '╭', 'tr': '╮', 'bl': '╰', 'br': '╯', 'h': '─', 'v': '│'},
        }
        c = chars.get(style, chars['double'])
        inner_width = width - 2
        
        lines = []
        lines.append(f"{c['tl']}{c['h'] * inner_width}{c['tr']}")
        
        # Handle multi-line text
        for line in text.split('\n'):
            padded = line.center(inner_width)
            lines.append(f"{c['v']}{padded}{c['v']}")
        
        lines.append(f"{c['bl']}{c['h'] * inner_width}{c['br']}")
        return '\n'.join(lines)
    
    def _create_progress_bar(self, current: int, total: int, width: int = 30) -> str:
        """Create a visual progress bar."""
        if total == 0:
            return ""
        progress = current / total
        filled = int(width * progress)
        empty = width - filled
        
        bar = Colors.colorize('█' * filled, Colors.GREEN) + \
              Colors.colorize('░' * empty, Colors.BRIGHT_BLACK)
        percentage = Colors.colorize(f"{progress * 100:5.1f}%", Colors.CYAN, Colors.BOLD)
        return f"[{bar}] {percentage}"
    
    def print_banner(self, mini: bool = False):
        """Print the project banner."""
        banner = MINI_BANNER if mini else BANNER
        print(Colors.colorize(banner, Colors.CYAN, Colors.BOLD))
    
    def print_header(self, title: str, subtitle: str = None):
        """Print a styled header section."""
        print()
        header = self._create_box(title, style='double')
        print(Colors.colorize(header, Colors.BRIGHT_CYAN, Colors.BOLD))
        if subtitle:
            print(Colors.colorize(f"  {self.ICONS['info']} {subtitle}", Colors.DIM))
        print()
    
    def set_phase(self, phase: str, total_steps: int, description: str = None):
        """Start a new phase with timing."""
        self.phase_name = phase
        self.total_steps = total_steps
        self.current_step = 0
        self.phase_start_time = time.time()
        self.step_times = []
        
        # Phase header
        print()
        print(Colors.colorize('━' * 75, Colors.BRIGHT_BLACK))
        
        # Phase title with icon
        icon = self.ICONS['arrow']
        phase_text = f"  {icon}  PHASE: {phase}"
        print(Colors.colorize(phase_text, Colors.BRIGHT_YELLOW, Colors.BOLD))
        
        if description:
            print(Colors.colorize(f"     {description}", Colors.WHITE, Colors.DIM))
        
        print(Colors.colorize('━' * 75, Colors.BRIGHT_BLACK))
        print()
    
    def step(self, message: str, detail: str = None):
        """Log a step with progress tracking."""
        step_start = time.time()
        self.current_step += 1
        
        # Progress bar
        progress_bar = self._create_progress_bar(self.current_step, self.total_steps)
        
        # Step indicator
        step_num = Colors.colorize(f"[{self.current_step}/{self.total_steps}]", Colors.CYAN, Colors.BOLD)
        step_msg = Colors.colorize(message, Colors.WHITE, Colors.BOLD)
        
        print(f"  {step_num} {step_msg}")
        print(f"       {progress_bar}")
        
        if detail:
            print(Colors.colorize(f"       {self.ICONS['info']} {detail}", Colors.DIM))
        
        return step_start
    
    def step_complete(self, start_time: float, result: str):
        """Mark a step as complete with timing."""
        elapsed = time.time() - start_time
        self.step_times.append(elapsed)
        
        time_str = Colors.colorize(f"({self._format_time(elapsed)})", Colors.BRIGHT_BLACK)
        result_colored = Colors.colorize(f"       {self.ICONS['success']} {result}", Colors.GREEN)
        print(f"{result_colored} {time_str}")
        print()
    
    def info(self, message: str):
        """Log info message."""
        icon = Colors.colorize(self.ICONS['info'], Colors.BLUE)
        print(f"       {icon} {message}")
    
    def success(self, message: str):
        """Log success message."""
        icon = Colors.colorize(self.ICONS['success'], Colors.GREEN, Colors.BOLD)
        msg = Colors.colorize(message, Colors.GREEN)
        print(f"       {icon} {msg}")
    
    def warning(self, message: str):
        """Log warning message."""
        icon = Colors.colorize(self.ICONS['warning'], Colors.YELLOW, Colors.BOLD)
        msg = Colors.colorize(message, Colors.YELLOW)
        print(f"       {icon} {msg}")
    
    def error(self, message: str):
        """Log error message."""
        icon = Colors.colorize(self.ICONS['error'], Colors.RED, Colors.BOLD)
        msg = Colors.colorize(message, Colors.RED)
        print(f"       {icon} {msg}")
    
    def phase_complete(self, additional_info: dict = None):
        """Complete a phase with summary."""
        phase_elapsed = time.time() - self.phase_start_time
        
        print()
        print(Colors.colorize('  ┌' + '─' * 71 + '┐', Colors.GREEN))
        
        # Phase complete message
        complete_msg = f"  │  {self.ICONS['check']} {self.phase_name} Complete"
        complete_msg = complete_msg.ljust(74) + '│'
        print(Colors.colorize(complete_msg, Colors.GREEN, Colors.BOLD))
        
        # Timing info
        time_msg = f"  │    Time: {self._format_time(phase_elapsed)}"
        time_msg = time_msg.ljust(74) + '│'
        print(Colors.colorize(time_msg, Colors.GREEN))
        
        # Additional info
        if additional_info:
            for key, value in additional_info.items():
                info_msg = f"  │    {key}: {value}"
                info_msg = info_msg.ljust(74) + '│'
                print(Colors.colorize(info_msg, Colors.GREEN))
        
        print(Colors.colorize('  └' + '─' * 71 + '┘', Colors.GREEN))
        print()
        
        # Store result for final summary
        self.phase_results.append({
            'phase': self.phase_name,
            'time': phase_elapsed,
            'steps': self.current_step,
            'info': additional_info
        })
    
    def print_summary_table(self, title: str, data: list, columns: list):
        """Print a formatted summary table."""
        # Calculate column widths
        widths = [len(col) for col in columns]
        for row in data:
            for i, cell in enumerate(row):
                widths[i] = max(widths[i], len(str(cell)))
        
        total_width = sum(widths) + len(columns) * 3 + 1
        
        # Header
        print()
        print(Colors.colorize(f"  {title}", Colors.BRIGHT_WHITE, Colors.BOLD))
        print(Colors.colorize('  ┌' + '─' * (total_width - 2) + '┐', Colors.BRIGHT_BLACK))
        
        # Column headers
        header = '  │'
        for col, width in zip(columns, widths):
            header += f' {col.center(width)} │'
        print(Colors.colorize(header, Colors.BRIGHT_CYAN, Colors.BOLD))
        
        # Separator
        sep = '  ├'
        for width in widths:
            sep += '─' * (width + 2) + '┼'
        sep = sep[:-1] + '┤'
        print(Colors.colorize(sep, Colors.BRIGHT_BLACK))
        
        # Data rows
        for row in data:
            row_str = '  │'
            for cell, width in zip(row, widths):
                row_str += f' {str(cell).ljust(width)} │'
            print(Colors.colorize(row_str, Colors.WHITE))
        
        # Footer
        print(Colors.colorize('  └' + '─' * (total_width - 2) + '┘', Colors.BRIGHT_BLACK))
    
    def print_final_summary(self, output_paths: dict = None):
        """Print a comprehensive final summary."""
        if self.total_start_time:
            total_elapsed = time.time() - self.total_start_time
        else:
            total_elapsed = sum(r['time'] for r in self.phase_results)
        
        print()
        print(Colors.colorize('═' * 75, Colors.BRIGHT_GREEN))
        
        # Success banner
        success_box = """
    ╔═════════════════════════════════════════════════════════════════╗
    ║                                                                 ║
    ║              ██████╗  ██████╗ ███╗   ██╗███████╗██╗             ║
    ║              ██╔══██╗██╔═══██╗████╗  ██║██╔════╝██║             ║
    ║              ██║  ██║██║   ██║██╔██╗ ██║█████╗  ██║             ║
    ║              ██║  ██║██║   ██║██║╚██╗██║██╔══╝  ╚═╝             ║
    ║              ██████╔╝╚██████╔╝██║ ╚████║███████╗██╗             ║
    ║              ╚═════╝  ╚═════╝ ╚═╝  ╚═══╝╚══════╝╚═╝             ║
    ║                                                                 ║
    ║            Pipeline execution completed successfully!           ║
    ║                                                                 ║
    ╚═════════════════════════════════════════════════════════════════╝
"""
        print(Colors.colorize(success_box, Colors.GREEN, Colors.BOLD))
        
        # Phase summary table
        if self.phase_results:
            table_data = []
            for result in self.phase_results:
                table_data.append([
                    result['phase'],
                    str(result['steps']),
                    self._format_time(result['time']),
                    '✓ Success'
                ])
            
            self.print_summary_table(
                'Execution Summary',
                table_data,
                ['Phase', 'Steps', 'Time', 'Status']
            )
        
        # Output paths
        if output_paths:
            print()
            print(Colors.colorize('  Output Locations:', Colors.BRIGHT_WHITE, Colors.BOLD))
            print(Colors.colorize('  ─' * 36, Colors.BRIGHT_BLACK))
            for name, path in output_paths.items():
                icon = self.ICONS['bullet']
                print(Colors.colorize(f"    {icon} {name}:", Colors.CYAN), end=' ')
                print(Colors.colorize(str(path), Colors.WHITE))
        
        # Total time
        print()
        total_box = f"""
  ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
  ┃  Total Execution Time: {self._format_time(total_elapsed).center(47)}  ┃
  ┃  Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S').center(58)} ┃
  ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
"""
        print(Colors.colorize(total_box, Colors.BRIGHT_GREEN, Colors.BOLD))
        
        # Next steps hint
        print(Colors.colorize('  💡 Next Steps:', Colors.BRIGHT_YELLOW, Colors.BOLD))
        print(Colors.colorize('     • Launch the dashboard: python main.py --app', Colors.WHITE))
        print(Colors.colorize('     • View figures in: data/processed/figures/', Colors.WHITE))
        print()
        print(Colors.colorize('═' * 75, Colors.BRIGHT_GREEN))
        print()


# Global console instance
console = Console()

# Legacy logger alias for backward compatibility
log = console


# =============================================================================
# PIPELINE: Core Data Flow
# =============================================================================
def run_pipeline(k_override=None):
    """
    Run the core data pipeline: Ingestion → Cleaning → Encoding → Clustering.
    
    Parameters:
    -----------
    k_override : int, optional
        If provided, use this k value instead of auto-detection.
    
    This is the main data processing workflow that transforms raw CSV files
    into clustered, analysis-ready datasets.
    """
    console.set_phase("DATA PIPELINE", 6, "Core data processing workflow")
    
    # Suppress verbose output from submodules for clean pipeline display
    # We use a null buffer to capture and discard module print statements
    devnull = io.StringIO()
    
    # Ensure output directories exist
    os.makedirs(PROCESSED_DATA_DIR, exist_ok=True)
    os.makedirs(FIGURES_DIR, exist_ok=True)
    os.makedirs(ARTIFACTS_DIR, exist_ok=True)
    os.makedirs(MODELS_DIR, exist_ok=True)
    
    output_dir = str(PROCESSED_DATA_DIR)
    data_dir = str(RAW_DATA_DIR)  # Raw CSVs are in data/raw/
    
    pipeline_stats = {}
    
    # Step 1: Data Ingestion
    step_time = console.step("Data Ingestion", "Loading and unifying raw CSV files")
    from preprocessing.data_ingestion import create_unified_dataset
    unified_path = os.path.join(output_dir, 'unified_raw_dataset.csv')
    with redirect_stdout(devnull):
        df_raw = create_unified_dataset(data_dir, unified_path)
    
    if df_raw.empty:
        console.error("No data loaded. Check your CSV files in project root.")
        return None
    pipeline_stats['raw_isolates'] = len(df_raw)
    console.step_complete(step_time, f"Unified dataset: {len(df_raw):,} isolates")
    
    # Step 2: Data Cleaning
    step_time = console.step("Data Cleaning", "Applying quality thresholds")
    from preprocessing.data_cleaning import clean_dataset, generate_cleaning_report
    with redirect_stdout(devnull):
        df_clean, cleaning_report = clean_dataset(
            df_raw,
            min_antibiotic_coverage=MIN_ANTIBIOTIC_COVERAGE,
            max_isolate_missing=MAX_ISOLATE_MISSING
        )
    clean_path = os.path.join(output_dir, 'cleaned_dataset.csv')
    df_clean.to_csv(clean_path, index=False)
    
    report_path = os.path.join(output_dir, 'cleaning_report.txt')
    generate_cleaning_report(cleaning_report, report_path)
    pipeline_stats['clean_isolates'] = len(df_clean)
    console.step_complete(step_time, f"Cleaned dataset: {len(df_clean):,} isolates retained")
    
    # Step 3: Resistance Encoding
    step_time = console.step("Resistance Encoding", "Converting S/I/R to numeric values")
    from preprocessing.resistance_encoding import create_encoded_dataset
    with redirect_stdout(devnull):
        df_encoded, encoding_info = create_encoded_dataset(df_clean)
    encoded_path = os.path.join(output_dir, 'encoded_dataset.csv')
    df_encoded.to_csv(encoded_path, index=False)
    pipeline_stats['antibiotics'] = len(encoding_info['encoded_columns'])
    console.step_complete(step_time, f"Encoded {len(encoding_info['encoded_columns'])} antibiotics")
    
    # Step 4: Feature Engineering
    step_time = console.step("Feature Engineering", "Preparing analysis-ready dataset")
    from preprocessing.feature_engineering import prepare_analysis_ready_dataset
    encoded_cols = encoding_info['encoded_columns']
    analysis_path = os.path.join(output_dir, 'analysis_ready_dataset.csv')
    with redirect_stdout(devnull):
        df_analysis, feature_matrix, metadata, feature_info = prepare_analysis_ready_dataset(
            df_encoded, encoded_cols, output_path=analysis_path
        )
    pipeline_stats['features'] = feature_matrix.shape[1]
    console.step_complete(step_time, f"Feature matrix: {feature_matrix.shape[0]:,} × {feature_matrix.shape[1]} features")
    
    # Step 5: Optimal k Selection
    step_time = console.step("Cluster Validation", "Finding optimal k using elbow + silhouette method")
    from validate_clustering import find_optimal_k, create_validation_plot
    import pandas as pd
    
    feature_cols = [c for c in df_analysis.columns if c.endswith('_encoded')]
    k_result = find_optimal_k(
        df_analysis, 
        feature_cols, 
        k_range=range(2, 11),
        min_k=3, 
        max_k=8,
        method='combined'
    )
    
    # Allow k override from CLI
    if k_override is not None:
        console.info(f"CLI override: k={k_override} (auto-detected: k={k_result['optimal_k']})")
        optimal_k = k_override
    else:
        optimal_k = k_result['optimal_k']
    pipeline_stats['optimal_k'] = optimal_k
    pipeline_stats['silhouette'] = f"{k_result['silhouette_score']:.4f}"
    console.step_complete(step_time, f"Optimal k={optimal_k} (Silhouette: {k_result['silhouette_score']:.4f})")
    
    # Create validation plot
    results_for_plot = pd.DataFrame({
        'k': list(k_result['all_silhouette_scores'].keys()),
        'silhouette_score': list(k_result['all_silhouette_scores'].values()),
        'wcss': list(k_result['all_wcss_scores'].values())
    })
    create_validation_plot(
        results_for_plot,
        os.path.join(str(FIGURES_DIR), 'cluster_validation.png'),
        selected_k=optimal_k,
        elbow_k=k_result['elbow_k']
    )
    
    # Step 6: Hierarchical Clustering
    step_time = console.step("Hierarchical Clustering", "Generating cluster assignments")
    from clustering.hierarchical_clustering import run_clustering_pipeline
    
    df_clustered, linkage_matrix, clustering_info = run_clustering_pipeline(
        df_analysis, 
        feature_cols, 
        n_clusters=optimal_k,
        perform_robustness=True,
        output_dir=str(ARTIFACTS_DIR)
    )
    
    clustered_path = os.path.join(output_dir, 'clustered_dataset.csv')
    df_clustered.to_csv(clustered_path, index=False)
    pipeline_stats['clusters'] = optimal_k
    console.step_complete(step_time, f"Clustered dataset saved with {optimal_k} clusters")
    
    # Phase complete with summary
    console.phase_complete({
        'Isolates Processed': f"{pipeline_stats['clean_isolates']:,}",
        'Antibiotics Encoded': pipeline_stats['antibiotics'],
        'Clusters Identified': pipeline_stats['clusters']
    })
    
    return df_clustered, linkage_matrix, feature_cols, clustering_info


# =============================================================================
# VALIDATE: Run Validation Scripts
# =============================================================================
def run_validate():
    """
    Execute validation and analysis scripts from scripts/ directory.
    
    Runs:
    - validate_clustering.py: Cluster quality assessment
    - coresistance_analysis.py: Co-resistance network analysis
    - antibiotic_clustering.py: Antibiotic clustering by co-resistance patterns
    """
    console.set_phase("VALIDATION SCRIPTS", 3, "Quality assessment and network analysis")
    
    # Suppress verbose output from submodules
    devnull = io.StringIO()
    
    validation_results = {'cluster': False, 'coresistance': False, 'antibiotic': False}
    
    # Step 1: Cluster Validation
    step_time = console.step("Cluster Validation", "Assessing cluster quality metrics")
    try:
        from validate_clustering import main as validate_main
        with redirect_stdout(devnull):
            validate_main()
        validation_results['cluster'] = True
        console.step_complete(step_time, "Cluster validation complete")
    except Exception as e:
        console.error(f"Cluster validation failed: {e}")
    
    # Step 2: Co-Resistance Analysis
    step_time = console.step("Co-Resistance Analysis", "Building resistance network")
    try:
        from coresistance_analysis import main as coresistance_main
        with redirect_stdout(devnull):
            coresistance_main()
        validation_results['coresistance'] = True
        console.step_complete(step_time, "Co-resistance network analysis complete")
    except Exception as e:
        console.error(f"Co-resistance analysis failed: {e}")
    
    # Step 3: Antibiotic Clustering
    step_time = console.step("Antibiotic Clustering", "Clustering antibiotics by resistance patterns")
    try:
        from antibiotic_clustering import main as antibiotic_clustering_main
        with redirect_stdout(devnull):
            antibiotic_clustering_main(args=[])
        validation_results['antibiotic'] = True
        console.step_complete(step_time, "Antibiotic clustering complete")
    except Exception as e:
        console.error(f"Antibiotic clustering failed: {e}")
    
    # Phase complete
    success_count = sum(validation_results.values())
    console.phase_complete({'Analyses Completed': f"{success_count}/3"})


# =============================================================================
# SENSITIVITY: Run Methodology Sensitivity Analysis
# =============================================================================
def run_sensitivity():
    """
    Execute sensitivity analysis scripts to validate methodology robustness.
    
    Runs:
    - encoding_sensitivity.py: Tests impact of different resistance encoding schemes
    - threshold_sensitivity.py: Tests impact of data cleaning thresholds
    """
    console.set_phase("SENSITIVITY ANALYSIS", 2, "Methodological robustness checks")
    
    sensitivity_results = {'encoding': False, 'threshold': False}
    
    # Step 1: Encoding Sensitivity
    step_time = console.step("Encoding Sensitivity", "Testing simplified vs linear vs exponential encoding")
    try:
        from scripts.encoding_sensitivity import main as encoding_main
        encoding_main(args=[])
        sensitivity_results['encoding'] = True
        console.step_complete(step_time, "Encoding sensitivity analysis complete")
    except Exception as e:
        console.error(f"Encoding sensitivity analysis failed: {e}")
    
    # Step 2: Threshold Sensitivity
    step_time = console.step("Threshold Sensitivity", "Testing data cleaning threshold impact")
    try:
        from scripts.threshold_sensitivity import main as threshold_main
        threshold_main(args=[])
        sensitivity_results['threshold'] = True
        console.step_complete(step_time, "Threshold sensitivity analysis complete")
    except Exception as e:
        console.error(f"Threshold sensitivity analysis failed: {e}")
    
    # Phase complete
    success_count = sum(sensitivity_results.values())
    console.phase_complete({'Validations Completed': f"{success_count}/2"})


# =============================================================================
# ANALYZE: Run Analysis Modules
# =============================================================================
def run_analyze():
    """
    Run the analysis modules from src/analysis/.
    
    Runs:
    - Regional and environmental distribution analysis
    - Integration and synthesis analysis
    - Species discrimination (supervised learning)
    """
    console.set_phase("ANALYSIS MODULES", 3, "Statistical and ML-based analysis")
    
    # Suppress verbose output from submodules
    devnull = io.StringIO()
    
    import pandas as pd
    
    # Load clustered dataset
    clustered_path = PROCESSED_DATA_DIR / 'clustered_dataset.csv'
    if not clustered_path.exists():
        console.error(f"Clustered dataset not found at {clustered_path}")
        console.info("Run --pipeline first to generate the clustered dataset.")
        return
    
    df_clustered = pd.read_csv(clustered_path)
    feature_cols = [c for c in df_clustered.columns if c.endswith('_encoded')]
    
    analysis_results = {'regional': False, 'integration': False, 'species': False}
    
    # Step 1: Regional & Environmental Analysis
    step_time = console.step("Regional & Environmental Analysis", "Analyzing geographic distribution patterns")
    try:
        from analysis.regional_environmental import run_regional_environmental_analysis
        with redirect_stdout(devnull):
            regional_results = run_regional_environmental_analysis(
                df_clustered, feature_cols, str(FIGURES_DIR)
            )
        analysis_results['regional'] = True
        console.step_complete(step_time, "Regional/environmental analysis complete")
    except Exception as e:
        console.error(f"Regional analysis failed: {e}")
    
    # Step 2: Integration & Synthesis
    step_time = console.step("Integration & Synthesis", "Cross-domain pattern synthesis")
    try:
        from analysis.integration_synthesis import run_integration_synthesis
        with redirect_stdout(devnull):
            integration_results = run_integration_synthesis(
                df_clustered, feature_cols, supervised_results=None
            )
        analysis_results['integration'] = True
        console.step_complete(step_time, "Integration/synthesis complete")
    except Exception as e:
        console.error(f"Integration synthesis failed: {e}")
    
    # Step 3: Supervised Learning (Species Discrimination)
    step_time = console.step("Species Discrimination", "Training species classifier with ML")
    try:
        from supervised.supervised_learning import run_species_discrimination, save_model
        
        if 'ISOLATE_ID' in df_clustered.columns and df_clustered['ISOLATE_ID'].nunique() > 1:
            with redirect_stdout(devnull):
                species_results = run_species_discrimination(df_clustered, feature_cols)
            
            # Save model
            with redirect_stdout(devnull):
                save_model(
                    species_results['best_model']['model_object'],
                    species_results['scaler'],
                    species_results['label_encoder'],
                    str(MODELS_DIR / 'species_classifier.joblib'),
                    imputer=species_results.get('imputer'),
                    preprocessing_info=species_results.get('preprocessing_info')
                )
            analysis_results['species'] = True
            console.step_complete(step_time, "Species classifier trained and saved")
        else:
            console.warning("Insufficient species diversity for discrimination analysis")
    except Exception as e:
        console.error(f"Species discrimination failed: {e}")
    
    # Phase complete
    success_count = sum(analysis_results.values())
    console.phase_complete({'Analyses Completed': f"{success_count}/3"})


# =============================================================================
# VIZ: Regenerate Visualizations
# =============================================================================
def run_viz():
    """
    Regenerate all visualizations from the clustered dataset.
    
    Outputs are saved to data/processed/figures/.
    """
    console.set_phase("VISUALIZATION", 1, "Generating publication-ready figures")
    
    # Suppress verbose output from submodules
    devnull = io.StringIO()
    
    import pandas as pd
    import pickle
    
    # Load clustered dataset
    clustered_path = PROCESSED_DATA_DIR / 'clustered_dataset.csv'
    if not clustered_path.exists():
        console.error(f"Clustered dataset not found at {clustered_path}")
        console.info("Run --pipeline first to generate the clustered dataset.")
        return
    
    # Load linkage matrix
    linkage_path = ARTIFACTS_DIR / 'linkage_matrix.pkl'
    clustering_info_path = ARTIFACTS_DIR / 'clustering_info.pkl'
    
    if not linkage_path.exists():
        console.error(f"Linkage matrix not found at {linkage_path}")
        console.info("Run --pipeline first to generate clustering artifacts.")
        return
    
    step_time = console.step("Generate Visualizations", "Creating dendrograms, heatmaps, and charts")
    
    df_clustered = pd.read_csv(clustered_path)
    feature_cols = [c for c in df_clustered.columns if c.endswith('_encoded')]
    
    with open(linkage_path, 'rb') as f:
        linkage_matrix = pickle.load(f)
    
    clustering_info = None
    if clustering_info_path.exists():
        with open(clustering_info_path, 'rb') as f:
            clustering_info = pickle.load(f)
    
    from visualization.visualization import generate_all_visualizations
    with redirect_stdout(devnull):
        generate_all_visualizations(
            df_clustered, 
            feature_cols, 
            linkage_matrix, 
            str(FIGURES_DIR),
            clustering_info
        )
    
    console.step_complete(step_time, f"Visualizations saved to {FIGURES_DIR}")
    console.phase_complete({'Output Directory': str(FIGURES_DIR)})


# =============================================================================
# APP: Launch Streamlit Dashboard
# =============================================================================
def run_app():
    """
    Launch the Streamlit interactive dashboard.
    """
    console.print_banner(mini=True)
    console.set_phase("STREAMLIT DASHBOARD", 1, "Interactive data exploration")
    
    step_time = console.step("Launch Dashboard", "Starting Streamlit server")
    
    app_path = PROJECT_ROOT / 'app' / 'streamlit_app.py'
    if not app_path.exists():
        console.error(f"Streamlit app not found at {app_path}")
        return
    
    # Use virtual environment's streamlit if available
    venv_streamlit = PROJECT_ROOT / 'venv' / 'Scripts' / 'streamlit.exe'
    if venv_streamlit.exists():
        streamlit_cmd = str(venv_streamlit)
    else:
        streamlit_cmd = "streamlit"  # Fall back to system PATH
    
    console.info(f"Command: {streamlit_cmd} run {app_path}")
    print()
    print(Colors.colorize('  ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓', Colors.CYAN))
    print(Colors.colorize('  ┃  🌐 Dashboard launching... Open your browser to view.                ┃', Colors.CYAN, Colors.BOLD))
    print(Colors.colorize('  ┃  📍 Default URL: http://localhost:8501                               ┃', Colors.CYAN, Colors.BOLD))
    print(Colors.colorize('  ┃  ⏹  Press Ctrl+C to stop the dashboard                               ┃', Colors.CYAN))
    print(Colors.colorize('  ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛', Colors.CYAN))
    print()
    
    subprocess.run(
        [streamlit_cmd, "run", str(app_path)],
        cwd=str(PROJECT_ROOT)
    )
    
    console.phase_complete()


# =============================================================================
# ALL: Run Everything
# =============================================================================
def run_all(k_override=None):
    """
    Run the complete pipeline in sequence:
    Pipeline → Validate → Analyze → Visualize
    
    Parameters:
    -----------
    k_override : int, optional
        If provided, use this k value instead of auto-detection.
    
    Note: Does not start the app automatically.
    """
    # Print banner and initialize timing
    console.print_banner()
    console.total_start_time = time.time()
    console.phase_results = []  # Reset phase results
    
    print(Colors.colorize(f"  Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", Colors.BRIGHT_BLACK))
    print()
    
    # 1. Core Pipeline
    result = run_pipeline(k_override=k_override)
    if result is None:
        console.error("Pipeline failed. Stopping execution.")
        return
    
    # 2. Validation Scripts
    run_validate()
    
    # 3. Sensitivity Analysis
    run_sensitivity()
    
    # 4. Analysis Modules
    run_analyze()
    
    # 5. Visualizations
    run_viz()
    
    # Final Summary
    console.print_final_summary(output_paths={
        'Processed Data': PROCESSED_DATA_DIR,
        'Figures': FIGURES_DIR,
        'Models': MODELS_DIR,
        'Artifacts': ARTIFACTS_DIR
    })


# =============================================================================
# CLI ARGUMENT PARSER
# =============================================================================
def create_parser():
    """Create the argument parser for the CLI."""
    parser = argparse.ArgumentParser(
        description="AMR Thesis Project - Central Pipeline Orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py --pipeline    # Run data processing pipeline
  python main.py --validate    # Run validation scripts
  python main.py --analyze     # Run analysis modules
  python main.py --viz         # Regenerate visualizations
  python main.py --app         # Launch Streamlit dashboard
  python main.py --all         # Run everything (except app)

For more information, see the project README.
        """
    )
    
    parser.add_argument(
        '--pipeline',
        action='store_true',
        help='Run core data pipeline: Ingestion → Cleaning → Encoding → Clustering'
    )
    
    parser.add_argument(
        '--validate',
        action='store_true',
        help='Run validation scripts (validate_clustering.py, coresistance_analysis.py)'
    )
    
    parser.add_argument(
        '--sensitivity',
        action='store_true',
        help='Run sensitivity analysis (encoding_sensitivity.py, threshold_sensitivity.py)'
    )
    
    parser.add_argument(
        '--analyze',
        action='store_true',
        help='Run analysis modules (regional_environmental, integration_synthesis, supervised)'
    )
    
    parser.add_argument(
        '--viz',
        action='store_true',
        help='Regenerate all visualizations from clustered dataset'
    )
    
    parser.add_argument(
        '--app',
        action='store_true',
        help='Launch the Streamlit interactive dashboard'
    )
    
    parser.add_argument(
        '--all',
        action='store_true',
        help='Run full pipeline: Pipeline → Validate → Analyze → Viz'
    )
    
    parser.add_argument(
        '--k',
        type=int,
        default=None,
        help='Override cluster count k (e.g., --k 4). Only affects --pipeline and --all.'
    )
    
    return parser


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================
def print_styled_help():
    """Print a beautiful styled help menu."""
    console.print_banner(mini=True)
    
    print(Colors.colorize('  Usage:', Colors.BRIGHT_WHITE, Colors.BOLD))
    print(Colors.colorize('  ──────────────────────────────────────────────────────────────────', Colors.BRIGHT_BLACK))
    print()
    
    commands = [
        ('--pipeline', 'Run core data pipeline', 'Ingestion → Cleaning → Encoding → Clustering'),
        ('--validate', 'Run validation scripts', 'Cluster quality, co-resistance analysis'),
        ('--sensitivity', 'Run sensitivity analysis', 'Methodology robustness checks'),
        ('--analyze', 'Run analysis modules', 'Regional, integration, supervised learning'),
        ('--viz', 'Regenerate visualizations', 'Dendrograms, heatmaps, distribution charts'),
        ('--app', 'Launch Streamlit dashboard', 'Interactive data exploration'),
        ('--all', 'Run full pipeline', 'Pipeline → Validate → Sensitivity → Analyze → Viz'),
    ]
    
    for cmd, desc, detail in commands:
        cmd_colored = Colors.colorize(f'    {cmd:<12}', Colors.CYAN, Colors.BOLD)
        desc_colored = Colors.colorize(desc, Colors.WHITE)
        detail_colored = Colors.colorize(f'({detail})', Colors.BRIGHT_BLACK)
        print(f'{cmd_colored} {desc_colored}')
        print(f'                 {detail_colored}')
        print()
    
    print(Colors.colorize('  Options:', Colors.BRIGHT_WHITE, Colors.BOLD))
    print(Colors.colorize('  ──────────────────────────────────────────────────────────────────', Colors.BRIGHT_BLACK))
    print()
    opt_colored = Colors.colorize('    --k N', Colors.CYAN, Colors.BOLD)
    opt_desc = Colors.colorize('Override cluster count', Colors.WHITE)
    opt_detail = Colors.colorize('(e.g., --k 4)', Colors.BRIGHT_BLACK)
    print(f'{opt_colored}       {opt_desc} {opt_detail}')
    print()
    
    print(Colors.colorize('  Examples:', Colors.BRIGHT_WHITE, Colors.BOLD))
    print(Colors.colorize('  ──────────────────────────────────────────────────────────────────', Colors.BRIGHT_BLACK))
    examples = [
        'python main.py --pipeline          # Run data processing',
        'python main.py --all               # Run everything',
        'python main.py --app               # Launch dashboard',
        'python main.py --pipeline --k 5    # Pipeline with k=5',
    ]
    for ex in examples:
        print(Colors.colorize(f'    {ex}', Colors.BRIGHT_BLACK))
    print()
    
    print(Colors.colorize('  ⚠ No command specified. Use one of the flags above.', Colors.YELLOW, Colors.BOLD))
    print()


def main():
    """Main entry point for the AMR pipeline orchestrator."""
    parser = create_parser()
    args = parser.parse_args()
    
    # Check if no arguments provided
    if not any([args.pipeline, args.validate, args.sensitivity, args.analyze, args.viz, args.app, args.all]):
        print_styled_help()
        return
    
    # Print mini banner for individual commands (not --all which prints full banner)
    if not args.all:
        console.print_banner(mini=True)
    
    # Execute requested commands
    if args.all:
        run_all(k_override=args.k)
    else:
        if args.pipeline:
            run_pipeline(k_override=args.k)
        if args.validate:
            run_validate()
        if args.sensitivity:
            run_sensitivity()
        if args.analyze:
            run_analyze()
        if args.viz:
            run_viz()
        if args.app:
            run_app()


if __name__ == "__main__":
    main()

