"""
Console Utilities Module for AMR Thesis Project
================================================

Provides consistent, colorful, and informative terminal output across all scripts.
This module is the SINGLE SOURCE OF TRUTH for console styling.

Usage:
    from utils.console import console, Colors
    
    console.header("My Section Title")
    console.step(1, 5, "Processing data")
    console.success("Done!")
"""

import sys
import time
from datetime import datetime


# =============================================================================
# TERMINAL COLOR & STYLING
# =============================================================================
class Colors:
    """ANSI color codes for terminal output with Windows fallback."""
    
    ENABLED = True
    
    @classmethod
    def init(cls):
        """Initialize color support, especially for Windows."""
        if sys.platform == 'win32':
            try:
                # Enable UTF-8 for Windows console
                import os
                os.system('')  # Enable ANSI on Windows 10+
                
                import ctypes
                kernel32 = ctypes.windll.kernel32
                # Enable virtual terminal processing
                kernel32.SetConsoleMode(kernel32.GetStdHandle(-11), 7)
                
                # Also try to set UTF-8 code page
                kernel32.SetConsoleOutputCP(65001)
            except Exception:
                cls.ENABLED = False
    
    @classmethod
    def c(cls, text: str, *codes) -> str:
        """Apply color codes to text if colors are enabled."""
        if not cls.ENABLED:
            return text
        code_str = ';'.join(str(c) for c in codes)
        return f"\033[{code_str}m{text}\033[0m"
    
    # Colors
    BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = 30, 31, 32, 33, 34, 35, 36, 37
    
    # Bright colors
    BRIGHT_BLACK, BRIGHT_RED, BRIGHT_GREEN, BRIGHT_YELLOW = 90, 91, 92, 93
    BRIGHT_BLUE, BRIGHT_MAGENTA, BRIGHT_CYAN, BRIGHT_WHITE = 94, 95, 96, 97
    
    # Styles
    BOLD, DIM, ITALIC, UNDERLINE = 1, 2, 3, 4


# Initialize colors on import
Colors.init()


# =============================================================================
# CONSOLE OUTPUT CLASS
# =============================================================================
class Console:
    """Enhanced console output with colors, formatting, and timing."""
    
    # Unicode symbols
    ICONS = {
        'success': '✓', 'error': '✗', 'warning': '⚠', 'info': '→',
        'bullet': '•', 'arrow': '▶', 'check': '✔', 'star': '★',
    }
    
    def __init__(self):
        self.start_time = None
        self.section_start = None
    
    def _format_time(self, seconds: float) -> str:
        """Format time duration."""
        if seconds < 60:
            return f"{seconds:.1f}s"
        mins = int(seconds // 60)
        secs = seconds % 60
        return f"{mins}m {secs:.1f}s"
    
    # =========================================================================
    # SECTION HEADERS
    # =========================================================================
    def header(self, title: str, subtitle: str = None):
        """Print a styled section header."""
        print()
        print(Colors.c('━' * 70, Colors.BRIGHT_BLACK))
        print(Colors.c(f"  {self.ICONS['arrow']}  {title}", Colors.BRIGHT_YELLOW, Colors.BOLD))
        if subtitle:
            print(Colors.c(f"     {subtitle}", Colors.WHITE, Colors.DIM))
        print(Colors.c('━' * 70, Colors.BRIGHT_BLACK))
        print()
        self.section_start = time.time()
    
    def subheader(self, title: str):
        """Print a smaller subheader."""
        print()
        print(Colors.c(f"  ── {title} ──", Colors.CYAN, Colors.BOLD))
        print()
    
    def separator(self, char: str = '─', width: int = 60):
        """Print a separator line."""
        print(Colors.c(f"  {char * width}", Colors.BRIGHT_BLACK))
    
    # =========================================================================
    # STEP TRACKING
    # =========================================================================
    def step(self, current: int, total: int, message: str, detail: str = None):
        """Log a step with numbering."""
        step_num = Colors.c(f"[{current}/{total}]", Colors.CYAN, Colors.BOLD)
        step_msg = Colors.c(message, Colors.WHITE, Colors.BOLD)
        print(f"  {step_num} {step_msg}")
        if detail:
            print(Colors.c(f"       {self.ICONS['info']} {detail}", Colors.DIM))
    
    def progress(self, current: int, total: int, width: int = 30):
        """Print a progress bar."""
        if total == 0:
            return
        pct = current / total
        filled = int(width * pct)
        bar = Colors.c('█' * filled, Colors.GREEN) + Colors.c('░' * (width - filled), Colors.BRIGHT_BLACK)
        pct_str = Colors.c(f"{pct * 100:5.1f}%", Colors.CYAN)
        print(f"       [{bar}] {pct_str}")
    
    # =========================================================================
    # STATUS MESSAGES
    # =========================================================================
    def info(self, message: str):
        """Log info message."""
        icon = Colors.c(self.ICONS['info'], Colors.BLUE)
        print(f"       {icon} {message}")
    
    def success(self, message: str):
        """Log success message."""
        icon = Colors.c(self.ICONS['success'], Colors.GREEN, Colors.BOLD)
        msg = Colors.c(message, Colors.GREEN)
        print(f"       {icon} {msg}")
    
    def warning(self, message: str):
        """Log warning message."""
        icon = Colors.c(self.ICONS['warning'], Colors.YELLOW, Colors.BOLD)
        msg = Colors.c(message, Colors.YELLOW)
        print(f"       {icon} {msg}")
    
    def error(self, message: str):
        """Log error message."""
        icon = Colors.c(self.ICONS['error'], Colors.RED, Colors.BOLD)
        msg = Colors.c(message, Colors.RED)
        print(f"       {icon} {msg}")
    
    def bullet(self, message: str, indent: int = 2):
        """Print a bullet point."""
        spaces = "  " * indent
        icon = Colors.c(self.ICONS['bullet'], Colors.CYAN)
        print(f"{spaces}{icon} {message}")
    
    # =========================================================================
    # DATA DISPLAY
    # =========================================================================
    def table_header(self, columns: list, widths: list = None):
        """Print a table header row."""
        if widths is None:
            widths = [max(15, len(c) + 2) for c in columns]
        
        # Header
        header = "  "
        for col, w in zip(columns, widths):
            header += f"{col:<{w}}"
        print(Colors.c(header, Colors.BRIGHT_CYAN, Colors.BOLD))
        
        # Separator
        sep = "  " + "─" * sum(widths)
        print(Colors.c(sep, Colors.BRIGHT_BLACK))
        
        return widths
    
    def table_row(self, values: list, widths: list, highlight: bool = False):
        """Print a table data row."""
        row = "  "
        for val, w in zip(values, widths):
            row += f"{str(val):<{w}}"
        if highlight:
            print(Colors.c(row, Colors.GREEN, Colors.BOLD))
        else:
            print(row)
    
    def kv(self, key: str, value, indent: int = 2):
        """Print a key-value pair."""
        spaces = "  " * indent
        key_str = Colors.c(f"{key}:", Colors.CYAN)
        print(f"{spaces}{key_str} {value}")
    
    # =========================================================================
    # SECTION COMPLETION
    # =========================================================================
    def complete(self, message: str = None, stats: dict = None):
        """Mark a section as complete with optional stats."""
        elapsed = time.time() - self.section_start if self.section_start else 0
        
        print()
        print(Colors.c('  ┌' + '─' * 66 + '┐', Colors.GREEN))
        
        msg = message or "Complete"
        line = f"  │  {self.ICONS['check']} {msg}"
        print(Colors.c(line.ljust(69) + '│', Colors.GREEN, Colors.BOLD))
        
        time_line = f"  │    Time: {self._format_time(elapsed)}"
        print(Colors.c(time_line.ljust(69) + '│', Colors.GREEN))
        
        if stats:
            for key, value in stats.items():
                stat_line = f"  │    {key}: {value}"
                print(Colors.c(stat_line.ljust(69) + '│', Colors.GREEN))
        
        print(Colors.c('  └' + '─' * 66 + '┘', Colors.GREEN))
        print()
    
    def done(self, message: str):
        """Simple completion message without box."""
        icon = Colors.c(self.ICONS['check'], Colors.GREEN, Colors.BOLD)
        msg = Colors.c(message, Colors.GREEN)
        print(f"       {icon} {msg}")


# Global console instance
console = Console()
