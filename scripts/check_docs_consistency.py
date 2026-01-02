#!/usr/bin/env python3
"""
Documentation Consistency Checker

Automated script to verify consistency across documentation files.
Run before commits or as part of CI/CD pipeline.

Usage:
    python check_docs_consistency.py
    python check_docs_consistency.py --fix  # Auto-fix simple issues (future feature)
"""

import re
import sys
from pathlib import Path
from typing import List, Dict, Tuple

# ANSI color codes for output
class Colors:
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    RESET = '\033[0m'
    BOLD = '\033[1m'

def print_header(text: str):
    """Print formatted section header"""
    print(f"\n{Colors.BOLD}{Colors.BLUE}{'='*70}{Colors.RESET}")
    print(f"{Colors.BOLD}{Colors.BLUE}{text}{Colors.RESET}")
    print(f"{Colors.BOLD}{Colors.BLUE}{'='*70}{Colors.RESET}\n")

def print_success(text: str):
    """Print success message"""
    print(f"{Colors.GREEN}✓ {text}{Colors.RESET}")

def print_warning(text: str):
    """Print warning message"""
    print(f"{Colors.YELLOW}⚠ {text}{Colors.RESET}")

def print_error(text: str):
    """Print error message"""
    print(f"{Colors.RED}✗ {text}{Colors.RESET}")

class DocumentationChecker:
    def __init__(self, root_dir: Path):
        self.root_dir = root_dir
        self.errors: List[str] = []
        self.warnings: List[str] = []
        self.checks_passed = 0
        self.checks_failed = 0

    def check_k_values(self) -> None:
        """Check for hardcoded k=5 references (should be k=4 or dynamic)"""
        print_header("Checking for Hardcoded k=5 References")
        
        files_to_check = [
            "README.md",
            "USER_MANUAL.md",
            "docs/TECHNICAL_REFERENCE.md",
            "docs/limitations.md",
        ]
        
        k5_pattern = re.compile(r'\bk\s*=\s*5\b')
        five_clusters_pattern = re.compile(r'\b5\s+clusters?\b', re.IGNORECASE)
        
        # Acceptable contexts for k=5 (comparison tables, metric recommendations)
        acceptable_contexts = [
            'Recommended k',  # In comparison tables showing what each metric recommends
            'Davies-Bouldin',  # DB Index recommends k=5
            'DEPRECATED',  # Historical documents
            'historical',
            'original',
        ]
        
        found_issues = False
        
        for file_path_str in files_to_check:
            file_path = self.root_dir / file_path_str
            if not file_path.exists():
                self.warnings.append(f"File not found: {file_path_str}")
                continue
            
            content = file_path.read_text(encoding='utf-8')
            lines = content.split('\n')
            
            for i, line in enumerate(lines, 1):
                if k5_pattern.search(line):
                    # Check if in acceptable context
                    context_ok = any(ctx.lower() in line.lower() for ctx in acceptable_contexts)
                    
                    # Also check surrounding lines for table context
                    start_idx = max(0, i-3)
                    end_idx = min(len(lines), i+2)
                    surrounding = ' '.join(lines[start_idx:end_idx])
                    context_ok = context_ok or any(ctx.lower() in surrounding.lower() for ctx in acceptable_contexts)
                    
                    if not context_ok:
                        self.errors.append(
                            f"{file_path_str}:{i} - Found 'k=5' (should be k=4 or dynamic)"
                        )
                        print_error(f"{file_path_str}:{i} - Found 'k=5'")
                        found_issues = True
                
                if five_clusters_pattern.search(line):
                    # Check if it's in acceptable context (e.g., historical reference)
                    if 'DEPRECATED' not in line and 'historical' not in line.lower():
                        self.warnings.append(
                            f"{file_path_str}:{i} - Found '5 clusters' (verify context)"
                        )
                        print_warning(f"{file_path_str}:{i} - Found '5 clusters'")
                        found_issues = True
        
        if not found_issues:
            print_success("No hardcoded k=5 references found")
            self.checks_passed += 1
        else:
            self.checks_failed += 1

    def check_silhouette_scores(self) -> None:
        """Check for old silhouette score (0.488 for k=5, should be 0.466 for k=4)"""
        print_header("Checking Silhouette Scores")
        
        files_to_check = [
            "docs/TECHNICAL_REFERENCE.md",
            "docs/limitations.md",
        ]
        
        old_score_pattern = re.compile(r'\b0\.488\b')
        new_score_pattern = re.compile(r'\b0\.466\b')
        
        found_old_score = False
        found_new_score = False
        
        for file_path_str in files_to_check:
            file_path = self.root_dir / file_path_str
            if not file_path.exists():
                continue
            
            content = file_path.read_text(encoding='utf-8')
            lines = content.split('\n')
            
            for i, line in enumerate(lines, 1):
                if old_score_pattern.search(line):
                    # Check if in deprecated context
                    if 'DEPRECATED' not in content[max(0, i-10):i+10]:
                        self.errors.append(
                            f"{file_path_str}:{i} - Found old silhouette score 0.488 (k=5)"
                        )
                        print_error(f"{file_path_str}:{i} - Old score 0.488 found")
                        found_old_score = True
                
                if new_score_pattern.search(line):
                    found_new_score = True
        
        if not found_old_score and found_new_score:
            print_success("Silhouette scores are current (0.466 for k=4)")
            self.checks_passed += 1
        elif not found_old_score and not found_new_score:
            self.warnings.append("No silhouette scores found in checked files")
            print_warning("No silhouette scores found")
        else:
            self.checks_failed += 1

    def check_encoding_consistency(self) -> None:
        """Check that encoding scheme (S=0, I=1, R=2) is consistent"""
        print_header("Checking Encoding Scheme Consistency")
        
        files_to_check = [
            "README.md",
            "USER_MANUAL.md",
            "docs/TECHNICAL_REFERENCE.md",
            "docs/archive/methods/REF_preprocessing_technical_notes.md",
            "docs/archive/methods/REF_clustering_technical_notes.md",
        ]
        
        encoding_pattern = re.compile(r'S\s*=\s*0.*I\s*=\s*1.*R\s*=\s*2', re.IGNORECASE)
        
        files_with_encoding = []
        
        for file_path_str in files_to_check:
            file_path = self.root_dir / file_path_str
            if not file_path.exists():
                continue
            
            content = file_path.read_text(encoding='utf-8')
            
            if encoding_pattern.search(content):
                files_with_encoding.append(file_path_str)
        
        if len(files_with_encoding) >= 3:  # Expect at least in 3 files
            print_success(f"Encoding scheme found consistently in {len(files_with_encoding)} files")
            self.checks_passed += 1
        else:
            self.warnings.append(f"Encoding scheme found in only {len(files_with_encoding)} files")
            print_warning(f"Encoding found in {len(files_with_encoding)} files (expected ≥3)")

    def check_mdr_definition(self) -> None:
        """Check MDR definition consistency (≥3 classes)"""
        print_header("Checking MDR Definition Consistency")
        
        files_to_check = [
            "README.md",
            "docs/TECHNICAL_REFERENCE.md",
            "docs/limitations.md",
            "docs/archive/methods/REF_preprocessing_technical_notes.md",
        ]
        
        mdr_pattern = re.compile(r'≥\s*3.*class', re.IGNORECASE)
        
        files_with_mdr = []
        
        for file_path_str in files_to_check:
            file_path = self.root_dir / file_path_str
            if not file_path.exists():
                continue
            
            content = file_path.read_text(encoding='utf-8')
            
            if mdr_pattern.search(content):
                files_with_mdr.append(file_path_str)
        
        if len(files_with_mdr) >= 2:
            print_success(f"MDR definition (≥3 classes) found in {len(files_with_mdr)} files")
            self.checks_passed += 1
        else:
            self.warnings.append("MDR definition not found consistently")
            print_warning("MDR definition found in too few files")

    def check_cross_references(self) -> None:
        """Check that cross-reference links exist"""
        print_header("Checking Cross-References")
        
        main_docs = [
            "README.md",
            "USER_MANUAL.md",
            "docs/TECHNICAL_REFERENCE.md",
        ]
        
        # Pattern for markdown links
        link_pattern = re.compile(r'\[([^\]]+)\]\(([^\)]+)\)')
        
        broken_links = []
        
        for doc_path_str in main_docs:
            doc_path = self.root_dir / doc_path_str
            if not doc_path.exists():
                continue
            
            content = doc_path.read_text(encoding='utf-8')
            doc_dir = doc_path.parent
            
            for match in link_pattern.finditer(content):
                link_text, link_target = match.groups()
                
                # Skip external links and anchors
                if link_target.startswith(('http://', 'https://', '#', 'mailto:')):
                    continue
                
                # Remove anchor from file path
                file_path = link_target.split('#')[0]
                if not file_path:
                    continue
                
                # Resolve path relative to document
                target_path = (doc_dir / file_path).resolve()
                
                if not target_path.exists():
                    broken_links.append(f"{doc_path_str}: Link to '{link_target}' (target not found)")
                    print_error(f"{doc_path_str}: Broken link to '{link_target}'")
        
        if not broken_links:
            print_success("All cross-reference links are valid")
            self.checks_passed += 1
        else:
            self.errors.extend(broken_links)
            self.checks_failed += 1

    def check_file_existence(self) -> None:
        """Check that all mentioned output files exist in documentation"""
        print_header("Checking Referenced Output Files")
        
        # Expected important files
        critical_files = [
            "data/processed/analysis_ready_dataset.csv",
            "data/processed/clustered_dataset.csv",
            "data/processed/figures/cluster_validation.png",
            "data/processed/figures/dendrogram_highres.png",
        ]
        
        missing_files = []
        
        for file_path_str in critical_files:
            file_path = self.root_dir / file_path_str
            if not file_path.exists():
                missing_files.append(file_path_str)
                print_warning(f"Referenced file not found: {file_path_str}")
        
        if not missing_files:
            print_success("All critical output files exist")
            self.checks_passed += 1
        else:
            self.warnings.extend([f"Missing file: {f}" for f in missing_files])
            print_warning(f"{len(missing_files)} referenced files not found (may need pipeline run)")

    def run_all_checks(self) -> bool:
        """Run all consistency checks"""
        print(f"\n{Colors.BOLD}Documentation Consistency Checker{Colors.RESET}")
        print(f"Root directory: {self.root_dir}\n")
        
        self.check_k_values()
        self.check_silhouette_scores()
        self.check_encoding_consistency()
        self.check_mdr_definition()
        self.check_cross_references()
        self.check_file_existence()
        
        # Print summary
        print_header("Summary")
        
        print(f"{Colors.BOLD}Checks Passed:{Colors.RESET} {Colors.GREEN}{self.checks_passed}{Colors.RESET}")
        print(f"{Colors.BOLD}Checks Failed:{Colors.RESET} {Colors.RED}{self.checks_failed}{Colors.RESET}")
        print(f"{Colors.BOLD}Warnings:{Colors.RESET} {Colors.YELLOW}{len(self.warnings)}{Colors.RESET}")
        print(f"{Colors.BOLD}Errors:{Colors.RESET} {Colors.RED}{len(self.errors)}{Colors.RESET}")
        
        if self.errors:
            print(f"\n{Colors.RED}{Colors.BOLD}❌ ERRORS FOUND:{Colors.RESET}")
            for error in self.errors:
                print(f"  {Colors.RED}• {error}{Colors.RESET}")
        
        if self.warnings:
            print(f"\n{Colors.YELLOW}{Colors.BOLD}⚠️  WARNINGS:{Colors.RESET}")
            for warning in self.warnings:
                print(f"  {Colors.YELLOW}• {warning}{Colors.RESET}")
        
        if not self.errors and not self.warnings:
            print(f"\n{Colors.GREEN}{Colors.BOLD}✅ All checks passed! Documentation is consistent.{Colors.RESET}\n")
            return True
        elif not self.errors:
            print(f"\n{Colors.YELLOW}{Colors.BOLD}⚠️  Checks passed with warnings. Review above.{Colors.RESET}\n")
            return True
        else:
            print(f"\n{Colors.RED}{Colors.BOLD}❌ Documentation consistency check FAILED. Fix errors above.{Colors.RESET}\n")
            return False

def main():
    """Main entry point"""
    # Get project root (assumes script is in project root)
    root_dir = Path(__file__).parent.absolute()
    
    # If script is in a subdirectory, adjust
    if (root_dir / "amr_thesis_project_code").exists():
        root_dir = root_dir / "amr_thesis_project_code"
    
    checker = DocumentationChecker(root_dir)
    success = checker.run_all_checks()
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
