"""
Architecture System Generator for NodeSkill
Generates complete NodeJS/TypeScript project architecture recommendations by combining
multi-domain search results with reasoning rules.
"""

import csv
import json
import os
from pathlib import Path
from datetime import datetime

from core import (
    DATA_DIR, load_csv, search_domain, multi_domain_search,
    BM25, build_search_text, detect_domain
)

REASONING_FILE = DATA_DIR / "nodejs-reasoning.csv"


def load_reasoning(query):
    """Load and match reasoning rules from nodejs-reasoning.csv."""
    rows = load_csv(REASONING_FILE)
    if not rows:
        return None

    search_cols = ["Project_Category", "Recommended_Architecture", "Package_Priority"]
    corpus = [build_search_text(row, search_cols) for row in rows]

    bm25 = BM25()
    bm25.fit(corpus)
    scores = bm25.score(query)

    ranked = sorted(enumerate(scores), key=lambda x: x[1], reverse=True)

    for idx, score in ranked[:1]:
        if score > 0:
            return rows[idx]

    return rows[0] if rows else None


def apply_reasoning(reasoning, search_results):
    """Apply reasoning rules to search results to select best combination."""
    if not reasoning:
        return search_results

    selected = {}

    # Extract decision rules
    rules_str = reasoning.get("Decision_Rules", "{}")
    try:
        rules = json.loads(rules_str.replace("'", '"'))
    except (json.JSONDecodeError, ValueError):
        rules = {}

    selected["architecture"] = reasoning.get("Recommended_Architecture", "Layered Architecture")
    selected["packages"] = reasoning.get("Package_Priority", "")
    selected["error_strategy"] = reasoning.get("Error_Strategy", "")
    selected["concurrency_strategy"] = reasoning.get("Concurrency_Strategy", "")
    selected["testing_strategy"] = reasoning.get("Testing_Strategy", "")
    selected["anti_patterns"] = reasoning.get("Anti_Patterns", "")
    selected["decision_rules"] = rules
    selected["severity"] = reasoning.get("Severity", "High")

    # Merge search results
    for domain, results in search_results.items():
        selected[f"{domain}_results"] = results

    return selected


def format_arch_system(selected, project_name, format_type="box"):
    """Format architecture recommendation as ASCII box or markdown."""
    if format_type == "markdown":
        return _format_markdown(selected, project_name)
    return _format_box(selected, project_name)


def _format_box(selected, project_name):
    """Format as ASCII box (default)."""
    W = 90
    lines = []

    def hr():
        lines.append("+" + "-" * (W - 2) + "+")

    def row(text=""):
        text = str(text)
        if len(text) > W - 6:
            text = text[:W - 9] + "..."
        lines.append("|  " + text.ljust(W - 4) + "|")

    def section(title):
        row()
        row(f"{title}:")

    hr()
    row(f"TARGET: {project_name} - RECOMMENDED NODEJS ARCHITECTURE")
    hr()

    # Architecture
    section("ARCHITECTURE")
    arch = selected.get("architecture", "Clean / Layered Architecture")
    row(f"   {arch}")

    # Packages
    section("CORE PACKAGES")
    packages = selected.get("packages", "")
    if packages:
        for pkg in packages.split(","):
            pkg = pkg.strip()
            if pkg:
                row(f"   {pkg}")

    # Also include package search results
    pkg_results = selected.get("package_results", [])
    if pkg_results:
        for r in pkg_results[:5]:
            use_case = r.get("Use Case", "")
            primary = r.get("Primary Package", "")
            if use_case and primary:
                row(f"   {use_case}: {primary}")

    # Error Strategy
    section("ERROR STRATEGY")
    error_strategy = selected.get("error_strategy", "")
    if error_strategy:
        for line in error_strategy.split(";"):
            line = line.strip()
            if line:
                row(f"   {line}")

    error_results = selected.get("error_results", [])
    if error_results:
        for r in error_results[:3]:
            pattern = r.get("Pattern Name", "")
            do = r.get("Do", "")
            if pattern:
                row(f"   {pattern}: {do[:60]}")

    # Patterns
    section("PATTERNS")
    pattern_results = selected.get("pattern_results", [])
    if pattern_results:
        for r in pattern_results[:4]:
            name = r.get("Pattern Name", "")
            desc = r.get("Description", "")
            if name:
                row(f"   {name}: {desc[:60]}")

    # Project Structure
    section("PROJECT STRUCTURE")
    arch_results = selected.get("arch_results", [])
    if arch_results:
        structure = arch_results[0].get("Structure", "")
        if structure:
            for line in structure.split("/"):
                line = line.strip()
                if line:
                    row(f"   {line}")
    else:
        # Default structure
        for d in [
            "src/", "   index.ts", "   app.ts", "   config/",
            "   routes/", "   controllers/", "   services/",
            "   models/", "   middlewares/", "   utils/",
            "package.json", "tsconfig.json", "Dockerfile"
        ]:
            row(f"   {d}")

    # Concurrency
    section("CONCURRENCY / ASYNC")
    conc_strategy = selected.get("async_results", [])
    if conc_strategy:
        for r in conc_strategy[:3]:
            name = r.get("Pattern Name", "")
            if name:
                row(f"   Pattern: {name}")

    # Testing
    section("TESTING")
    test_strategy = selected.get("testing_strategy", "")
    if test_strategy:
        for line in test_strategy.split(";"):
            line = line.strip()
            if line:
                row(f"   {line}")

    # Makefile / package.json scripts
    section("NPM SCRIPTS ESSENTIALS")
    for cmd in [
        "build:    tsc or tsup",
        "start:    node dist/index.js",
        "dev:      tsx watch src/index.ts",
        "test:     vitest or jest",
        "lint:     eslint . --ext .ts",
    ]:
        row(f"   {cmd}")

    # Dockerfile
    section("DOCKERFILE (multi-stage)")
    row("   Stage 1: node:20-alpine -> build (npm install, npm run build)")
    row("   Stage 2: node:20-alpine -> production (copy dist, node dist/index.js)")

    # Anti-patterns
    section("AVOID (Anti-patterns)")
    anti = selected.get("anti_patterns", "")
    if anti:
        row(f"   {anti[:W - 8]}")

    row()
    hr()

    return "\n".join(lines)


def _format_markdown(selected, project_name):
    """Format as markdown."""
    lines = []
    lines.append(f"# {project_name} - NodeJS Architecture Recommendation")
    lines.append("")
    lines.append(f"*Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}*")
    lines.append("")

    # Architecture
    arch = selected.get("architecture", "Clean / Layered Architecture")
    lines.append(f"## Architecture: {arch}")
    lines.append("")

    # Packages
    lines.append("## Core Packages")
    lines.append("")
    packages = selected.get("packages", "")
    if packages:
        for pkg in packages.split(","):
            pkg = pkg.strip()
            if pkg:
                lines.append(f"- {pkg}")

    pkg_results = selected.get("package_results", [])
    if pkg_results:
        lines.append("")
        lines.append("| Use Case | Package |")
        lines.append("|----------|---------|")
        for r in pkg_results[:5]:
            use_case = r.get("Use Case", "")
            primary = r.get("Primary Package", "")
            if use_case and primary:
                lines.append(f"| {use_case} | `{primary}` |")
    lines.append("")

    # Error Strategy
    lines.append("## Error Strategy")
    lines.append("")
    error_strategy = selected.get("error_strategy", "")
    if error_strategy:
        for line in error_strategy.split(";"):
            line = line.strip()
            if line:
                lines.append(f"- {line}")
    lines.append("")

    # Patterns
    lines.append("## Recommended Patterns")
    lines.append("")
    pattern_results = selected.get("pattern_results", [])
    if pattern_results:
        for r in pattern_results[:4]:
            name = r.get("Pattern Name", "")
            desc = r.get("Description", "")
            if name:
                lines.append(f"### {name}")
                lines.append(f"{desc}")
                code = r.get("Code Example", "")
                if code:
                    lines.append("")
                    lines.append("```typescript")
                    lines.append(code)
                    lines.append("```")
                lines.append("")

    # Project Structure
    lines.append("## Project Structure")
    lines.append("")
    lines.append("```")
    arch_results = selected.get("arch_results", [])
    if arch_results:
        structure = arch_results[0].get("Structure", "")
        if structure:
            lines.append(structure)
    else:
        lines.append("src/")
        lines.append("   index.ts")
        lines.append("   app.ts")
        lines.append("   config/")
        lines.append("   controllers/")
        lines.append("   services/")
        lines.append("   middlewares/")
        lines.append("package.json")
        lines.append("tsconfig.json")
        lines.append("Dockerfile")
    lines.append("```")
    lines.append("")

    # Concurrency
    lines.append("## Async Strategy")
    lines.append("")
    conc_strategy = selected.get("concurrency_strategy", "")
    if conc_strategy:
        for line in conc_strategy.split(";"):
            line = line.strip()
            if line:
                lines.append(f"- {line}")
    lines.append("")

    # Testing
    lines.append("## Testing Strategy")
    lines.append("")
    test_strategy = selected.get("testing_strategy", "")
    if test_strategy:
        for line in test_strategy.split(";"):
            line = line.strip()
            if line:
                lines.append(f"- {line}")
    lines.append("")

    # Anti-patterns
    lines.append("## Anti-patterns to Avoid")
    lines.append("")
    anti = selected.get("anti_patterns", "")
    if anti:
        for item in anti.split(","):
            item = item.strip()
            if item:
                lines.append(f"- {item}")
    lines.append("")

    return "\n".join(lines)


def persist_arch_system(output, project_name, page=1, output_dir=None):
    """Save architecture recommendation to file."""
    if output_dir is None:
        output_dir = Path.cwd()
    else:
        output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    safe_name = project_name.lower().replace(" ", "-").replace("/", "-")
    if page > 1:
        filename = f"{safe_name}-arch-p{page}.md"
    else:
        filename = f"{safe_name}-arch.md"

    filepath = output_dir / filename
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(output)

    return str(filepath)


def generate_arch_system(query, project_name="MyProject", format_type="box", persist=False, page=1, output_dir=None):
    """Main entry point for Architecture System Generator.

    1. Search multiple domains in parallel
    2. Load reasoning rules from nodejs-reasoning.csv
    3. Apply reasoning to select best combination
    4. Format output (ASCII box or Markdown)
    5. Persist if requested
    """
    # 1. Search multiple domains
    search_results = {}
    for domain, max_r in [("package", 3), ("pattern", 3), ("error", 2), ("arch", 2), ("async", 2), ("perf", 2)]:
        results = search_domain(query, domain, max_results=max_r)
        if results:
            search_results[domain] = results

    # 2. Load reasoning rules
    reasoning = load_reasoning(query)

    # 3. Apply reasoning
    selected = apply_reasoning(reasoning, search_results)

    # 4. Format output
    if persist:
        format_type = "markdown"
    output = format_arch_system(selected, project_name, format_type)

    # 5. Persist if requested
    saved_path = None
    if persist:
        saved_path = persist_arch_system(output, project_name, page, output_dir)

    return output, saved_path
