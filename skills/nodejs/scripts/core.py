"""
BM25 + Regex Hybrid Search Engine for NodeSkill
Adapted from GoSkill search engine for Node.js/TypeScript domains.
"""

import csv
import math
import os
import re
from pathlib import Path

# ─── CSV Configuration ───────────────────────────────────────────────

DATA_DIR = Path(__file__).parent.parent / "data"

CSV_CONFIG = {
    "pattern": {
        "file": "patterns.csv",
        "search_cols": ["Pattern Name", "Category", "Keywords", "Description", "When To Use"],
        "output_cols": ["Pattern Name", "Category", "Keywords", "Description", "When To Use", "When Not To Use", "Code Example", "Anti-Pattern", "Complexity", "Common Packages", "Docs URL"]
    },
    "package": {
        "file": "packages.csv",
        "search_cols": ["Use Case", "Keywords", "Primary Package", "Key Considerations"],
        "output_cols": ["Use Case", "Keywords", "Primary Package", "Secondary Packages", "Architecture Pattern", "Error Strategy", "Key Considerations"]
    },
    "error": {
        "file": "error-handling.csv",
        "search_cols": ["Pattern Name", "Keywords", "Context", "Do"],
        "output_cols": ["Pattern Name", "Keywords", "Context", "Do", "Dont", "Code Good", "Code Bad", "Severity", "Package"]
    },
    "perf": {
        "file": "performance.csv",
        "search_cols": ["Category", "Issue", "Keywords", "Description"],
        "output_cols": ["Category", "Issue", "Description", "Do", "Dont", "Code Good", "Code Bad", "Impact", "Measurement"]
    },
    "test": {
        "file": "testing.csv",
        "search_cols": ["Pattern Name", "Category", "Keywords", "Description"],
        "output_cols": ["Pattern Name", "Category", "Description", "Setup", "Code Example", "Common Mistake", "Package", "Docs URL"]
    },
    "async": {
        "file": "async.csv",
        "search_cols": ["Pattern Name", "Keywords", "Description"],
        "output_cols": ["Pattern Name", "Description", "Do", "Dont", "Code Good", "Code Bad", "Severity", "Docs URL"]
    },
    "typescript": {
        "file": "typescript.csv",
        "search_cols": ["Guideline", "Category", "Keywords", "Description"],
        "output_cols": ["Guideline", "Category", "Description", "Do", "Dont", "Code Good", "Code Bad", "Severity", "Docs URL"]
    },
    "arch": {
        "file": "architecture.csv",
        "search_cols": ["Architecture", "Keywords", "Description", "When To Use"],
        "output_cols": ["Architecture", "Description", "Structure", "When To Use", "When Not To Use", "Example Packages", "Key Principles"]
    },
    "idiom": {
        "file": "idioms.csv",
        "search_cols": ["Idiom", "Category", "Keywords", "Description"],
        "output_cols": ["Idiom", "Category", "Description", "Do", "Dont", "Code Good", "Code Bad", "Docs URL"]
    },
    "anti": {
        "file": "anti-patterns.csv",
        "search_cols": ["Anti-Pattern", "Category", "Keywords", "What People Do Wrong"],
        "output_cols": ["Anti-Pattern", "Category", "What People Do Wrong", "Why Its Bad", "Better Approach", "Code Bad", "Code Good", "Severity", "Detection"]
    },
    "tooling": {
        "file": "tooling.csv",
        "search_cols": ["Tool", "Category", "Keywords", "Description"],
        "output_cols": ["Tool", "Category", "Description", "Command", "When To Use", "Output"]
    },
    "module": {
        "file": "modules.csv",
        "search_cols": ["Topic", "Keywords", "Description"],
        "output_cols": ["Topic", "Description", "Do", "Dont", "Example", "Severity", "Docs URL"]
    }
}

STACK_CONFIG = {
    "web-express": "stacks/web-express.csv",
    "web-fastify": "stacks/web-fastify.csv",
    "web-nestjs": "stacks/web-nestjs.csv",
    "web-hono": "stacks/web-hono.csv",
    "web-nextjs": "stacks/web-nextjs.csv",
    "web-remix": "stacks/web-remix.csv",
    "orm-prisma": "stacks/orm-prisma.csv",
    "orm-drizzle": "stacks/orm-drizzle.csv",
    "graphql": "stacks/graphql.csv",
    "database": "stacks/database.csv",
    "auth": "stacks/auth.csv",
    "realtime": "stacks/realtime.csv",
}

STACK_SEARCH_COLS = ["Category", "Guideline", "Description", "Do"]
STACK_OUTPUT_COLS = ["Category", "Guideline", "Description", "Do", "Dont", "Code Good", "Code Bad", "Severity", "Docs URL"]

# ─── Domain Auto-Detection Keywords ──────────────────────────────────

DOMAIN_KEYWORDS = {
    "pattern": ["pattern", "factory", "builder", "singleton", "decorator", "adapter", "facade", "proxy", "strategy", "observer", "middleware", "pipeline", "di container", "dependency injection", "command pattern"],
    "package": ["package", "library", "npm", "install", "recommend", "which package", "framework", "express", "fastify", "prisma", "zod", "what library", "alternative"],
    "error": ["error", "try", "catch", "throw", "exception", "error handling", "custom error", "validation error", "error boundary", "unhandled rejection", "uncaught exception"],
    "perf": ["performance", "fast", "slow", "optimize", "memory leak", "event loop", "profiling", "benchmark", "bundle", "tree shaking", "caching", "clinic", "heap snapshot"],
    "test": ["test", "vitest", "jest", "mock", "e2e", "playwright", "supertest", "coverage", "snapshot", "fixture", "testing", "describe", "it.each"],
    "async": ["async", "await", "promise", "stream", "event loop", "worker thread", "cluster", "abort", "queue", "event emitter", "observable", "callback", "readable", "writable", "transform"],
    "typescript": ["typescript", "type", "interface", "generic", "utility type", "discriminated union", "type guard", "branded", "zod", "satisfies", "infer", "keyof", "mapped type"],
    "arch": ["architecture", "structure", "organize", "layout", "project", "hexagonal", "clean", "modular", "monorepo", "serverless", "microservice", "directory", "feature-sliced"],
    "idiom": ["idiom", "best practice", "convention", "style", "proper way", "node way", "barrel", "esm", "cjs", "nullish", "optional chaining", "as const"],
    "anti": ["anti-pattern", "bad practice", "mistake", "wrong", "avoid", "dont", "smell", "code smell", "gotcha", "common mistake", "callback hell"],
    "tooling": ["tool", "eslint", "prettier", "tsc", "tsx", "node", "npm", "pnpm", "bun", "nvm", "esbuild", "turbo", "vitest", "nodemon"],
    "module": ["module", "package.json", "lockfile", "semver", "workspace", "monorepo", "esm", "cjs", "tree shaking", "peer dep", "exports field"]
}


# ─── BM25 Search Engine ─────────────────────────────────────────────

class BM25:
    """BM25 ranking algorithm for text search."""

    def __init__(self, k1=1.5, b=0.75):
        self.k1 = k1
        self.b = b
        self.corpus = []
        self.doc_len = []
        self.avgdl = 0
        self.doc_freqs = {}
        self.idf = {}
        self.N = 0

    @staticmethod
    def tokenize(text):
        """Tokenize text into lowercase terms."""
        if not text:
            return []
        text = str(text).lower()
        text = re.sub(r'[^\w\s./\-]', ' ', text)
        tokens = text.split()
        return [t for t in tokens if len(t) > 1]

    def fit(self, corpus):
        """Fit BM25 model on a corpus of documents."""
        self.corpus = corpus
        self.N = len(corpus)
        self.doc_len = []
        nd = {}

        for doc in corpus:
            tokens = self.tokenize(doc)
            self.doc_len.append(len(tokens))
            seen = set()
            for token in tokens:
                if token not in seen:
                    nd[token] = nd.get(token, 0) + 1
                    seen.add(token)

        self.avgdl = sum(self.doc_len) / self.N if self.N > 0 else 0
        self.doc_freqs = nd

        # Calculate IDF
        for word, freq in nd.items():
            self.idf[word] = math.log((self.N - freq + 0.5) / (freq + 0.5) + 1)

    def score(self, query):
        """Score all documents against a query."""
        query_tokens = self.tokenize(query)
        scores = [0.0] * self.N

        for token in query_tokens:
            if token not in self.idf:
                continue
            idf = self.idf[token]
            for idx, doc in enumerate(self.corpus):
                doc_tokens = self.tokenize(doc)
                tf = doc_tokens.count(token)
                dl = self.doc_len[idx]
                numerator = tf * (self.k1 + 1)
                denominator = tf + self.k1 * (1 - self.b + self.b * dl / self.avgdl) if self.avgdl > 0 else tf + self.k1
                scores[idx] += idf * numerator / denominator

        return scores


# ─── CSV Loading ─────────────────────────────────────────────────────

def load_csv(filepath):
    """Load a CSV file and return list of row dicts."""
    rows = []
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                rows.append(row)
    except FileNotFoundError:
        pass
    return rows


def build_search_text(row, search_cols):
    """Concatenate search columns into a single searchable string."""
    parts = []
    for col in search_cols:
        val = row.get(col, "")
        if val:
            parts.append(str(val))
    return " ".join(parts)


# ─── Search Functions ────────────────────────────────────────────────

def search_domain(query, domain, max_results=5):
    """Search a specific domain CSV using BM25."""
    if domain not in CSV_CONFIG:
        return []

    config = CSV_CONFIG[domain]
    filepath = DATA_DIR / config["file"]
    rows = load_csv(filepath)

    if not rows:
        return []

    # Build corpus from search columns
    search_cols = config["search_cols"]
    corpus = [build_search_text(row, search_cols) for row in rows]

    # BM25 search
    bm25 = BM25()
    bm25.fit(corpus)
    scores = bm25.score(query)

    # Rank and filter
    ranked = sorted(enumerate(scores), key=lambda x: x[1], reverse=True)
    results = []
    output_cols = config["output_cols"]

    for idx, score in ranked[:max_results]:
        if score > 0:
            row = rows[idx]
            result = {col: row.get(col, "") for col in output_cols if col in row}
            result["_score"] = round(score, 4)
            results.append(result)

    return results


def search_stack(query, stack, max_results=5):
    """Search a stack-specific CSV using BM25."""
    if stack not in STACK_CONFIG:
        return []

    filepath = DATA_DIR / STACK_CONFIG[stack]
    rows = load_csv(filepath)

    if not rows:
        return []

    corpus = [build_search_text(row, STACK_SEARCH_COLS) for row in rows]

    bm25 = BM25()
    bm25.fit(corpus)
    scores = bm25.score(query)

    ranked = sorted(enumerate(scores), key=lambda x: x[1], reverse=True)
    results = []

    for idx, score in ranked[:max_results]:
        if score > 0:
            row = rows[idx]
            result = {col: row.get(col, "") for col in STACK_OUTPUT_COLS if col in row}
            result["_score"] = round(score, 4)
            results.append(result)

    return results


def detect_domain(query):
    """Auto-detect the most relevant domain based on query keywords."""
    query_lower = query.lower()
    scores = {}

    for domain, keywords in DOMAIN_KEYWORDS.items():
        score = 0
        for kw in keywords:
            if kw in query_lower:
                score += len(kw)  # Longer keyword matches score higher
        if score > 0:
            scores[domain] = score

    if not scores:
        return None

    return max(scores, key=scores.get)


def multi_domain_search(query, domains=None, max_results=3):
    """Search across multiple domains."""
    if domains is None:
        domains = ["package", "pattern", "error", "arch", "async", "perf"]

    all_results = {}
    for domain in domains:
        results = search_domain(query, domain, max_results=max_results)
        if results:
            all_results[domain] = results

    return all_results


# ─── Output Formatting ───────────────────────────────────────────────

def format_results(results, domain=None, format_type="text"):
    """Format search results for display."""
    if not results:
        return "No results found."

    if format_type == "markdown":
        return _format_markdown(results, domain)
    return _format_text(results, domain)


def _format_text(results, domain=None):
    """Format results as plain text."""
    lines = []
    header = f"=== {domain.upper()} SEARCH RESULTS ===" if domain else "=== SEARCH RESULTS ==="
    lines.append(header)
    lines.append("")

    for i, result in enumerate(results, 1):
        score = result.pop("_score", 0)
        lines.append(f"--- Result {i} (score: {score}) ---")
        for key, value in result.items():
            if value:
                # Truncate long values
                display_val = str(value)
                if len(display_val) > 200:
                    display_val = display_val[:200] + "..."
                lines.append(f"  {key}: {display_val}")
        lines.append("")

    return "\n".join(lines)


def _format_markdown(results, domain=None):
    """Format results as markdown."""
    lines = []
    header = f"## {domain.upper()} Search Results" if domain else "## Search Results"
    lines.append(header)
    lines.append("")

    for i, result in enumerate(results, 1):
        score = result.pop("_score", 0)
        # Use first meaningful field as title
        title = list(result.values())[0] if result else f"Result {i}"
        lines.append(f"### {i}. {title}")
        lines.append(f"*Relevance: {score}*")
        lines.append("")

        for key, value in result.items():
            if value and key != list(result.keys())[0]:
                lines.append(f"**{key}:** {value}")
                lines.append("")

        lines.append("---")
        lines.append("")

    return "\n".join(lines)
