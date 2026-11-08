# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///

"""Builds an interactive HTML gallery for 100 elementary reaction mechanisms."""

from __future__ import annotations

import argparse
import json
import os
import sys

_SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPTS_DIR not in sys.path:
  sys.path.insert(0, _SCRIPTS_DIR)

import reactions_data


def build_gallery(
    base_dir: str = ".", img_dir: str = "", out_html: str = ""
) -> str:
  base_dir = os.path.abspath(base_dir)
  if not img_dir:
    img_dir = os.path.join(base_dir, "images")
  if not out_html:
    out_html = os.path.join(base_dir, "gallery.html")
  rxns = reactions_data.REACTIONS

  # Prepare data payload for HTML
  cards_data = []
  categories = sorted(list(set(r["category"] for r in rxns)))

  # Image srcs must point at the directory the caller actually passed, made
  # relative to the output file so the gallery stays portable.
  html_dir = os.path.dirname(os.path.abspath(out_html))
  try:
    img_rel_base = os.path.relpath(os.path.abspath(img_dir), html_dir)
  except ValueError:
    # Different drives/roots: fall back to an absolute src.
    img_rel_base = os.path.abspath(img_dir)

  for r in rxns:
    rid = r["id"]
    rel_img_path = os.path.join(img_rel_base, f"{rid}.png").replace(os.sep, "/")

    top_c = (
        r["cond"]
        .get("top", "")
        .replace(r"$\text{", "")
        .replace(r"}$", "")
        .replace(r"^\circ", "°")
        .replace(r"\text{", "")
        .replace("}", "")
        .replace("$", "")
    )
    bot_c = (
        r["cond"]
        .get("bottom", "")
        .replace(r"$\text{", "")
        .replace(r"}$", "")
        .replace(r"\text{", "")
        .replace("}", "")
        .replace("$", "")
    )
    cond_str = f"{top_c} | {bot_c}" if bot_c else top_c

    cards_data.append({
        "id": rid,
        "title": r["title"],
        "category": r["category"],
        "desc": r["desc"],
        "cond": cond_str,
        "img": rel_img_path,
        "r1_name": r["r1"]["name"],
        "r1_smi": r["r1"]["smi"],
        "r2_name": r["r2"]["name"] if r["r2"] else None,
        "r2_smi": r["r2"]["smi"] if r["r2"] else None,
        "p1_name": r["p1"]["name"],
        "p1_smi": r["p1"]["smi"],
        "p2_name": r["p2"]["name"] if r["p2"] else None,
        "p2_smi": r["p2"]["smi"] if r["p2"] else None,
    })

  html_content = f"""<!DOCTYPE html>
    <html lang="en">
    <head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>The 100 Elementary Reaction Mechanisms Library</title>
    <style>
      :root {{
        --bg-primary: #0f172a;
        --bg-secondary: #1e293b;
        --bg-card: #1e293b;
        --border-color: #334155;
        --text-primary: #f8fafc;
        --text-secondary: #94a3b8;
        --accent-color: #38bdf8;
        --accent-hover: #0ea5e9;
        --tag-bg: #0369a1;
        --card-hover-border: #38bdf8;
      }}
    
      * {{ box-sizing: border-box; margin: 0; padding: 0; }}
      body {{
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
        background-color: var(--bg-primary);
        color: var(--text-primary);
        line-height: 1.5;
        padding: 24px;
      }}
    
      header {{
        max-width: 1400px;
        margin: 0 auto 32px auto;
        border-bottom: 1px solid var(--border-color);
        padding-bottom: 24px;
      }}
    
      h1 {{
        font-size: 28px;
        font-weight: 700;
        color: #ffffff;
        margin-bottom: 8px;
        display: flex;
        align-items: center;
        gap: 12px;
      }}
    
      .subtitle {{
        font-size: 15px;
        color: var(--text-secondary);
        margin-bottom: 20px;
      }}
    
      .stats-bar {{
        display: flex;
        gap: 16px;
        flex-wrap: wrap;
        margin-bottom: 20px;
      }}
    
      .stat-pill {{
        background: var(--bg-secondary);
        border: 1px solid var(--border-color);
        padding: 6px 14px;
        border-radius: 9999px;
        font-size: 13px;
        font-weight: 500;
        color: var(--accent-color);
      }}
    
      .controls {{
        display: flex;
        flex-direction: column;
        gap: 16px;
        max-width: 1400px;
        margin: 0 auto 28px auto;
      }}
    
      .search-box {{
        width: 100%;
        padding: 12px 18px;
        background: var(--bg-secondary);
        border: 1px solid var(--border-color);
        border-radius: 8px;
        color: #fff;
        font-size: 15px;
        outline: none;
        transition: border-color 0.2s;
      }}
      .search-box:focus {{ border-color: var(--accent-color); }}
    
      .category-tabs {{
        display: flex;
        gap: 8px;
        flex-wrap: wrap;
      }}
    
      .tab-btn {{
        background: var(--bg-secondary);
        border: 1px solid var(--border-color);
        color: var(--text-secondary);
        padding: 6px 14px;
        border-radius: 6px;
        font-size: 13px;
        cursor: pointer;
        transition: all 0.2s;
      }}
    
      .tab-btn:hover, .tab-btn.active {{
        background: var(--tag-bg);
        color: #ffffff;
        border-color: var(--accent-color);
      }}
    
      .gallery-grid {{
        max-width: 1400px;
        margin: 0 auto;
        display: grid;
        grid-template-columns: repeat(auto-fill, minmax(420px, 1fr));
        gap: 24px;
      }}
    
      .rxn-card {{
        background: var(--bg-card);
        border: 1px solid var(--border-color);
        border-radius: 12px;
        overflow: hidden;
        display: flex;
        flex-direction: column;
        transition: transform 0.2s, border-color 0.2s, box-shadow 0.2s;
      }}
    
      .rxn-card:hover {{
        border-color: var(--card-hover-border);
        transform: translateY(-2px);
        box-shadow: 0 10px 25px -5px rgba(0, 0, 0, 0.4);
      }}
    
      .card-header {{
        padding: 14px 18px;
        border-bottom: 1px solid var(--border-color);
        display: flex;
        justify-content: space-between;
        align-items: center;
      }}
    
      .card-id {{
        font-size: 11px;
        font-weight: 700;
        text-transform: uppercase;
        letter-spacing: 0.5px;
        background: rgba(56, 189, 248, 0.15);
        color: var(--accent-color);
        padding: 3px 8px;
        border-radius: 4px;
      }}
    
      .card-cat {{
        font-size: 12px;
        color: var(--text-secondary);
      }}
    
      .card-body {{
        padding: 16px 18px;
        display: flex;
        flex-direction: column;
        flex-grow: 1;
      }}
    
      .card-title {{
        font-size: 16px;
        font-weight: 600;
        color: #ffffff;
        margin-bottom: 8px;
        line-height: 1.3;
      }}
    
      .card-cond {{
        display: inline-block;
        font-size: 12px;
        color: #38bdf8;
        background: #0f172a;
        padding: 3px 8px;
        border-radius: 4px;
        margin-bottom: 12px;
        border: 1px solid #1e293b;
        align-self: flex-start;
      }}
    
      .card-img-container {{
        background: #ffffff;
        border-radius: 8px;
        padding: 12px;
        margin-bottom: 14px;
        cursor: zoom-in;
        display: flex;
        justify-content: center;
        align-items: center;
        min-height: 140px;
        border: 1px solid var(--border-color);
      }}
    
      .card-img {{
        max-width: 100%;
        max-height: 160px;
        object-fit: contain;
      }}
    
      .card-desc {{
        font-size: 13px;
        color: var(--text-secondary);
        line-height: 1.4;
        margin-top: auto;
      }}
    
      /* Modal Zoom View */
      .modal-overlay {{
        display: none;
        position: fixed;
        top: 0; left: 0; right: 0; bottom: 0;
        background: rgba(0, 0, 0, 0.85);
        z-index: 1000;
        justify-content: center;
        align-items: center;
        padding: 24px;
      }}
    
      .modal-overlay.active {{ display: flex; }}
    
      .modal-content {{
        background: var(--bg-card);
        border: 1px solid var(--border-color);
        border-radius: 14px;
        max-width: 1000px;
        width: 100%;
        max-height: 90vh;
        overflow-y: auto;
        padding: 24px;
        position: relative;
        box-shadow: 0 25px 50px -12px rgba(0, 0, 0, 0.7);
      }}
    
      .modal-close {{
        position: absolute;
        top: 18px; right: 20px;
        background: transparent;
        border: none;
        color: var(--text-secondary);
        font-size: 26px;
        cursor: pointer;
      }}
      .modal-close:hover {{ color: #ffffff; }}
    
      .modal-img-container {{
        background: #ffffff;
        border-radius: 8px;
        padding: 20px;
        margin: 18px 0;
        display: flex;
        justify-content: center;
        align-items: center;
      }}
    
      .modal-img {{
        max-width: 100%;
        max-height: 480px;
        object-fit: contain;
      }}
    
      .modal-meta {{
        font-size: 14px;
        color: var(--text-secondary);
        line-height: 1.6;
      }}
    
      .smiles-tag {{
        font-family: monospace;
        font-size: 12px;
        background: #0f172a;
        padding: 2px 6px;
        border-radius: 4px;
        color: #38bdf8;
      }}
    
    </style>
    </head>
    <body>
    
    <header>
      <h1>The 100 Elementary Reaction Mechanisms Library</h1>
      <div class="subtitle">Curated organic synthesis transformations rendered with RDKit 2D depictions, authentic ChemDraw-style curved electron-pushing arrows, and verified collision-free clearances.</div>
      <div class="stats-bar">
        <div class="stat-pill">100 Total Mechanisms</div>
        <div class="stat-pill">10 Reaction Classes</div>
        <div class="stat-pill">100% RDKit Depictions</div>
        <div class="stat-pill">Slender Stealth Arrowheads</div>
        <div class="stat-pill">Zero Bond Overlap</div>
      </div>
    </header>
    
    <div class="controls">
      <input type="text" id="search" class="search-box" placeholder="Search by reaction name, category, reagents, or mechanism keyword...">
      <div class="category-tabs" id="categoryTabs">
        <button class="tab-btn active" data-cat="all">All (100)</button>
        {"".join(f'<button class="tab-btn" data-cat="{c}">{c} (10)</button>' for c in categories)}
      </div>
    </div>
    
    <div class="gallery-grid" id="galleryGrid"></div>
    
    <!-- Modal Zoom -->
    <div class="modal-overlay" id="modalOverlay">
      <div class="modal-content">
        <button class="modal-close" id="modalClose">&times;</button>
        <div class="card-id" id="modalId"></div>
        <h2 id="modalTitle" style="margin-top: 10px; color:#fff; font-size:22px;"></h2>
        <div id="modalCat" style="color:var(--text-secondary); font-size:13px; margin-bottom:12px;"></div>
        <div id="modalCond" class="card-cond"></div>
        
        <div class="modal-img-container">
          <img id="modalImg" class="modal-img" src="" alt="Mechanism Diagram">
        </div>
    
        <div class="modal-meta" id="modalDesc"></div>
        <div style="margin-top: 16px; border-top: 1px solid var(--border-color); padding-top: 14px; font-size: 13px;">
          <div style="margin-bottom: 6px;"><strong>Reactants:</strong> <span id="modalReactants"></span></div>
          <div><strong>Products:</strong> <span id="modalProducts"></span></div>
        </div>
      </div>
    </div>
    
    <script>
      const data = {json.dumps(cards_data)};
      const grid = document.getElementById('galleryGrid');
      const searchInput = document.getElementById('search');
      const tabs = document.querySelectorAll('.tab-btn');
    
      const modalOverlay = document.getElementById('modalOverlay');
      const modalClose = document.getElementById('modalClose');
      const modalId = document.getElementById('modalId');
      const modalTitle = document.getElementById('modalTitle');
      const modalCat = document.getElementById('modalCat');
      const modalCond = document.getElementById('modalCond');
      const modalImg = document.getElementById('modalImg');
      const modalDesc = document.getElementById('modalDesc');
      const modalReactants = document.getElementById('modalReactants');
      const modalProducts = document.getElementById('modalProducts');
    
      let currentCategory = 'all';
      let currentQuery = '';
    
      function renderCards() {{
        grid.innerHTML = '';
        const filtered = data.filter(item => {{
          const matchCat = currentCategory === 'all' || item.category === currentCategory;
          const q = currentQuery.toLowerCase();
          const matchQuery = !q || 
            item.title.toLowerCase().includes(q) || 
            item.category.toLowerCase().includes(q) || 
            item.desc.toLowerCase().includes(q) || 
            item.cond.toLowerCase().includes(q) ||
            item.r1_name.toLowerCase().includes(q) ||
            item.p1_name.toLowerCase().includes(q);
          return matchCat && matchQuery;
        }});
    
        filtered.forEach(item => {{
          const card = document.createElement('div');
          card.className = 'rxn-card';
          card.innerHTML = `
            <div class="card-header">
              <span class="card-id">${{item.id}}</span>
              <span class="card-cat">${{item.category}}</span>
            </div>
            <div class="card-body">
              <div class="card-title">${{item.title}}</div>
              <div class="card-cond">${{item.cond}}</div>
              <div class="card-img-container" onclick="openModal('${{item.id}}')">
                <img class="card-img" src="${{item.img}}" alt="${{item.title}}">
              </div>
              <div class="card-desc">${{item.desc}}</div>
            </div>
          `;
          grid.appendChild(card);
        }});
      }}
    
      function openModal(id) {{
        const item = data.find(x => x.id === id);
        if (!item) return;
        modalId.textContent = item.id;
        modalTitle.textContent = item.title;
        modalCat.textContent = item.category;
        modalCond.textContent = item.cond;
        modalImg.src = item.img;
        modalDesc.textContent = item.desc;
        
        let rStr = `${{item.r1_name}} <span class="smiles-tag">${{item.r1_smi}}</span>`;
        if (item.r2_name) rStr += ` + ${{item.r2_name}} <span class="smiles-tag">${{item.r2_smi}}</span>`;
        modalReactants.innerHTML = rStr;
    
        let pStr = `${{item.p1_name}} <span class="smiles-tag">${{item.p1_smi}}</span>`;
        if (item.p2_name) pStr += ` + ${{item.p2_name}} <span class="smiles-tag">${{item.p2_smi}}</span>`;
        modalProducts.innerHTML = pStr;
    
        modalOverlay.classList.add('active');
      }}
    
      modalClose.addEventListener('click', () => modalOverlay.classList.remove('active'));
      modalOverlay.addEventListener('click', (e) => {{
        if (e.target === modalOverlay) modalOverlay.classList.remove('active');
      }});
    
      searchInput.addEventListener('input', (e) => {{
        currentQuery = e.target.value;
        renderCards();
      }});
    
      tabs.forEach(tab => {{
        tab.addEventListener('click', () => {{
          tabs.forEach(t => t.classList.remove('active'));
          tab.classList.add('active');
          currentCategory = tab.dataset.cat;
          renderCards();
        }});
      }});
    
      renderCards();
    </script>
    
    </body>
    </html>
    """

  os.makedirs(os.path.dirname(os.path.abspath(out_html)), exist_ok=True)
  with open(out_html, "w", encoding="utf-8") as f:
    f.write(html_content)

  print(
      f"Gallery written successfully to {out_html} (size: {len(html_content)}"
      " bytes)"
  )
  return out_html


def main():
  parser = argparse.ArgumentParser(
      description="Build reaction mechanism HTML gallery."
  )
  parser.add_argument(
      "--base_dir",
      default=".",
      help="Base directory containing images/ and outputs",
  )
  parser.add_argument(
      "--img_dir",
      default="",
      help="Images directory (defaults to base_dir/images)",
  )
  parser.add_argument(
      "--output",
      default="",
      help="Output HTML file path (defaults to base_dir/gallery.html)",
  )
  args = parser.parse_args()
  build_gallery(
      base_dir=args.base_dir, img_dir=args.img_dir, out_html=args.output
  )


if __name__ == "__main__":
  main()
