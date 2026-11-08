import os
import shutil
from pathlib import Path

WORKSPACE = Path(r"c:\Users\ntruo\Desktop\Project\IT.Skills\AI Coding Skills Kit")
DEST_REPO = WORKSPACE / "ai-agent-skills-backup"
CLAUDE_SKILLS = Path(r"C:\Users\ntruo\.claude\skills")
GEMINI_PLUGINS = Path(r"C:\Users\ntruo\.gemini\config\plugins")

def ignore_patterns(path, names):
    return [n for n in names if n in ('.git', 'node_modules', '.cache', '__pycache__', '.venv')]

print("=== Starting Comprehensive Backup Sync ===")

# Ensure directories exist
(DEST_REPO / "skills").mkdir(parents=True, exist_ok=True)
(DEST_REPO / "plugins").mkdir(parents=True, exist_ok=True)
(DEST_REPO / "rules").mkdir(parents=True, exist_ok=True)
(DEST_REPO / "agents").mkdir(parents=True, exist_ok=True)
(DEST_REPO / "scripts").mkdir(parents=True, exist_ok=True)

# 1. Sync skills from ~/.claude/skills
claude_count = 0
if CLAUDE_SKILLS.exists():
    for item in CLAUDE_SKILLS.iterdir():
        if item.is_dir() and item.name not in ('.git', 'cache'):
            target = DEST_REPO / "skills" / item.name
            if target.exists():
                shutil.rmtree(target)
            shutil.copytree(item, target, ignore=ignore_patterns)
            claude_count += 1
print(f"Synced {claude_count} skills from ~/.claude/skills")

# 2. Sync skills from ~/.gemini/config/plugins/*/skills
gemini_count = 0
if GEMINI_PLUGINS.exists():
    for plugin_dir in GEMINI_PLUGINS.iterdir():
        if plugin_dir.is_dir():
            skills_sub = plugin_dir / "skills"
            if skills_sub.exists() and skills_sub.is_dir():
                for skill_folder in skills_sub.iterdir():
                    if skill_folder.is_dir() and skill_folder.name not in ('.git', 'cache'):
                        target = DEST_REPO / "skills" / skill_folder.name
                        if target.exists():
                            shutil.rmtree(target)
                        shutil.copytree(skill_folder, target, ignore=ignore_patterns)
                        gemini_count += 1
print(f"Synced skills from ~/.gemini/config/plugins (processed {gemini_count} skill instances)")

total_skills = len([d for d in (DEST_REPO / "skills").iterdir() if d.is_dir()])
print(f"--> Total unique skills in backup now: {total_skills}")

# 3. Sync all plugins to DEST_REPO/plugins
plugin_count = 0
for plugin_dir in GEMINI_PLUGINS.iterdir():
    if plugin_dir.is_dir() and plugin_dir.name not in ('.git', 'cache'):
        target = DEST_REPO / "plugins" / plugin_dir.name
        if target.exists():
            shutil.rmtree(target)
        shutil.copytree(plugin_dir, target, ignore=ignore_patterns)
        plugin_count += 1
print(f"--> Synced {plugin_count} complete plugins into plugins/")

# 4. Sync Agents from maestro-growth-and-fullstack to agents/
maestro_agents = GEMINI_PLUGINS / "maestro-growth-and-fullstack" / "agents"
if maestro_agents.exists():
    for agent_file in maestro_agents.glob("*.md"):
        shutil.copy2(agent_file, DEST_REPO / "agents" / agent_file.name)
    agent_count = len(list((DEST_REPO / "agents").glob("*.md")))
    print(f"--> Synced {agent_count} agent definition files to agents/")

# 5. Sync Rules
src_agents = GEMINI_PLUGINS / "antigravity-kit-plugin" / "rules" / "AGENTS.md"
if src_agents.exists():
    shutil.copy2(src_agents, DEST_REPO / "rules" / "AGENTS.md")
    shutil.copy2(src_agents, DEST_REPO / "AGENTS.md")
    print("--> Updated AGENTS.md in root and rules/")

maestro_gemini = GEMINI_PLUGINS / "maestro-growth-and-fullstack" / "rules" / "GEMINI.md"
if maestro_gemini.exists():
    shutil.copy2(maestro_gemini, DEST_REPO / "rules" / "GEMINI.md")
    print("--> Updated GEMINI.md in rules/")

# Ensure SOUL.md is in both places
if (DEST_REPO / "rules" / "SOUL.md").exists() and not (DEST_REPO / "SOUL.md").exists():
    shutil.copy2(DEST_REPO / "rules" / "SOUL.md", DEST_REPO / "SOUL.md")
elif (DEST_REPO / "SOUL.md").exists() and not (DEST_REPO / "rules" / "SOUL.md").exists():
    shutil.copy2(DEST_REPO / "SOUL.md", DEST_REPO / "rules" / "SOUL.md")

print("=== Backup Synchronization Completed Successfully ===")
