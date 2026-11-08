---
name: solve.exec
description: Executive-level problem analysis with maximum detail, all frameworks, and mental models for leadership and decision-makers.
---

# /solve.exec - Problem Solving (Executive Depth)

## Problem Solving (Executive Depth)

1. Read the skill instructions:
// turbo
```
cat .agents/skills/problem-solving-pro/SKILL.md
```

2. Run executive analysis:
// turbo
```
cd .agents/skills/problem-solving-pro/scripts && python search.py "$ARGUMENTS" --plan --depth executive -f markdown
```

3. If user mentions "save", "persist", "step-by-step", "workspace":
// turbo
```
cd .agents/skills/problem-solving-pro/scripts && python search.py "$ARGUMENTS" --plan --depth executive --persist --step-docs -p "$PROJECT_NAME" -f markdown
```

4. Present the output, then append:

```
---
🎯 **Next Steps:**
| Command | Description |
|---------|-------------|
| `/decide.exec` | Executive-level decision from this analysis |
| Add "save step-by-step" | Create full markdown workspace |
| `/decide` | Compare options at standard depth |
```
