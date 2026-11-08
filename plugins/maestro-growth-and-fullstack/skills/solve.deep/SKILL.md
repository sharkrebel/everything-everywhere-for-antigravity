---
name: solve.deep
description: Deep problem analysis with multiple frameworks, alternatives, and detailed mental models for complex or high-stakes problems.
---

# /solve.deep - Problem Solving (Deep Depth)

## Problem Solving (Deep Depth)

1. Read the skill instructions:
// turbo
```
cat .agents/skills/problem-solving-pro/SKILL.md
```

2. Run deep analysis:
// turbo
```
cd .agents/skills/problem-solving-pro/scripts && python search.py "$ARGUMENTS" --plan --depth deep -f markdown
```

3. If user mentions "save", "persist", "step-by-step", "workspace":
// turbo
```
cd .agents/skills/problem-solving-pro/scripts && python search.py "$ARGUMENTS" --plan --depth deep --persist --step-docs -p "$PROJECT_NAME" -f markdown
```

4. Present the output, then append:

```
---
🎯 **Next Steps:**
| Command | Description |
|---------|-------------|
| `/solve.exec` | Executive summary for leadership/stakeholders |
| `/decide.deep` | Detailed comparison of options from this analysis |
| Add "save step-by-step" | Create markdown workspace for each step |
```
