---
name: decide
description: Make a decision with standard analysis depth when comparing alternatives, weighing trade-offs, or evaluating pros and cons.
---

# /decide - Decision Making (Standard Depth)

## Decision Making (Standard Depth)

1. Read the skill instructions:
// turbo
```
cat .agents/skills/make-decision/SKILL.md
```

2. Run the analysis:
// turbo
```
cd .agents/skills/make-decision/scripts && python search.py "$ARGUMENTS" --plan -f markdown
```

3. If user mentions "save", "persist", "step-by-step", "workspace":
// turbo
```
cd .agents/skills/make-decision/scripts && python search.py "$ARGUMENTS" --plan --persist --step-docs -p "$PROJECT_NAME" -f markdown
```

4. Present the output, then append:

```
---
🎯 **Next Steps:**
| Command | Description |
|---------|-------------|
| `/decide.deep` | More detailed comparison with additional criteria & frameworks |
| `/decide.exec` | Executive-level analysis for board/leadership |
| `/solve` | Deep-dive into the problem before deciding |
```
