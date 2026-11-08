---
name: solve.quick
description: Quick problem scan that delivers fast, essential insights and root-cause overview when the user needs a rapid answer.
---

# /solve.quick - Problem Solving (Quick Depth)

## Problem Solving (Quick Depth)

1. Read the skill instructions:
// turbo
```
cat .agents/skills/problem-solving-pro/SKILL.md
```

2. Run quick analysis:
// turbo
```
cd .agents/skills/problem-solving-pro/scripts && python search.py "$ARGUMENTS" --plan --depth quick -f markdown
```

3. Present the output, then append:

```
---
🎯 **Next Steps:**
| Command | Description |
|---------|-------------|
| `/solve` | Full standard analysis |
| `/solve.deep` | Deep dive with alternatives & mental models |
| `/decide.quick` | Quick decision from this analysis |
```
