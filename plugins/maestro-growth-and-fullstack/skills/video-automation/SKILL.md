---
name: video-automation
description: Automated video production using Remotion, combining Antigravity aesthetics with official Remotion best practices.
allowed-tools: Read, RunCommand, Write
---

# Video Automation - The Full Remotion Suite

> This skill combines the **[remotion-best-practices](../remotion-best-practices/SKILL.md)** skill with the Antigravity design system, enabling production-grade video ads.

---

## 📚 Required Reading

For detailed technical implementation, refer to the **`remotion-best-practices`** skill which contains 16 comprehensive rule files:

| Core Rules | Description |
|------------|-------------|
| `animations.md` | Frame-based animations with `useCurrentFrame()` |
| `timing.md` | Springs, easing, interpolation curves |
| `sequencing.md` | `<Sequence>` and `<Series>` patterns |
| `transitions.md` | Scene transitions with `@remotion/transitions` |
| `audio.md` | Audio import, trimming, volume control |
| `videos.md` | Video embedding, looping, speed control |
| `display-captions.md` | TikTok-style captions with word highlighting |
| `charts.md` | Animated bar/pie charts |
| `3d.md` | Three.js integration with `@remotion/three` |
| `fonts.md` | Google Fonts & Vietnamese font support |

---

## 1. Antigravity Core Aesthetics (The Look)

Apply these design tokens to all Remotion components:

```tsx
const ANTIGRAVITY_TOKENS = {
  background: '#0a0a0c',
  primary: '#00f2ff',    // Neon Cyan
  secondary: '#ff00e5',  // Neon Pink
  glass: 'rgba(255, 255, 255, 0.05)',
  fonts: ['Inter', 'Be Vietnam Pro'],
};
```

**Key styling patterns:**
- **Glassmorphism**: `backdropFilter: 'blur(20px)'` with neon borders
- **Staggered Motion**: Map array indices to delays (`i * 10`)
- **Micro-animations**: Use `Math.sin(frame/10)` for pulsing effects

---

## 2. Usage Patterns

### High-Fidelity Landing Page Cloning
1. **Research**: Use `browser_subagent` to extract brand colors, images, headlines
2. **Implementation**: Use `interpolate`, `spring`, `Sequence` to animate sections
3. **Aspect Ratios**: Always offer 16:9 (YouTube) and 9:16 (TikTok/Reels)

### TikTok-Style Captions
```tsx
// Refer to display-captions.md in remotion-best-practices
import { createTikTokStyleCaptions } from '@remotion/captions';
```

---

## 3. Quick Commands

For non-tech users, just 3 commands:

| Command | Purpose |
|---------|---------|
| `npx remotion preview` | Preview in browser |
| `npx remotion render <id> out.mp4` | Export final video |
| `--port=300x` | Run multiple previews |

---

## 4. Troubleshooting

### FFmpeg Installation (macOS)
```bash
brew install ffmpeg
```

### Headless Chrome Error
```bash
npx remotion render <id> --browser-executable="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
```

---

## Credits

- **Technical Rules**: Based on [remotion-dev/skills](https://github.com/remotion-dev/skills)
- **Design System**: Antigravity Marketing Kit aesthetics
