---
name: campaign
description: Create a comprehensive marketing campaign from strategy and content creation to launch, optimization, HTML presentation, and automated video production.
---

# /campaign - Create Marketing Campaign

Follow this workflow to create a marketing campaign.

## Step 1: Campaign Brief

Define the campaign basics:

1. **Goal**: What are we trying to achieve?
   - Awareness
   - Lead generation
   - Sales/conversions
   - Product launch
   - Brand building

2. **Target Audience**: Who are we reaching?
   - Demographics
   - Psychographics
   - Pain points
   - Where they hang out

3. **Budget & Timeline**: Resources and deadlines

4. **Success Metrics**: How will we measure success?

## Step 2: Strategy Development

1. **Positioning**: How do we want to be perceived?
2. **Messaging**: Key messages and value proposition
3. **Channels**: Which channels to use
   - Paid (Google, Meta, LinkedIn)
   - Organic (SEO, social, email)
   - Earned (PR, influencers)
4. **Creative direction**: Visual and copy themes

## Step 3: Content Creation

Create campaign assets:

1. **Landing pages**: Design and copy
2. **Ad creative**: Images, videos, copy
3. **Email sequences**: Nurture flows
4. **Social content**: Platform-specific posts
5. **Blog/SEO content**: Supporting content

## Step 4: Setup & Launch

1. **Tracking**: UTMs, pixels, conversion events
2. **Audience setup**: Targeting configuration
3. **A/B tests**: Variations to test
4. **Launch checklist**: Pre-flight checks
5. **Go live**: Execute launch

## Step 5: Optimize

1. **Monitor**: Daily/weekly check-ins
2. **Analyze**: Performance vs goals
3. **Iterate**: Adjust based on data
4. **Report**: Share learnings

## Output

- Campaign brief document
- Content calendar
- Creative assets
- Tracking setup
- Results report

## Step 6: Presentation Builder (Optional)

Ask the user: **"Bạn có muốn tạo file thuyết trình HTML cho campaign này không?"**

If YES, follow these steps:

### 6.1 Build HTML Presentation

Use the `ui-ux-pro-max` skill to create a stunning HTML presentation:

1. Create a single-file HTML presentation with:
   - Modern dark theme with gradient accents
   - Slide-based navigation (arrow keys + click)
   - Progress bar
   - Touch/swipe support for mobile
   - All CSS inline (no external dependencies)
   - **Font selection**: Use Vietnamese-friendly fonts like `Inter`, `Be Vietnam Pro`, or `Roboto` if the content is in Vietnamese to avoid broken characters.

2. Content slides should include:
   - Title slide with campaign name
   - Goals & KPIs
   - Target audience
   - Key messaging
   - Marketing channels
   - Timeline
   - Content plan
   - Success metrics
   - CTA / Next steps

3. Save to `docs/campaign-presentation.html`

### 6.2 Preview Locally

Open the presentation in browser:

```bash
# macOS
open docs/campaign-presentation.html

# Linux
xdg-open docs/campaign-presentation.html

# Windows
start docs/campaign-presentation.html
```

### 6.3 Deploy to Internet (Optional)

Ask the user: **"Bạn có muốn đẩy lên internet để chia sẻ không?"**

If YES, use nport to create a public URL:

```bash
# Install nport if not installed
npm install -g nport

# Start a local server and expose it
npx serve docs -l 3000 &
nport 3000 -s <campaign-name>
```

This will create a public URL like: `https://<campaign-name>.nport.link`

**Note**: The tunnel remains active for ~4 hours. Perfect for sharing presentations in meetings!

## Step 7: Automated Video Production (Optional)

Ask the user: **"Bạn có muốn tạo video ad tự động cho campaign này không? (TikTok/Reels/Meta)"**
## Bước 7: Sản xuất Video Tự động (Tùy chọn) ✨
Dành cho người dùng không rành kỹ thuật, Agent sẽ lo toàn bộ khâu code.

### 7.1 Cách yêu cầu Agent
Mầy chỉ cần đưa yêu cầu bằng ngôn ngữ tự nhiên. Ví dụ:
- *"Clone cái landing page này thành video Facebook: [URL]"*
- *"Tạo một video TikTok giới thiệu tính năng X, phong cách neon, có hoa rơi."*
- *"Làm một bản 16:9 và một bản 9:16 cho chiến dịch Tết."*

### 7.2 Xem Preview (Cực dễ)
Sau khi Agent code xong, mầy chỉ cần mở link Localhost mà Agent cung cấp (thường là `http://localhost:3000`).
- Nhấn **Space** để xem.
- Kéo thanh Timeline để soi từng frame.

### 7.3 Xuất Video (Render)
Khi đã ưng ý, mầy chỉ cần bảo Agent: *"Xuất mp4 cho tao"*. 
- Agent sẽ tự động kiểm tra FFmpeg.
- Nếu thiếu, Agent sẽ hướng dẫn mầy chạy: `brew install ffmpeg`.
- Sau đó Agent sẽ tự chạy lệnh render và mở file video lên cho mầy.

**Composition ID mầy cần biết:**
- Video Facebook (16:9): `JupViecTetAd`
- Video TikTok/Reels (9:16): `JupViecTetAd916` hoặc `JupViecTikTokAd`

## Skills Used

- social-media-expert
- email-marketing
- analytics-marketing
- ui-ux-pro-max
- video-automation
