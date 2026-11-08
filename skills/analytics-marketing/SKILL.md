---
name: analytics-marketing
description: Marketing analytics, attribution models, KPI tracking, dashboards, and data-driven decision making.
allowed-tools: Read, Glob, Grep
---

# Marketing Analytics

> Transform data into actionable marketing insights.

---

## 1. Analytics Framework

| Layer | Focus | Tools |
|-------|-------|-------|
| **Collection** | Tracking, tagging | GA4, GTM, Segment |
| **Storage** | Data warehouse | BigQuery, Snowflake |
| **Analysis** | Insights, patterns | Looker, Tableau |
| **Action** | Optimization | A/B tools, personalization |

---

## 2. Key Marketing Metrics

### Acquisition Metrics

| Metric | Formula | Benchmark |
|--------|---------|-----------|
| **CAC** | Total acquisition cost / New customers | Industry varies |
| **CPL** | Ad spend / Leads generated | $5-$50 B2B |
| **CPC** | Ad spend / Clicks | $1-$5 avg |
| **CPM** | (Ad spend / Impressions) × 1000 | $5-$15 avg |

### Engagement Metrics

| Metric | Formula | Good |
|--------|---------|------|
| **Bounce rate** | Single-page sessions / Total sessions | < 40% |
| **Pages/session** | Pageviews / Sessions | > 2 |
| **Avg session duration** | Total duration / Sessions | > 2 min |
| **Engagement rate** | Engaged sessions / Total sessions | > 50% |

### Conversion Metrics

| Metric | Formula | Good |
|--------|---------|------|
| **Conversion rate** | Conversions / Visitors | 2-5% |
| **Lead-to-customer** | Customers / Leads | 10-20% |
| **Cart abandonment** | Abandonments / Carts created | < 70% |
| **Checkout abandonment** | Abandonments / Checkouts started | < 50% |

### Revenue Metrics

| Metric | Formula |
|--------|---------|
| **LTV** | ARPU × Customer lifetime |
| **LTV:CAC** | LTV / CAC (target > 3:1) |
| **ROAS** | Revenue / Ad spend |
| **MER** | Revenue / Total marketing spend |

---

## 3. Attribution Models

### Common Models

| Model | How It Works | Best For |
|-------|--------------|----------|
| **Last click** | 100% to final touchpoint | Direct response |
| **First click** | 100% to first touchpoint | Awareness |
| **Linear** | Equal across all touchpoints | Long cycles |
| **Time decay** | More weight to recent | B2B |
| **Position-based** | 40-20-40 first/mid/last | Balanced |
| **Data-driven** | ML-based allocation | High volume |

### Multi-Touch Attribution

| Stage | Touchpoint Credit |
|-------|-------------------|
| **Awareness** | First-touch heavy |
| **Consideration** | Mid-touch matters |
| **Decision** | Last-touch heavy |

---

## 4. Google Analytics 4 (GA4)

### Key Differences from UA

| UA | GA4 |
|----|-----|
| Sessions | Events |
| Bounce rate | Engagement rate |
| Goals | Conversions |
| Views | Data streams |

### Essential Events

| Event | Purpose |
|-------|---------|
| **page_view** | Page tracking |
| **purchase** | Revenue tracking |
| **generate_lead** | Lead capture |
| **sign_up** | Registration |
| **begin_checkout** | Funnel start |
| **add_to_cart** | E-commerce |

### Custom Dimensions

| Dimension | Use Case |
|-----------|----------|
| **user_id** | Cross-device |
| **membership_tier** | Segmentation |
| **content_group** | Content analysis |
| **experiment_id** | A/B testing |

---

## 5. Dashboard Design

### Marketing Dashboard Structure

| Section | Metrics |
|---------|---------|
| **Overview** | Revenue, sessions, conversions |
| **Acquisition** | Traffic sources, CAC, new users |
| **Engagement** | Bounce rate, time on site |
| **Conversion** | Conversion rate, funnel |
| **Revenue** | ROAS, LTV, AOV |

### Visualization Best Practices

| Metric Type | Best Visual |
|-------------|-------------|
| **Trend over time** | Line chart |
| **Comparison** | Bar chart |
| **Composition** | Pie/stacked bar |
| **Distribution** | Histogram |
| **Relationship** | Scatter plot |
| **KPIs** | Big numbers |

---

## 6. Funnel Analysis

### Standard Marketing Funnel

| Stage | Metrics |
|-------|---------|
| **Awareness** | Impressions, reach |
| **Interest** | Clicks, site visits |
| **Consideration** | Lead form starts, cart adds |
| **Intent** | Form completions, checkouts |
| **Purchase** | Transactions, revenue |
| **Loyalty** | Repeat purchases, referrals |

### Funnel Metrics

| Metric | Formula |
|--------|---------|
| **Stage conversion** | Next stage / Current stage |
| **Overall conversion** | Final stage / First stage |
| **Drop-off rate** | 1 - Stage conversion |
| **Funnel velocity** | Avg time between stages |

---

## 7. Cohort Analysis

### Cohort Types

| Type | Example |
|------|---------|
| **Acquisition** | Users by signup month |
| **Behavioral** | Users by first action |
| **Demographic** | Users by segment |

### Retention Cohort Table

| Cohort | Week 0 | Week 1 | Week 2 | Week 3 |
|--------|--------|--------|--------|--------|
| Jan W1 | 100% | 45% | 30% | 25% |
| Jan W2 | 100% | 42% | 28% | 22% |
| Jan W3 | 100% | 48% | 32% | 27% |

---

## 8. Marketing Mix Modeling (MMM)

### Purpose

- Measure channel effectiveness
- Optimize budget allocation
- Account for external factors

### Variables

| Type | Examples |
|------|----------|
| **Media** | TV, digital, OOH |
| **Base** | Brand, distribution |
| **External** | Seasonality, economy |
| **Competition** | Competitor activity |

---

## 9. A/B Test Analysis

### Statistical Significance

| Sample | Min Detectable Effect |
|--------|----------------------|
| 1,000 | 20%+ |
| 10,000 | 5-10% |
| 100,000 | 2-5% |

### Test Metrics

| Metric | Purpose |
|--------|---------|
| **Conversion rate** | Primary metric |
| **Revenue/user** | Revenue impact |
| **Confidence** | Statistical validity |
| **P-value** | < 0.05 for significance |

---

## 10. Reporting Cadence

### Report Types

| Report | Frequency | Audience |
|--------|-----------|----------|
| **Real-time** | Live | Operations |
| **Daily** | Every day | Marketing team |
| **Weekly** | Every week | Marketing leads |
| **Monthly** | Every month | Leadership |
| **Quarterly** | Every quarter | Executive |

### Report Structure

| Section | Content |
|---------|---------|
| **Summary** | Key highlights, changes |
| **Performance** | Metrics vs goals |
| **Insights** | Why it happened |
| **Actions** | What to do next |

---

> **Remember:** Data without action is just noise. Always tie analytics to decisions and outcomes.
