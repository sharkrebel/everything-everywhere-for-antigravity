# 🚀 AI Agent Skills Hub & Backup Repository

Kho lưu trữ chuẩn hóa toàn bộ kỹ năng (**Skills**), cấu hình quy tắc (**AGENTS.md**), bản sắc hoạt động (**SOUL.md**), và hệ thống tác nhân đa vai trò (**Agents**) cho các AI Coding Agent hàng đầu: **Google Antigravity (Gemini CLI)** và **Claude Code**.

Được thiết kế theo cấu trúc module tương thích chuẩn mã nguồn mở (tương tự `anthropics/skills` và `fission-ai/openspec`), tích hợp **449+ skills chuyên biệt**, **20 tác nhân chuyên gia (Agents)**, và **13 plugins hoàn chỉnh**, tự động kích hoạt ngầm qua ngôn ngữ tự nhiên tiếng Việt.

---

## 📂 Cấu Trúc Thư Mục (Repository Layout)

```text
ai-agent-skills-backup/
├── AGENTS.md                  # Bản quy tắc toàn cục & bộ phân loại ý định tiếng Việt (Root)
├── SOUL.md                    # Triết lý kiến trúc, linh hồn & chuẩn mực Senior AI Agent (Root)
├── .gitignore                 # Bỏ qua cache, file rác và log tạm
├── README.md                  # Hướng dẫn chi tiết & danh mục kỹ năng
│
├── rules/                     # Thư mục chứa quy tắc hệ thống
│   ├── AGENTS.md              # Global Rules & Natural Intent Classifier (64+ triggers)
│   ├── SOUL.md                # Senior Staff Engineer Persona & Engineering Tenets
│   └── GEMINI.md              # Maestro Orchestrator v4.0 Rules & Protocols
│
├── agents/                    # Hệ thống 20 tác nhân chuyên gia đa vai trò (Maestro Agents)
│   ├── orchestrator.md        # Nhạc trưởng điều phối toàn bộ workflow
│   ├── project-planner.md     # Lập kế hoạch, chia nhỏ roadmap và quản lý scope
│   ├── frontend-specialist.md # Chuyên gia Frontend (React, Vue, Next.js, Tailwind, Remotion)
│   ├── backend-specialist.md  # Chuyên gia Backend (NestJS, Node.js, Python, Prisma)
│   ├── database-architect.md  # Kiến trúc sư CSDL & Indexing
│   ├── devops-engineer.md     # Quản trị Docker, Linux, CI/CD Pipeline
│   ├── debugger.md            # Chuyên gia điều tra và fix bug tận gốc
│   ├── test-engineer.md       # Kỹ sư kiểm thử định hướng hành vi (smart-test)
│   ├── growth-specialist.md   # Tăng trưởng, phễu chuyển đổi & Growth Hacking
│   ├── marketing-strategist.md# Chiến lược tiếp thị, phân tích đối thủ & Campaign
│   ├── content-creator.md     # Sáng tạo nội dung, copywriting & viral hooks
│   ├── seo-specialist.md      # Tối ưu SEO on-page, programmatic SEO
│   ├── security-auditor.md    # Rà soát bảo mật mã nguồn
│   ├── penetration-tester.md  # Red Team & kiểm thử xâm nhập
│   └── ...                    # Mobile, Game, Analytics, Documentation, Performance
│
├── scripts/                   # Bộ công cụ tự động hóa triển khai
│   ├── install-to-claude.ps1      # Cài đặt / đồng bộ toàn bộ skills sang ~/.claude/
│   ├── install-to-antigravity.ps1 # Cài đặt / đồng bộ plugins sang ~/.gemini/config/plugins/
│   └── sync_backup.py             # Script tự động đồng bộ hóa kho sao lưu
│
├── skills/                    # Kho 449+ kỹ năng chuẩn hóa (Mỗi thư mục chứa SKILL.md)
│   ├── hallmark/              # [MỚI] Thiết kế web cao cấp, chống AI slop
│   ├── app-builder/           # [MỚI] Khởi tạo dự án ứng dụng từ đầu
│   ├── solve.deep/ & decide.deep/ # [MỚI] Giải quyết vấn đề sâu & ra quyết định đa khung
│   ├── conversion-optimization/ # [MỚI] Tối ưu hóa tỷ lệ chuyển đổi (CRO)
│   ├── ad-creative-variations/# [MỚI] Viết bài quảng cáo đa kênh
│   ├── copywriting/           # [MỚI] Viết bài bán hàng & landing page copy
│   ├── competitor-teardown/   # [MỚI] Phân tích đối thủ cạnh tranh chuyên sâu
│   ├── signup-flow-cro/       # [MỚI] Tối ưu luồng đăng ký & onboarding
│   ├── campaign/              # [MỚI] Lên kế hoạch chiến dịch marketing tổng thể
│   ├── programmatic-seo/      # [MỚI] Tối ưu SEO tự động hàng loạt
│   ├── nestjs-expert/         # [MỚI] Kiến trúc backend NestJS chuyên sâu
│   ├── vue-expert/            # [MỚI] Lập trình Vue 3 Composition API & Pinia
│   ├── prisma-expert/         # [MỚI] Thiết kế schema & tối ưu truy vấn Prisma
│   ├── remotion-best-practices/# [MỚI] Dựng video lập trình bằng React & Remotion
│   ├── docker-expert/         # [MỚI] Tối ưu Dockerfile & Docker Compose
│   ├── server-management/     # [MỚI] Quản trị Linux server & bash scripting
│   ├── red-team-tactics/      # [MỚI] Pentest bảo mật & kiểm thử xâm nhập
│   ├── smart-test/            # Kỹ nghệ kiểm thử hành vi & tối ưu CI (Blast-radius lanes)
│   ├── subagent-driven-development/ # Điều phối và ủy quyền agent con tự động
│   ├── dispatching-parallel-agents/ # Chạy song song nhiều agent độc lập
│   ├── systematic-debugging/  # Điều tra bug từ gốc rễ theo phương pháp 4 pha
│   ├── verification-before-completion/ # Nghiệm thu chứng cứ trước khi kết thúc task
│   ├── writing-plans/ & executing-plans/ # Lập roadmap & thực thi kế hoạch chia nhỏ
│   ├── using-git-worktrees/   # Quản lý nhiều nhánh song song bằng git worktree
│   ├── impeccable/            # Audit & trau chuốt UI đỉnh cao
│   ├── design-taste-frontend/ # Chống AI slop, nâng tầm thẩm mỹ frontend
│   ├── minimalist-ui/         # Giao diện tối giản phong cách báo chí
│   ├── industrial-brutalist-ui/ # Giao diện thô mộc, kỹ thuật
│   ├── full-output-enforcement/ # Chống viết tắt, bắt buộc code trọn vẹn
│   ├── plan-ceo-review/       # Đánh giá ý tưởng sản phẩm từ góc nhìn Founder/CEO
│   ├── plan-eng-review/       # Đánh giá kiến trúc kỹ thuật (Engineering Manager)
│   ├── qa/ & ship/            # Kiểm thử giao diện & đóng gói phát hành
│   ├── code-review-skill/     # Review mã nguồn cho hơn 25+ ngôn ngữ
│   ├── claude-mem/            # Bộ nhớ dài hạn cross-session
│   ├── strategic-compact/     # Quản lý và nén ngữ cảnh chủ động
│   ├── humanizer/             # Khử giọng văn máy móc của AI, viết tự nhiên
│   ├── ui-ux-pro-max/         # Thư viện giao diện, bảng màu, typography 50+ styles
│   ├── clean-code/            # Quy chuẩn viết code tinh gọn, chuẩn SOLID
│   ├── database-design/       # Thiết kế CSDL, indexing, quan hệ thực thể
│   ├── openspec-*/            # 12 kỹ năng Spec-Driven Change Management
│   └── ...                    # Cloudflare, Vercel, WordPress, Khoa học, Y sinh...
│
└── plugins/                   # Đóng gói 13 Plugins hoàn chỉnh cho Google Antigravity
    ├── maestro-growth-and-fullstack/ # [MỚI] Bộ suite 90 skills & 20 agent roles
    ├── antigravity-kit-plugin/
    ├── developer-power-skills/
    ├── elite-workflows-and-taste/
    ├── openspec-plugin/
    ├── social-media-skills/
    ├── chrome-devtools-plugin/
    ├── android-cli-plugin/
    ├── modern-web-guidance-plugin/
    ├── science/
    ├── scientific-agent-skills/
    ├── firebase/
    └── google-antigravity-sdk/
```

---

## 🎯 Danh Mục Kỹ Năng Tinh Hoa Nổi Bật

### 1. Maestro Orchestrator v4.0 & Đa Tác Nhân (Multi-Agent Suite - MỚI)
- **20 Tác nhân chuyên môn hóa**: Từ Nhạc trưởng (`orchestrator`), Kiến trúc sư backend (`backend-specialist`), Kỹ sư DevOps (`devops-engineer`), Chuyên gia bảo mật (`penetration-tester`) đến Chuyên gia tăng trưởng (`growth-specialist`).
- **Phân rã & Giải quyết sâu**: `solve.deep` (phân tích đa tầng), `decide.deep` (ma trận đánh giá phương án đa khung), `orchestrate` (điều phối luồng đa nhiệm vụ).

### 2. Tăng Trưởng, Tiếp Thị & Tối Ưu Chuyển Đổi (Growth & CRO - MỚI)
- **`conversion-optimization`**: Tối ưu toàn diện tỷ lệ chuyển đổi cho website, landing page và funnel.
- **`signup-flow-cro` & `onboarding-cro`**: Giảm thiểu tỷ lệ drop-off trong quy trình đăng ký và hướng dẫn người dùng mới.
- **`ad-creative-variations` & `copywriting`**: Soạn thảo thông điệp quảng cáo đa kênh và trang bán hàng chốt đơn không dùng từ ngữ sáo rỗng.
- **`competitor-teardown`**: Bóc tách chiến lược sản phẩm, định giá và kênh chuyển đổi của đối thủ cạnh tranh.
- **`programmatic-seo`**: Xây dựng kiến trúc nội dung và tạo hàng loạt trang landing page chuẩn SEO.

### 3. Fullstack Frameworks & Lập Trình Video (Modern Fullstack - MỚI)
- **`nestjs-expert`**: Kiến trúc microservices, modular design, dependency injection và tối ưu API chuẩn NestJS.
- **`vue-expert`**: Lập trình Vue 3 hiện đại với Composition API, Pinia, Vue Router và tối ưu reactivity.
- **`prisma-expert`**: Thiết kế schema CSDL quan hệ, quan hệ n-n, migrations an toàn và tối ưu truy vấn n+1.
- **`remotion-best-practices`**: Dựng video động, infographic chuyển động trực tiếp bằng mã nguồn React & Remotion.
- **`app-builder` & `hallmark`**: Tạo lập dự án hoàn chỉnh từ số 0 với gu thẩm mỹ chuẩn cao cấp, chống AI slop.

### 4. Quản Trị Hệ Thống, DevOps & Bảo Mật (DevOps & Red Team - MỚI)
- **`docker-expert`**: Tối ưu Dockerfile đa tầng (multi-stage builds), giảm kích thước image và cấu hình Docker Compose.
- **`server-management` & `bash-linux`**: Quản trị máy chủ Linux, tự động hóa tác vụ hệ thống bằng bash scripts.
- **`red-team-tactics`**: Rà soát lỗ hổng OWASP, kiểm thử xâm nhập và mô phỏng tấn công an ninh mạng.

### 5. Điều Phối Đa Agent & Quản Trị Nhánh Song Song (Multi-Agent & Worktrees)
- **`subagent-driven-development`**: Phân rã tác vụ lớn thành chuỗi nhiệm vụ độc lập, tự động ủy quyền và giám sát agent con.
- **`dispatching-parallel-agents`**: Chạy đồng thời nhiều agent trên các module khác nhau, tổng hợp kết quả tự động mà không nghẽn context.
- **`using-git-worktrees`**: Tạo và quản lý nhiều workspace độc lập trên cùng một repository, giúp chạy thử nghiệm hoặc code đa tính năng song song mà không xung đột git branch.

### 6. Kỹ Nghệ Debug Chuyên Sâu, Kiểm Thử & Nghiệm Thu Bằng Chứng
- **`systematic-debugging`**: Quy trình 4 pha điều tra bug từ gốc rễ: Cô lập biến ➔ Tái hiện lỗi tối thiểu ➔ Phân tích giả thuyết ➔ Vá lỗi có kiểm chứng.
- **`smart-test`**: Kiểm thử định hướng hành vi (Behavior-First), phân làn vùng ảnh hưởng (Blast-Radius Lanes) và loại bỏ hoàn toàn test rác.
- **`verification-before-completion`**: Chặn đứng thói quen tự mãn của AI bằng checklist bằng chứng bắt buộc trước khi tuyên bố hoàn thành task.

### 7. Thẩm Mỹ Giao Diện & Chống "AI Slop"
- **`impeccable` & `design-taste-frontend`**: Khắc phục các giao diện AI nhạt nhẽo (generic gradient, card vô hồn), áp dụng ngôn ngữ thiết kế độc bản (Minimalist, Industrial Brutalist).
- **`full-output-enforcement`**: Loại bỏ vĩnh viễn thói quen viết tắt `// TODO: tự code tiếp` của AI.

---

## ⚡ Hướng Dẫn Cài Đặt & Sử Dụng

### Dành cho Claude Code:
Chạy script PowerShell để sao chép toàn bộ skills, agents và thiết lập rules:
```powershell
.\scripts\install-to-claude.ps1
```

### Dành cho Google Antigravity (Gemini CLI):
Chạy script PowerShell để kích hoạt toàn bộ 13 Plugins và Intent Classifier:
```powershell
.\scripts\install-to-antigravity.ps1
```

---

## 🗣️ Điều Khiển Tự Nhiên Bằng Tiếng Việt

Nhờ cơ chế nhận diện ý định tại `AGENTS.md`, bạn không cần nhớ tên file hay câu lệnh tiếng Anh:
- *"Thiết kế web không slop / redesign giao diện đẳng cấp"* ➔ Tự động gọi `hallmark`.
- *"Khởi tạo app từ đầu / tạo dự án mới"* ➔ Tự động gọi `app-builder`.
- *"Giải quyết vấn đề sâu / phân tích bài toán khó"* ➔ Tự động gọi `solve.deep`.
- *"Đưa ra quyết định đa khung / so sánh phương án sâu"* ➔ Tự động gọi `decide.deep`.
- *"Tối ưu tỷ lệ chuyển đổi / tăng conversion"* ➔ Tự động gọi `conversion-optimization`.
- *"Viết bài quảng cáo đa kênh / ad creative"* ➔ Tự động gọi `ad-creative-variations`.
- *"Viết bài bán hàng / landing page copy"* ➔ Tự động gọi `copywriting`.
- *"Phân tích đối thủ / competitor teardown"* ➔ Tự động gọi `competitor-teardown`.
- *"Tối ưu luồng đăng ký / giảm drop-off"* ➔ Tự động gọi `signup-flow-cro`.
- *"Lên plan chiến dịch marketing tổng thể"* ➔ Tự động gọi `campaign`.
- *"Tối ưu SEO programmatic hàng loạt"* ➔ Tự động gọi `programmatic-seo`.
- *"Kiến trúc NestJS / viết code backend NestJS"* ➔ Tự động gọi `nestjs-expert`.
- *"Lập trình Vue 3 / Composition API"* ➔ Tự động gọi `vue-expert`.
- *"Thiết kế schema Prisma / tối ưu query Prisma"* ➔ Tự động gọi `prisma-expert`.
- *"Tạo video bằng code / dựng video React Remotion"* ➔ Tự động gọi `remotion-best-practices`.
- *"Cấu hình Docker / tối ưu Dockerfile"* ➔ Tự động gọi `docker-expert`.
- *"Quản trị server Linux / lệnh bash"* ➔ Tự động gọi `server-management`.
- *"Red Team / kiểm thử xâm nhập / pentest"* ➔ Tự động gọi `red-team-tactics`.
- *"Brainstorm tính năng mới / thảo luận ý tưởng"* ➔ Tự động gọi `brainstorming`.
- *"Debug lỗi từ gốc rễ / tìm nguyên nhân bug khó"* ➔ Tự động gọi `systematic-debugging`.
- *"Ủy quyền agent con / chia việc đa agent"* ➔ Tự động gọi `subagent-driven-development`.
- *"Chạy song song nhiều agent"* ➔ Tự động gọi `dispatching-parallel-agents`.
- *"Dùng git worktree / làm nhiều branch song song"* ➔ Tự động gọi `using-git-worktrees`.
- *"Kiểm tra trước khi hoàn tất / nghiệm thu bằng chứng"* ➔ Tự động gọi `verification-before-completion`.
- *"Viết unit test cho tính năng này / kiểm thử hành vi"* ➔ Tự động gọi `smart-test`.
- *"Trau chuốt lại UI trang này cho thật đẳng cấp"* ➔ Tự động gọi `impeccable`.
- *"Nén ngữ cảnh lại để tránh tràn token"* ➔ Tự động gọi `strategic-compact`.
- *"Viết code đầy đủ, không để placeholder"* ➔ Tự động gọi `full-output-enforcement`.
- *"Sửa bài viết này cho tự nhiên, khử mùi AI"* ➔ Tự động gọi `humanizer`.

---
*Kho lưu trữ được cập nhật ngày 08/11/2026 với 449+ kỹ năng chuẩn hóa và 20 tác nhân đa vai trò.*
