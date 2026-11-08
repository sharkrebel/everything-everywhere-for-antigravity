# Global Workspace Rules & Intent Classifier

> Áp dụng toàn cục cho mọi dự án, hỗ trợ tự nhiên tiếng Việt và tiếng Anh.

## 1. Nguyên tắc nhận diện ngôn ngữ tiếng Việt tự nhiên
- Luôn thấu hiểu mọi yêu cầu bằng tiếng Việt (bao gồm khẩu ngữ, từ lóng kỹ thuật, cách diễn đạt đời thường).
- Không yêu cầu người dùng phải gõ đúng tên skill hay cú pháp tiếng Anh.
- Tự động map ý định người dùng tới các bộ skills tương ứng:
  * "viết code sạch / gọn lại / dễ đọc" -> clean-code
  * "thiết kế bảng / database / csdl" -> database-design
  * "giao diện đẹp / làm lại ui / bảng màu" -> ui-ux-pro-max
  * "vẽ biểu đồ / làm chart / dashboard" -> lieflat-charts
  * "làm cho nó người hơn / khử giọng ai / nhân tính hơn" -> humanizer
  * "tạo file word / xuất word" -> docx
  * "làm slide / thuyết trình" -> pptx
  * "tính toán / xuất excel" -> xlsx
  * "soi lỗi pr / review pr" -> bmad-os-review-pr
  * "bị lỗi / tìm nguyên nhân / tại sao bị vậy" -> problem-solving-pro
  * "phân vân / nên chọn cái nào / so sánh" -> make-decision
  * "lập spec / viết đặc tả / tạo đề xuất thay đổi" -> openspec-propose
  * "thảo luận ý tưởng / khảo sát giải pháp" -> openspec-explore
  * "triển khai spec / code theo spec" -> openspec-apply-change
  * "nghiệm thu spec / kiểm tra khớp spec" -> openspec-verify-change
  * "lưu trữ spec / đóng spec" -> openspec-archive-change
  * "nghiên cứu khoa học / ý tưởng nghiên cứu / brainstorm đề tài" -> scientific-brainstorming
  * "phản biện khoa học / đánh giá bằng chứng / tư duy phản biện" -> scientific-critical-thinking
  * "viết bài báo khoa học / viết paper / bản thảo nghiên cứu" -> scientific-writing
  * "tổng quan tài liệu / tìm paper / tổng thuật y sinh" -> literature-review
  * "thiết kế thí nghiệm / phương pháp nghiên cứu" -> experimental-design
  * "phân tích thống kê khoa học / kiểm định số liệu" -> statistical-analysis
  * "thuật toán mạng xã hội / content từng nền tảng / viral threads" -> platform-content-strategy
  * "trau chuốt ui / làm đẹp giao diện đỉnh cao / thiết kế đẳng cấp / audit ui" -> impeccable
  * "chống slop / thẩm mỹ giao diện / gu thiết kế / landing page đẹp" -> design-taste-frontend
  * "giao diện tối giản / phong cách báo chí / minimalist" -> minimalist-ui
  * "giao diện brutalist / thô mộc / phong cách kỹ thuật" -> industrial-brutalist-ui
  * "viết code đầy đủ / không viết tắt / không bỏ sót code" -> full-output-enforcement
  * "review code đa ngôn ngữ / kiểm tra pr chi tiết" -> code-review-skill
  * "founder review / góc nhìn ceo / đánh giá ý tưởng lớn" -> plan-ceo-review
  * "eng manager review / đánh giá kiến trúc kỹ thuật" -> plan-eng-review
  * "tìm bug / test giao diện / qa hệ thống" -> qa
  * "ship code / ra mắt tính năng / đóng gói release" -> ship
  * "bàn giao phiên bản / handoff công việc" -> handoff
  * "nén context / quản lý ngữ cảnh chủ động" -> strategic-compact
  * "viết unit test / tạo test / kiểm thử hành vi" -> smart-test
  * "audit test / tìm test rác / xóa test vô nghĩa" -> smart-test
  * "tăng tốc CI / blast radius test / ci chậm" -> smart-test
  * "thiết kế logo / tạo logo / vẽ logo / làm thương hiệu" -> logo-design
  * "logo cho app / logo startup / brand identity / nhận diện thương hiệu" -> logo-design
  * "brainstorm tính năng / thảo luận ý tưởng / socratic design" -> brainstorming
  * "debug hệ thống / tìm nguyên nhân gốc rễ bug / fix bug khó" -> systematic-debugging
  * "ủy quyền agent con / chia việc cho subagent / đa agent" -> subagent-driven-development
  * "chạy song song agent / điều phối agent song song" -> dispatching-parallel-agents
  * "kiểm tra trước khi hoàn tất / nghiệm thu bằng chứng / checklist hoàn thành" -> verification-before-completion
  * "dùng git worktree / làm nhiều branch song song" -> using-git-worktrees
  * "lập kế hoạch triển khai / chia nhỏ task / lên roadmap code" -> writing-plans
  * "thực thi plan / code theo kế hoạch từng bước" -> executing-plans
  * "pattern frontend / chuẩn react component / custom hook sạch" -> frontend-patterns
  * "pattern backend / kiến trúc backend / cấu trúc api chuẩn" -> backend-patterns
  * "đánh giá output ai / eval harness / benchmark mô hình" -> eval-harness
  * "vòng lặp kiểm chứng / test-fix tự động lặp lại" -> verification-loop
  * "bộ quy chuẩn thương hiệu / brand guideline" -> brand-guidelines
  * "bảng màu theme / thiết kế theme giao diện" -> theme-factory
  * "tạo web artifact / giao diện html tương tác tức thì" -> web-artifacts-builder
  * "thiết kế web không slop / hallmark / web cao cấp / redesign giao diện" -> hallmark
  * "khởi tạo app từ đầu / app builder / tạo dự án mới" -> app-builder
  * "giải quyết vấn đề sâu / solve deep / phân tích bài toán khó" -> solve.deep
  * "đưa ra quyết định đa khung / so sánh phương án sâu" -> decide.deep
  * "tối ưu tỷ lệ chuyển đổi / cro trang web / tăng conversion" -> conversion-optimization
  * "viết bài quảng cáo / ad creative / mẫu quảng cáo đa kênh" -> ad-creative-variations
  * "viết bài bán hàng / copywriting chốt đơn / landing page copy" -> copywriting
  * "nghiên cứu đối thủ / phân tích đối thủ / competitor teardown" -> competitor-teardown
  * "tối ưu luồng đăng ký / signup flow cro / giảm drop-off" -> signup-flow-cro
  * "chiến dịch marketing tổng thể / lên plan campaign" -> campaign
  * "tối ưu seo programmatic / seo hàng loạt" -> programmatic-seo
  * "kiến trúc nestjs / backend nestjs / viết code nestjs" -> nestjs-expert
  * "lập trình vue / vue 3 composition api / vue router" -> vue-expert
  * "thiết kế schema prisma / query prisma / tối ưu prisma" -> prisma-expert
  * "tạo video bằng code / dựng video react remotion" -> remotion-best-practices
  * "cấu hình docker / tối ưu dockerfile / docker compose" -> docker-expert
  * "quản trị server linux / lệnh bash linux / server management" -> server-management
  * "red team / kiểm thử xâm nhập / pentest bảo mật" -> red-team-tactics

## 2. Tiêu chuẩn phản hồi
- Giữ phong cách ngắn gọn, súc tích, giải quyết thẳng vào vấn đề.
- Ưu tiên hiển thị giải pháp hoặc code mẫu trực tiếp, không vòng vo giải thích thừa.
